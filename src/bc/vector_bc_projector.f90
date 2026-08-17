! Copyright (c) 2026, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Implements boundary condition projectors for vector fields. Two types concrete
!! types are provided: `segregated_vector_bc_projector_t` for simple
!! component-wise resolution, and `coupled_vector_bc_projector_t` for resolving
!! mixed boundary conditions.
module vector_bc_projector
  use bc, only : bc_t, BC_DIRICHLET
  use bc_list, only : bc_list_t
  use mixed_bc, only : mixed_bc_t
  use mask, only : mask_t
  use coefs, only : coef_t
  use dofmap, only : dofmap_t
  use field, only : field_t
  use field_list, only : field_list_t
  use fld_file, only : fld_file_t
  use hex, only : edge_nodes, edge_faces, node_faces
  use htable, only : htable_i4_t
  use logger, only : LOG_SIZE
  use matrix, only : matrix_t
  use scratch_registry, only : neko_scratch_registry
  use math, only : cfill_mask, masked_scatter_copy, rzero, cfill
  use gs_ops, only : GS_OP_ADD, GS_OP_MIN
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use tuple, only : tuple_i4_t
  use utils, only : neko_error, nonlinear_index, linear_index
  use device, only : device_get_ptr, device_memcpy, device_map, device_unmap, &
       HOST_TO_DEVICE, DEVICE_TO_HOST, glb_cmd_queue
  use device_math, only : device_cfill_mask
  use device_coupled_vector_bc_projector, only : &
       device_coupled_vector_bc_projector_apply
  use scalar_bc_projector, only : scalar_bc_projector_t
  use operators, only : rotate_cyc
  use, intrinsic :: iso_c_binding, only : c_ptr, c_null_ptr, c_associated
  implicit none
  private

  public :: vector_bc_projector_components

  !> Abstract type for resolving vector boundary conditions.
  type, public, abstract :: vector_bc_projector_t
   contains
     !> Constructor.
     procedure(vector_bc_projector_init_intrf), pass(this), deferred :: init
     !> Destructor.
     procedure(vector_bc_projector_free_intrf), pass(this), deferred :: free
     !> Finalize the vector boundary-condition projector.
     procedure(vector_bc_projector_finalize_intrf), pass(this), deferred :: &
          finalize
     !> Apply homogeneous resolved vector boundary constraints.
     procedure(vector_bc_projector_apply_intrf), pass(this), deferred :: apply
     !> Mark a boundary condition in the vector boundary-condition projector.
     procedure(vector_bc_projector_mark_bc_intrf), pass(this), deferred :: &
          mark_bc
     !> Mark a boundary condition on one Cartesian component.
     procedure, pass(this) :: mark_bc_component => &
          vector_bc_projector_mark_bc_component
     !> Mark a list of boundary conditions in the projector.
     procedure(vector_bc_projector_mark_bc_list_intrf), pass(this), deferred :: &
          mark_bc_list
     !> Mark a list of boundary conditions on one Cartesian component.
     procedure, pass(this) :: mark_bc_list_component => &
          vector_bc_projector_mark_bc_list_component
     generic :: mark => mark_bc, mark_bc_component, mark_bc_list, &
          mark_bc_list_component
  end type vector_bc_projector_t

  !> A projector for vector fields that acts component-wise.
  !! @details Stores one `scalar_bc_projector_t` per global component and applies
  !! vector boundary conditions independently to `x`, `y` and `z`. This is
  !! suitable only for fully component-wise constraints. Like the underlying
  !! scalar projectors, this does not perform any real global resolution, we just
  !! feed it the Dirichlet dofs component-wise. Therefore, `finalize` is no-op.
  type, public, extends(vector_bc_projector_t) :: segregated_vector_bc_projector_t
     !> Projector for the x-component.
     type(scalar_bc_projector_t) :: x
     !> Projector for the y-component.
     type(scalar_bc_projector_t) :: y
     !> Projector for the z-component.
     type(scalar_bc_projector_t) :: z
   contains
     !> Constructor, no-op.
     procedure, pass(this) :: init => segregated_vector_bc_projector_init
     !> Destructor.
     procedure, pass(this) :: free => segregated_vector_bc_projector_free
     !> No-op finalize routine.
     procedure, pass(this) :: finalize => segregated_vector_bc_projector_finalize
     !> Mark a boundary condition.
     procedure, pass(this) :: mark_bc => segregated_vector_bc_projector_mark_bc
     !> Mark a boundary condition on one Cartesian component.
     procedure, pass(this) :: mark_bc_component => &
          segregated_vector_bc_projector_mark_bc_component
     !> Mark a list of boundary conditions.
     procedure, pass(this) :: mark_bc_list => &
          segregated_vector_bc_projector_mark_bc_list
     !> Mark a list of boundary conditions on one Cartesian component.
     procedure, pass(this) :: mark_bc_list_component => &
          segregated_vector_bc_projector_mark_bc_list_component
     !> Zero-out constrained dofs using the underlying component projectors.
     procedure, pass(this) :: apply => segregated_vector_bc_projector_apply
  end type segregated_vector_bc_projector_t

  !> A coupled projector for vector fields, suitable for mixed boundary
  !! conditions.
  !! @details This projector first accumulates the marked vector boundary
  !! conditions, then resolves them onto the global velocity dofs in
  !! `finalize()`. The resolved support is split into two disjoint parts:
  !! plain Dirichlet nodes, stored in `dirichlet_dof_mask`, and genuinely
  !! mixed nodes, stored in `mixed_dof_mask`. For the mixed nodes, the projector
  !! also builds a local orthonormal basis `(n, t1, t2)` together with the
  !! resolved local constraint flags `constraint_n`, `constraint_t1`, and
  !! `constraint_t2`. The `apply()` routine then projects the vector field to
  !! this local basis, zeroes the constrained components, and reconstructs the
  !! Cartesian vector.
  !! Additionally, the projector propagates the mixed masks and basis data to
  !! the boundary conditions of class `mixed_bc_t`, which is necessary for them
  !! to work correctly.
  type, public, extends(vector_bc_projector_t) :: coupled_vector_bc_projector_t
     !> DOFs that are fully constrained in Cartesian space.
     type(mask_t) :: dirichlet_dof_mask
     !> DOFs that require basis-aware mixed treatment.
     type(mask_t) :: mixed_dof_mask
     !> Boundary conditions queued for resolution during `finalize()`.
     type(bc_list_t), private :: bcs
     !> SEM coefficients.
     type(coef_t), pointer, private :: coef => null()
     !> Degree-of-freedom map.
     type(dofmap_t), pointer, private :: dof => null()
     !> Per-boundary-DOF node classification using the `BC_*` constants.
     real(kind=rp), allocatable :: node_type(:)
     !> Per-boundary-face classification using the `BC_*` constants.
     real(kind=rp), allocatable :: face_type(:,:)
     !> Linear indices of boundary DOFs.
     integer, allocatable :: boundary_dof(:)
     !> Hash table from a linear DOF index to `boundary_dof`.
     type(htable_i4_t) :: boundary_idx
     !> Reference-space coordinates of the element corner nodes.
     integer, allocatable :: node_rst(:,:)
     !> Reference-space coordinates of the edge-midpoint nodes.
     integer, allocatable :: edge_mid_rst(:,:)
     !> Linear indices of the element corner nodes.
     integer, allocatable :: node_linear_idx(:)
     !> Local normal-component constraint flags on `mixed_dof_mask`.
     integer, allocatable :: constraint_n(:)
     !> Local first-tangent constraint flags on `mixed_dof_mask`.
     integer, allocatable :: constraint_t1(:)
     !> Local second-tangent constraint flags on `mixed_dof_mask`.
     integer, allocatable :: constraint_t2(:)
     !> Local normal basis vectors for the mixed-node support.
     type(matrix_t) :: n
     !> Local first tangential basis vectors for the mixed-node support.
     type(matrix_t) :: t1
     !> Local second tangential basis vectors for the mixed-node support.
     type(matrix_t) :: t2

     !> Device mirror of `constraint_n`.
     type(c_ptr) :: constraint_n_d = c_null_ptr
     !> Device mirror of `constraint_t1`.
     type(c_ptr) :: constraint_t1_d = c_null_ptr
     !> Device mirror of `constraint_t2`.
     type(c_ptr) :: constraint_t2_d = c_null_ptr
   contains
     !> Constructor.
     procedure, pass(this) :: init => coupled_vector_bc_projector_init
     !> Destructor.
     procedure, pass(this) :: free => coupled_vector_bc_projector_free
     !> Mark a boundary condition for later resolution.
     procedure, pass(this) :: mark_bc => coupled_vector_bc_projector_mark_bc
     !> Mark a list of boundary conditions for later resolution.
     procedure, pass(this) :: mark_bc_list => &
          coupled_vector_bc_projector_mark_bc_list
     !> Resolve the queued boundary conditions into masks and basis data.
     procedure, pass(this) :: finalize => coupled_vector_bc_projector_finalize
     !> Apply the resolved homogeneous mixed constraints to a vector field.
     procedure, pass(this) :: apply => coupled_vector_bc_projector_apply
     !> Write diagnostic output for the resolved masks and basis.
     procedure, pass(this) :: debug_output => &
          coupled_vector_bc_projector_debug_output
     !> Write debug output of normal component on mixed BC nodes.
     procedure, pass(this) :: debug_output_normal_component => &
          coupled_vector_bc_projector_debug_output_normal_component
     !> Clear the resolved mask-side state.
     procedure, pass(this), private :: clear_masks => &
          coupled_vector_bc_projector_clear_masks
     !> Rebuild the Dirichlet and mixed-node masks from the queued BCs.
     procedure, pass(this), private :: rebuild_masks => &
          coupled_vector_bc_projector_rebuild_masks
     !> Clear the resolved basis-side state.
     procedure, pass(this), private :: clear_basis => &
          coupled_vector_bc_projector_clear_basis
     !> Rebuild the local basis on the resolved mixed-node support.
     procedure, pass(this), private :: rebuild_basis => &
          coupled_vector_bc_projector_rebuild_basis
  end type coupled_vector_bc_projector_t

  abstract interface
     !> Free the vector boundary-condition projector.
     subroutine vector_bc_projector_free_intrf(this)
       import :: vector_bc_projector_t
       class(vector_bc_projector_t), intent(inout) :: this
     end subroutine vector_bc_projector_free_intrf
  end interface

  abstract interface
     !> Initialize the vector boundary-condition projector.
     !! @param[in] coef SEM coefficients defining the dof layout.
     subroutine vector_bc_projector_init_intrf(this, coef)
       import :: vector_bc_projector_t, coef_t
       class(vector_bc_projector_t), intent(inout) :: this
       type(coef_t), target, intent(in) :: coef
     end subroutine vector_bc_projector_init_intrf
  end interface

  abstract interface
     !> Finalize the vector boundary-condition projector.
     subroutine vector_bc_projector_finalize_intrf(this, rebuild_mask)
       import :: vector_bc_projector_t
       class(vector_bc_projector_t), intent(inout) :: this
       logical, intent(in) :: rebuild_mask
     end subroutine vector_bc_projector_finalize_intrf
  end interface

  abstract interface
     !> Apply the resolved vector boundary constraints.
     !! @param[inout] x x-component field values.
     !! @param[inout] y y-component field values.
     !! @param[inout] z z-component field values.
     !! @param[in] n Number of entries in each component array.
     !! @param[inout] strm Optional backend stream used on device backends.
     subroutine vector_bc_projector_apply_intrf(this, x, y, z, n, strm)
       import :: vector_bc_projector_t, rp, c_ptr
       class(vector_bc_projector_t), intent(in) :: this
       integer, intent(in) :: n
       real(kind=rp), intent(inout) :: x(n)
       real(kind=rp), intent(inout) :: y(n)
       real(kind=rp), intent(inout) :: z(n)
       type(c_ptr), intent(inout), optional :: strm
     end subroutine vector_bc_projector_apply_intrf
  end interface

  abstract interface
     !> Mark a boundary condition in the vector boundary-condition projector.
     !! @param[inout] bc Boundary condition to register.
     subroutine vector_bc_projector_mark_bc_intrf(this, bc)
       import :: vector_bc_projector_t, bc_t
       class(vector_bc_projector_t), intent(inout) :: this
       class(bc_t), intent(inout), target :: bc
     end subroutine vector_bc_projector_mark_bc_intrf
  end interface

  abstract interface
     !> Mark a list of boundary conditions in the projector.
     !! @param[in] bclst Boundary-condition list to register.
     subroutine vector_bc_projector_mark_bc_list_intrf(this, bclst)
       import :: vector_bc_projector_t, bc_list_t
       class(vector_bc_projector_t), intent(inout) :: this
       type(bc_list_t), intent(in) :: bclst
     end subroutine vector_bc_projector_mark_bc_list_intrf
  end interface

contains

  !> Report unsupported component-specific marking.
  subroutine vector_bc_projector_mark_bc_component(this, bc, &
       component)
    class(vector_bc_projector_t), intent(inout) :: this
    class(bc_t), intent(inout), target :: bc
    character(len=1), intent(in) :: component

    call neko_error("Component-specific marking is only supported by " // &
         "segregated vector BC projectors.")
  end subroutine vector_bc_projector_mark_bc_component

  !> Report unsupported component-specific list marking.
  subroutine vector_bc_projector_mark_bc_list_component(this, &
       bclst, component)
    class(vector_bc_projector_t), intent(inout) :: this
    type(bc_list_t), intent(in) :: bclst
    character(len=1), intent(in) :: component

    call neko_error("Component-specific marking is only supported by " // &
         "segregated vector BC projectors.")
  end subroutine vector_bc_projector_mark_bc_list_component

  !> Destructor
  subroutine segregated_vector_bc_projector_free(this)
    class(segregated_vector_bc_projector_t), intent(inout) :: this
    call this%x%free()
    call this%y%free()
    call this%z%free()
  end subroutine segregated_vector_bc_projector_free

  !> Constructor.
  !! @param[in] coef SEM coefficients defining the dof layout.
  subroutine segregated_vector_bc_projector_init(this, coef)
    class(segregated_vector_bc_projector_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef

    ! Just ensure we get a fresh object
    call this%free()
  end subroutine segregated_vector_bc_projector_init

  !> Finalize the segregated vector boundary-condition projector.
  subroutine segregated_vector_bc_projector_finalize(this, rebuild_mask)
    class(segregated_vector_bc_projector_t), intent(inout) :: this
    logical, intent(in) :: rebuild_mask
  end subroutine segregated_vector_bc_projector_finalize

  !> Mark a boundary condition in the segregated projector.
  !! @details Fully constrained vector boundary conditions are applied to all
  !! three component projectors.
  !! @param[inout] bc Boundary condition to register.
  subroutine segregated_vector_bc_projector_mark_bc(this, bc)
    class(segregated_vector_bc_projector_t), intent(inout) :: this
    class(bc_t), intent(inout), target :: bc

    if (bc%bc_type .ne. BC_DIRICHLET) then
       call neko_error("Segregated vector BC projector only accepts " // &
            "Dirichlet boundary conditions.")
    end if

    call this%x%mark_bc(bc)
    call this%y%mark_bc(bc)
    call this%z%mark_bc(bc)
  end subroutine segregated_vector_bc_projector_mark_bc

  !> Mark a boundary condition on one Cartesian component.
  !! @param[inout] bc Boundary condition to register.
  !! @param[in] component Component selector.
  subroutine segregated_vector_bc_projector_mark_bc_component(this, bc, &
       component)
    class(segregated_vector_bc_projector_t), intent(inout) :: this
    class(bc_t), intent(inout), target :: bc
    character(len=1), intent(in) :: component

    if (bc%bc_type .ne. BC_DIRICHLET) then
       call neko_error("Segregated vector BC projector only accepts " // &
            "Dirichlet boundary conditions.")
    end if

    select case (component)
    case ('x')
       call this%x%mark_bc(bc)
    case ('y')
       call this%y%mark_bc(bc)
    case ('z')
       call this%z%mark_bc(bc)
    case default
       call neko_error("Invalid component for segregated vector BC " // &
            "projector mark.")
    end select
  end subroutine segregated_vector_bc_projector_mark_bc_component

  !> Mark a list of boundary conditions in the segregated projector.
  !! @param[in] bclst Boundary-condition list to register.
  subroutine segregated_vector_bc_projector_mark_bc_list(this, bclst)
    class(segregated_vector_bc_projector_t), intent(inout) :: this
    type(bc_list_t), intent(in) :: bclst
    class(bc_t), pointer :: bc_i
    integer :: i

    do i = 1, bclst%size()
       bc_i => bclst%get(i)
       call this%mark(bc_i)
    end do
  end subroutine segregated_vector_bc_projector_mark_bc_list

  !> Mark a list of boundary conditions on one Cartesian component.
  !! @param[in] bclst Boundary-condition list to register.
  !! @param[in] component Component selector.
  subroutine segregated_vector_bc_projector_mark_bc_list_component(this, &
       bclst, component)
    class(segregated_vector_bc_projector_t), intent(inout) :: this
    type(bc_list_t), intent(in) :: bclst
    character(len=1), intent(in) :: component
    class(bc_t), pointer :: bc_i
    integer :: i

    do i = 1, bclst%size()
       bc_i => bclst%get(i)
       call this%mark(bc_i, component)
    end do
  end subroutine segregated_vector_bc_projector_mark_bc_list_component

  !> Apply the segregated vector boundary constraints.
  !! @param[inout] x x-component field values.
  !! @param[inout] y y-component field values.
  !! @param[inout] z z-component field values.
  !! @param[in] n Number of entries in each component array.
  !! @param[inout] strm Optional backend stream/queue used on device backends.
  subroutine segregated_vector_bc_projector_apply(this, x, y, z, n, strm)
    class(segregated_vector_bc_projector_t), intent(in) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: x(n)
    real(kind=rp), intent(inout) :: y(n)
    real(kind=rp), intent(inout) :: z(n)
    type(c_ptr), intent(inout), optional :: strm

    call this%x%apply(x, n, strm = strm)
    call this%y%apply(y, n, strm = strm)
    call this%z%apply(z, n, strm = strm)
  end subroutine segregated_vector_bc_projector_apply

  !> Access the component scalar projectors from a segregated vector projector.
  !! @param[inout] x Pointer to the x-component scalar projector.
  !! @param[inout] y Pointer to the y-component scalar projector.
  !! @param[inout] z Pointer to the z-component scalar projector.
  subroutine vector_bc_projector_components(this, x, y, z)
    class(vector_bc_projector_t), target, intent(inout) :: this
    type(scalar_bc_projector_t), pointer, intent(inout) :: x
    type(scalar_bc_projector_t), pointer, intent(inout) :: y
    type(scalar_bc_projector_t), pointer, intent(inout) :: z

    select type (this)
    type is (segregated_vector_bc_projector_t)
       x => this%x
       y => this%y
       z => this%z
    class default
       call neko_error("Component access is only available for " // &
            "segregated vector BC projectors. You have likely forgotten to " // &
            "select a coupled linear solver for velocity in the fluid " // &
            "configuration.")
    end select
  end subroutine vector_bc_projector_components


  !
  ! Coupled projector TBPs
  !

  !> Destructor.
  subroutine coupled_vector_bc_projector_free(this)
    class(coupled_vector_bc_projector_t), intent(inout) :: this

    call this%bcs%free()
    if (allocated(this%node_rst)) deallocate(this%node_rst)
    if (allocated(this%edge_mid_rst)) deallocate(this%edge_mid_rst)
    if (allocated(this%node_linear_idx)) deallocate(this%node_linear_idx)
    if (allocated(this%face_type)) deallocate(this%face_type)
    call this%clear_masks()
    call this%clear_basis()
    this%constraint_n_d = c_null_ptr
    this%constraint_t1_d = c_null_ptr
    this%constraint_t2_d = c_null_ptr
    nullify(this%coef)
    nullify(this%dof)
  end subroutine coupled_vector_bc_projector_free

  !> Constructor.
  !! @param[in] coef SEM coefficients.
  !! @details This resets any previous resolved state, stores pointers to the
  !! coefficient and dofmap objects, and precomputes reference-element lookup
  !! tables used later to reconstruct mixed-node normals on edges and corners.
  subroutine coupled_vector_bc_projector_init(this, coef)
    class(coupled_vector_bc_projector_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    integer :: lx, ly, lz
    integer :: mid_i, mid_j, mid_k
    integer :: nface

    call this%free()
    call this%bcs%init()

    this%coef => coef
    this%dof => coef%dof

    lx = coef%Xh%lx
    ly = coef%Xh%ly
    lz = coef%Xh%lz
    mid_i = (lx + 1) / 2
    mid_j = (ly + 1) / 2
    mid_k = (lz + 1) / 2

    allocate(this%node_rst(3, 8))
    allocate(this%edge_mid_rst(3, 12))
    allocate(this%node_linear_idx(8))
    nface = 2 * coef%msh%gdim
    allocate(this%face_type(nface, coef%msh%nelv))
    this%face_type = 5.0_rp

    this%node_rst(:,1) = [1, 1, 1]
    this%node_rst(:,2) = [lx, 1, 1]
    this%node_rst(:,3) = [1, ly, 1]
    this%node_rst(:,4) = [lx, ly, 1]
    this%node_rst(:,5) = [1, 1, lz]
    this%node_rst(:,6) = [lx, 1, lz]
    this%node_rst(:,7) = [1, ly, lz]
    this%node_rst(:,8) = [lx, ly, lz]

    this%edge_mid_rst(:,1) = [mid_i, 1, 1]
    this%edge_mid_rst(:,2) = [mid_i, ly, 1]
    this%edge_mid_rst(:,3) = [mid_i, 1, lz]
    this%edge_mid_rst(:,4) = [mid_i, ly, lz]
    this%edge_mid_rst(:,5) = [1, mid_j, 1]
    this%edge_mid_rst(:,6) = [lx, mid_j, 1]
    this%edge_mid_rst(:,7) = [1, mid_j, lz]
    this%edge_mid_rst(:,8) = [lx, mid_j, lz]
    this%edge_mid_rst(:,9) = [1, 1, mid_k]
    this%edge_mid_rst(:,10) = [lx, 1, mid_k]
    this%edge_mid_rst(:,11) = [1, ly, mid_k]
    this%edge_mid_rst(:,12) = [lx, ly, mid_k]

    this%node_linear_idx(1) = linear_index(1, 1, 1, 1, lx, ly, lz)
    this%node_linear_idx(2) = linear_index(lx, 1, 1, 1, lx, ly, lz)
    this%node_linear_idx(3) = linear_index(1, ly, 1, 1, lx, ly, lz)
    this%node_linear_idx(4) = linear_index(lx, ly, 1, 1, lx, ly, lz)
    this%node_linear_idx(5) = linear_index(1, 1, lz, 1, lx, ly, lz)
    this%node_linear_idx(6) = linear_index(lx, 1, lz, 1, lx, ly, lz)
    this%node_linear_idx(7) = linear_index(1, ly, lz, 1, lx, ly, lz)
    this%node_linear_idx(8) = linear_index(lx, ly, lz, 1, lx, ly, lz)
  end subroutine coupled_vector_bc_projector_init

  !> Register one vector boundary condition in the coupled projector.
  !! @param[inout] bc Boundary condition to queue for resolution.
  subroutine coupled_vector_bc_projector_mark_bc(this, bc)
    class(coupled_vector_bc_projector_t), intent(inout) :: this
    class(bc_t), intent(inout), target :: bc

    if (.not. associated(this%coef)) then
       call neko_error("Coupled vector BC projector must be initialized " // &
            "before mark().")
    end if

    call this%bcs%append(bc)
  end subroutine coupled_vector_bc_projector_mark_bc

  !> Register a list of vector boundary conditions in the coupled projector.
  !! @param[in] bclst Boundary-condition list to queue.
  subroutine coupled_vector_bc_projector_mark_bc_list(this, bclst)
    class(coupled_vector_bc_projector_t), intent(inout) :: this
    type(bc_list_t), intent(in) :: bclst
    class(bc_t), pointer :: bc_i
    integer :: i

    do i = 1, bclst%size()
       bc_i => bclst%get(i)
       call this%mark(bc_i)
    end do
  end subroutine coupled_vector_bc_projector_mark_bc_list

  !> Clear the resolved mask-side state of the coupled projector.
  !! @details Releases the resolved Dirichlet and mixed masks, compact
  !! boundary-node caches, and the per-mixed-node constraint flags together
  !! with their device mirrors.
  subroutine coupled_vector_bc_projector_clear_masks(this)
    class(coupled_vector_bc_projector_t), intent(inout) :: this

    call this%dirichlet_dof_mask%free()
    call this%mixed_dof_mask%free()
    if (allocated(this%boundary_dof)) deallocate(this%boundary_dof)
    if (allocated(this%node_type)) deallocate(this%node_type)
    call this%boundary_idx%free()

    if (allocated(this%constraint_n)) then
       if (NEKO_BCKND_DEVICE .eq. 1 .and. &
            c_associated(this%constraint_n_d)) then
          call device_unmap(this%constraint_n, this%constraint_n_d)
       end if
       deallocate(this%constraint_n)
    end if
    if (allocated(this%constraint_t1)) then
       if (NEKO_BCKND_DEVICE .eq. 1 .and. &
            c_associated(this%constraint_t1_d)) then
          call device_unmap(this%constraint_t1, this%constraint_t1_d)
       end if
       deallocate(this%constraint_t1)
    end if
    if (allocated(this%constraint_t2)) then
       if (NEKO_BCKND_DEVICE .eq. 1 .and. &
            c_associated(this%constraint_t2_d)) then
          call device_unmap(this%constraint_t2, this%constraint_t2_d)
       end if
       deallocate(this%constraint_t2)
    end if
    this%constraint_n_d = c_null_ptr
    this%constraint_t1_d = c_null_ptr
    this%constraint_t2_d = c_null_ptr
  end subroutine coupled_vector_bc_projector_clear_masks

  !> Clear the resolved basis-side state of the coupled projector.
  subroutine coupled_vector_bc_projector_clear_basis(this)
    class(coupled_vector_bc_projector_t), intent(inout) :: this

    call this%n%free()
    call this%t1%free()
    call this%t2%free()
  end subroutine coupled_vector_bc_projector_clear_basis

  !> Finalize the coupled projector by resolving the accumulated BC list.
  !! @details This routine builds the dof masks for Dirichlet and mixed nodes,
  !! and computes the local basis for the mixed ones.
  subroutine coupled_vector_bc_projector_finalize(this, rebuild_mask)
    class(coupled_vector_bc_projector_t), intent(inout) :: this
    logical, intent(in) :: rebuild_mask
    if (this%bcs%size() .eq. 0) return

    if (rebuild_mask) then
       call this%rebuild_masks()
    end if

    call this%rebuild_basis()
  end subroutine coupled_vector_bc_projector_finalize

  !> Rebuild the resolved masks and local constraint flags.
  !! @details This routine reduces the queued boundary conditions to a nodal
  !! classification, constructs the global Dirichlet and mixed-node masks,
  !! computes BC-local resolved supports for `mixed_bc_t` instances, and fills
  !! the local constraint flags `constraint_n`, `constraint_t1`, and
  !! `constraint_t2` on the mixed-node mask.
  subroutine coupled_vector_bc_projector_rebuild_masks(this)
    class(coupled_vector_bc_projector_t), intent(inout) :: this
    type(field_t), pointer :: boundary_mask_field
    type(field_t), pointer :: node_type_field
    type(tuple_i4_t), pointer :: marked_faces(:)
    type(tuple_i4_t) :: marked_face
    integer, allocatable :: dirichlet_mask_values(:)
    integer, allocatable :: mixed_mask_values(:)
    integer, allocatable :: resolved_mask_values(:)
    integer :: scratch_idx(2)
    integer :: boundary_size
    integer :: compact_node_type_idx
    integer :: boundary_dof_key
    integer :: i, j, k, dof_size, m
    integer :: dirichlet_mask_size, mixed_mask_size, resolved_mask_size
    integer :: facet, el
    real(kind=rp) :: bc_type
    class(bc_t), pointer :: bc

    call this%clear_masks()

    call neko_scratch_registry%request_field(boundary_mask_field, &
         scratch_idx(1), .true.)
    call neko_scratch_registry%request_field(node_type_field, &
         scratch_idx(2), .true.)

    dof_size = this%dof%size()
    this%face_type = 5.0_rp

    ! Build a mask of all dofs on the boundary.
    ! Scratch clearing follows the active backend, whereas the resolution
    ! below is deliberately assembled on the host. Clear the host copy
    ! explicitly before setting its marked entries.
    call rzero(boundary_mask_field%x, dof_size)
    do i = 1, this%bcs%size()
       bc => this%bcs%get(i)

       if (.not. allocated(bc%msk)) then
          call neko_error("Attempting to finalize coupled projector " // &
               "unfinalized BC.")
       end if

       ! Mask all the dofs touched by this BC. Since %msk is propagated to all
       ! local dofs via gather-scatter, boundary_mask_field will contain all
       ! local nodes on the boundary, including those elements that don't touch
       ! it with a face.
       ! Note that bc%msk stores its length in slot 0 and is therefore passed
       ! with the `_0` masked-wrapper convention.
       call cfill_mask(boundary_mask_field%x, 1.0_rp, dof_size, &
            bc%msk(1:bc%msk(0)), &
            bc%msk(0))
    end do

    ! Set priority values (see the BC_* constants) for constraint assignment.
    ! Mimics the procedure in Nek5000 directly.
    ! The values are chosen so that a min reduction applies the
    ! highest-priority constraint.
    ! 5 -> unconstrained
    ! 3 -> tangentially constrained
    ! 2 -> normally constrained
    ! 0 -> fully constrained

    ! Fill the field to not mess up gather-scatter reduction later.
    call cfill(node_type_field%x, 5.0_rp, dof_size)

    do i = 1, this%bcs%size()
       bc => this%bcs%get(i)
       bc_type = bc%bc_type

       ! Store the type on each boundary face touched by this BC.
       ! This is the compact analogue of Nek's face-resident HFMASK field:
       ! one scalar type value per local (facet, element) pair.
       marked_faces => bc%marked_facet%array()
       do j = 1, bc%marked_facet%size()
          marked_face = marked_faces(j)
          facet = marked_face%x(1)
          el = marked_face%x(2)
          this%face_type(facet, el) = bc_type
       end do

       ! Note that facet_node_msk is used, so constraints are only directly
       ! applied to elements that touch the boundary at this point.
       do j = 1, bc%facet_node_msk(0)
          m = bc%facet_node_msk(j)
          ! The min here ensures that the highest-priority constraint is kept
          ! within a single element.
          node_type_field%x(m,1,1,1) = min(bc_type, node_type_field%x(m,1,1,1))
       end do
    end do

    ! Propagate constraints to all local dofs via gather-scatter.
    ! Ensures the highest-priority constraint is kept across element
    ! boundaries.
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call node_type_field%copy_from(HOST_TO_DEVICE, .true.)
    end if
    call this%coef%gs_h%op(node_type_field, GS_OP_MIN)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call node_type_field%copy_from(DEVICE_TO_HOST, .true.)
    end if

    ! Build compact nodal type cache for all boundary dofs.
    ! First pass to count the size of the boundary dof set.

    boundary_size = 0
    !$omp parallel do reduction(+:boundary_size)
    do i = 1, dof_size
       if (boundary_mask_field%x(i,1,1,1) .gt. 0.5_rp) then
          boundary_size = boundary_size + 1
       end if
    end do
    !$omp end parallel do

    ! Linear indices of boundary dofs
    allocate(this%boundary_dof(boundary_size))
    ! Node BC type per boundary dof
    allocate(this%node_type(boundary_size))
    ! Mapping from global dof index to compact boundary dof index and type.
    call this%boundary_idx%init(boundary_size, compact_node_type_idx)

    boundary_size = 0
    do i = 1, dof_size
       if (boundary_mask_field%x(i,1,1,1) .lt. 0.5_rp) cycle

       boundary_size = boundary_size + 1
       this%boundary_dof(boundary_size) = i
       this%node_type(boundary_size) = node_type_field%x(i,1,1,1)
       boundary_dof_key = i
       compact_node_type_idx = boundary_size
       call this%boundary_idx%set(boundary_dof_key, compact_node_type_idx)
    end do

    ! For mixed BCs, build a resolved subset of the original bc%msk support.
    ! A dof survives in the resolved mask only if the globally reduced type
    ! still matches the type semantics of that BC. Nodes that were touched by
    ! the BC originally, but whose meaning changed after shared-node reduction,
    ! are therefore dropped here.
    do i = 1, this%bcs%size()
       bc => this%bcs%get(i)

       select type (bc)
       class is (mixed_bc_t)
          call bc%resolved_msk%free()
          call bc%n%free()
          call bc%t1%free()
          call bc%t2%free()

          bc_type = bc%bc_type

          ! First pass to count the size of the resolved mask.
          resolved_mask_size = 0
          do j = 1, bc%msk(0)
             k = bc%msk(j)
             ! Compare face and node bc type
             if (abs(node_type_field%x(k,1,1,1) - bc_type) .lt. 1.0e-6_rp) then
                resolved_mask_size = resolved_mask_size + 1
             end if
          end do

          allocate(resolved_mask_values(resolved_mask_size))

          ! Fill in the mask values
          resolved_mask_size = 0
          do j = 1, bc%msk(0)
             k = bc%msk(j)
             if (abs(node_type_field%x(k,1,1,1) - bc_type) .lt. 1.0e-6_rp) then
                resolved_mask_size = resolved_mask_size + 1
                resolved_mask_values(resolved_mask_size) = k
             end if
          end do

          call bc%resolved_msk%init(resolved_mask_values, resolved_mask_size)
          deallocate(resolved_mask_values)
       end select
    end do

    ! Partition the resolved boundary dofs into the fully constrained subset
    ! and the mixed subset. Only the latter needs a local basis.
    dirichlet_mask_size = 0
    mixed_mask_size = 0

    !$omp parallel do reduction(+:dirichlet_mask_size,mixed_mask_size)
    do i = 1, dof_size
       ! Internal node
       if (boundary_mask_field%x(i,1,1,1) .lt. 0.5_rp) cycle

       if (node_type_field%x(i,1,1,1) .lt. 1.9_rp) then
          dirichlet_mask_size = dirichlet_mask_size + 1
       else if (node_type_field%x(i,1,1,1) .gt. 1.9_rp .and. &
            node_type_field%x(i,1,1,1) .lt. 3.9_rp) then
          mixed_mask_size = mixed_mask_size + 1
       end if
    end do
    !$omp end parallel do

    allocate(dirichlet_mask_values(dirichlet_mask_size))
    allocate(mixed_mask_values(mixed_mask_size))

    ! We reuse the variables as counters in the loop below, which actually
    ! fills the masks, now that we know the size.
    dirichlet_mask_size = 0
    mixed_mask_size = 0
    do i = 1, dof_size
       if (boundary_mask_field%x(i,1,1,1) .lt. 0.5_rp) cycle

       if (node_type_field%x(i,1,1,1) .lt. 1.9_rp) then
          dirichlet_mask_size = dirichlet_mask_size + 1
          dirichlet_mask_values(dirichlet_mask_size) = i
       else if (node_type_field%x(i,1,1,1) .gt. 1.9_rp .and. &
            node_type_field%x(i,1,1,1) .lt. 3.9_rp) then
          mixed_mask_size = mixed_mask_size + 1
          mixed_mask_values(mixed_mask_size) = i
       end if
    end do

    ! Initialize both masks. The temporary arrays remain valid actual
    ! arguments when their allocated extent is zero.
    call this%dirichlet_dof_mask%init(dirichlet_mask_values, &
         dirichlet_mask_size)
    call this%mixed_dof_mask%init(mixed_mask_values, mixed_mask_size)

    ! Allocate mixed-node constraints and fill them from the reduced type.
    allocate(this%constraint_n(mixed_mask_size))
    allocate(this%constraint_t1(mixed_mask_size))
    allocate(this%constraint_t2(mixed_mask_size))
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%constraint_n, this%constraint_n_d, &
            size(this%constraint_n))
       call device_map(this%constraint_t1, this%constraint_t1_d, &
            size(this%constraint_t1))
       call device_map(this%constraint_t2, this%constraint_t2_d, &
            size(this%constraint_t2))
    end if

    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp do
    do i = 1, mixed_mask_size
       j = mixed_mask_values(i)

       if (node_type_field%x(j,1,1,1) .lt. 1.9_rp) then
          this%constraint_n(i) = 1
          this%constraint_t1(i) = 1
          this%constraint_t2(i) = 1
       else if (node_type_field%x(j,1,1,1) .gt. 1.9_rp .and. &
            node_type_field%x(j,1,1,1) .lt. 2.9_rp) then
          this%constraint_n(i) = 1
          this%constraint_t1(i) = 0
          this%constraint_t2(i) = 0
       else if (node_type_field%x(j,1,1,1) .gt. 2.9_rp .and. &
            node_type_field%x(j,1,1,1) .lt. 3.9_rp) then
          this%constraint_n(i) = 0
          this%constraint_t1(i) = 1
          this%constraint_t2(i) = 1
       else if (node_type_field%x(j,1,1,1) .gt. 3.9_rp) then
          this%constraint_n(i) = 0
          this%constraint_t1(i) = 0
          this%constraint_t2(i) = 0
       end if
    end do
    !$omp end do

    deallocate(dirichlet_mask_values)
    deallocate(mixed_mask_values)
    call neko_scratch_registry%relinquish_field(scratch_idx)
  end subroutine coupled_vector_bc_projector_rebuild_masks

  !> Rebuild the local mixed-node basis.
  !! @details Using the resolved mixed-node mask and constraint metadata built
  !! by `rebuild_masks()`, this routine reconstructs consistent normals on
  !! face, edge, and corner nodes and then builds an orthonormal local basis
  !! `(n, t1, t2)` for each mixed node. It then propagates the basis data to
  !! all marked mixed boundary conditions.
  subroutine coupled_vector_bc_projector_rebuild_basis(this)
    class(coupled_vector_bc_projector_t), intent(inout) :: this
    type(field_t), pointer :: normal_x_field
    type(field_t), pointer :: normal_y_field
    type(field_t), pointer :: normal_z_field
    class(bc_t), pointer :: bc
    integer, pointer :: mixed_dof_values(:)
    integer, allocatable :: dof_to_mixed_idx(:)
    integer :: scratch_idx(3)
    integer :: node_type_lookup_status, compact_node_type_idx
    integer :: i, j, k, dof_size, m
    integer :: idx(4), facet, el, edge, node, ii, p
    integer :: rst(3), rst1(3), rst2(3), step_rst(3)
    integer :: edge_len, edge_idx, node_idx
    real(kind=rp) :: normal(3), t1_vec(3), t2_vec(3), len, bc_type
    real(kind=rp), parameter :: normal_tol = 100.0_rp * epsilon(1.0_rp)
    character(len=LOG_SIZE) :: error_msg

    call this%clear_basis()

    call neko_scratch_registry%request_field(normal_x_field, scratch_idx(1), &
         .true.)
    call neko_scratch_registry%request_field(normal_y_field, scratch_idx(2), &
         .true.)
    call neko_scratch_registry%request_field(normal_z_field, scratch_idx(3), &
         .true.)

    dof_size = this%dof%size()
    m = this%mixed_dof_mask%size()
    mixed_dof_values => this%mixed_dof_mask%get()

    call this%n%init(3, m)
    call this%t1%init(3, m)
    call this%t2%init(3, m)

    ! A mapping between the field linear index of a mixed node into its index
    ! in the mixed_dof_mask.
    allocate(dof_to_mixed_idx(dof_size))
    dof_to_mixed_idx = 0

    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp do
    do i = 1, m
       dof_to_mixed_idx(mixed_dof_values(i)) = i
    end do
    !$omp end do

    ! this%face_type stores the face-based bc_type values, and this%node_type
    ! stores them node-wise, after propagation with min reduction.
    ! Note that the propagation means that some nodes may have a different,
    ! higher-priority class than their owning face!
    ! The algorithm for constructing normals below will make use of both
    ! classifications when looking at edges and corners. We will really only
    ! care about types 2 and 3, i.e. mixed bcs. The key question will be
    ! whether a given face should contribute its normal to the edge and corner
    ! dofs. The idea is that the face only contributes its normal if its type
    ! is the same as that of the node.

    ! Set normals at unambiguous boundary dofs based on face normals.
    call rzero(normal_x_field%x, dof_size)
    call rzero(normal_y_field%x, dof_size)
    call rzero(normal_z_field%x, dof_size)

    this%n = 0.0_rp
    this%t1 = 0.0_rp
    this%t2 = 0.0_rp

    do i = 1, this%bcs%size()
       bc => this%bcs%get(i)

       ! First pass: seed the local normal field on all nodes that lie
       ! on directly marked mixed faces.
       ! Since several faces will own edge and corner nodes, the normals there
       ! will be overwritten in arbitrary order, but we don't care because we
       ! will reset those later and treat them specially.
       do j = 1, bc%facet_node_msk(0)
          ! Global linear index of the node on which to set the normal.
          k = bc%facet_node_msk(j)

          ! Grab face and ijke indices to address face_type and get_normal.
          facet = bc%facet(j)
          idx = nonlinear_index(k, this%coef%Xh%lx, this%coef%Xh%ly, &
               this%coef%Xh%lz)

          if (this%face_type(facet, idx(4)) .lt. 1.9_rp .or. &
               this%face_type(facet, idx(4)) .gt. 3.1_rp) cycle

          normal = this%coef%get_normal(idx(1), idx(2), idx(3), idx(4), &
               facet)
          normal_x_field%x(k,1,1,1) = normal(1)
          normal_y_field%x(k,1,1,1) = normal(2)
          normal_z_field%x(k,1,1,1) = normal(3)
       end do
    end do

    !write(*,*) "Seeded normals at directly marked mixed nodes in " // &
    !     "coupled vector BC projector."

    ! We now treat the special edges and conrners. Everything is done locally
    ! per element, using reference element address tables found in hex.f90
    ! and inside this type.

    ! Mixed edge interiors are rebuilt from the normals of the adjacent
    ! faces whose local face type matches the reduced nodal type.
    ! This is the central point: if the adjacent face is a different type,
    ! which by construction can only be a lower-priority type, then it
    ! should not contribute its normal.
    !
    ! Consider the following 2D example. In 2D an edge becomes a node in the
    ! corner of the element. Look at the node marked with X. After the nodal
    ! type is propagated, it will have type 2---the highest-priority of the
    ! adjacent. So, only the face with type 2 in El 2 will contribute to the
    ! normal. This is a rather extreme example, but it illustrates well what
    ! can happen.
    !
    ! ---------
    ! | El 1  |
    ! |     3 |
    ! |       |
    ! |   3   |   3
    ! --------X---------
    !         | El 2  |
    !         |       |
    !         | 2     |
    !         |       |
    !         --------0

    ! We loop over the edges of all elements, so we catch those that touch
    ! the boundary with an edge or a corner but not a face. Note that we will
    ! only treat the interior nodes of the edge here. The endpoints, i.e.
    ! the corner nodes are handled in the next loop.
    if (m .gt. 0) then
       do el = 1, this%coef%msh%nelv
          do edge = 1, size(edge_nodes, 2)
             ! Representitive rst index in the middle of an edge.
             rst = this%edge_mid_rst(:, edge)

             ! Global linear index of the midpoint node.
             edge_idx = linear_index(rst(1), rst(2), rst(3), el, &
                  this%coef%Xh%lx, this%coef%Xh%ly, this%coef%Xh%lz)

             ! Get the BC type.
             node_type_lookup_status = this%boundary_idx%get( &
                  edge_idx, compact_node_type_idx)

             ! If we did not find the node in the lookup, it means it is not
             ! a boundary node and we can skip.
             if (node_type_lookup_status .ne. 0) cycle
             bc_type = abs(this%node_type(compact_node_type_idx))

             ! If this is not a mixed bc edge, just leave it alone.
             if (bc_type .lt. 1.9_rp .or. bc_type .gt. 3.1_rp) cycle

             ! Recall that "node" in the lookup table names refer to element
             ! corners. This is to stay consistent with the hex_t notation.

             ! Get edge endpoints index triples. For example,
             ! (1, 1, 1) and (lx, 1, 1)
             rst1 = this%node_rst(:, edge_nodes(1, edge))
             rst2 = this%node_rst(:, edge_nodes(2, edge))

             ! Compute number of gll nodes on the edge, so lx, ly, or lz.
             ! Which currently in Neko is one and the same.
             edge_len = maxval(abs(rst2 - rst1)) + 1

             ! The running index direction along the edge in rst-space.
             ! This is a triple, but only one component is nonzero.
             step_rst = 0
             do ii = 1, 3
                if (rst2(ii) .gt. rst1(ii)) then
                   step_rst(ii) = 1
                else if (rst2(ii) .lt. rst1(ii)) then
                   step_rst(ii) = -1
                end if
             end do

             ! Loop over interior edge nodes and reset the normals.
             do p = 2, edge_len - 1
                rst = rst1 + (p - 1) * step_rst
                k = linear_index(rst(1), rst(2), rst(3), el, &
                     this%coef%Xh%lx, this%coef%Xh%ly, this%coef%Xh%lz)
                normal_x_field%x(k,1,1,1) = 0.0_rp
                normal_y_field%x(k,1,1,1) = 0.0_rp
                normal_z_field%x(k,1,1,1) = 0.0_rp
             end do

             ! Loop over the faces adjacent to this edge and add the normals
             ! if the BC type matches between the face and the edge.
             do ii = 1, size(edge_faces, 1)
                facet = edge_faces(ii, edge)

                ! Skip if the BC type is not the same.
                if (abs(bc_type - this%face_type(facet, el)) .gt. 1.0e-6_rp) then
                   cycle
                end if

                ! Loop over the interior edge nodes again and add the normals.
                do p = 2, edge_len - 1
                   rst = rst1 + (p - 1) * step_rst
                   k = linear_index(rst(1), rst(2), rst(3), el, &
                        this%coef%Xh%lx, this%coef%Xh%ly, this%coef%Xh%lz)
                   normal = this%coef%get_normal(rst(1), rst(2), rst(3), &
                        el, facet)
                   normal_x_field%x(k,1,1,1) = &
                        normal_x_field%x(k,1,1,1) + normal(1)
                   normal_y_field%x(k,1,1,1) = &
                        normal_y_field%x(k,1,1,1) + normal(2)
                   normal_z_field%x(k,1,1,1) = &
                        normal_z_field%x(k,1,1,1) + normal(3)
                end do
             end do
          end do
       end do

       !write(*,*) "Finished reconstructing normals at mixed edges in " // &
       !     "coupled vector BC projector."

       ! Mixed corner node normals are rebuilt from the adjacent faces whose
       ! local face type matches the reduced nodal type at that node.
       do el = 1, this%coef%msh%nelv
          do node = 1, size(this%node_linear_idx)
             rst = this%node_rst(:, node)
             node_idx = linear_index(rst(1), rst(2), rst(3), el, &
                  this%coef%Xh%lx, this%coef%Xh%ly, this%coef%Xh%lz)

             node_type_lookup_status = this%boundary_idx%get( &
                  node_idx, compact_node_type_idx)
             if (node_type_lookup_status .ne. 0) cycle
             bc_type = abs(this%node_type(compact_node_type_idx))

             ! Ignore if the BC type is not a mixed one.
             if (bc_type .lt. 1.9_rp .or. bc_type .gt. 3.1_rp) cycle

             ! Kill the normal to start fresh.
             normal_x_field%x(node_idx,1,1,1) = 0.0_rp
             normal_y_field%x(node_idx,1,1,1) = 0.0_rp
             normal_z_field%x(node_idx,1,1,1) = 0.0_rp

             ! Note, 3 faces share a corner node in 3D.
             do ii = 1, this%coef%msh%gdim
                facet = node_faces(ii, node)

                ! Check type agreement
                if (abs(bc_type - this%face_type(facet, el)) .gt. 1.0e-6_rp) then
                   cycle
                end if

                ! Add the normal.
                normal = this%coef%get_normal(rst(1), rst(2), rst(3), &
                     el, facet)
                normal_x_field%x(node_idx,1,1,1) = &
                     normal_x_field%x(node_idx,1,1,1) + normal(1)
                normal_y_field%x(node_idx,1,1,1) = &
                     normal_y_field%x(node_idx,1,1,1) + normal(2)
                normal_z_field%x(node_idx,1,1,1) = &
                     normal_z_field%x(node_idx,1,1,1) + normal(3)
             end do
          end do
       end do
    end if

    ! We are done element-wise. We now need global consistency at shared nodes.
    !
    ! For cyclic meshes, however, periodic counterparts may not be co-planar.
    ! If we directly gather Cartesian normal components, we mix vectors that are
    ! expressed in different local frames across the periodic map, which can
    ! tilt the reconstructed normal.
    !
    ! We therefore follow the same cyclic treatment used in other vector
    ! operators in Neko:
    ! 1) Rotate cyclic-marked nodes into the cyclic-normal/tangential frame.
    ! 2) Perform the global GS_OP_ADD while all contributions use that frame.
    ! 3) Rotate cyclic-marked nodes back to Cartesian.
    !
    ! This keeps MPI/shared-node averaging intact while preventing frame-mixing
    ! on cyclic couplings.
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call normal_x_field%copy_from(HOST_TO_DEVICE, .false.)
       call normal_y_field%copy_from(HOST_TO_DEVICE, .false.)
       call normal_z_field%copy_from(HOST_TO_DEVICE, .true.)
    end if

    if (this%coef%cyclic) then
       call rotate_cyc(normal_x_field%x, normal_y_field%x, normal_z_field%x, &
            1, this%coef)
    end if

    call this%coef%gs_h%op(normal_x_field, GS_OP_ADD)
    call this%coef%gs_h%op(normal_y_field, GS_OP_ADD)
    call this%coef%gs_h%op(normal_z_field, GS_OP_ADD)

    if (this%coef%cyclic) then
       call rotate_cyc(normal_x_field%x, normal_y_field%x, normal_z_field%x, &
            0, this%coef)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call normal_x_field%copy_from(DEVICE_TO_HOST, .false.)
       call normal_y_field%copy_from(DEVICE_TO_HOST, .false.)
       call normal_z_field%copy_from(DEVICE_TO_HOST, .true.)
    end if

    ! Normalize normals and build tangential directions for all mixed nodes.
    ! Populate this into the basis components in the type.

    !OCL NORECURRENCE, NOVREC, NOALIAS
    !DIR$ CONCURRENT
    !DIR$ IVDEP
    !GCC$ ivdep
    !$omp do
    do i = 1, m
       j = mixed_dof_values(i)

       ! Normalize the normal
       normal(1) = normal_x_field%x(j,1,1,1)
       normal(2) = normal_y_field%x(j,1,1,1)
       normal(3) = normal_z_field%x(j,1,1,1)
       len = sqrt(sum(normal**2))
       if (len .le. normal_tol) then
          write(error_msg, '(A,I0,A,ES13.6,A)') &
               "Coupled vector BC projector could not construct a normal " // &
               "at local DOF ", j, " (norm = ", len, ")."
          call neko_error(error_msg)
       end if

       this%n%x(:,i) = normal / len

       ! Select the first tangent direction by crossing the normal with a
       ! coordinate axis. Use y near the z-axis, where the usual z-axis
       ! construction becomes ill-conditioned.
       if (abs(this%n%x(3,i)) .gt. 0.999_rp) then
          t1_vec = [ this%n%x(3,i), 0.0_rp, -this%n%x(1,i) ]
       else
          t1_vec = [ -this%n%x(2,i), this%n%x(1,i), 0.0_rp ]
       end if
       len = sqrt(sum(t1_vec**2))
       if (len .gt. 0.0_rp) then
          this%t1%x(:,i) = t1_vec / len
       end if

       ! Get t2 as a cross product of n and t1.
       t2_vec(1) = this%n%x(2,i) * this%t1%x(3,i) - &
            this%n%x(3,i) * this%t1%x(2,i)
       t2_vec(2) = this%n%x(3,i) * this%t1%x(1,i) - &
            this%n%x(1,i) * this%t1%x(3,i)
       t2_vec(3) = this%n%x(1,i) * this%t1%x(2,i) - &
            this%n%x(2,i) * this%t1%x(1,i)
       len = sqrt(sum(t2_vec**2))
       if (len .gt. 0.0_rp) then
          this%t2%x(:,i) = t2_vec / len
       end if
    end do
    !$omp end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%constraint_n, this%constraint_n_d, &
            size(this%constraint_n), HOST_TO_DEVICE, sync = .true.)
       call device_memcpy(this%constraint_t1, this%constraint_t1_d, &
            size(this%constraint_t1), HOST_TO_DEVICE, sync = .true.)
       call device_memcpy(this%constraint_t2, this%constraint_t2_d, &
            size(this%constraint_t2), HOST_TO_DEVICE, sync = .true.)
       call this%n%copy_from(HOST_TO_DEVICE, .true.)
       call this%t1%copy_from(HOST_TO_DEVICE, .true.)
       call this%t2%copy_from(HOST_TO_DEVICE, .true.)
    end if

    ! Transfer the final mixed-node basis into each mixed BC on its
    ! resolved support, so strong application on the physical field can use
    ! BC-local data rather than the global projector internals.
    do i = 1, this%bcs%size()
       bc => this%bcs%get(i)

       select type (bc)
       class is (mixed_bc_t)
          m = bc%resolved_msk%size()
          call bc%n%init(3, m)
          call bc%t1%init(3, m)
          call bc%t2%init(3, m)

          do j = 1, m
             k = bc%resolved_msk%get(j)
             p = dof_to_mixed_idx(k)

             if (p .eq. 0) then
                call neko_error("Mixed BC resolved_msk entry missing from " // &
                     "the coupled projector mixed basis.")
             end if

             bc%n%x(:,j) = this%n%x(:,p)
             bc%t1%x(:,j) = this%t1%x(:,p)
             bc%t2%x(:,j) = this%t2%x(:,p)
          end do

          if (NEKO_BCKND_DEVICE .eq. 1) then
             call bc%n%copy_from(HOST_TO_DEVICE, .true.)
             call bc%t1%copy_from(HOST_TO_DEVICE, .true.)
             call bc%t2%copy_from(HOST_TO_DEVICE, .true.)
          end if
       end select
    end do
    call neko_scratch_registry%relinquish_field(scratch_idx)
    if (allocated(dof_to_mixed_idx)) deallocate(dof_to_mixed_idx)
  end subroutine coupled_vector_bc_projector_rebuild_basis

  !> Apply homogeneous boundary constraints in the local basis.
  !! @param[inout] x x-component field values.
  !! @param[inout] y y-component field values.
  !! @param[inout] z z-component field values.
  !! @param[in] n Number of entries in each component array.
  !! @param[inout] strm Optional backend stream/queue used on device backends.
  subroutine coupled_vector_bc_projector_apply(this, x, y, z, n, strm)
    class(coupled_vector_bc_projector_t), intent(in) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: x(n)
    real(kind=rp), intent(inout) :: y(n)
    real(kind=rp), intent(inout) :: z(n)
    type(c_ptr), intent(inout), optional :: strm

    integer, pointer :: dirichlet_msk(:)
    integer, pointer :: mixed_msk(:)
    integer :: i, j, m
    real(kind=rp) :: u(3), uloc(3)
    type(c_ptr) :: x_d, y_d, z_d
    type(c_ptr):: strm_

    if (present(strm)) then
       strm_ = strm
    else
       strm_ = glb_cmd_queue
    end if

    ! Fully constrained nodes do not need the local basis. They are simply
    ! zeroed in Cartesian space before the mixed-node pass.
    m = this%dirichlet_dof_mask%size()
    if (m .gt. 0) then
       if (NEKO_BCKND_DEVICE .eq. 1) then


          x_d = device_get_ptr(x)
          y_d = device_get_ptr(y)
          z_d = device_get_ptr(z)
          call device_cfill_mask(x_d, 0.0_rp, n, &
               this%dirichlet_dof_mask%get_d(), m, strm = strm_)
          call device_cfill_mask(y_d, 0.0_rp, n, &
               this%dirichlet_dof_mask%get_d(), m, strm = strm_)
          call device_cfill_mask(z_d, 0.0_rp, n, &
               this%dirichlet_dof_mask%get_d(), m, strm = strm_)
       else
          dirichlet_msk => this%dirichlet_dof_mask%get()
          call cfill_mask(x, 0.0_rp, n, dirichlet_msk, m)
          call cfill_mask(y, 0.0_rp, n, dirichlet_msk, m)
          call cfill_mask(z, 0.0_rp, n, dirichlet_msk, m)
       end if
    end if

    m = this%mixed_dof_mask%size()

    if (m .gt. 0) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          x_d = device_get_ptr(x)
          y_d = device_get_ptr(y)
          z_d = device_get_ptr(z)
          call device_coupled_vector_bc_projector_apply( &
               this%mixed_dof_mask%get_d(), x_d, y_d, z_d, &
               this%constraint_n_d, this%constraint_t1_d, &
               this%constraint_t2_d, this%n%x_d, this%t1%x_d, &
               this%t2%x_d, m, strm_)
       else
          mixed_msk => this%mixed_dof_mask%get()

          do i = 1, m
             j = mixed_msk(i)

             u(1) = x(j)
             u(2) = y(j)
             u(3) = z(j)

             uloc(1) = u(1) * this%n%x(1,i) + u(2) * this%n%x(2,i) + &
                  u(3) * this%n%x(3,i)
             uloc(2) = u(1) * this%t1%x(1,i) + u(2) * this%t1%x(2,i) + &
                  u(3) * this%t1%x(3,i)
             uloc(3) = u(1) * this%t2%x(1,i) + u(2) * this%t2%x(2,i) + &
                  u(3) * this%t2%x(3,i)

             if (this%constraint_n(i) .ne. 0) uloc(1) = 0.0_rp
             if (this%constraint_t1(i) .ne. 0) uloc(2) = 0.0_rp
             if (this%constraint_t2(i) .ne. 0) uloc(3) = 0.0_rp

             u = uloc(1) * this%n%x(:,i) + uloc(2) * this%t1%x(:,i) + &
                  uloc(3) * this%t2%x(:,i)

             x(j) = u(1)
             y(j) = u(2)
             z(j) = u(3)
          end do
       end if
    end if
  end subroutine coupled_vector_bc_projector_apply

  !> Write fields showing the coupled projector mask and basis.
  !! @param[in] field_name Optional base name for the output file. The `.fld`
  !! suffix is appended automatically.
  subroutine coupled_vector_bc_projector_debug_output(this, field_name)
    use device_math, only : device_cfill, device_cfill_mask
    class(coupled_vector_bc_projector_t), intent(inout) :: this
    character(len=*), intent(in), optional :: field_name
    type(field_t), pointer :: mask_field
    type(field_t), pointer :: nx_field, ny_field, nz_field
    type(field_list_t) :: basis_fields
    type(fld_file_t) :: basis_file
    integer :: scratch_idx(4)
    integer, pointer :: mixed_mask_values(:)
    integer :: dof_size, mixed_mask_size
    character(len=:), allocatable :: field_name_

    if (present(field_name)) then
       field_name_ = trim(field_name)
    else
       field_name_ = 'bc_projector'
    end if

    call neko_scratch_registry%request_field(mask_field, scratch_idx(1), .true.)
    call neko_scratch_registry%request_field(nx_field, scratch_idx(2), .true.)
    call neko_scratch_registry%request_field(ny_field, scratch_idx(3), .true.)
    call neko_scratch_registry%request_field(nz_field, scratch_idx(4), .true.)

    dof_size = this%dof%size()
    mixed_mask_size = this%mixed_dof_mask%size()

    call rzero(mask_field%x, dof_size)
    call rzero(nx_field%x, dof_size)
    call rzero(ny_field%x, dof_size)
    call rzero(nz_field%x, dof_size)

    if (this%mixed_dof_mask%is_set()) then
       mixed_mask_values => this%mixed_dof_mask%get()
       call masked_scatter_copy(nx_field%x(:,1,1,1), this%n%x(1,:), &
            mixed_mask_values, dof_size, mixed_mask_size)
       call masked_scatter_copy(ny_field%x(:,1,1,1), this%n%x(2,:), &
            mixed_mask_values, dof_size, mixed_mask_size)
       call masked_scatter_copy(nz_field%x(:,1,1,1), this%n%x(3,:), &
            mixed_mask_values, dof_size, mixed_mask_size)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(mask_field%x_d, 5.0_rp, dof_size)

       if (this%dirichlet_dof_mask%is_set()) then
          call device_cfill_mask(mask_field%x_d, 1.0_rp, dof_size, &
               this%dirichlet_dof_mask%get_d(), &
               this%dirichlet_dof_mask%size())
       end if

       if (this%mixed_dof_mask%is_set()) then
          call device_cfill_mask(mask_field%x_d, 2.0_rp, dof_size, &
               this%mixed_dof_mask%get_d(), this%mixed_dof_mask%size())
       end if

       call device_memcpy(mask_field%x, mask_field%x_d, dof_size, &
            DEVICE_TO_HOST, &
            sync = .true.)
    else
       call cfill(mask_field%x, 5.0_rp, dof_size)

       if (this%dirichlet_dof_mask%is_set()) then
          call cfill_mask(mask_field%x, 1.0_rp, dof_size, &
               this%dirichlet_dof_mask%get(), this%dirichlet_dof_mask%size())
       end if

       if (this%mixed_dof_mask%is_set()) then
          call cfill_mask(mask_field%x, 2.0_rp, dof_size, &
               this%mixed_dof_mask%get(), this%mixed_dof_mask%size())
       end if
    end if

    call basis_fields%init(4)
    call basis_fields%assign(1, mask_field)
    call basis_fields%assign(2, nx_field)
    call basis_fields%assign(3, ny_field)
    call basis_fields%assign(4, nz_field)

    call basis_file%init(field_name_ // '.fld')
    call basis_file%write(basis_fields)
    call basis_fields%free()

    call neko_scratch_registry%relinquish_field(scratch_idx)
  end subroutine coupled_vector_bc_projector_debug_output

  !> Write scalar fields with the normal component `u.n` on mixed BC nodes.
  !! @details The output is zero outside `mixed_dof_mask`. The normal vectors
  !! come from (1) the resolved basis stored in this projector and
  !! (2) the normals stored in `coef_t`. No area scaling is applied.
  !! @param[in] x x-component of the vector field.
  !! @param[in] y y-component of the vector field.
  !! @param[in] z z-component of the vector field.
  !! @param[in] n Number of entries in x, y, z.
  !! @param[in] field_name Optional base name for the output file. The `.fld`
  !! suffix is appended automatically.
  subroutine coupled_vector_bc_projector_debug_output_normal_component( &
       this, x, y, z, n, field_name)
    class(coupled_vector_bc_projector_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: x(n)
    real(kind=rp), intent(in) :: y(n)
    real(kind=rp), intent(in) :: z(n)
    character(len=*), intent(in), optional :: field_name
    type(field_t), pointer :: normal_component_field
    type(field_t), pointer :: normal_component_coef_field
    type(field_list_t) :: output_fields
    type(fld_file_t) :: output_file
    type(field_t), pointer :: projector_nx_field, projector_ny_field, projector_nz_field
    type(field_t), pointer :: coef_nx_field, coef_ny_field, coef_nz_field
    type(field_list_t) :: normals_fields
    type(fld_file_t) :: normals_file
    integer :: scratch_idx(8)
    integer, pointer :: mixed_mask_values(:)
    integer, allocatable :: mixed_lut(:)
    integer :: dof_size, mixed_mask_size
    integer :: i, j, m, k, facet
    integer :: idx(4)
    real(kind=rp) :: coef_normal(3)
    class(bc_t), pointer :: bc
    character(len=:), allocatable :: field_name_

    if (present(field_name)) then
       field_name_ = trim(field_name)
    else
       field_name_ = 'bc_projector_normal_component'
    end if

    call neko_scratch_registry%request_field(normal_component_field, &
         scratch_idx(1), .true.)
    call neko_scratch_registry%request_field(normal_component_coef_field, &
         scratch_idx(2), .true.)
    call neko_scratch_registry%request_field(projector_nx_field, &
         scratch_idx(3), .true.)
    call neko_scratch_registry%request_field(projector_ny_field, &
         scratch_idx(4), .true.)
    call neko_scratch_registry%request_field(projector_nz_field, &
         scratch_idx(5), .true.)
    call neko_scratch_registry%request_field(coef_nx_field, &
         scratch_idx(6), .true.)
    call neko_scratch_registry%request_field(coef_ny_field, &
         scratch_idx(7), .true.)
    call neko_scratch_registry%request_field(coef_nz_field, &
         scratch_idx(8), .true.)

    dof_size = this%dof%size()
    mixed_mask_size = this%mixed_dof_mask%size()
    call rzero(normal_component_field%x, dof_size)
    call rzero(normal_component_coef_field%x, dof_size)
    call rzero(projector_nx_field%x, dof_size)
    call rzero(projector_ny_field%x, dof_size)
    call rzero(projector_nz_field%x, dof_size)
    call rzero(coef_nx_field%x, dof_size)
    call rzero(coef_ny_field%x, dof_size)
    call rzero(coef_nz_field%x, dof_size)

    allocate(mixed_lut(dof_size))
    mixed_lut = 0
    if (this%mixed_dof_mask%is_set()) then
       mixed_mask_values => this%mixed_dof_mask%get()
       do i = 1, mixed_mask_size
          mixed_lut(mixed_mask_values(i)) = i
       end do
    end if

    do i = 1, this%bcs%size()
       bc => this%bcs%get(i)
       do m = 1, bc%facet_node_msk(0)
          k = bc%facet_node_msk(m)
          facet = bc%facet(m)
          idx = nonlinear_index(k, this%coef%Xh%lx, this%coef%Xh%ly, &
               this%coef%Xh%lz)

          coef_normal = this%coef%get_normal(idx(1), idx(2), idx(3), idx(4), &
               facet)
          coef_nx_field%x(k,1,1,1) = coef_normal(1)
          coef_ny_field%x(k,1,1,1) = coef_normal(2)
          coef_nz_field%x(k,1,1,1) = coef_normal(3)
          normal_component_coef_field%x(k,1,1,1) = &
               x(k) * coef_normal(1) + y(k) * coef_normal(2) + z(k) * coef_normal(3)

          j = mixed_lut(k)
          if (j .gt. 0) then
             projector_nx_field%x(k,1,1,1) = this%n%x(1,j)
             projector_ny_field%x(k,1,1,1) = this%n%x(2,j)
             projector_nz_field%x(k,1,1,1) = this%n%x(3,j)
             normal_component_field%x(k,1,1,1) = &
                  x(k) * this%n%x(1,j) + y(k) * this%n%x(2,j) + z(k) * this%n%x(3,j)
          end if
       end do
    end do

    deallocate(mixed_lut)

    call output_fields%init(2)
    call output_fields%assign(1, normal_component_field)
    call output_fields%assign(2, normal_component_coef_field)
    call output_file%init(field_name_ // '.fld')
    call output_file%write(output_fields)
    call output_fields%free()

    call normals_fields%init(3)
    call normals_fields%assign(1, projector_nx_field)
    call normals_fields%assign(2, projector_ny_field)
    call normals_fields%assign(3, projector_nz_field)
    call normals_file%init(field_name_ // '_projector_normals.fld')
    call normals_file%write(normals_fields)
    call normals_fields%free()

    call normals_fields%init(3)
    call normals_fields%assign(1, coef_nx_field)
    call normals_fields%assign(2, coef_ny_field)
    call normals_fields%assign(3, coef_nz_field)
    call normals_file%init(field_name_ // '_coef_normals.fld')
    call normals_file%write(normals_fields)
    call normals_fields%free()

    call neko_scratch_registry%relinquish_field(scratch_idx)
  end subroutine coupled_vector_bc_projector_debug_output_normal_component

end module vector_bc_projector
