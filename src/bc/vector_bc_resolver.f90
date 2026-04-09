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
!> Vector boundary-condition resolvers.
module vector_bc_resolver
  use bc, only : bc_t
  use bc_list, only : bc_list_t
  use mixed_bc, only : mixed_bc_t
  use mask, only : mask_t
  use coefs, only : coef_t
  use dofmap, only : dofmap_t
  use field, only : field_t
  use field_list, only : field_list_t
  use fld_file, only : fld_file_t
  use hex, only : edge_nodes, edge_faces, node_faces
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
  use device_coupled_vector_bc_resolver, only : &
       device_coupled_vector_bc_resolver_apply
  use scalar_bc_resolver, only : scalar_bc_resolver_t
  use, intrinsic :: iso_c_binding, only : c_ptr, c_null_ptr, c_associated
  implicit none
  private

  public ::  vector_bc_resolver_components

  !> Abstract type for resolving vector boundary conditions.
  type, public, abstract :: vector_bc_resolver_t
   contains
     procedure(vector_bc_resolver_init_intrf), pass(this), deferred :: init
     procedure(vector_bc_resolver_free_intrf), pass(this), deferred :: free
     procedure(vector_bc_resolver_finalize_intrf), pass(this), deferred :: &
          finalize
     procedure(vector_bc_resolver_apply_intrf), pass(this), deferred :: apply
     procedure(vector_bc_resolver_mark_bc_intrf), pass(this), deferred :: &
          mark_bc
     procedure(vector_bc_resolver_mark_bc_list_intrf), pass(this), deferred :: &
          mark_bc_list
     generic :: mark => mark_bc, mark_bc_list
  end type vector_bc_resolver_t

  !> A resolver for vector fields that acts component-wise.
  type, public, extends(vector_bc_resolver_t) :: segregated_vector_bc_resolver_t
     type(scalar_bc_resolver_t) :: x
     type(scalar_bc_resolver_t) :: y
     type(scalar_bc_resolver_t) :: z
   contains
     procedure, pass(this) :: init => segregated_vector_bc_resolver_init
     procedure, pass(this) :: free => segregated_vector_bc_resolver_free
     procedure, pass(this) :: finalize => segregated_vector_bc_resolver_finalize
     procedure, pass(this) :: mark_bc => segregated_vector_bc_resolver_mark_bc
     procedure, pass(this) :: mark_bc_list => &
          segregated_vector_bc_resolver_mark_bc_list
     procedure, pass(this) :: apply => segregated_vector_bc_resolver_apply
  end type segregated_vector_bc_resolver_t

  !> A coupled resolver for vector fields, suitable for mixed boundary
  !! conditions.
  type, public, extends(vector_bc_resolver_t) :: coupled_vector_bc_resolver_t
     type(mask_t) :: dirichlet_dof_mask
     type(mask_t) :: mixed_dof_mask
     type(bc_list_t), private :: bcs
     type(coef_t), pointer, private :: coef => null()
     type(dofmap_t), pointer, private :: dof => null()
     integer, allocatable :: node_rst(:,:)
     integer, allocatable :: edge_mid_rst(:,:)
     integer, allocatable :: node_linear_idx(:)
     real(kind=rp), allocatable :: face_class(:,:)
     integer, allocatable :: constraint_n(:)
     integer, allocatable :: constraint_t1(:)
     integer, allocatable :: constraint_t2(:)
     real(kind=rp), allocatable :: n(:,:)
     real(kind=rp), allocatable :: t1(:,:)
     real(kind=rp), allocatable :: t2(:,:)

     type(c_ptr) :: constraint_n_d = c_null_ptr
     type(c_ptr) :: constraint_t1_d = c_null_ptr
     type(c_ptr) :: constraint_t2_d = c_null_ptr
     type(c_ptr) :: n_d = c_null_ptr
     type(c_ptr) :: t1_d = c_null_ptr
     type(c_ptr) :: t2_d = c_null_ptr
   contains
     procedure, pass(this) :: free => coupled_vector_bc_resolver_free
     procedure, pass(this) :: init => coupled_vector_bc_resolver_init
     procedure, pass(this) :: mark_bc => coupled_vector_bc_resolver_mark_bc
     procedure, pass(this) :: mark_bc_list => &
          coupled_vector_bc_resolver_mark_bc_list
     procedure, pass(this) :: finalize => coupled_vector_bc_resolver_finalize
     procedure, pass(this) :: apply => coupled_vector_bc_resolver_apply
  end type coupled_vector_bc_resolver_t

  abstract interface
     subroutine vector_bc_resolver_free_intrf(this)
       import :: vector_bc_resolver_t
       class(vector_bc_resolver_t), intent(inout) :: this
     end subroutine vector_bc_resolver_free_intrf
  end interface

  abstract interface
     subroutine vector_bc_resolver_init_intrf(this, coef)
       import :: vector_bc_resolver_t, coef_t
       class(vector_bc_resolver_t), intent(inout) :: this
       type(coef_t), target, intent(in) :: coef
     end subroutine vector_bc_resolver_init_intrf
  end interface

  abstract interface
     subroutine vector_bc_resolver_finalize_intrf(this)
       import :: vector_bc_resolver_t
       class(vector_bc_resolver_t), intent(inout) :: this
     end subroutine vector_bc_resolver_finalize_intrf
  end interface

  abstract interface
     subroutine vector_bc_resolver_apply_intrf(this, x, y, z, n, strm)
       import :: vector_bc_resolver_t, rp, c_ptr
       class(vector_bc_resolver_t), intent(in) :: this
       integer, intent(in) :: n
       real(kind=rp), intent(inout) :: x(n)
       real(kind=rp), intent(inout) :: y(n)
       real(kind=rp), intent(inout) :: z(n)
       type(c_ptr), intent(inout), optional :: strm
     end subroutine vector_bc_resolver_apply_intrf
  end interface

  abstract interface
     subroutine vector_bc_resolver_mark_bc_intrf(this, bc, component)
       import :: vector_bc_resolver_t, bc_t
       class(vector_bc_resolver_t), intent(inout) :: this
       class(bc_t), intent(inout), target :: bc
       character(len=1), optional, intent(in) :: component
     end subroutine vector_bc_resolver_mark_bc_intrf
  end interface

  abstract interface
     subroutine vector_bc_resolver_mark_bc_list_intrf(this, bclst, component)
       import :: vector_bc_resolver_t, bc_list_t
       class(vector_bc_resolver_t), intent(inout) :: this
       type(bc_list_t), intent(in) :: bclst
       character(len=1), optional, intent(in) :: component
     end subroutine vector_bc_resolver_mark_bc_list_intrf
  end interface

contains

  subroutine segregated_vector_bc_resolver_free(this)
    class(segregated_vector_bc_resolver_t), intent(inout) :: this
    call this%x%free()
    call this%y%free()
    call this%z%free()
  end subroutine segregated_vector_bc_resolver_free

  subroutine segregated_vector_bc_resolver_init(this, coef)
    class(segregated_vector_bc_resolver_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
  end subroutine segregated_vector_bc_resolver_init

  subroutine segregated_vector_bc_resolver_finalize(this)
    class(segregated_vector_bc_resolver_t), intent(inout) :: this
  end subroutine segregated_vector_bc_resolver_finalize

  subroutine segregated_vector_bc_resolver_mark_bc(this, bc, component)
    class(segregated_vector_bc_resolver_t), intent(inout) :: this
    class(bc_t), intent(inout), target :: bc
    character(len=1), optional, intent(in) :: component

    if (.not. all(bc%constraints)) then
       call neko_error("Segregated vector BC resolver only accepts " // &
            "BCs with constraints = (.true., .true., .true.).")
    end if

    if (.not. present(component)) then
       call this%x%mark_bc(bc)
       call this%y%mark_bc(bc)
       call this%z%mark_bc(bc)
    else
       select case (component)
       case ('x')
          call this%x%mark_bc(bc)
       case ('y')
          call this%y%mark_bc(bc)
       case ('z')
          call this%z%mark_bc(bc)
       case default
          call neko_error("Invalid component for segregated vector BC " // &
               "resolver mark.")
       end select
    end if
  end subroutine segregated_vector_bc_resolver_mark_bc

  subroutine segregated_vector_bc_resolver_mark_bc_list(this, bclst, component)
    class(segregated_vector_bc_resolver_t), intent(inout) :: this
    type(bc_list_t), intent(in) :: bclst
    class(bc_t), pointer :: bc_i
    character(len=1), optional, intent(in) :: component
    integer :: i

    do i = 1, bclst%size()
       bc_i => bclst%get(i)
       call this%mark_bc(bc_i, component)
    end do
  end subroutine segregated_vector_bc_resolver_mark_bc_list

  subroutine segregated_vector_bc_resolver_apply(this, x, y, z, n, strm)
    class(segregated_vector_bc_resolver_t), intent(in) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: x(n)
    real(kind=rp), intent(inout) :: y(n)
    real(kind=rp), intent(inout) :: z(n)
    type(c_ptr), intent(inout), optional :: strm

    call this%x%apply(x, n, strm = strm)
    call this%y%apply(y, n, strm = strm)
    call this%z%apply(z, n, strm = strm)
  end subroutine segregated_vector_bc_resolver_apply

  subroutine vector_bc_resolver_components(this, x, y, z)
    class(vector_bc_resolver_t), target, intent(inout) :: this
    type(scalar_bc_resolver_t), pointer, intent(inout) :: x
    type(scalar_bc_resolver_t), pointer, intent(inout) :: y
    type(scalar_bc_resolver_t), pointer, intent(inout) :: z

    select type (this)
    type is (segregated_vector_bc_resolver_t)
       x => this%x
       y => this%y
       z => this%z
    class default
       call neko_error("Component access is only available for " // &
            "segregated vector BC resolvers. You have likely forgotten to " // &
            "select a coupled linear solver for velocity in the fluid " // &
            "configuration.")
    end select
  end subroutine vector_bc_resolver_components


  !
  ! Coupled resolver TBPs
  !

  subroutine coupled_vector_bc_resolver_free(this)
    class(coupled_vector_bc_resolver_t), intent(inout) :: this

    call this%dirichlet_dof_mask%free()
    call this%mixed_dof_mask%free()
    call this%bcs%free()

    if (allocated(this%node_rst)) deallocate(this%node_rst)
    if (allocated(this%edge_mid_rst)) deallocate(this%edge_mid_rst)
    if (allocated(this%node_linear_idx)) deallocate(this%node_linear_idx)
    if (allocated(this%face_class)) deallocate(this%face_class)

    if (allocated(this%constraint_n)) then
       if (NEKO_BCKND_DEVICE .eq. 1 .and. c_associated(this%constraint_n_d)) then
          call device_unmap(this%constraint_n, this%constraint_n_d)
       end if
       deallocate(this%constraint_n)
    end if
    if (allocated(this%constraint_t1)) then
       if (NEKO_BCKND_DEVICE .eq. 1 .and. c_associated(this%constraint_t1_d)) then
          call device_unmap(this%constraint_t1, this%constraint_t1_d)
       end if
       deallocate(this%constraint_t1)
    end if
    if (allocated(this%constraint_t2)) then
       if (NEKO_BCKND_DEVICE .eq. 1 .and. c_associated(this%constraint_t2_d)) then
          call device_unmap(this%constraint_t2, this%constraint_t2_d)
       end if
       deallocate(this%constraint_t2)
    end if
    if (allocated(this%n)) then
       if (NEKO_BCKND_DEVICE .eq. 1 .and. c_associated(this%n_d)) then
          call device_unmap(this%n, this%n_d)
       end if
       deallocate(this%n)
    end if
    if (allocated(this%t1)) then
       if (NEKO_BCKND_DEVICE .eq. 1 .and. c_associated(this%t1_d)) then
          call device_unmap(this%t1, this%t1_d)
       end if
       deallocate(this%t1)
    end if
    if (allocated(this%t2)) then
       if (NEKO_BCKND_DEVICE .eq. 1 .and. c_associated(this%t2_d)) then
          call device_unmap(this%t2, this%t2_d)
       end if
       deallocate(this%t2)
    end if

    this%constraint_n_d = c_null_ptr
    this%constraint_t1_d = c_null_ptr
    this%constraint_t2_d = c_null_ptr
    this%n_d = c_null_ptr
    this%t1_d = c_null_ptr
    this%t2_d = c_null_ptr
    nullify(this%coef)
    nullify(this%dof)
  end subroutine coupled_vector_bc_resolver_free

  subroutine coupled_vector_bc_resolver_init(this, coef)
    class(coupled_vector_bc_resolver_t), intent(inout) :: this
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
    allocate(this%face_class(nface, coef%msh%nelv))
    this%face_class = 5.0_rp

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
  end subroutine coupled_vector_bc_resolver_init

  subroutine coupled_vector_bc_resolver_mark_bc(this, bc, component)
    class(coupled_vector_bc_resolver_t), intent(inout) :: this
    class(bc_t), intent(inout), target :: bc
    character(len=1), optional, intent(in) :: component

    if (.not. associated(this%coef)) then
       call neko_error("Coupled vector BC resolver must be initialized " // &
            "before mark().")
    end if

    call this%bcs%append(bc)
  end subroutine coupled_vector_bc_resolver_mark_bc

  subroutine coupled_vector_bc_resolver_mark_bc_list(this, bclst, component)
    class(coupled_vector_bc_resolver_t), intent(inout) :: this
    type(bc_list_t), intent(in) :: bclst
    character(len=1), optional, intent(in) :: component
    class(bc_t), pointer :: bc_i
    integer :: i

    do i = 1, bclst%size()
       bc_i => bclst%get(i)
       call this%mark_bc(bc_i, component)
    end do
  end subroutine coupled_vector_bc_resolver_mark_bc_list
  !> Finalize the coupled resolver by resolving the accumulated BC list.
  !! @details This routine builds the dof masks for Dirichlet and mixed nodes,
  !! and computes the local basis for the mixed ones.
  subroutine coupled_vector_bc_resolver_finalize(this)
    class(coupled_vector_bc_resolver_t), intent(inout) :: this
    type(field_t), pointer :: work1
    type(field_t), pointer :: node_class
    type(field_t), pointer :: work3
    type(field_t), pointer :: work4
    type(tuple_i4_t), pointer :: marked_faces(:)
    type(tuple_i4_t) :: marked_face
    integer :: scratch_idx(4)
    integer, allocatable :: dirichlet_mask_values(:)
    integer, allocatable :: mixed_mask_values(:)
    integer, allocatable :: resolved_mask_values(:)
    integer, allocatable :: dof_to_mixed_idx(:)
    integer :: i, j, k, dof_size, m
    integer :: dirichlet_mask_size, mixed_mask_size, resolved_mask_size
    integer :: idx(4), facet, el, edge, node, ii, p
    integer :: rst(3), rst1(3), rst2(3), step_rst(3)
    integer :: edge_len, edge_idx, node_idx
    real(kind=rp) :: normal(3), t1_vec(3), t2_vec(3), len, prio
    class(bc_t), pointer :: bc

    call this%dirichlet_dof_mask%free()
    call this%mixed_dof_mask%free()

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
    if (allocated(this%n)) then
       if (NEKO_BCKND_DEVICE .eq. 1 .and. c_associated(this%n_d)) then
          call device_unmap(this%n, this%n_d)
       end if
       deallocate(this%n)
    end if
    if (allocated(this%t1)) then
       if (NEKO_BCKND_DEVICE .eq. 1 .and. c_associated(this%t1_d)) then
          call device_unmap(this%t1, this%t1_d)
       end if
       deallocate(this%t1)
    end if
    if (allocated(this%t2)) then
       if (NEKO_BCKND_DEVICE .eq. 1 .and. c_associated(this%t2_d)) then
          call device_unmap(this%t2, this%t2_d)
       end if
       deallocate(this%t2)
    end if

    if (this%bcs%size() .eq. 0) return

    call neko_scratch_registry%request_field(work1, scratch_idx(1), .true.)
    call neko_scratch_registry%request_field(node_class, scratch_idx(2), .true.)
    call neko_scratch_registry%request_field(work3, scratch_idx(3), .true.)
    call neko_scratch_registry%request_field(work4, scratch_idx(4), .true.)

    dof_size = this%dof%size()
    this%face_class = 5.0_rp

    ! Build a mask of all dofs on the boundary.
    do i = 1, this%bcs%size()
       bc => this%bcs%get(i)
       write(*,*) "Processing BC ", i, " called ", bc%name, &
            " and contraints, ", bc%constraints

       if (.not. allocated(bc%msk)) then
          call neko_error("Attempting to finalize coupled resolver " // &
               "unfinalized BC.")
       end if

       ! Mask all the dofs touched by this BC. Since %msk is propagated to all
       ! local dofs via gather-scatter, work1 will contain all local nodes on
       ! the boundary, including those elements that don't touch it with a face.
       ! Note that bc%msk stores its length in slot 0 and is therefore passed
       ! with the `_0` masked-wrapper convention.
       call cfill_mask(work1%x, 1.0_rp, dof_size, bc%msk(1:bc%msk(0)), &
            bc%msk(0))
    end do


    ! Set priority values (class) for constraint assignment. Mimics the
    ! procedure in Nek5000 directly.
    ! The values are chosen in a way that a min reduction applies the most
    ! restrictive constraint.
    ! 5 -> unconstrained
    ! 3 -> tangentially constrained
    ! 2 -> normally constrained
    ! 0 -> fully constrained
    ! Fill the field to not mess up gather-scatter reduction later.
    call cfill(node_class%x, 5.0_rp, dof_size)
    do i = 1, this%bcs%size()
       bc => this%bcs%get(i)


       prio = bc%bc_type

       ! Store the class on each boundary face touched by this BC.
       ! This is the compact analogue of Nek's face-resident HFMASK field:
       ! one scalar class value per local (facet, element) pair.
       marked_faces => bc%marked_facet%array()
       do j = 1, bc%marked_facet%size()
          marked_face = marked_faces(j)
          facet = marked_face%x(1)
          el = marked_face%x(2)
          this%face_class(facet, el) = prio
       end do

       ! Note that the face_msk is used, so constraints are only directly
       ! applied to elements that touch the boundary at this point.
       ! The min here ensures that the most restricive constraint is kept
       ! within a single element.
       do j = 1, bc%facet_msk(0)
          m = bc%facet_msk(j)
          node_class%x(m,1,1,1) = min(prio, node_class%x(m,1,1,1))
       end do
    end do

    ! Propagate constraints to all local dofs via gather-scatter.
    ! Ensures most restrictive constraint is kept across element boundaries.
    call this%coef%gs_h%op(node_class, GS_OP_MIN)

    ! For mixed BCs, build a resolved subset of the original bc%msk support.
    ! A dof survives in the resolved mask only if the globally reduced class
    ! still matches the class semantics of that BC. Nodes that were touched by
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

          prio = bc%bc_type

          resolved_mask_size = 0
          do j = 1, bc%msk(0)
             k = bc%msk(j)
             if (abs(node_class%x(k,1,1,1) - prio) .lt. 1.0e-6_rp) then
                resolved_mask_size = resolved_mask_size + 1
             end if
          end do

          allocate(resolved_mask_values(resolved_mask_size))
          resolved_mask_size = 0
          do j = 1, bc%msk(0)
             k = bc%msk(j)
             if (abs(node_class%x(k,1,1,1) - prio) .lt. 1.0e-6_rp) then
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
    do i = 1, dof_size
       if (work1%x(i,1,1,1) .lt. 0.5_rp) cycle

       if (node_class%x(i,1,1,1) .lt. 1.9_rp) then
          dirichlet_mask_size = dirichlet_mask_size + 1
       else if (node_class%x(i,1,1,1) .gt. 1.9_rp .and. &
            node_class%x(i,1,1,1) .lt. 3.9_rp) then
          mixed_mask_size = mixed_mask_size + 1
       end if
    end do


    if (dirichlet_mask_size .gt. 0) then
       allocate(dirichlet_mask_values(dirichlet_mask_size))
    end if
    if (mixed_mask_size .gt. 0) then
       allocate(mixed_mask_values(mixed_mask_size))
    end if

    ! We reuse the variables as counters in the loop below, which actually
    ! fills the masks, now that we know the size.
    dirichlet_mask_size = 0
    mixed_mask_size = 0
    do i = 1, dof_size
       if (work1%x(i,1,1,1) .lt. 0.5_rp) cycle

       if (node_class%x(i,1,1,1) .lt. 1.9_rp) then
          dirichlet_mask_size = dirichlet_mask_size + 1
          dirichlet_mask_values(dirichlet_mask_size) = i
       else if (node_class%x(i,1,1,1) .gt. 1.9_rp .and. &
            node_class%x(i,1,1,1) .lt. 3.9_rp) then
          mixed_mask_size = mixed_mask_size + 1
          mixed_mask_values(mixed_mask_size) = i
       end if
    end do

    ! We init even if the size can be zero. Should be OK.
    call this%dirichlet_dof_mask%init(dirichlet_mask_values, &
         dirichlet_mask_size)
    call this%mixed_dof_mask%init(mixed_mask_values, mixed_mask_size)

    ! A mapping between the field linear index of a mixed node into its index
    ! in the mixed_dof_mask.
    allocate(dof_to_mixed_idx(dof_size))
    dof_to_mixed_idx = 0
    do i = 1, mixed_mask_size
       dof_to_mixed_idx(mixed_mask_values(i)) = i
    end do

    ! Allocate mixed-node constraints and fill them from the reduced class.
    allocate(this%constraint_n(mixed_mask_size))
    allocate(this%constraint_t1(mixed_mask_size))
    allocate(this%constraint_t2(mixed_mask_size))
    allocate(this%n(3, mixed_mask_size))
    allocate(this%t1(3, mixed_mask_size))
    allocate(this%t2(3, mixed_mask_size))
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%constraint_n, this%constraint_n_d, &
            size(this%constraint_n))
       call device_map(this%constraint_t1, this%constraint_t1_d, &
            size(this%constraint_t1))
       call device_map(this%constraint_t2, this%constraint_t2_d, &
            size(this%constraint_t2))
       call device_map(this%n, this%n_d, size(this%n))
       call device_map(this%t1, this%t1_d, size(this%t1))
       call device_map(this%t2, this%t2_d, size(this%t2))
    end if


    do i = 1, mixed_mask_size
       j = mixed_mask_values(i)

       if (node_class%x(j,1,1,1) .lt. 1.9_rp) then
          this%constraint_n(i) = 1
          this%constraint_t1(i) = 1
          this%constraint_t2(i) = 1
       else if (node_class%x(j,1,1,1) .gt. 1.9_rp .and. &
            node_class%x(j,1,1,1) .lt. 2.9_rp) then
          this%constraint_n(i) = 1
          this%constraint_t1(i) = 0
          this%constraint_t2(i) = 0
       else if (node_class%x(j,1,1,1) .gt. 2.9_rp .and. &
            node_class%x(j,1,1,1) .lt. 3.9_rp) then
          this%constraint_n(i) = 0
          this%constraint_t1(i) = 1
          this%constraint_t2(i) = 1
       else if (node_class%x(j,1,1,1) .gt. 3.9_rp) then
          this%constraint_n(i) = 0
          this%constraint_t1(i) = 0
          this%constraint_t2(i) = 0
       end if
    end do

    write(*,*) "Filled mixed node constraints in coupled vector BC resolver."

    ! At this point, face_class stores the face-based bc class values, and
    ! node_class stores them node-wise, after propagation with min reduction.
    ! Note that the propagation means that some nodes may have a different,
    ! more restrictive class than owning face!
    ! The algorithm for constructing normals below will make use of both
    ! classifications when looking at edges and corners. We will really only
    ! care about classes 2 and 3, i.e. mixed bcs. The key question will be
    ! whether a given face should contribute its normal to the edge and corner
    ! dofs.


    ! Set normals at unambiguous boundary dofs based on face normals.
    associate(nx => work1, ny => work3, nz => work4)
      call rzero(nx%x, dof_size)
      call rzero(ny%x, dof_size)
      call rzero(nz%x, dof_size)

      this%n = 0.0_rp
      this%t1 = 0.0_rp
      this%t2 = 0.0_rp

      ! First pass: seed the local normal field on all nodes that lie
      ! on directly marked mixed faces. This is the compact analogue of the
      ! SETCSYS face sweep in Nek5000 before edge and corner reconstruction.
      ! Since several faces will own edge and corner nodes, the normals there
      ! will be overwritten in arbitrary order, but we don't care because we
      ! will reset those later and treat them specially.
      do i = 1, this%bcs%size()
         bc => this%bcs%get(i)

         do j = 1, bc%facet_msk(0)
            ! Global linear index of the node on which to set the normal.
            k = bc%facet_msk(j)

            ! Grab face and ijke indices to address face_class and get_normal.
            facet = bc%facet(j)
            idx = nonlinear_index(k, this%coef%Xh%lx, this%coef%Xh%ly, &
                 this%coef%Xh%lz)

            if (this%face_class(facet, idx(4)) .lt. 1.9_rp .or. &
                 this%face_class(facet, idx(4)) .gt. 3.1_rp) cycle

            normal = this%coef%get_normal(idx(1), idx(2), idx(3), idx(4), &
                 facet)
            nx%x(k,1,1,1) = normal(1)
            ny%x(k,1,1,1) = normal(2)
            nz%x(k,1,1,1) = normal(3)
         end do
      end do

      write(*,*) "Seeded normals at directly marked mixed nodes in " // &
           "coupled vector BC resolver."

      ! We now treat the special edges and conrners. Everything is done locally
      ! per element, using reference element address tables found in hex.f90
      ! and inside this type.

      ! Mixed edge interiors are rebuilt from the normals of the adjacent
      ! faces whose local face class matches the reduced nodal class.
      ! This is the central point: if the adjacent face is a different class,
      ! which by construction can only be a less restrictive class, then it
      ! should not contribute its normal.
      !
      ! Consider the following 2D example. In 2D an edge becomes a node in the
      ! corner of the element. Look at the node marked with X. After the nodal
      ! class is propagated, it will have class 2---the most restrictive of the
      ! adjacent. So, only the face with class 2 in El 2 will contribute to the
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
      if (mixed_mask_size .gt. 0) then
         do el = 1, this%coef%msh%nelv
            do edge = 1, size(edge_nodes, 2)

               ! Representitive rst index in the middle of an edge.
               rst = this%edge_mid_rst(:, edge)

               ! Global linear index of the midpoint node.
               edge_idx = linear_index(rst(1), rst(2), rst(3), el, &
                    this%coef%Xh%lx, this%coef%Xh%ly, this%coef%Xh%lz)

               ! Get the prio class.
               prio = abs(node_class%x(edge_idx,1,1,1))

               ! If this is not a mixed bc edge, just leave it alone.
               if (prio .lt. 1.9_rp .or. prio .gt. 3.1_rp) cycle

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
               ! This is a tripe, but only one component is nonzero.
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
                  nx%x(k,1,1,1) = 0.0_rp
                  ny%x(k,1,1,1) = 0.0_rp
                  nz%x(k,1,1,1) = 0.0_rp
               end do

               ! Loop over the faces adjacent to this edge and add the normals
               ! if the bc prior class matches between the face and the edge.
               do ii = 1, size(edge_faces, 1)
                  facet = edge_faces(ii, edge)

                  ! Skip if prio class is not the same.
                  if (abs(prio - this%face_class(facet, el)) .gt. &
                       1.0e-6_rp) then
                     cycle
                  end if

                  ! Loop over the interior edge nodes again and add the normals.
                  do p = 2, edge_len - 1
                     rst = rst1 + (p - 1) * step_rst
                     k = linear_index(rst(1), rst(2), rst(3), el, &
                          this%coef%Xh%lx, this%coef%Xh%ly, this%coef%Xh%lz)
                     normal = this%coef%get_normal(rst(1), rst(2), rst(3), &
                          el, facet)
                     nx%x(k,1,1,1) = nx%x(k,1,1,1) + normal(1)
                     ny%x(k,1,1,1) = ny%x(k,1,1,1) + normal(2)
                     nz%x(k,1,1,1) = nz%x(k,1,1,1) + normal(3)
                  end do
               end do
            end do
         end do

         write(*,*) "Finished reconstructing normals at mixed edges in " // &
              "coupled vector BC resolver."

         ! Mixed corner node normals are rebuilt from the adjacent faces whose
         ! local face class matches the reduced nodal class at that node.
         do el = 1, this%coef%msh%nelv
            do node = 1, size(this%node_linear_idx)
               rst = this%node_rst(:, node)
               node_idx = linear_index(rst(1), rst(2), rst(3), el, &
                    this%coef%Xh%lx, this%coef%Xh%ly, this%coef%Xh%lz)
               prio = abs(node_class%x(node_idx,1,1,1))

               ! Ignore if the bc class is not a mixed one.
               if (prio .lt. 1.9_rp .or. prio .gt. 3.1_rp) cycle

               ! Kill the normal to start fresh.
               nx%x(node_idx,1,1,1) = 0.0_rp
               ny%x(node_idx,1,1,1) = 0.0_rp
               nz%x(node_idx,1,1,1) = 0.0_rp

               ! Note, 3 faces share a corner node in 3D.
               do ii = 1, this%coef%msh%gdim
                  facet = node_faces(ii, node)

                  ! Check class agreement
                  if (abs(prio - this%face_class(facet, el)) .gt. &
                       1.0e-6_rp) then
                     cycle
                  end if

                  ! Add the normal.
                  normal = this%coef%get_normal(rst(1), rst(2), rst(3), &
                       el, facet)
                  nx%x(node_idx,1,1,1) = nx%x(node_idx,1,1,1) + normal(1)
                  ny%x(node_idx,1,1,1) = ny%x(node_idx,1,1,1) + normal(2)
                  nz%x(node_idx,1,1,1) = nz%x(node_idx,1,1,1) + normal(3)
               end do
            end do
         end do
      end if

      write(*,*) "Finished reconstructing normals at mixed corners in " // &
           "coupled vector BC resolver."

      ! We are done element-wise. Now we can just sum the normals across nodes
      ! shared by multiple elements.
      call this%coef%gs_h%op(nx, GS_OP_ADD)
      call this%coef%gs_h%op(ny, GS_OP_ADD)
      call this%coef%gs_h%op(nz, GS_OP_ADD)

      do i = 1, mixed_mask_size
         j = mixed_mask_values(i)

         ! Normalize the normal
         normal(1) = nx%x(j,1,1,1)
         normal(2) = ny%x(j,1,1,1)
         normal(3) = nz%x(j,1,1,1)
         len = sqrt(sum(normal**2))
         if (len .le. 0.0_rp) cycle

         this%n(:,i) = normal / len

         ! This selects the first tangent direction.
         ! Generally, we pick a perpendicular to the normal in the x-y plane,
         ! but if the normal is almost aligned with z, that formula degenerates
         ! to 0, so we just choose x as the first tangent direction.
         if (abs(this%n(3,i)) .gt. 0.999_rp) then
            this%t1(:,i) = [ 1.0_rp, 0.0_rp, 0.0_rp ]
         else
            t1_vec = [ -this%n(2,i), this%n(1,i), 0.0_rp ]
            len = sqrt(sum(t1_vec**2))
            if (len .gt. 0.0_rp) then
               this%t1(:,i) = t1_vec / len
            end if
         end if

         ! Get t2 as a cross product of n and t1.
         t2_vec(1) = this%n(2,i) * this%t1(3,i) - &
              this%n(3,i) * this%t1(2,i)
         t2_vec(2) = this%n(3,i) * this%t1(1,i) - &
              this%n(1,i) * this%t1(3,i)
         t2_vec(3) = this%n(1,i) * this%t1(2,i) - &
              this%n(2,i) * this%t1(1,i)
         len = sqrt(sum(t2_vec**2))
         if (len .gt. 0.0_rp) then
            this%t2(:,i) = t2_vec / len
         end if
      end do
    end associate

    write(*,*) "Finished building local basis for coupled vector BC resolver."

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%constraint_n, this%constraint_n_d, &
            size(this%constraint_n), HOST_TO_DEVICE, sync = .true.)
       call device_memcpy(this%constraint_t1, this%constraint_t1_d, &
            size(this%constraint_t1), HOST_TO_DEVICE, sync = .true.)
       call device_memcpy(this%constraint_t2, this%constraint_t2_d, &
            size(this%constraint_t2), HOST_TO_DEVICE, sync = .true.)
       call device_memcpy(this%n, this%n_d, size(this%n), HOST_TO_DEVICE, &
            sync = .true.)
       call device_memcpy(this%t1, this%t1_d, size(this%t1), HOST_TO_DEVICE, &
            sync = .true.)
       call device_memcpy(this%t2, this%t2_d, size(this%t2), HOST_TO_DEVICE, &
            sync = .true.)
    end if

    ! Transfer the final mixed-node basis into each mixed BC on its
    ! resolved support, so strong application on the physical field can use
    ! BC-local data rather than the global resolver internals.
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
                     "the coupled resolver mixed basis.")
             end if

             bc%n%x(:,j) = this%n(:,p)
             bc%t1%x(:,j) = this%t1(:,p)
             bc%t2%x(:,j) = this%t2(:,p)
          end do

          if (NEKO_BCKND_DEVICE .eq. 1) then
             call bc%n%copy_from(HOST_TO_DEVICE, .true.)
             call bc%t1%copy_from(HOST_TO_DEVICE, .true.)
             call bc%t2%copy_from(HOST_TO_DEVICE, .true.)
          end if
       end select
    end do

    block
      use device_math, only : device_cfill, device_cfill_mask
      type(field_list_t) :: basis_fields
      type(fld_file_t) :: basis_file

      call rzero(work1%x, dof_size)
      call rzero(work3%x, dof_size)
      call rzero(work4%x, dof_size)

      if (this%mixed_dof_mask%is_set()) then
         call masked_scatter_copy(work1%x(:,1,1,1), this%n(1,:), &
              mixed_mask_values, dof_size, mixed_mask_size)
         call masked_scatter_copy(work3%x(:,1,1,1), this%n(2,:), &
              mixed_mask_values, dof_size, mixed_mask_size)
         call masked_scatter_copy(work4%x(:,1,1,1), this%n(3,:), &
              mixed_mask_values, dof_size, mixed_mask_size)
      end if

      call basis_fields%init(3)
      call basis_fields%assign(1, work1)
      call basis_fields%assign(2, work3)
      call basis_fields%assign(3, work4)


      call basis_file%init('bc_resolver_basis.fld')
      call basis_file%write(basis_fields)
      call basis_fields%free()

      if (NEKO_BCKND_DEVICE .eq. 1) then
         call device_cfill(work1%x_d, 5.0_rp, dof_size)

         if (this%dirichlet_dof_mask%is_set()) then
            call device_cfill_mask(work1%x_d, 1.0_rp, dof_size, &
                 this%dirichlet_dof_mask%get_d(), &
                 this%dirichlet_dof_mask%size())
         end if

         if (this%mixed_dof_mask%is_set()) then
            call device_cfill_mask(work1%x_d, 2.0_rp, dof_size, &
                 this%mixed_dof_mask%get_d(), &
                 this%mixed_dof_mask%size())
         end if

         call device_memcpy(work1%x, work1%x_d, dof_size, DEVICE_TO_HOST, &
              sync = .true.)

      else
         call cfill(work1%x, 5.0_rp, dof_size)
         call cfill_mask(work1%x, 1.0_rp, dof_size, &
              this%dirichlet_dof_mask%get(), &
              this%dirichlet_dof_mask%size())
         call cfill_mask(work1%x, 2.0_rp, dof_size, this%mixed_dof_mask%get(), &
              this%mixed_dof_mask%size())
      end if

      call basis_file%init('bc_resolver_mask.fld')
      call basis_file%write(work1)
    end block

    call neko_scratch_registry%relinquish_field(scratch_idx)

    if (allocated(dirichlet_mask_values)) deallocate(dirichlet_mask_values)
    if (allocated(mixed_mask_values)) deallocate(mixed_mask_values)
    if (allocated(dof_to_mixed_idx)) deallocate(dof_to_mixed_idx)
  end subroutine coupled_vector_bc_resolver_finalize

  !> Apply the coupled vector boundary constraints in the local basis.
  subroutine coupled_vector_bc_resolver_apply(this, x, y, z, n, strm)
    class(coupled_vector_bc_resolver_t), intent(in) :: this
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
          call device_coupled_vector_bc_resolver_apply( &
               this%mixed_dof_mask%get_d(), x_d, y_d, z_d, &
               this%constraint_n_d, this%constraint_t1_d, &
               this%constraint_t2_d, this%n_d, this%t1_d, &
               this%t2_d, m, strm_)
       else
          mixed_msk => this%mixed_dof_mask%get()

          do i = 1, m
             j = mixed_msk(i)

             u(1) = x(j)
             u(2) = y(j)
             u(3) = z(j)

             uloc(1) = u(1) * this%n(1,i) + u(2) * this%n(2,i) + &
                  u(3) * this%n(3,i)
             uloc(2) = u(1) * this%t1(1,i) + u(2) * this%t1(2,i) + &
                  u(3) * this%t1(3,i)
             uloc(3) = u(1) * this%t2(1,i) + u(2) * this%t2(2,i) + &
                  u(3) * this%t2(3,i)

             if (this%constraint_n(i) .ne. 0) uloc(1) = 0.0_rp
             if (this%constraint_t1(i) .ne. 0) uloc(2) = 0.0_rp
             if (this%constraint_t2(i) .ne. 0) uloc(3) = 0.0_rp

             u = uloc(1) * this%n(:,i) + uloc(2) * this%t1(:,i) + &
                  uloc(3) * this%t2(:,i)

             x(j) = u(1)
             y(j) = u(2)
             z(j) = u(3)
          end do
       end if
    end if
  end subroutine coupled_vector_bc_resolver_apply

end module vector_bc_resolver
