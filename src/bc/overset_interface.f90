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
!> Defines overset interface scalar boundary conditions
module overset_interface
  use comm, only : NEKO_GLOBAL_COMM
  use neko_config, only : NEKO_BCKND_DEVICE
  use registry, only : neko_registry
  use num_types, only : rp
  use coefs, only : coef_t
  use global_interpolation, only : global_interpolation_t, &
       global_interpolation_settings_t
  use mask, only : mask_t
  use bc, only : bc_t, BC_TYPES
  use field_list, only : field_list_t
  use math, only : masked_copy_0, copy
  use device_math, only : device_masked_copy_0, device_copy
  use vector, only : vector_t
  use vector_math, only : vector_masked_gather_copy, vector_masked_scatter_copy
  use device, only : DEVICE_TO_HOST
  use field_dirichlet, only : field_dirichlet_t
  use utils, only : neko_error, nonlinear_index, linear_index
  use stack, only : stack_i4_t
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use field, only : field_t
  use, intrinsic :: iso_c_binding, only : c_ptr
  use time_state, only : time_state_t
  implicit none
  private

  !> Overset interface BC for a scalar field.
  type, public, extends(bc_t) :: overset_interface_t
     !> Underlying scalar field Dirichlet bc.
     type(field_dirichlet_t) :: bc_s
     !> Single-field list for compatibility with field-based update patterns.
     type(field_list_t) :: field_list
     !> Name of the scalar field to interpolate from the registry.
     character(len=:), allocatable :: field_name
     !> Interpolator
     type(global_interpolation_t) :: interface_interpolator
     !> Mask for overset interface points.
     type(mask_t) :: interface_dof_mask
     type(mask_t) :: domain_element_mask
     !> Vectors holding dof coordinates for cases where mesh moves.
     type(vector_t) :: x_dof, y_dof, z_dof
     type(vector_t) :: x_interface_dof, y_interface_dof, z_interface_dof
     !> Interpolated scalar values on the interface.
     type(vector_t) :: s_interface
     !> Interpolation settings.
     type(global_interpolation_settings_t) :: interpolation_settings
   contains
     !> Constructor.
     procedure, pass(this) :: init => overset_interface_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          overset_interface_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => overset_interface_free
     !> Finalize.
     procedure, pass(this) :: finalize => overset_interface_finalize
     !> Apply scalar by performing a masked copy.
     procedure, pass(this) :: apply_scalar => overset_interface_apply_scalar
     !> (No-op) Apply vector.
     procedure, pass(this) :: apply_vector => overset_interface_apply_vector
     !> (No-op) Apply vector (device).
     procedure, pass(this) :: apply_vector_dev => overset_interface_apply_vector_dev
     !> Apply scalar (device).
     procedure, pass(this) :: apply_scalar_dev => overset_interface_apply_scalar_dev
     procedure, pass(this) :: update => overset_interface_update

     !> Build domain masks for the overset interface.
     procedure, pass(this), private :: build_masks_ => build_masks_
     !> Gather the dofs at the interface.
     procedure, pass(this), private :: gather_interface_dofs_ => gather_interface_dofs_
     !> Set up the interpolator.
     procedure, pass(this), private :: setup_interpolator_ => setup_interpolator_
  end type overset_interface_t

contains

  !> Constructor
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine overset_interface_init(this, coef, json)
    class(overset_interface_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: field_name
    real(kind=rp) :: tol, pad

    call json_get(json, "field_name", field_name)
    call json_get_or_default(json, "interpolation.tolerance", tol, -1.0_rp)
    call json_get_or_default(json, "interpolation.padding", pad, -1.0_rp)
    call this%init_from_components(coef, field_name, tol, pad)
    if (allocated(field_name)) deallocate(field_name)

  end subroutine overset_interface_init

  !> Constructor from components
  !! @param[in] coef The SEM coefficients.
  subroutine overset_interface_init_from_components(this, coef, field_name, &
       tol, pad)
    class(overset_interface_t), intent(inout), target :: this
    type(coef_t), intent(in) :: coef
    character(len=*), intent(in) :: field_name
    real(kind=rp), intent(in), optional :: tol, pad

    call this%init_base(coef)
    this%bc_type = BC_TYPES%DIRICHLET

    if (present(tol)) then
       if (tol .gt. 0.0_rp) then
          this%interpolation_settings%tolerance = tol
       end if
    end if

    if (present(pad)) then
       if (pad .gt. 0.0_rp) then
          this%interpolation_settings%padding = pad
       end if
    end if

    this%field_name = field_name

    call this%bc_s%init_from_components(coef, this%field_name)
    call this%field_list%init(1)
    call this%field_list%assign_to_field(1, this%bc_s%field_bc)

    call this%x_dof%init(this%dof%size(), 'x')
    call this%y_dof%init(this%dof%size(), 'y')
    call this%z_dof%init(this%dof%size(), 'z')

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%x_dof%x_d, this%dof%x_d, this%dof%size())
       call device_copy(this%y_dof%x_d, this%dof%y_d, this%dof%size())
       call device_copy(this%z_dof%x_d, this%dof%z_d, this%dof%size())

       call this%x_dof%copy_from(DEVICE_TO_HOST, sync = .false.)
       call this%y_dof%copy_from(DEVICE_TO_HOST, sync = .false.)
       call this%z_dof%copy_from(DEVICE_TO_HOST, sync = .true.)
    else
       call copy(this%x_dof%x, this%dof%x, this%dof%size())
       call copy(this%y_dof%x, this%dof%y, this%dof%size())
       call copy(this%z_dof%x, this%dof%z, this%dof%size())
    end if

  end subroutine overset_interface_init_from_components

  !> Destructor.
  subroutine overset_interface_free(this)
    class(overset_interface_t), target, intent(inout) :: this

    call this%bc_s%free()
    call this%field_list%free()

    call this%x_dof%free()
    call this%y_dof%free()
    call this%z_dof%free()

    call this%x_interface_dof%free()
    call this%y_interface_dof%free()
    call this%z_interface_dof%free()
    call this%s_interface%free()

    if (allocated(this%field_name)) then
       deallocate(this%field_name)
    end if

    call this%interface_interpolator%free()

    call this%interface_dof_mask%free()
    call this%domain_element_mask%free()

    call this%free_base()
  end subroutine overset_interface_free

  !> Apply scalar.
  !! @param x Field onto which to copy the values.
  !! @param n Size of the array `x`.
  !! @param time The current time state.
  subroutine overset_interface_apply_scalar(this, x, n, time, strong)
    class(overset_interface_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (strong_) then
       if (.not. this%updated) then
          call this%update(time)
          this%updated = .true.
       end if

       call masked_copy_0(x, this%bc_s%field_bc%x, this%msk, n, this%msk(0))
    end if

  end subroutine overset_interface_apply_scalar

  !> Apply scalar (device).
  !! @param x_d Device pointer to the field onto which to copy the values.
  !! @param time The current time state.
  !! @param strm Device stream
  subroutine overset_interface_apply_scalar_dev(this, x_d, time, strong, strm)
    class(overset_interface_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (strong_) then
       if (.not. this%updated) then
          call this%update(time)
          this%updated = .true.
       end if

       if (this%msk(0) .gt. 0) then
          call device_masked_copy_0(x_d, this%bc_s%field_bc%x_d, this%bc_s%msk_d, &
               this%bc_s%dof%size(), this%msk(0), strm)
       end if
    end if

  end subroutine overset_interface_apply_scalar_dev

  !> (No-op) Apply vector.
  subroutine overset_interface_apply_vector(this, x, y, z, n, time, strong)
    class(overset_interface_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong

    call neko_error("overset_interface cannot apply vector BCs.&
    & Use overset_interface_vector instead!")

  end subroutine overset_interface_apply_vector

  !> (No-op) Apply vector (device).
  subroutine overset_interface_apply_vector_dev(this, x_d, y_d, z_d, time, &
       strong, strm)
    class(overset_interface_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(c_ptr), intent(inout) :: y_d
    type(c_ptr), intent(inout) :: z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm

    call neko_error("overset_interface cannot apply vector BCs.&
    & Use overset_interface_vector instead!")

  end subroutine overset_interface_apply_vector_dev

  !> Finalize by building the mask arrays and preparing interpolation data.
  subroutine overset_interface_finalize(this)
    class(overset_interface_t), target, intent(inout) :: this

    call this%finalize_base()

    call this%bc_s%mark_facets(this%marked_facet)
    call this%bc_s%finalize()

    call this%build_masks_()

    call this%x_interface_dof%init(this%interface_dof_mask%size(), 'x_interface')
    call this%y_interface_dof%init(this%interface_dof_mask%size(), 'y_interface')
    call this%z_interface_dof%init(this%interface_dof_mask%size(), 'z_interface')
    call this%gather_interface_dofs_()

    call this%setup_interpolator_()

    call this%s_interface%init(this%interface_dof_mask%size(), 's_interface')

  end subroutine overset_interface_finalize

  !> Update values at the overset interface.
  subroutine overset_interface_update(this, time)
    class(overset_interface_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: s

    s => neko_registry%get_field(trim(this%field_name))

    call this%interface_interpolator%evaluate_masked(this%s_interface%x, s%x, &
         this%domain_element_mask, .false.)

    call vector_masked_scatter_copy(this%bc_s%field_bc%x(:,1,1,1), this%s_interface, &
         this%interface_dof_mask, this%bc_s%dof%size())

  end subroutine overset_interface_update

  !===================
  ! Helper subroutines
  !===================

  !> Build masks.
  subroutine build_masks_(this)
    class(overset_interface_t), intent(inout) :: this
    type(mask_t) :: temp_mask
    logical, allocatable :: found(:)
    integer :: i, j, k, e, nelems
    integer :: lx, ly, lz
    integer :: nonlinear_idx(4), linear_idx
    type(stack_i4_t) :: idx_stack

    call this%interface_dof_mask%init(this%msk(1:this%msk(0)), this%msk(0))

    lx = this%Xh%lx
    ly = this%Xh%ly
    lz = this%Xh%lz

    allocate(found(this%msh%nelv))
    found = .false.

    do i = 1, this%msk(0)
       linear_idx = this%msk(i)
       nonlinear_idx = nonlinear_index(linear_idx, lx, ly, lz)
       found(nonlinear_idx(4)) = .true.
    end do

    nelems = 0
    call idx_stack%init()
    do e = 1, this%msh%nelv
       if (found(e)) then
          nelems = nelems + 1
          do k = 1, this%Xh%lz
             do j = 1, this%Xh%ly
                do i = 1, this%Xh%lx
                   linear_idx = linear_index(i, j, k, e, lx, ly, lz)
                   call idx_stack%push(linear_idx)
                end do
             end do
          end do
       end if
    end do

    deallocate(found)

    call temp_mask%init(idx_stack%array(), idx_stack%size())
    call idx_stack%free()

    call this%domain_element_mask%invert_mask(temp_mask, this%dof%size())
    call temp_mask%free()

  end subroutine build_masks_

  !> Gather interface dofs.
  subroutine gather_interface_dofs_(this)
    class(overset_interface_t), intent(inout) :: this

    call vector_masked_gather_copy(this%x_interface_dof, this%dof%x(:,1,1,1), &
         this%interface_dof_mask, this%dof%size())
    call vector_masked_gather_copy(this%y_interface_dof, this%dof%y(:,1,1,1), &
         this%interface_dof_mask, this%dof%size())
    call vector_masked_gather_copy(this%z_interface_dof, this%dof%z(:,1,1,1), &
         this%interface_dof_mask, this%dof%size())

    call this%x_interface_dof%copy_from(DEVICE_TO_HOST, sync = .false.)
    call this%y_interface_dof%copy_from(DEVICE_TO_HOST, sync = .false.)
    call this%z_interface_dof%copy_from(DEVICE_TO_HOST, sync = .true.)

  end subroutine gather_interface_dofs_

  !> Set up the global interpolator.
  subroutine setup_interpolator_(this)
    class(overset_interface_t), intent(inout) :: this

    call this%interface_interpolator%init(this%dof, &
         comm = NEKO_GLOBAL_COMM, &
         tol = this%interpolation_settings%tolerance, &
         pad = this%interpolation_settings%padding, &
         mask = this%domain_element_mask)

    call this%interface_interpolator%find_points(this%x_interface_dof%x, &
         this%y_interface_dof%x, this%z_interface_dof%x, &
         this%x_interface_dof%size())

  end subroutine setup_interpolator_

end module overset_interface
