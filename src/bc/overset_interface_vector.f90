! Copyright (c) 2020-2026, The Neko Authors
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
!> Defines overset interface vector boundary conditions
module overset_interface_vector
  use comm, only : NEKO_GLOBAL_COMM
  use neko_config, only: NEKO_BCKND_DEVICE
  use registry, only : neko_registry
  use num_types, only : rp
  use coefs, only : coef_t
  use dirichlet, only : dirichlet_t
  use global_interpolation, only : global_interpolation_t
  use mask, only : mask_t
  use dofmap, only : dofmap_t
  use bc, only : bc_t
  use bc_list, only : bc_list_t
  use utils, only : split_string
  use field, only : field_t
  use field_list, only : field_list_t
  use math, only : masked_copy_0
  use device_math, only : device_masked_copy_0, device_copy
  use dofmap, only : dofmap_t
  use vector, only : vector_t
  use vector_math, only : vector_masked_gather_copy, vector_masked_scatter_copy
  use math, only : copy
  use device, only : DEVICE_TO_HOST, HOST_TO_DEVICE
  use vector_math, only : vector_copy
  use field_dirichlet, only : field_dirichlet_t, field_dirichlet_update
  use utils, only : neko_error, nonlinear_index, linear_index
  use stack, only: stack_i4_t
  use json_module, only : json_file
  use json_utils, only : json_get_subdict_or_empty
  use field_list, only : field_list_t
  use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t
  use time_state, only : time_state_t
  implicit none
  private

  !> Extension of the user defined dirichlet condition `field_dirichlet`
  ! for the application on a vector field.
  type, public, extends(bc_t) :: overset_interface_vector_t
     ! The bc for the first compoent.
     type(field_dirichlet_t) :: bc_u
     ! The bc for the second compoent.
     type(field_dirichlet_t) :: bc_v
     ! The bc for the third compoent.
     type(field_dirichlet_t) :: bc_w
     !> A field list to store the bcs for passing to various subroutines.
     type(field_list_t) :: field_list
     !> Interpolator
     type(global_interpolation_t) :: interface_interpolator
     !> Mask for the overset interface points.
     type(mask_t) :: interface_dof_mask
     type(mask_t) :: domain_element_mask
     !> Vectors holding the dof coordinates for cases where mesh moves
     type(vector_t) :: x_dof, y_dof, z_dof
     type(vector_t) :: x_interface_dof, y_interface_dof, z_interface_dof
     type(vector_t) :: u_interface, v_interface, w_interface
     !> Interpolation settings.
     type(json_file) :: interpolation_settings
   contains
     !> Constructor.
     procedure, pass(this) :: init => overset_interface_vector_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          overset_interface_vector_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => overset_interface_vector_free
     !> Finalize.
     procedure, pass(this) :: finalize => overset_interface_vector_finalize
     !> Apply scalar by performing a masked copy.
     procedure, pass(this) :: apply_scalar => &
          overset_interface_vector_apply_scalar
     !> (No-op) Apply vector.
     procedure, pass(this) :: apply_vector => &
          overset_interface_vector_apply_vector
     !> (No-op) Apply vector (device).
     procedure, pass(this) :: apply_vector_dev => &
          overset_interface_vector_apply_vector_dev
     !> Apply scalar (device).
     procedure, pass(this) :: apply_scalar_dev => &
          overset_interface_vector_apply_scalar_dev
     procedure, pass(this) :: update => overset_interface_update
  end type overset_interface_vector_t

contains

  !> Constructor
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine overset_interface_vector_init(this, coef, json)
    class(overset_interface_vector_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) ::json

    call this%init_from_components(coef)

    !> Store the interpolation settings
    call json_get_subdict_or_empty(json, "interpolation", &
             this%interpolation_settings)

  end subroutine overset_interface_vector_init

  !> Constructor from components
  !! @param[in] coef The SEM coefficients.
  subroutine overset_interface_vector_init_from_components(this, coef)
    class(overset_interface_vector_t), intent(inout), target :: this
    type(coef_t), intent(in) :: coef

    !> This initializes coef, dof, msh, and Xh pointers
    call this%init_base(coef)

    call this%bc_u%init_from_components(coef, "u")
    call this%bc_v%init_from_components(coef, "v")
    call this%bc_w%init_from_components(coef, "w")

    call this%field_list%init(3)
    call this%field_list%assign_to_field(1, this%bc_u%field_bc)
    call this%field_list%assign_to_field(2, this%bc_v%field_bc)
    call this%field_list%assign_to_field(3, this%bc_w%field_bc)

    !> init coord vectors
    call this%x_dof%init(this%dof%size(), 'x')
    call this%y_dof%init(this%dof%size(), 'y')
    call this%z_dof%init(this%dof%size(), 'z')

    !> Copy the dof coordinates to the internal vectors
    !> keep this to later check if the coordinates have changed.
    if (NEKO_BCKND_DEVICE .eq. 1) then
       ! copy
       call device_copy(this%x_dof%x_d, this%dof%x_d, this%dof%size())
       call device_copy(this%y_dof%x_d, this%dof%y_d, this%dof%size())
       call device_copy(this%z_dof%x_d, this%dof%z_d, this%dof%size())
       ! synchronize
       call this%x_dof%copy_from(DEVICE_TO_HOST, sync = .false.)
       call this%y_dof%copy_from(DEVICE_TO_HOST, sync = .false.)
       call this%z_dof%copy_from(DEVICE_TO_HOST, sync = .true.)
    else
       call copy(this%x_dof%x, this%dof%x, this%dof%size())
       call copy(this%y_dof%x, this%dof%y, this%dof%size())
       call copy(this%z_dof%x, this%dof%z, this%dof%size())
    end if

  end subroutine overset_interface_vector_init_from_components

  !> Destructor. Currently unused as is, all field_dirichlet attributes
  !! are freed in `fluid_scheme_incompressible::free`.
  subroutine overset_interface_vector_free(this)
    class(overset_interface_vector_t), target, intent(inout) :: this

    call this%bc_u%free()
    call this%bc_v%free()
    call this%bc_w%free()

    call this%field_list%free()

    call this%x_dof%free()
    call this%y_dof%free()
    call this%z_dof%free()

    call this%x_interface_dof%free()
    call this%y_interface_dof%free()
    call this%z_interface_dof%free()
    call this%u_interface%free()
    call this%v_interface%free()
    call this%w_interface%free()
    call this%interface_interpolator%free()
    call this%interface_dof_mask%free()
    call this%domain_element_mask%free()
    call this%free_base()

    !if (associated(this%update_)) then
    !   nullify(this%update_)
    !end if
  end subroutine overset_interface_vector_free

  !> No-op apply scalar.
  !! @param x Field onto which to copy the values (e.g. u,v,w,p or s).
  !! @param n Size of the array `x`.
  !! @param t Current time state.
  subroutine overset_interface_vector_apply_scalar(this, x, n, time, strong)
    class(overset_interface_vector_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong

    call neko_error("overset_interface_vector cannot apply scalar BCs.&
    & Use overset_interface_vector::apply_vector instead!")

  end subroutine overset_interface_vector_apply_scalar

  !> No-op apply scalar (device).
  !! @param x_d Device pointer to the field onto which to copy the values.
  !! @param time The current time state.
  subroutine overset_interface_vector_apply_scalar_dev(this, x_d, time, &
       strong, strm)
    class(overset_interface_vector_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm

    call neko_error("overset_interface_vector cannot apply scalar BCs.&
    & Use overset_interface_vector::apply_vector instead!")

  end subroutine overset_interface_vector_apply_scalar_dev

  !> Apply the boundary condition to a vector field.
  !! @param x x-component of the field onto which to apply the values.
  !! @param y y-component of the field onto which to apply the values.
  !! @param z z-component of the field onto which to apply the values.
  !! @param n Size of the `x`, `y` and `z` arrays.
  !! @param time The current time state.
  subroutine overset_interface_vector_apply_vector(this, x, y, z, n, time, &
       strong)
    class(overset_interface_vector_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (strong_) then

       ! We can send any of the 3 bcs we have as argument, since they are all
       ! the same boundary.
       if (.not. this%updated) then
          call this%update(time)
          this%updated = .true.
       end if

       !> TODO: Change the mask to mask_t
       call masked_copy_0(x, this%bc_u%field_bc%x, this%msk, n, this%msk(0))
       call masked_copy_0(y, this%bc_v%field_bc%x, this%msk, n, this%msk(0))
       call masked_copy_0(z, this%bc_w%field_bc%x, this%msk, n, this%msk(0))
    end if

  end subroutine overset_interface_vector_apply_vector

  !> Apply the boundary condition to a vector field on the device.
  !! @param x x-component of the field onto which to apply the values.
  !! @param y y-component of the field onto which to apply the values.
  !! @param z z-component of the field onto which to apply the values.
  !! @param time The current time state.
  !! @param strm Device stream
  subroutine overset_interface_vector_apply_vector_dev(this, x_d, y_d, z_d, &
       time, strong, strm)
    class(overset_interface_vector_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(c_ptr), intent(inout) :: y_d
    type(c_ptr), intent(inout) :: z_d
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
          call device_masked_copy_0(x_d, this%bc_u%field_bc%x_d, &
               this%bc_u%msk_d, this%bc_u%dof%size(), this%msk(0), strm) ! adperez: change the masks used here
          call device_masked_copy_0(y_d, this%bc_v%field_bc%x_d, &
               this%bc_v%msk_d, this%bc_v%dof%size(), this%msk(0), strm)
          call device_masked_copy_0(z_d, this%bc_w%field_bc%x_d, &
               this%bc_w%msk_d, this%bc_w%dof%size(), this%msk(0), strm)
       end if
    end if

  end subroutine overset_interface_vector_apply_vector_dev

  !> Finalize by building the mask arrays and propagating to underlying bcs.
  subroutine overset_interface_vector_finalize(this, only_facets)
    class(overset_interface_vector_t), target, intent(inout) :: this
    logical, optional, intent(in) :: only_facets
    logical :: only_facets_

    if (present(only_facets)) then
       only_facets_ = only_facets
    else
       only_facets_ = .false.
    end if

    !> From field_dirichlet_vector_t
    call this%finalize_base(only_facets_)

    call this%bc_u%mark_facets(this%marked_facet)
    call this%bc_v%mark_facets(this%marked_facet)
    call this%bc_w%mark_facets(this%marked_facet)

    call this%bc_u%finalize(only_facets_)
    call this%bc_v%finalize(only_facets_)
    call this%bc_w%finalize(only_facets_)

    !> Build heper masks
    call build_masks_(this)

    !> Gather the interface boundary points
    call this%x_interface_dof%init(this%interface_dof_mask%size(), 'x_interface')
    call this%y_interface_dof%init(this%interface_dof_mask%size(), 'y_interface')
    call this%z_interface_dof%init(this%interface_dof_mask%size(), 'z_interface')
    call gather_interface_dofs_(this)

    !> Initialize the interpolator and find the points
    call setup_interpolator_(this)

    !> Keep a vector list that holds the values of interface fields
    call this%u_interface%init(this%interface_dof_mask%size(), 'u_interface')
    call this%v_interface%init(this%interface_dof_mask%size(), 'v_interface')
    call this%w_interface%init(this%interface_dof_mask%size(), 'w_interface')


  end subroutine overset_interface_vector_finalize

  !> Update the values at the interfaces
  subroutine overset_interface_update(this, time)
    class(overset_interface_vector_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, v, w


    !> Update in sub-step 1 should be an extrapolation of the boundary values
    ! not implemented for now
    ! if (substep .eq. 1) then
    !   call this%extrapolate()
    ! end if

    !> At some point check if the coordintes have changed. If so, find points again
    ! not implemented for now

    !> For more substep than 1, then we just interpolate
    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")

    !> Interpolate the values
    call this%interface_interpolator%evaluate_masked(this%u_interface%x, u%x, this%domain_element_mask, .false.)
    call this%interface_interpolator%evaluate_masked(this%v_interface%x, v%x, this%domain_element_mask, .false.)
    call this%interface_interpolator%evaluate_masked(this%w_interface%x, w%x, this%domain_element_mask, .false.)

    !> Scatter them to the bc fields
    call vector_masked_scatter_copy(this%bc_u%field_bc%x(:,1,1,1), this%u_interface, &
         this%interface_dof_mask, this%bc_u%dof%size())
    call vector_masked_scatter_copy(this%bc_v%field_bc%x(:,1,1,1), this%v_interface, &
         this%interface_dof_mask, this%bc_v%dof%size())
    call vector_masked_scatter_copy(this%bc_w%field_bc%x(:,1,1,1), this%w_interface, &
         this%interface_dof_mask, this%bc_w%dof%size())


  end subroutine overset_interface_update

  !===================
  ! Helper subroutines
  !===================

  !> Build masks
  subroutine build_masks_(this)
    class(overset_interface_vector_t), intent(inout) :: this
    type(mask_t) :: temp_mask
    logical, allocatable :: found(:)
    integer :: i, j, k, e, new_size, nelems
    integer :: lx, ly, lz
    integer :: nonlinear_idx(4), linear_idx
    type(stack_i4_t) :: stack

    !> Create a new mask_t based on the 0-mask for the overset boundary points
    call this%interface_dof_mask%init(this%msk(1:this%msk(0)), this%msk(0))

    !> Create a mask that marks the full elements
    !! containing overset boundary points
    lx = this%Xh%lx
    ly = this%Xh%ly
    lz = this%Xh%lz
    allocate(found(this%msh%nelv))
    found = .false.
    !! Find sem elements that contain the boundary points
    do i = 1, this%msk(0)
       linear_idx = this%msk(i)
       nonlinear_idx = nonlinear_index(linear_idx, lx, ly, lz)
       found(nonlinear_idx(4)) = .true.
    end do
    !! fill the stack containing the gll indices
    nelems = 0
    call stack%init()
    do e = 1, this%msh%nelv
       if (found(e)) then
          nelems = nelems + 1
          do k = 1, this%Xh%lz
             do j = 1, this%Xh%ly
                do i = 1, this%Xh%lx
                   linear_idx = linear_index(i, j, k, e, lx, ly, lz)
                   call stack%push(linear_idx)
                end do
             end do
          end do
       end if
    end do
    deallocate(found)
    call temp_mask%init(stack%array(), stack%size())
    call stack%free()

    !> Create a domain mask that exclude the boundary elements.
    call this%domain_element_mask%invert_mask(temp_mask, this%dof%size())

    call temp_mask%free()

  end subroutine build_masks_

  !> Gather interface dofs
  subroutine gather_interface_dofs_(this)
    class(overset_interface_vector_t), intent(inout) :: this

    !> Gather the coordinates to the vectors
    call vector_masked_gather_copy(this%x_interface_dof, &
         this%dof%x(:,1,1,1), &
         this%interface_dof_mask, &
         this%dof%size())
    call vector_masked_gather_copy(this%y_interface_dof, &
         this%dof%y(:,1,1,1), &
         this%interface_dof_mask, &
         this%dof%size())
    call vector_masked_gather_copy(this%z_interface_dof, &
         this%dof%z(:,1,1,1), &
         this%interface_dof_mask, &
         this%dof%size())
    !> synchronize if on device
    call this%x_interface_dof%copy_from(DEVICE_TO_HOST, sync = .false.)
    call this%y_interface_dof%copy_from(DEVICE_TO_HOST, sync = .false.)
    call this%z_interface_dof%copy_from(DEVICE_TO_HOST, sync = .true.)

  end subroutine gather_interface_dofs_


  !> Set up the interpolator
  subroutine setup_interpolator_(this)
    class(overset_interface_vector_t), intent(inout) :: this

    !> Initialize the interpolator with the mask and the dof coords
    call this%interface_interpolator%init(this%dof, this%interpolation_settings, &
         NEKO_GLOBAL_COMM, &
         mask=this%domain_element_mask)

    !> Find the interface points on the global domain
    call this%interface_interpolator%find_points(this%x_interface_dof%x, &
         this%y_interface_dof%x, &
         this%z_interface_dof%x, &
         this%x_interface_dof%size())

  end subroutine setup_interpolator_

end module overset_interface_vector
