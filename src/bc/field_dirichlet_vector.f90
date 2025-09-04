! Copyright (c) 2020-2025, The Neko Authors
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
!> Defines inflow dirichlet conditions
module field_dirichlet_vector
  use num_types, only: rp
  use coefs, only: coef_t
  use dirichlet, only: dirichlet_t
  use bc, only: bc_t
  use bc_list, only : bc_list_t
  use utils, only: split_string
  use field, only : field_t
  use field_list, only : field_list_t
  use math, only: masked_copy_0
  use device_math, only: device_masked_copy_0
  use dofmap, only : dofmap_t
  use field_dirichlet, only: field_dirichlet_t, field_dirichlet_update
  use utils, only: neko_error
  use json_module, only : json_file
  use field_list, only : field_list_t
  use, intrinsic :: iso_c_binding, only : c_ptr, c_size_t
  use time_state, only : time_state_t
  implicit none
  private

  !> Extension of the user defined dirichlet condition `field_dirichlet`
  ! for the application on a vector field.
  type, public, extends(bc_t) :: field_dirichlet_vector_t
     ! The bc for the first compoent.
     type(field_dirichlet_t) :: bc_u
     ! The bc for the second compoent.
     type(field_dirichlet_t) :: bc_v
     ! The bc for the third compoent.
     type(field_dirichlet_t) :: bc_w
     !> A field list to store the bcs for passing to various subroutines.
     type(field_list_t) :: field_list
     !> Function pointer to the user routine performing the update of the values
     !! of the boundary fields.
     procedure(field_dirichlet_update), nopass, pointer :: update => null()
   contains
     !> Constructor.
     procedure, pass(this) :: init => field_dirichlet_vector_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          field_dirichlet_vector_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => field_dirichlet_vector_free
     !> Finalize.
     procedure, pass(this) :: finalize => field_dirichlet_vector_finalize
     !> Apply scalar by performing a masked copy.
     procedure, pass(this) :: apply_scalar => &
          field_dirichlet_vector_apply_scalar
     !> (No-op) Apply vector.
     procedure, pass(this) :: apply_vector => &
          field_dirichlet_vector_apply_vector
     !> (No-op) Apply vector (device).
     procedure, pass(this) :: apply_vector_dev => &
          field_dirichlet_vector_apply_vector_dev
     !> Apply scalar (device).
     procedure, pass(this) :: apply_scalar_dev => &
          field_dirichlet_vector_apply_scalar_dev
  end type field_dirichlet_vector_t

contains

  !> Constructor
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine field_dirichlet_vector_init(this, coef, json)
    class(field_dirichlet_vector_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) ::json

    call this%init_from_components(coef)

  end subroutine field_dirichlet_vector_init

  !> Constructor from components
  !! @param[in] coef The SEM coefficients.
  subroutine field_dirichlet_vector_init_from_components(this, coef)
    class(field_dirichlet_vector_t), intent(inout), target :: this
    type(coef_t), intent(in) :: coef

    call this%init_base(coef)

    call this%bc_u%init_from_components(coef, "u")
    call this%bc_v%init_from_components(coef, "v")
    call this%bc_w%init_from_components(coef, "w")

    ! TODO set to u v w values

    call this%field_list%init(3)
    call this%field_list%assign_to_field(1, this%bc_u%field_bc)
    call this%field_list%assign_to_field(2, this%bc_v%field_bc)
    call this%field_list%assign_to_field(3, this%bc_w%field_bc)

  end subroutine field_dirichlet_vector_init_from_components

  !> Destructor. Currently unused as is, all field_dirichlet attributes
  !! are freed in `fluid_scheme_incompressible::free`.
  subroutine field_dirichlet_vector_free(this)
    class(field_dirichlet_vector_t), target, intent(inout) :: this

    call this%bc_u%free()
    call this%bc_v%free()
    call this%bc_w%free()

    call this%field_list%free()

    if (associated(this%update)) then
       nullify(this%update)
    end if
  end subroutine field_dirichlet_vector_free

  !> No-op apply scalar.
  !! @param x Field onto which to copy the values (e.g. u,v,w,p or s).
  !! @param n Size of the array `x`.
  !! @param t Current time state.
  subroutine field_dirichlet_vector_apply_scalar(this, x, n, time, strong)
    class(field_dirichlet_vector_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong

    call neko_error("field_dirichlet_vector cannot apply scalar BCs.&
    & Use field_dirichlet instead!")

  end subroutine field_dirichlet_vector_apply_scalar

  !> No-op apply scalar (device).
  !! @param x_d Device pointer to the field onto which to copy the values.
  !! @param time The current time state.
  subroutine field_dirichlet_vector_apply_scalar_dev(this, x_d, time, &
       strong, strm)
    class(field_dirichlet_vector_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm

    call neko_error("field_dirichlet_vector cannot apply scalar BCs.&
    & Use field_dirichlet instead!")

  end subroutine field_dirichlet_vector_apply_scalar_dev

  !> Apply the boundary condition to a vector field.
  !! @param x x-component of the field onto which to apply the values.
  !! @param y y-component of the field onto which to apply the values.
  !! @param z z-component of the field onto which to apply the values.
  !! @param n Size of the `x`, `y` and `z` arrays.
  !! @param time The current time state.
  subroutine field_dirichlet_vector_apply_vector(this, x, y, z, n, time, &
       strong)
    class(field_dirichlet_vector_t), intent(inout) :: this
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
          call this%update(this%field_list, this%bc_u, time)
          this%updated = .true.
       end if

       call masked_copy_0(x, this%bc_u%field_bc%x, this%msk, n, this%msk(0))
       call masked_copy_0(y, this%bc_v%field_bc%x, this%msk, n, this%msk(0))
       call masked_copy_0(z, this%bc_w%field_bc%x, this%msk, n, this%msk(0))
    end if

  end subroutine field_dirichlet_vector_apply_vector

  !> Apply the boundary condition to a vector field on the device.
  !! @param x x-component of the field onto which to apply the values.
  !! @param y y-component of the field onto which to apply the values.
  !! @param z z-component of the field onto which to apply the values.
  !! @param time The current time state.
  !! @param strm Device stream
  subroutine field_dirichlet_vector_apply_vector_dev(this, x_d, y_d, z_d, &
       time, strong, strm)
    class(field_dirichlet_vector_t), intent(inout), target :: this
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
          call this%update(this%field_list, this%bc_u, time)
          this%updated = .true.
       end if

       if (this%msk(0) .gt. 0) then
          call device_masked_copy_0(x_d, this%bc_u%field_bc%x_d, this%bc_u%msk_d,&
               this%bc_u%dof%size(), this%msk(0), strm)
          call device_masked_copy_0(y_d, this%bc_v%field_bc%x_d, this%bc_v%msk_d,&
               this%bc_v%dof%size(), this%msk(0), strm)
          call device_masked_copy_0(z_d, this%bc_w%field_bc%x_d, this%bc_w%msk_d,&
               this%bc_w%dof%size(), this%msk(0), strm)
       end if
    end if

  end subroutine field_dirichlet_vector_apply_vector_dev

  !> Finalize by building the mask arrays and propagating to underlying bcs.
  subroutine field_dirichlet_vector_finalize(this, only_facets)
    class(field_dirichlet_vector_t), target, intent(inout) :: this
    logical, optional, intent(in) :: only_facets
    logical :: only_facets_

    if (present(only_facets)) then
       only_facets_ = only_facets
    else
       only_facets_ = .false.
    end if

    call this%finalize_base(only_facets_)

    call this%bc_u%mark_facets(this%marked_facet)
    call this%bc_v%mark_facets(this%marked_facet)
    call this%bc_w%mark_facets(this%marked_facet)

    call this%bc_u%finalize(only_facets_)
    call this%bc_v%finalize(only_facets_)
    call this%bc_w%finalize(only_facets_)

  end subroutine field_dirichlet_vector_finalize

end module field_dirichlet_vector
