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
!> Defines user neumann condition for a scalar field.
module field_neumann
  use num_types, only : rp
  use coefs, only : coef_t
  use bc, only : bc_t
  use field, only : field_t
  use field_list, only : field_list_t
  use vector, only : vector_t
  use utils, only : neko_error, nonlinear_index
  use json_module, only : json_file
  use json_utils, only : json_get
  use math, only : masked_gather_copy_0 
  use device_math, only : device_masked_gather_copy_0
  use device_neumann, only : device_neumann_apply_scalar
  use neko_config, only : NEKO_BCKND_DEVICE
  use, intrinsic :: iso_c_binding, only : c_ptr
  use time_state, only : time_state_t
  implicit none
  private

  !> User defined neumann condition, for which the user can work
  !! with an entire field.
  !! The type stores a separate dummy field `field_bc`, which is passed
  !! to the user routine and can be populated with arbitrary values. The
  !! boundary condition then gathers these values at the bc mask locations
  !! into a compact flux vector and applies the weak neumann contribution.
  type, public, extends(bc_t) :: field_neumann_t
     !> A dummy field which can be manipulated by the user to set flux values.
     type(field_t) :: field_bc
     !> A field list, which just stores `field_bc`, for convenience.
     type(field_list_t) :: field_list
     !> Flux values gathered on boundary-mask ordering.
     type(vector_t) :: flux
     !> Function pointer to the user routine performing the update of the flux.
     procedure(field_neumann_update), nopass, pointer :: update => null()
   contains
     !> Constructor.
     procedure, pass(this) :: init => field_neumann_init
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => field_neumann_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => field_neumann_free
     !> Finalize.
     procedure, pass(this) :: finalize => field_neumann_finalize
     !> Apply scalar by adding weak neumann contribution.
     procedure, pass(this) :: apply_scalar => field_neumann_apply_scalar
     !> (No-op) Apply vector.
     procedure, pass(this) :: apply_vector => field_neumann_apply_vector
     !> (No-op) Apply vector (device).
     procedure, pass(this) :: apply_vector_dev => field_neumann_apply_vector_dev
     !> Apply scalar (device).
     procedure, pass(this) :: apply_scalar_dev => field_neumann_apply_scalar_dev
     !> Gather flux values at masked points.
     procedure, pass(this), private :: gather_flux => field_neumann_gather_flux

  end type field_neumann_t

  !> Abstract interface defining a neumann condition on a list of fields.
  abstract interface
     subroutine field_neumann_update(fields, bc, time)
       import field_list_t, field_neumann_t, time_state_t
       type(field_list_t), intent(inout) :: fields
       type(field_neumann_t), intent(in) :: bc
       type(time_state_t), intent(in) :: time
     end subroutine field_neumann_update
  end interface

  public :: field_neumann_update

contains
  !> Constructor.
  !! @param[in] coef The SEM coefficients.
  !! @param[inout] json The JSON object configuring the boundary condition.
  subroutine field_neumann_init(this, coef, json)
    class(field_neumann_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: field_name

    call json_get(json, "field_name", field_name)
    call this%init_from_components(coef, field_name)

  end subroutine field_neumann_init

  !> Constructor from components.
  !! @param[in] coef The SEM coefficients.
  !! @param[in] field_name Name of the dummy field.
  subroutine field_neumann_init_from_components(this, coef, field_name)
    class(field_neumann_t), intent(inout), target :: this
    type(coef_t), intent(in) :: coef
    character(len=*), intent(in) :: field_name

    call this%init_base(coef)
    this%strong = .false.

    call this%field_bc%init(this%dof, field_name)
    call this%field_list%init(1)
    call this%field_list%assign_to_field(1, this%field_bc)

  end subroutine field_neumann_init_from_components

  !> Destructor.
  subroutine field_neumann_free(this)
    class(field_neumann_t), target, intent(inout) :: this

    call this%free_base
    call this%field_bc%free()
    call this%field_list%free()
    call this%flux%free()

    if (associated(this%update)) then
       this%update => null()
    end if

  end subroutine field_neumann_free

  !> Gather field-defined values into compact boundary flux storage.
  subroutine field_neumann_gather_flux(this)
    class(field_neumann_t), intent(inout) :: this
    integer :: i
    integer :: idx(4)

    if (this%msk(0) .gt. 0) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_masked_gather_copy_0(this%flux%x_d, this%field_bc%x_d, &
               this%msk_d, this%field_bc%dof%size(), this%msk(0))
       else
          call masked_gather_copy_0(this%flux%x, this%field_bc%x, &
               this%msk, this%field_bc%dof%size(), this%msk(0))
       end if
    end if

  end subroutine field_neumann_gather_flux

  !> Apply scalar by adding weak neumann contribution.
  !! @param x Field to which the weak neumann contribution is added.
  !! @param n Size of the array `x`.
  !! @param time The current time state.
  subroutine field_neumann_apply_scalar(this, x, n, time, strong)
    class(field_neumann_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    integer :: i, m, k, facet
    integer :: idx(4)
    logical :: strong_

    if (present(strong)) then
       strong_ = strong
    else
       strong_ = .true.
    end if

    if (.not. strong_) then

       if (.not. this%updated) then
          call this%update(this%field_list, this, time)
          call this%gather_flux()
          this%updated = .true.
       end if

       m = this%msk(0)
       !$omp parallel do private(k, facet, idx)
       do i = 1, m
          k = this%msk(i)
          facet = this%facet(i)
          idx = nonlinear_index(k, this%coef%Xh%lx, this%coef%Xh%lx, &
               this%coef%Xh%lx)
          select case (facet)
          case (1,2)
             x(k) = x(k) + &
                  this%flux%x(i) * &
                  this%coef%area(idx(2), idx(3), facet, idx(4))
          case (3,4)
             x(k) = x(k) + &
                  this%flux%x(i) * &
                  this%coef%area(idx(1), idx(3), facet, idx(4))
          case (5,6)
             x(k) = x(k) + &
                  this%flux%x(i) * &
                  this%coef%area(idx(1), idx(2), facet, idx(4))
          end select
       end do
       !$omp end parallel do
    end if

  end subroutine field_neumann_apply_scalar

  !> Apply scalar (device).
  !! @param x_d Device pointer to the field to update.
  !! @param time The current time state.
  !! @param strm Device stream.
  subroutine field_neumann_apply_scalar_dev(this, x_d, time, strong, strm)
    class(field_neumann_t), intent(inout), target :: this
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

    if (.not. strong_) then
       if (.not. this%updated) then
          call this%update(this%field_list, this, time)
          call this%gather_flux()
          this%updated = .true.
       end if

       if (this%msk(0) .gt. 0) then
          call device_neumann_apply_scalar(this%msk_d, this%facet_d, x_d, &
               this%flux%x_d, this%coef%area_d, this%coef%Xh%lx, &
               size(this%msk), strm)
       end if
    end if

  end subroutine field_neumann_apply_scalar_dev

  !> (No-op) Apply vector.
  subroutine field_neumann_apply_vector(this, x, y, z, n, time, strong)
    class(field_neumann_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    real(kind=rp), intent(inout), dimension(n) :: y
    real(kind=rp), intent(inout), dimension(n) :: z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong

    call neko_error("field_neumann cannot apply vector BCs.")

  end subroutine field_neumann_apply_vector

  !> (No-op) Apply vector (device).
  subroutine field_neumann_apply_vector_dev(this, x_d, y_d, z_d, time, &
       strong, strm)
    class(field_neumann_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(c_ptr), intent(inout) :: y_d
    type(c_ptr), intent(inout) :: z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm

    call neko_error("field_neumann cannot apply vector BCs.")

  end subroutine field_neumann_apply_vector_dev

  !> Finalize.
  subroutine field_neumann_finalize(this, only_facets)
    class(field_neumann_t), target, intent(inout) :: this
    logical, optional, intent(in) :: only_facets

    call this%finalize_base(.true.)
    call this%flux%init(this%msk(0))

  end subroutine field_neumann_finalize

end module field_neumann
