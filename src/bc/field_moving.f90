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

!> Defines moving boundary condition (ALE) - CPU Only (by now)
module field_moving
  use num_types, only : rp
  use coefs, only : coef_t
  use math, only: masked_copy_0
  use device_math, only: device_masked_copy_0
  use field, only : field_t
  use bc, only : bc_t
  use json_module, only : json_file
  use time_state, only : time_state_t
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr
  implicit none
  private

  !> Moving Boundary Condition Class
  type, public, extends(bc_t) :: field_moving_t
     type(field_t), pointer :: wx => null()
     type(field_t), pointer :: wy => null()
     type(field_t), pointer :: wz => null()
  contains
     procedure, pass(this) :: init => field_moving_init
     procedure, pass(this) :: free => field_moving_free
     procedure, pass(this) :: finalize => field_moving_finalize
     procedure, pass(this) :: apply_vector => field_moving_apply_vector
     procedure, pass(this) :: apply_vector_dev => field_moving_apply_vector_dev
   
     ! Todo: dummy for now
     procedure, pass(this) :: apply_scalar => field_moving_apply_scalar
     procedure, pass(this) :: apply_scalar_dev => field_moving_apply_scalar_dev
  end type field_moving_t

contains

  !> Constructor
  subroutine field_moving_init(this, coef, json)
    class(field_moving_t), intent(inout), target :: this
    type(coef_t), target, intent(in) :: coef
    type(json_file), intent(inout) :: json
    call this%init_base(coef)
  end subroutine field_moving_init

  !> Destructor
  subroutine field_moving_free(this)
    class(field_moving_t), target, intent(inout) :: this
    call this%free_base()
    nullify(this%wx)
    nullify(this%wy)
    nullify(this%wz)
  end subroutine field_moving_free

  !> Finalize
  subroutine field_moving_finalize(this, only_facets)
    class(field_moving_t), target, intent(inout) :: this
    logical, optional, intent(in) :: only_facets
    logical :: only_facets_

    if (present(only_facets)) then
       only_facets_ = only_facets
    else
       only_facets_ = .false.
    end if
    call this%finalize_base(only_facets_)
  end subroutine field_moving_finalize

  subroutine field_moving_apply_vector(this, x, y, z, n, time, strong)
    class(field_moving_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x, y, z
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    logical :: strong_

    strong_ = .true.
    if (present(strong)) strong_ = strong

    if (strong_) then
       if (.not. associated(this%wx) .or. .not. associated(this%wy) .or. &
       .not. associated(this%wz)) then
          call neko_error("field_moving error: mesh velocity pointers (wx/wy/wz) not linked!")
       end if

       call masked_copy_0(x, this%wx%x, this%msk, n, this%msk(0))
       call masked_copy_0(y, this%wy%x, this%msk, n, this%msk(0))
       call masked_copy_0(z, this%wz%x, this%msk, n, this%msk(0))
    end if

  end subroutine field_moving_apply_vector

  subroutine field_moving_apply_vector_dev(this, x_d, y_d, z_d, time, strong, strm)
    class(field_moving_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d, y_d, z_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
    logical :: strong_

    strong_ = .true.
    if (present(strong)) strong_ = strong

    if (strong_) then
       if (.not. associated(this%wx) .or. .not. associated(this%wy) .or. &
           .not. associated(this%wz)) then
             call neko_error("field_moving error: mesh velocity pointers not linked!")
       end if

       if (this%msk(0) .gt. 0) then
          call device_masked_copy_0(x_d, this%wx%x_d, this%msk_d, &
                                    this%wx%dof%size(), this%msk(0), strm)
          call device_masked_copy_0(y_d, this%wy%x_d, this%msk_d, &
                                    this%wy%dof%size(), this%msk(0), strm)
          call device_masked_copy_0(z_d, this%wz%x_d, this%msk_d, &
                                    this%wz%dof%size(), this%msk(0), strm)
       end if
    end if
  end subroutine field_moving_apply_vector_dev

  !-----------------------------------------------------------------------------
  ! Todo
  !-----------------------------------------------------------------------------
  
  subroutine field_moving_apply_scalar(this, x, n, time, strong)
    class(field_moving_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(inout), dimension(n) :: x
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    call neko_error("field_moving cannot apply scalar BCs")
  end subroutine field_moving_apply_scalar

  subroutine field_moving_apply_scalar_dev(this, x_d, time, strong, strm)
    class(field_moving_t), intent(inout), target :: this
    type(c_ptr), intent(inout) :: x_d
    type(time_state_t), intent(in), optional :: time
    logical, intent(in), optional :: strong
    type(c_ptr), intent(inout) :: strm
    call neko_error("field_moving cannot apply scalar BCs")
  end subroutine field_moving_apply_scalar_dev

end module field_moving