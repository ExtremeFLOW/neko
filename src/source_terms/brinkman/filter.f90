! Copyright (c) 2023, The Neko Authors
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
!
!> Filter to be applied to a scalar field

! (sorry about the naming convention of filter vs filters, I want the
! name "filters" for this eventually)
module filter
  use num_types, only : rp
  use json_module, only : json_file
  use coefs, only : coef_t
  use json_utils, only : json_get_or_default, json_get
  use field, only: field_t
  implicit none
  private

  !> Base abstract class for filter.
  type, abstract, public :: filter_t
  	  !> Coefficients for the SEM.
     type(coef_t), pointer :: coef => null() 
 
   contains
     !> Constructor for the filter_t class.
     procedure, pass(this) :: init_base => filter_init_base
     !> Destructor for the filter_t (base) class.
     procedure, pass(this) :: free_base => filter_free_base
     !> The common constructor using a JSON dictionary.
     procedure(filter_init), pass(this), deferred :: init
     !> Destructor.
     procedure(filter_free), pass(this), deferred :: free
     !> The main function to be executed during the run.
     procedure(filter_apply), pass(this), deferred :: apply
     !> MAYBE in the future we need an adjoint!
     ! procedure(filter_apply_adjoint), pass(this), deferred :: apply_adjoint
  end type filter_t




  abstract interface
     !> The common constructor using a JSON dictionary.
     !! @param json The JSON with properties.
     !! @param case The case_t object.
     subroutine filter_init(this, json, coef)
       import filter_t, json_file, coef_t
       class(filter_t), intent(inout) :: this
       type(json_file), intent(inout) :: json
       type(coef_t), intent(inout) :: coef
     end subroutine filter_init
  end interface

  abstract interface
     !> Destructor.
     subroutine filter_free(this)
       import filter_t
       class(filter_t), intent(inout) :: this
     end subroutine filter_free
  end interface

  abstract interface
     !> The application of the filter.
     !! @param F_out The output field
     !! @param F_in The input field
     subroutine filter_apply(this, F_out, F_in)
       import filter_t, field_t
       class(filter_t), intent(inout) :: this
       type(field_t), intent(in) ::  F_in
    	 type(field_t), intent(inout) ::  F_out
     end subroutine filter_apply
  end interface

contains
  !> Constructor for the `filter_t` (base) class.
  subroutine filter_init_base(this, json, coef)
    class(filter_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(coef_t), intent(inout), target :: coef
    character(len=:), allocatable :: compute_control, output_control
    real(kind=rp) :: compute_value, output_value
    integer :: order


    this%coef => coef
!    call json_get_or_default(json, "compute_control", compute_control, &
!                             "tsteps")
!    call json_get_or_default(json, "compute_value", compute_value, 1.0_rp)
!
!    ! We default to output whenever we execute
!    call json_get_or_default(json, "output_control", output_control, &
!                             compute_control)
!    call json_get_or_default(json, "output_value", output_value, &
!                             compute_value)
!
!
!    if (output_control == "global") then
!       call json_get(this%case%params, 'case.fluid.output_control', &
!                     output_control)
!       call json_get(this%case%params, 'case.fluid.output_value', &
!                     output_value)
!    end if
!
!    call json_get(json, "order", order)
!    this%order = order
!
!    call this%compute_controller%init(case%end_time, compute_control, &
!                                        compute_value)
!    call this%output_controller%init(case%end_time, output_control, &
!                                     output_value)

  end subroutine filter_init_base

  !> Destructor for the `filter_t` (base) class.
  subroutine filter_free_base(this)
    class(filter_t), intent(inout) :: this

    nullify(this%coef)
  end subroutine filter_free_base



end module filter
