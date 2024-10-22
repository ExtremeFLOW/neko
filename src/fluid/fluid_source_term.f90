
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
!> Implements the `fluid_source_term_t` type.
module fluid_source_term
  use fluid_user_source_term, only: fluid_user_source_term_t
  use source_term, only: source_term_t
  use source_term_handler, only: source_term_handler_t
  use field, only: field_t
  use field_list, only: field_list_t
  use coefs, only: coef_t
  use user_intf, only: user_t
  implicit none
  private

  !> Wrapper contaning and executing the fluid source terms.
  !! @details
  !! Exists mainly to keep the `fluid_scheme_t` type smaller and also as
  !! placeholder for future optimizations.
  type, public, extends(source_term_handler_t) :: fluid_source_term_t

   contains
     !> Constructor.
     procedure, pass(this) :: init => fluid_source_term_init
     !> Initialize the user source term.
     procedure, nopass :: init_user_source => fluid_init_user_source

  end type fluid_source_term_t

contains

  !> Constructor.
  subroutine fluid_source_term_init(this, f_x, f_y, f_z, coef, user)
    class(fluid_source_term_t), intent(inout) :: this
    type(field_t), pointer, intent(in) :: f_x, f_y, f_z
    type(coef_t), target, intent(inout) :: coef
    type(user_t), target, intent(in) :: user

    type(field_list_t) :: rhs_fields

    ! We package the fields for the source term to operate on in a field list.
    call rhs_fields%init(3)
    call rhs_fields%assign(1, f_x)
    call rhs_fields%assign(2, f_y)
    call rhs_fields%assign(3, f_z)

    call this%init_base(rhs_fields, coef, user)
  end subroutine fluid_source_term_init

  !> Initialize the user source term.
  !! @param source_term The allocatable source term to be initialized to a user.
  !! @param rhs_fields The field list with the 3 right-hand-side components.
  !! @param coef The SEM coefs.
  !! @param type The type of the user source term, "user_vector" or
  !! "user_poinwise".
  !! @param user The user type containing the user source term routines.
  subroutine fluid_init_user_source(source_term, rhs_fields, coef, type, user)
    class(source_term_t), allocatable, intent(inout) :: source_term
    type(field_list_t) :: rhs_fields
    type(coef_t), intent(inout) :: coef
    character(len=*) :: type
    type(user_t), intent(in) :: user

    allocate(fluid_user_source_term_t::source_term)

    select type (source_term)
      type is (fluid_user_source_term_t)
       call source_term%init_from_components(rhs_fields, coef, type, &
            user%fluid_user_f_vector, &
            user%fluid_user_f)
    end select
  end subroutine fluid_init_user_source

end module fluid_source_term
