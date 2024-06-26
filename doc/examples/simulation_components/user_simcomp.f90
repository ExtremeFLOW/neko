! Copyright (c) 2024, The Neko Authors
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

! Implements the `user_simcomp_t` type.
module simcomp_example
  use num_types, only : rp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use case, only : case_t
  implicit none
  private

  ! An empty user defined simulation component.
  ! This is a simple example of a user-defined simulation component.
  type, public, extends(simulation_component_t) :: user_simcomp_t

   contains
     ! Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => simcomp_test_init_from_json
     ! Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
       simcomp_test_init_from_attributes
     ! Destructor.
     procedure, pass(this) :: free => simcomp_test_free
     ! Compute the simcomp_test field.
     procedure, pass(this) :: compute_ => simcomp_test_compute
  end type user_simcomp_t

contains

  ! Constructor from json.
  subroutine simcomp_test_init_from_json(this, json, case)
    class(user_simcomp_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case

    call this%init_from_attributes()
    call this%init_base(json, case)

  end subroutine simcomp_test_init_from_json

  ! Actual constructor.
  subroutine simcomp_test_init_from_attributes(this)
    class(user_simcomp_t), intent(inout) :: this
  end subroutine simcomp_test_init_from_attributes

  ! Destructor.
  subroutine simcomp_test_free(this)
    class(user_simcomp_t), intent(inout) :: this
    call this%free_base()
  end subroutine simcomp_test_free

  ! Compute the simcomp_test field.
  subroutine simcomp_test_compute(this, t, tstep)
    class(user_simcomp_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

  end subroutine simcomp_test_compute

end module simcomp_example


! User module for the user defined simulation component
module user
  use neko
  use simcomp_example, only: user_simcomp_t
  implicit none

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%init_user_simcomp => user_simcomp
  end subroutine user_setup

  subroutine user_simcomp(params)
    type(json_file), intent(inout) :: params
    type(user_simcomp_t), allocatable :: my_simcomp
    type(json_file) :: simcomp_settings

    ! Allocate a simulation component
    allocate(my_simcomp)
    simcomp_settings = simulation_component_user_settings("my_comp", params)
    call neko_simcomps%add_user_simcomp(my_simcomp, simcomp_settings)

  end subroutine user_simcomp

end module user
