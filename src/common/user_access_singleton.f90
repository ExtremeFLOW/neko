! Copyright (c) 2025, The Neko Authors
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

!> User access singleton
!! Defines a singleton object available in the user file. Intended to allow
!! unrestricted access to the entire simulation case and enable all sorts of
!! hacking, while keeping neko proper clear.
!! The object is only initialized by makeneko, so trying to use it
!! outside user code will result in a segfault---intentionally.
module user_access_singleton
  use case, only : case_t
  implicit none
  private

  !> Helper type to give users global access to the simulation case.
  type, public :: user_access_t
     type(case_t), pointer :: case => null()
   contains
     !> Constructor.
     procedure, pass(this) :: init => user_access_init
     !> Destructor.
     procedure, pass(this) :: free => user_access_free
  end type user_access_t

  !> The singleton object.
  type(user_access_t), target, public :: neko_user_access

contains

  !> Constructor.
  subroutine user_access_init(this, case)
    class(user_access_t), intent(inout) :: this
    type(case_t), target :: case

    this%case => case
  end subroutine user_access_init

  !> Destructor.
  subroutine user_access_free(this)
    class(user_access_t), intent(inout) :: this

    if (associated(this%case)) this%case => null()
  end subroutine user_access_free

end module user_access_singleton
