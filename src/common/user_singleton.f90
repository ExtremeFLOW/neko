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

! Defines a singleton object available in the user file. Intended to allow
! unrestricted access to the entire simulation case and enable all sorts of
! hacking, while keeping neko proper clear.
! The neko_user object is only initialized by makeneko, so trying to use it
! outside user code will result in a segfault---intentionally.
module user_singleton
  use case, only : case_t
  use user_intf, only : user_t
  implicit none
  private

  !> Helper type to give users global access to the simulation case.
  type, public :: neko_user_t
     type(user_t), pointer :: user
     type(case_t), pointer :: case
   contains
     !> Constructor.
     procedure, pass(this) :: init => neko_user_init
     !> Destructor.
     procedure, pass(this) :: free => neko_user_free
  end type neko_user_t

  !> The singleton object.
  type(neko_user_t), target, public :: neko_user

contains

  !> Consturctor.
  subroutine neko_user_init(this, user, case)
    class(neko_user_t), intent(inout) :: this
    type(user_t), target :: user
    type(case_t), target :: case

    this%user => user
    this%case => case
  end subroutine neko_user_init

  !> Destructor.
  subroutine neko_user_free(this)
    class(neko_user_t), intent(inout) :: this

    if (associated(this%user)) this%user => null()
    if (associated(this%case)) this%case => null()
  end subroutine neko_user_free

end module user_singleton
