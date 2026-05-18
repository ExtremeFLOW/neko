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
!> Module containing device only array type.
module device_array
  use num_types, only : rp, sp, dp
  use device, only : device_alloc, device_free
  use device_math, only : device_rzero
  use, intrinsic :: iso_c_binding

  implicit none
  private

  type, public :: device_array_t
     type(c_ptr) :: x_d = C_NULL_PTR
     integer :: n = 0
   contains
     procedure, pass(this) :: init => device_array_init
     procedure, pass(this) :: free => device_array_free
     procedure, pass(this) :: size => device_array_size
     procedure, pass(this) :: is_allocated => device_array_is_allocated
  end type device_array_t

contains

  subroutine device_array_init(this, size)
    class(device_array_t), intent(inout) :: this
    integer, intent(in) :: size
    integer(c_size_t) :: c_size

    call this%free()
    this%n = size

    select case (rp)
    case (sp)
       c_size = size * int(4, c_size_t)
    case (dp)
       c_size = size * int(8, c_size_t)
    case default
       call neko_error('Unknown Fortran type')
    end select

    call device_alloc(this%x_d, c_size)
    call device_rzero(this%x_d, this%n)

  end subroutine device_array_init

  subroutine device_array_free(this)
    class(device_array_t), intent(inout) :: this

    this%n = 0
    if (c_associated(this%x_d)) call device_free(this%x_d)

  end subroutine device_array_free

  pure function device_array_size(this) result(n)
    class(device_array_t), intent(in) :: this
    integer :: n
    n = this%n
  end function device_array_size

  pure function device_array_is_allocated(this) result(is_alloc)
    class(device_array_t), intent(in) :: this
    logical :: is_alloc

    is_alloc = c_associated(this%x_d)

  end function device_array_is_allocated

end module device_array
