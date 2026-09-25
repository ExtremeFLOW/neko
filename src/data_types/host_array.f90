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
!> Module containing host-only array type.
module host_array
  use num_types, only : rp
  use math, only : rzero, copy

  implicit none
  private

  !> Host-only temporary array.
  type, public :: host_array_t
     !> Host values.
     real(kind=rp), allocatable :: x(:)
     !> Number of entries.
     integer, private :: n = 0
   contains
     !> Initialize the array.
     procedure, pass(this) :: init => host_array_init
     !> Free the array.
     procedure, pass(this) :: free => host_array_free
     !> Allocate the array
     procedure, private, pass(this) :: allocate => host_array_allocate
     !> Return the number of entries.
     procedure, pass(this) :: size => host_array_size
     !> Check whether storage is allocated.
     procedure, pass(this) :: is_allocated => host_array_is_allocated
     !> Assignment with deep-copy ownership semantics.
     procedure, pass(this) :: assign => host_array_assign

     !> Assignments
     generic :: assignment(=) => assign


  end type host_array_t

contains

  !> Initialize a host array of size `size`.
  !! @param this Host array object.
  !! @param size Number of entries to allocate.
  subroutine host_array_init(this, size)
    class(host_array_t), intent(inout) :: this
    integer, intent(in) :: size

    call this%free()
    call this%allocate(size)
    call rzero(this%x, this%n)

  end subroutine host_array_init

  !> Free a host array.
  !! @param this Host array object.
  subroutine host_array_free(this)
    class(host_array_t), intent(inout) :: this

    this%n = 0
    if (allocated(this%x)) deallocate(this%x)

  end subroutine host_array_free

  !> Allocate a host array (used by init and assignment).
  !! @param this host array object.
  !! @param size Number of entries to allocate.
  subroutine host_array_allocate(this, size)
    class(host_array_t), intent(inout) :: this
    integer, intent(in) :: size

    this%n = size
    allocate(this%x(size))

  end subroutine host_array_allocate

  !> Assignment with deep-copy ownership semantics.
  !! @param this Destination host array object.
  !! @param source Source host array object.
  subroutine host_array_assign(this, source)
    class(host_array_t), intent(inout) :: this
    class(host_array_t), intent(in) :: source

    if (.not. source%is_allocated()) then
       call this%free()
       return
    else if (source%size() .ne. this%size()) then
       call this%free()
       call this%allocate(source%size())
    end if

    call copy(this%x, source%x, this%n)

  end subroutine host_array_assign

  !> Return the number of entries in the host array.
  !! @param this Host array object.
  pure function host_array_size(this) result(n)
    class(host_array_t), intent(in) :: this
    integer :: n
    n = this%n
  end function host_array_size

  !> Check whether the host array is allocated.
  !! @param this Host array object.
  pure function host_array_is_allocated(this) result(is_alloc)
    class(host_array_t), intent(in) :: this
    logical :: is_alloc

    is_alloc = allocated(this%x)

  end function host_array_is_allocated

end module host_array
