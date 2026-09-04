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
!> Defines an array
module array
  use, intrinsic :: iso_c_binding, only : c_ptr, C_NULL_PTR
  use num_types, only : rp
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_map, device_unmap, &
       device_memcpy, device_sync, HOST_TO_DEVICE
  use math, only : cfill, copy
  use device_math, only : device_copy, device_cfill
  use utils, only : neko_error, NEKO_VARNAME_LEN
  implicit none
  private

  type, abstract, public :: array_t
     !> Name of the array
     character(len=NEKO_VARNAME_LEN) :: name = ""
     !> Array entries.
     real(kind=rp), allocatable :: data(:)
     !> Device pointer.
     type(c_ptr) :: x_d = C_NULL_PTR
     !> Size of array.
     integer, private :: n = 0
   contains
     !> All array types must implement the free method.
     procedure(array_free), pass(this), deferred :: free

     !> Initialise a array of size `n` and optional name.
     procedure, pass(this) :: init_base => array_init_base
     !> Deallocate a array.
     procedure, pass(this) :: free_base => array_free_base
     !> Copy data between host and device
     procedure, pass(this) :: copy_from => array_copy_from
     !> Returns the number of entries in the array.
     procedure, pass(this) :: size => array_size
     !> Determine if the array has been allocated.
     procedure, pass(this) :: is_allocated => array_is_allocated

     !> Assignments
     generic :: assignment(=) => array_assign_array, &
          array_assign_scalar, array_assign_real

     !> Assignment \f$ v = array type \f$
     procedure, pass(this) :: array_assign_array
     !> Assignment \f$ v = scalar real \f$.
     procedure, pass(this) :: array_assign_scalar
     !> Assignment \f$ v = array of reals \f$.
     procedure, pass(this) :: array_assign_real

     ! Private interfaces
     procedure, pass(this), private :: alloc => array_allocate

  end type array_t

  ! Define the abstract interfaces for the deferred procedures
  abstract interface
     subroutine array_free(this)
       import array_t
       class(array_t), intent(inout) :: this
     end subroutine array_free
  end interface

contains

  !> Initialise a array of size @a n.
  subroutine array_init_base(this, n, name)
    class(array_t), intent(inout), target :: this
    integer, intent(in) :: n
    character(len=*), intent(in), optional :: name

    call this%alloc(n)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       ! Zero the device side first: under zero-copy the device then
       ! faults the pages (device first touch), which gives contiguous
       ! physical mappings and thus better GPU TLB utilisation
       call device_cfill(this%x_d, 0.0_rp, n)
       call device_sync()
    end if
    call cfill(this%data, 0.0_rp, n)

    if (present(name)) then
       this%name = name
    end if

  end subroutine array_init_base

  !> Vector allocation without initialisation.
  subroutine array_allocate(this, n)
    class(array_t), intent(inout) :: this
    integer, intent(in) :: n

    if (this%n .eq. n) return
    call this%free_base()

    this%n = n
    allocate(this%data(n))
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%data, this%x_d, n)
    end if

  end subroutine array_allocate

  !> Deallocate a array.
  subroutine array_free_base(this)
    class(array_t), intent(inout) :: this

    if (allocated(this%data)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_unmap(this%data, this%x_d)
       end if
       deallocate(this%data)
    end if

    this%n = 0
    this%name = ""

  end subroutine array_free_base

  !> Return the number of entries in the array.
  pure function array_size(this) result(s)
    class(array_t), intent(in) :: this
    integer :: s
    s = this%n
  end function array_size

  !> Easy way to copy between host and device.
  !! @param this array to copy to/from device/host
  !! @memdir direction to copy (HOST_TO_DEVICE or DEVICE_TO_HOST)
  !! @sync whether the memcopy to be blocking or not
  subroutine array_copy_from(this, memdir, sync)
    class(array_t), intent(inout) :: this
    integer, intent(in) :: memdir
    logical, intent(in) :: sync

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%data, this%x_d, this%size(), memdir, sync)
    end if

  end subroutine array_copy_from

  !> Check the the object is currently allocated
  pure function array_is_allocated(this) result(is_allocated)
    class(array_t), intent(in) :: this
    logical :: is_allocated
    is_allocated = allocated(this%data)
  end function array_is_allocated

  !> Assignment \f$ this = w \f$.
  subroutine array_assign_array(this, w)
    class(array_t), intent(inout) :: this
    class(array_t), intent(in) :: w

    if (this%size() .ne. w%size()) then
       call neko_error('Error in array assignment: incompatible size')
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%x_d, w%x_d, this%size())
    else
       call copy(this%data, w%data, this%size())
    end if

  end subroutine array_assign_array

  !> Assignment \f$ this = s \f$.
  subroutine array_assign_scalar(this, s)
    class(array_t), intent(inout) :: this
    real(kind=rp), intent(in) :: s

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(this%x_d, s, this%size())
    else
       call cfill(this%data, s, this%size())
    end if

  end subroutine array_assign_scalar

  !> Assignment \f$ this = array \f$.
  subroutine array_assign_real(this, array)
    class(array_t), intent(inout) :: this
    real(kind=rp), intent(in) :: array(:)

    if (this%size() .ne. size(array)) then
       call neko_error('Error in array assignment: incompatible size')
    end if

    call copy(this%data, array, this%size())
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call this%copy_from(HOST_TO_DEVICE, .true.)
    end if

  end subroutine array_assign_real

end module array
