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
!
!> Object for handling masks in Neko.
module mask
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, &
       c_size_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: device_map, device_free, device_memcpy, &
       HOST_TO_DEVICE, DEVICE_TO_HOST, DEVICE_TO_DEVICE
  use device_math, only: device_cadd
  use utils, only: neko_error

  implicit none
  private

  !> Type for consistently handling masks in Neko.
  !! This type encapsulates the mask array and its associated device pointer.
  !!
  !! @note The mask is platform-agnostic, meaning it can be used on both CPU and
  !! GPU, without adjusting for 1-based or 0-based indexing.
  type, public :: mask_t
     private
     integer :: n_elements = 0 ! Number of elements in the mask
     integer, allocatable :: mask(:) ! The mask array
     type(c_ptr) :: mask_d = c_null_ptr ! Pointer to the device mask array
     logical :: is_set_ = .false. ! Flag to indicate if the mask is set

   contains
     ! Public interface for the mask type
     generic, public :: init => init_from_array, init_from_array_device, &
          init_from_mask
     procedure, public, pass(this) :: free => mask_free

     procedure, public, pass(this) :: size => mask_size
     procedure, public, pass(this) :: is_set => mask_is_set
     procedure, public, pass(this) :: get_d => mask_get_d

     generic, public :: set => mask_set, mask_set_d
     generic, public :: get => mask_get, mask_get_i

     ! Private procedures
     procedure, pass(this) :: allocate => mask_allocate
     procedure, pass(this) :: init_from_array
     procedure, pass(this) :: init_from_array_device
     procedure, pass(this) :: init_from_mask

     ! Setters
     procedure, pass(this) :: mask_set
     procedure, pass(this) :: mask_set_d

     ! Getters
     procedure, pass(this) :: mask_get
     procedure, pass(this) :: mask_get_i

  end type mask_t

contains

  !> Allocate the mask object.
  subroutine mask_allocate(this, n_elements)
    class(mask_t), intent(inout) :: this
    integer, intent(in) :: n_elements

    this%is_set_ = .false.
    if (n_elements .eq. this%n_elements) return
    call this%free()

    allocate(this%mask(n_elements))

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(this%mask, this%mask_d, n_elements)
    end if

    this%n_elements = n_elements
  end subroutine mask_allocate

  !> Free the mask object.
  subroutine mask_free(this)
    class(mask_t), intent(inout) :: this

    if (allocated(this%mask)) then
       deallocate(this%mask)
    end if

    if (c_associated(this%mask_d)) then
       call device_free(this%mask_d)
    end if

    this%n_elements = 0
    this%mask_d = c_null_ptr
    this%is_set_ = .false.
  end subroutine mask_free

  !> Initialize the mask from a 1-indexed host array.
  subroutine init_from_array(this, mask_array, n_elements)
    class(mask_t), intent(inout) :: this
    integer, intent(in) :: n_elements
    integer, intent(in) :: mask_array(n_elements)

    call this%allocate(n_elements)

    this%mask = mask_array
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%mask, this%mask_d, this%n_elements, &
            HOST_TO_DEVICE, sync = .true.)
       call device_cadd(this%mask_d, -1, this%n_elements)
    end if

    this%is_set_ = .true.
  end subroutine init_from_array

  !> Initialize the mask from a 0-indexed device array.
  subroutine init_from_array_device(this, mask_array_d, n_elements)
    class(mask_t), intent(inout) :: this
    integer, intent(in) :: n_elements
    type(c_ptr), intent(inout):: mask_array_d
    integer(kind=c_size_t) :: size_c

    size_c = n_elements * int(4, c_size_t)

    call this%allocate(n_elements)
    call device_memcpy(this%mask_d, mask_array_d, size_c, &
         DEVICE_TO_DEVICE, sync = .false.)
    call device_memcpy(this%mask, mask_array_d, n_elements, &
         DEVICE_TO_HOST, sync = .true.)
    this%mask = this%mask - 1 ! Adjust for 0-based indexing

    this%is_set_ = .true.
  end subroutine init_from_array_device

  !> Initialize the mask from another mask object.
  subroutine init_from_mask(this, other)
    class(mask_t), intent(inout) :: this
    class(mask_t), intent(inout) :: other
    integer(kind=c_size_t) :: size_c

    call this%allocate(other%n_elements)

    size_c = other%n_elements * int(4, c_size_t)

    this%mask = other%mask
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%mask_d, other%mask_d, size_c, &
            DEVICE_TO_DEVICE, sync = .true.)
    end if

    this%n_elements = other%n_elements
    this%is_set_ = other%is_set_
  end subroutine init_from_mask

  !> Get the size of the mask.
  pure function mask_size(this) result(n_elements)
    class(mask_t), intent(in) :: this
    integer :: n_elements

    n_elements = this%n_elements
  end function mask_size

  !> Check if the mask is set.
  pure function mask_is_set(this) result(is_set)
    class(mask_t), intent(in) :: this
    logical :: is_set

    is_set = this%is_set_
  end function mask_is_set

  !> Get the mask array.
  function mask_get(this) result(mask_array)
    class(mask_t), intent(in), target :: this
    integer, pointer :: mask_array(:)

    if (.not. this%is_set()) call neko_error("Mask is not set.")

    mask_array => this%mask
  end function mask_get

  !> Get the mask array.
  function mask_get_i(this, index) result(mask_value)
    class(mask_t), intent(in), target :: this
    integer, intent(in) :: index
    integer :: mask_value

    if (.not. this%is_set()) call neko_error("Mask is not set.")
    if (index < 1 .or. index > this%n_elements) then
       call neko_error("Index out of bounds in mask_get_i")
    end if

    mask_value = this%mask(index)
  end function mask_get_i

  !> Get the device pointer to the mask array.
  function mask_get_d(this) result(mask_array_d)
    class(mask_t), intent(in) :: this
    type(c_ptr) :: mask_array_d

    if (.not. this%is_set()) call neko_error("Mask is not set.")

    mask_array_d = this%mask_d
  end function mask_get_d

  !> Set the mask from a 1-indexed host array.
  subroutine mask_set(this, mask_array, n_elements)
    class(mask_t), intent(inout) :: this
    integer, intent(in) :: n_elements
    integer, intent(in) :: mask_array(n_elements)

    call this%allocate(n_elements)

    this%mask = mask_array
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(this%mask, this%mask_d, n_elements, &
            HOST_TO_DEVICE, sync = .true.)
       call device_cadd(this%mask_d, -1, n_elements)
    end if

    this%is_set_ = .true.
  end subroutine mask_set

  !> Set the mask from a 0-indexed device array.
  subroutine mask_set_d(this, mask_array_d, n_elements)
    class(mask_t), intent(inout) :: this
    integer, intent(in) :: n_elements
    type(c_ptr), intent(inout) :: mask_array_d
    integer(kind=c_size_t) :: size_c

    call this%allocate(n_elements)
    size_c = n_elements * int(4, c_size_t)

    call device_memcpy(this%mask_d, mask_array_d, size_c, &
         DEVICE_TO_DEVICE, sync = .false.)
    call device_memcpy(this%mask, mask_array_d, n_elements, &
         DEVICE_TO_HOST, sync = .true.)
    this%mask = this%mask - 1 ! Adjust for 0-based indexing

    this%is_set_ = .true.
  end subroutine mask_set_d


end module mask
