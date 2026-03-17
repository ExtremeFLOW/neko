! Copyright (c) 2022-2025, The Neko Authors
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
!> Defines a vector
module vector
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use device, only : device_map, device_unmap, &
       device_memcpy, device_sync, HOST_TO_DEVICE
  use math, only : cfill, copy
  use device_math, only : device_copy, device_cfill, device_cmult, &
       device_sub3, device_cmult2, device_add3, device_cadd2, device_col3, &
       device_col2, device_invcol3, device_cdiv2
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding
  implicit none
  private

  type, public :: vector_t
     !> Vector entries.
     real(kind=rp), allocatable :: x(:)
     character(len=80) :: name = "" !< Name of the vector
     !> Device pointer.
     type(c_ptr) :: x_d = C_NULL_PTR
     !> Size of vector.
     integer, private :: n = 0
   contains
     !> Initialise a vector of size `n`.
     procedure, pass(v) :: init => vector_init
     !> Deallocate a vector.
     procedure, pass(v) :: free => vector_free
     !> Copy data between host and device
     procedure, pass(v) :: copy_from => vector_copy_from
     !> Returns the number of entries in the vector.
     procedure, pass(v) :: size => vector_size
     !> Assignment \f$ v = w \f$
     procedure, pass(v) :: vector_assign_vector
     !> Assignment \f$ v = s \f$.
     procedure, pass(v) :: vector_assign_scalar
     !> Assignment \f$ v = array \f$.
     procedure, pass(v) :: vector_assign_array

     !> Assignments
     generic :: assignment(=) => vector_assign_vector, &
          vector_assign_scalar, vector_assign_array

     ! Private interfaces
     procedure, pass(a), private :: alloc => vector_allocate

  end type vector_t

  type, public :: vector_ptr_t
     type(vector_t), pointer :: ptr => null()
   contains
     !> Constructor. Just assigns the pointer
     procedure, pass(this) :: init => vector_ptr_init
     !> Destructor. Just nullifies the pointer.
     procedure, pass(this) :: free => vector_ptr_free
  end type vector_ptr_t

contains

  !> Initialise a vector of size @a n.
  subroutine vector_init(v, n, name)
    class(vector_t), intent(inout) :: v
    integer, intent(in) :: n
    character(len=*), intent(in), optional :: name

    call v%alloc(n)
    call cfill(v%x, 0.0_rp, n)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(v%x_d, 0.0_rp, n)
       call device_sync()
    end if

    if (present(name)) then
       v%name = name
    end if

  end subroutine vector_init

  !> Vector allocation without initialisation.
  subroutine vector_allocate(a, n)
    class(vector_t), intent(inout) :: a
    integer, intent(in) :: n


    if (a%n .eq. n) return
    call a%free()

    a%n = n
    allocate(a%x(n))
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(a%x, a%x_d, n)
    end if

  end subroutine vector_allocate

  !> Deallocate a vector.
  subroutine vector_free(v)
    class(vector_t), intent(inout) :: v

    if (allocated(v%x)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_unmap(v%x, v%x_d)
       end if
       deallocate(v%x)
    end if

    v%n = 0
    v%name = ""

  end subroutine vector_free

  !> Return the number of entries in the vector.
  pure function vector_size(v) result(s)
    class(vector_t), intent(in) :: v
    integer :: s
    s = v%n
  end function vector_size

  !> Easy way to copy between host and device.
  !! @param v vector to copy to/from device/host
  !! @memdir direction to copy (HOST_TO_DEVICE or DEVICE_TO_HOST)
  !! @sync whether the memcopy to be blocking or not
  subroutine vector_copy_from(v, memdir, sync)
    class(vector_t), intent(inout) :: v
    integer, intent(in) :: memdir
    logical, intent(in) :: sync

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(v%x, v%x_d, v%n, memdir, sync)
    end if

  end subroutine vector_copy_from


  !> Assignment \f$ v = w \f$.
  subroutine vector_assign_vector(v, w)
    class(vector_t), intent(inout) :: v
    type(vector_t), intent(in) :: w

    call v%alloc(w%n)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(v%x_d, w%x_d, v%n)
    else
       call copy(v%x, w%x, v%n)
    end if

    v%name = w%name

  end subroutine vector_assign_vector

  !> Assignment \f$ v = s \f$.
  subroutine vector_assign_scalar(v, s)
    class(vector_t), intent(inout) :: v
    real(kind=rp), intent(in) :: s

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(v%x_d, s, v%n)
    else
       call cfill(v%x, s, v%n)
    end if

  end subroutine vector_assign_scalar

  !> Assignment \f$ v = array \f$.
  subroutine vector_assign_array(v, array)
    class(vector_t), intent(inout) :: v
    real(kind=rp), intent(in) :: array(:)

    call v%alloc(size(array))
    v%x = array

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call v%copy_from(HOST_TO_DEVICE, .true.)
    end if

  end subroutine vector_assign_array

  ! ========================================================================== !
  ! vector pointer type subroutines

  subroutine vector_ptr_init(this, ptr)
    class(vector_ptr_t), intent(inout) :: this
    type(vector_t), target, intent(in) :: ptr

    call this%free()
    this%ptr => ptr
  end subroutine vector_ptr_init

  subroutine vector_ptr_free(this)
    class(vector_ptr_t), intent(inout) :: this

    if (associated(this%ptr)) then
       nullify(this%ptr)
    end if

  end subroutine vector_ptr_free

end module vector
