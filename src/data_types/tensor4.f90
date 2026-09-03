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
!> Defines a rank-4 tensor
module tensor4
  use neko_config, only : NEKO_BCKND_DEVICE
  use num_types, only : rp
  use device, only : device_map, device_unmap, device_memcpy, &
       device_sync
  use device_math, only : device_copy, device_cfill
  use utils, only : neko_error, NEKO_VARNAME_LEN
  use, intrinsic :: iso_c_binding
  implicit none
  private

  type, public :: tensor4_t
     real(kind=rp), allocatable :: x(:,:,:,:) !< Tensor entries.
     character(len=NEKO_VARNAME_LEN) :: name = "" !< Name of the tensor
     type(c_ptr) :: x_d = C_NULL_PTR !< Device pointer.
     integer, private :: n1 = 0 !< Size of the first dimension.
     integer, private :: n2 = 0 !< Size of the second dimension.
     integer, private :: n3 = 0 !< Size of the third dimension.
     integer, private :: n4 = 0 !< Size of the fourth dimension.
     integer, private :: n = 0 !< Total size n1*n2*n3*n4.
   contains
     !> Initialise a tensor of size `n1*n2*n3*n4`.
     procedure, pass(t), private :: init_dims => tensor4_init
     !> Initialise a tensor of size `n1*n2*n3*n4`.
     !! @note Declared as a generic (rather than a plain binding) so that
     !! types extending tensor4_t can override `init_dims` and have that
     !! override correctly replace this specific within the inherited
     !! `init` generic, instead of ambiguously merging with it.
     generic :: init => init_dims
     !> Deallocate a tensor.
     procedure, pass(t) :: free => tensor4_free
     !> Returns the number of entries in the tensor.
     procedure, pass(t) :: size => tensor4_size
     !> Copy between host and device
     procedure, pass(t) :: copy_from => tensor4_copy_from
     !> Returns the size of the first dimension.
     procedure, pass(t) :: get_n1 => tensor4_n1
     !> Returns the size of the second dimension.
     procedure, pass(t) :: get_n2 => tensor4_n2
     !> Returns the size of the third dimension.
     procedure, pass(t) :: get_n3 => tensor4_n3
     !> Returns the size of the fourth dimension.
     procedure, pass(t) :: get_n4 => tensor4_n4
     !> Assignment \f$ t = w \f$
     procedure, pass(t) :: tensor4_assign_tensor4
     !> Assignment \f$ t = s \f$.
     procedure, pass(t) :: tensor4_assign_scalar

     generic :: assignment(=) => tensor4_assign_tensor4, &
          tensor4_assign_scalar

     !> Allocate a tensor of size `n1*n2*n3*n4`.
     procedure, pass(t), private :: alloc => tensor4_allocate
  end type tensor4_t

  type, public :: tensor4_ptr_t
     type(tensor4_t), pointer :: ptr => null()
   contains
     !> Constructor. Just assigns the pointer
     procedure, pass(this) :: init => tensor4_ptr_init
     !> Destructor. Just nullifies the pointer.
     procedure, pass(this) :: free => tensor4_ptr_free
  end type tensor4_ptr_t

contains

  !> Initialise a tensor of size `n1*n2*n3*n4`.
  !! @param t Tensor to initialise.
  !! @param n1 Size of the first dimension.
  !! @param n2 Size of the second dimension.
  !! @param n3 Size of the third dimension.
  !! @param n4 Size of the fourth dimension.
  !! @param name Optional name of the tensor.
  subroutine tensor4_init(t, n1, n2, n3, n4, name)
    class(tensor4_t), intent(inout) :: t
    integer, intent(in) :: n1
    integer, intent(in) :: n2
    integer, intent(in) :: n3
    integer, intent(in) :: n4
    character(len=*), intent(in), optional :: name

    ! t%alloc zeroes the device side (and synchronizes) before any
    ! host-side touch: under zero-copy the device then faults the
    ! pages first (device first touch), which gives contiguous
    ! physical mappings and thus better GPU TLB utilisation; rewriting
    ! the zeros on the host afterwards is benign.
    call t%alloc(n1, n2, n3, n4)
    t%x = 0.0_rp

    if (present(name)) then
       t%name = name
    end if

  end subroutine tensor4_init

  !> Allocate a tensor of size `n1*n2*n3*n4`.
  !! @param t Tensor to allocate.
  !! @param n1 Size of the first dimension.
  !! @param n2 Size of the second dimension.
  !! @param n3 Size of the third dimension.
  !! @param n4 Size of the fourth dimension.
  subroutine tensor4_allocate(t, n1, n2, n3, n4)
    class(tensor4_t), intent(inout) :: t
    integer, intent(in) :: n1
    integer, intent(in) :: n2
    integer, intent(in) :: n3
    integer, intent(in) :: n4

    call t%free()

    allocate(t%x(n1, n2, n3, n4))
    t%n1 = n1
    t%n2 = n2
    t%n3 = n3
    t%n4 = n4
    t%n = n1*n2*n3*n4

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_map(t%x, t%x_d, t%n)
       call device_cfill(t%x_d, 0.0_rp, t%n)
       call device_sync()
    end if

  end subroutine tensor4_allocate

  !> Deallocate a tensor.
  !! @param t Tensor to deallocate.
  subroutine tensor4_free(t)
    class(tensor4_t), intent(inout) :: t

    if (allocated(t%x)) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_unmap(t%x, t%x_d)
       end if
       deallocate(t%x)
    end if

    t%n1 = 0
    t%n2 = 0
    t%n3 = 0
    t%n4 = 0
    t%n = 0
    t%name = ""

  end subroutine tensor4_free

  !> Returns the number of entries in the tensor.
  !! @param t Tensor to query.
  pure function tensor4_size(t) result(s)
    class(tensor4_t), intent(in) :: t
    integer :: s
    s = t%n
  end function tensor4_size

  !> Easy way to copy between host and device.
  !! @param t Tensor to copy to/from device/host.
  !! @param memdir Direction to copy (HOST_TO_DEVICE or DEVICE_TO_HOST).
  !! @param sync Whether the memcopy is to be blocking or not.
  subroutine tensor4_copy_from(t, memdir, sync)
    class(tensor4_t), intent(inout) :: t
    integer, intent(in) :: memdir
    logical, intent(in) :: sync

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(t%x, t%x_d, t%n, memdir, sync)
    end if

  end subroutine tensor4_copy_from

  !> Returns the size of the first dimension.
  !! @param t Tensor to query.
  pure function tensor4_n1(t) result(n1)
    class(tensor4_t), intent(in) :: t
    integer :: n1
    n1 = t%n1
  end function tensor4_n1

  !> Returns the size of the second dimension.
  !! @param t Tensor to query.
  pure function tensor4_n2(t) result(n2)
    class(tensor4_t), intent(in) :: t
    integer :: n2
    n2 = t%n2
  end function tensor4_n2

  !> Returns the size of the third dimension.
  !! @param t Tensor to query.
  pure function tensor4_n3(t) result(n3)
    class(tensor4_t), intent(in) :: t
    integer :: n3
    n3 = t%n3
  end function tensor4_n3

  !> Returns the size of the fourth dimension.
  !! @param t Tensor to query.
  pure function tensor4_n4(t) result(n4)
    class(tensor4_t), intent(in) :: t
    integer :: n4
    n4 = t%n4
  end function tensor4_n4

  !> Assignment \f$ t = w \f$
  !! @param t Tensor to assign to.
  !! @param w Tensor to assign from.
  subroutine tensor4_assign_tensor4(t, w)
    class(tensor4_t), intent(inout) :: t
    type(tensor4_t), intent(in) :: w

    if (allocated(t%x)) then
       call t%free()
    end if

    if (.not. allocated(t%x)) then

       t%n1 = w%n1
       t%n2 = w%n2
       t%n3 = w%n3
       t%n4 = w%n4
       t%n = w%n
       allocate(t%x(t%n1, t%n2, t%n3, t%n4))

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_map(t%x, t%x_d, t%n)
       end if

    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(t%x_d, w%x_d, t%n)
    else
       t%x = w%x
    end if

    t%name = w%name

  end subroutine tensor4_assign_tensor4

  !> Assignment \f$ t = s \f$.
  !! @param t Tensor to assign to.
  !! @param s Scalar to fill the tensor with.
  subroutine tensor4_assign_scalar(t, s)
    class(tensor4_t), intent(inout) :: t
    real(kind=rp), intent(in) :: s

    if (.not. allocated(t%x)) then
       call neko_error('tensor4 not allocated')
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(t%x_d, s, t%n)
    else
       t%x = s
    end if

  end subroutine tensor4_assign_scalar

  ! ========================================================================== !
  ! tensor4 pointer type subroutines

  !> Constructor. Just assigns the pointer.
  !! @param this Pointer wrapper to initialise.
  !! @param ptr Tensor to point to.
  subroutine tensor4_ptr_init(this, ptr)
    class(tensor4_ptr_t), intent(inout) :: this
    type(tensor4_t), target, intent(in) :: ptr

    call this%free()
    this%ptr => ptr
  end subroutine tensor4_ptr_init

  !> Destructor. Just nullifies the pointer.
  !! @param this Pointer wrapper to destroy.
  subroutine tensor4_ptr_free(this)
    class(tensor4_ptr_t), intent(inout) :: this

    if (associated(this%ptr)) then
       nullify(this%ptr)
    end if

  end subroutine tensor4_ptr_free

end module tensor4
