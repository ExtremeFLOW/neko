!> Device abstraction, common interface for various accelerators
module device
  use num_types
  use cuda_intf
  use hip_intf
  use htable
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

  integer, parameter :: HOST_TO_DEVICE = 1, DEVICE_TO_HOST = 2

  !> Table of host to device address mappings
  type(htable_cptr_t), private :: device_addrtbl
  
contains

  subroutine device_init
#if defined(HAVE_HIP) || defined(HAVE_CUDA)
    call device_addrtbl%init(64)
#endif
  end subroutine device_init

  subroutine device_finalize
#if defined(HAVE_HIP) || defined(HAVE_CUDA)
    call device_addrtbl%free()
#endif
  end subroutine device_finalize
  
  !> Allocate memory on the device
  subroutine device_alloc(x_d, s)
    type(c_ptr), intent(inout) :: x_d
    integer(c_size_t) :: s
#ifdef HAVE_HIP
    if (hipmalloc(x_d, s) .ne. HIP_SUCCESS) then
       call neko_error('Memory allocation on device failed')
    end if
#elif HAVE_CUDA
    if (cudamalloc(x_d, s) .ne. CUDA_SUCCESS) then
       call neko_error('Memory allocation on device failed')
    end if
#endif
  end subroutine device_alloc

  !> Deallocate memory on the device
  subroutine device_free(x_d)
    type(c_ptr), intent(inout) :: x_d
#ifdef HAVE_HIP
    if (hipfree(x_d) .ne. HIP_SUCCESS) then
       call neko_error('Memory deallocation on device failed')
    end if
#elif HAVE_CUDA
    if (cudafree(x_d) .ne. CUDA_SUCCESS) then
       call neko_error('Memory deallocation on device failed')
    end if
#endif
  end subroutine device_free

  !> Copy data between host and device
  subroutine device_memcpy(x, x_d, n, dir)
    integer, intent(in) :: n
    class(*), intent(inout), target :: x(..)
    type(c_ptr), intent(inout) :: x_d
    integer, intent(in), value :: dir
    type(c_ptr) :: ptr_h
    integer(c_size_t) :: s

    select rank(x)
    rank(1)
       select type(x)
       type is (integer)
          s = n * 4
          ptr_h = c_loc(x)
       type is (integer(8))       
          s = n * 8
          ptr_h = c_loc(x)
       type is (real)
          s = n * 4
          ptr_h = c_loc(x)
       type is (double precision)
          s = n * 8
          ptr_h = c_loc(x)
       class default
          call neko_error('Unknown Fortran type')
       end select
    rank(2)
       select type(x)
       type is (integer)
          s = n * 4
          ptr_h = c_loc(x)
       type is (integer(8))       
          s = n * 8
          ptr_h = c_loc(x)
       type is (real)
          s = n * 4
          ptr_h = c_loc(x)
       type is (double precision)
          s = n * 8
          ptr_h = c_loc(x)
       class default
          call neko_error('Unknown Fortran type')
       end select
    rank(3)
       select type(x)
       type is (integer)
          s = n * 4
          ptr_h = c_loc(x)
       type is (integer(8))       
          s = n * 8
          ptr_h = c_loc(x)
       type is (real)
          s = n * 4
          ptr_h = c_loc(x)
       type is (double precision)
          s = n * 8
          ptr_h = c_loc(x)
       class default
          call neko_error('Unknown Fortran type')
       end select
    rank(4)
       select type(x)
       type is (integer)
          s = n * 4
          ptr_h = c_loc(x)
       type is (integer(8))       
          s = n * 8
          ptr_h = c_loc(x)
       type is (real)
          s = n * 4
          ptr_h = c_loc(x)
       type is (double precision)
          s = n * 8
          ptr_h = c_loc(x)
       class default
          call neko_error('Unknown Fortran type')
       end select
    end select

#ifdef HAVE_HIP
    if (hipmemcpy(ptr_h, x_d, s, dir) .ne. HIP_SUCCESS) then
       if (dir .eq. HOST_TO_DEVICE) then
          call neko_error('Device memcpy (host-to-device) failed')
       else if (dir .eq. DEVICE_TO_HOST) then
          call neko_error('Device memcpy (device-to-host) failed')
       else
          call neko_error('Device memcpy failed (invalid direction')
       end if
    end if
#elif HAVE_CUDA
    if (cudamemcpy(ptr_h, x_d, s, dir) .ne. CUDA_SUCCESS) then
       if (dir .eq. HOST_TO_DEVICE) then
          call neko_error('Device memcpy (host-to-device) failed')
       else if (dir .eq. DEVICE_TO_HOST) then
          call neko_error('Device memcpy (device-to-host) failed')
       else
          call neko_error('Device memcpy failed (invalid direction')
       end if
    end if
#endif
    
  end subroutine device_memcpy

  !> Associate a Fortran array to a (allocated) device pointer
  subroutine device_associate(x, x_d, n)
    integer, intent(in) :: n
    type(*), intent(inout), target :: x(..)
    type(c_ptr), intent(inout) :: x_d
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d

    htbl_ptr_h%ptr = c_loc(x)
    htbl_ptr_d%ptr = x_d
    
    call device_addrtbl%set(htbl_ptr_h, htbl_ptr_d)

  end subroutine device_associate
  
  !> Map a Fortran array to a device (allocate and associate)
  subroutine device_map(x, x_d, n)
    integer, intent(in) :: n
    class(*), intent(inout), target :: x(..)
    type(c_ptr), intent(inout) :: x_d
    integer(c_size_t) :: s

    if (c_associated(x_d)) then
       call neko_error('Device pointer already associated')
    end if

    select rank(x)
    rank (1)
       select type(x)
       type is (integer)
          s = n * 4
       type is (integer(8))
          s = n * 8
       type is (real)
          s = n * 4
       type is (double precision)
          s = n * 8
       class default
          call neko_error('Unknown Fortran type')
       end select
    rank (2)
       select type(x)
       type is (integer)
          s = n * 4
       type is (integer(8))
          s = n * 8
       type is (real)
          s = n * 4
       type is (double precision)
          s = n * 8
       class default
          call neko_error('Unknown Fortran type')
       end select
    rank (3)
       select type(x)
       type is (integer)
          s = n * 4
       type is (integer(8))
          s = n * 8
       type is (real)
          s = n * 4
       type is (double precision)
          s = n * 8
       class default
          call neko_error('Unknown Fortran type')
       end select
    rank (4)
       select type(x)
       type is (integer)
          s = n * 4
       type is (integer(8))
          s = n * 8
       type is (real)
          s = n * 4
       type is (double precision)
          s = n * 8
       class default
          call neko_error('Unknown Fortran type')
       end select
    end select

    call device_alloc(x_d, s)
    call device_associate(x, x_d, n)

  end subroutine device_map
 
  !> Check if a Fortran array is assoicated with a device pointer
  function device_associated(x) result(assoc)
    type(*), intent(inout), target :: x(..)
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d
    logical :: assoc

    htbl_ptr_h%ptr = c_loc(x)

    if (device_addrtbl%get(htbl_ptr_h, htbl_ptr_d) .eq. 0) then
       assoc = .true.
    else
       assoc = .false.       
    end if
    
  end function device_associated
  
  !> Return the device pointer for an associated Fortran
  function device_get_ptr(x, n)
    integer, intent(in) :: n
    type(*), intent(in), target :: x(..)
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d
    type(c_ptr) :: device_get_ptr

    device_get_ptr = C_NULL_PTR
    htbl_ptr_h%ptr = c_loc(x)

    if (device_addrtbl%get(htbl_ptr_h, htbl_ptr_d) .eq. 0) then
       device_get_ptr = htbl_ptr_d%ptr
    else
       call neko_error('Array not associated with device')
    end if
  end function device_get_ptr
  
  !> Synchronize the device
  subroutine device_sync()
#ifdef HAVE_HIP
    if (hipdevicesynchronize() .ne. HIP_SUCCESS) then
       call neko_error('Error during device sync')
    end if
#elif HAVE_CUDA
    if (cudadevicesynchronize() .ne. CUDA_SUCCESS) then
       call neko_error('Error during device sync')
    end if
#endif
  end subroutine device_sync
  
end module device
