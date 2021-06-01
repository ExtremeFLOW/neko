module device
  use num_types
  use hip_intf
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

  integer, parameter :: HOST_TO_DEVICE = 1, DEVICE_TO_HOST = 2
  
contains

  !> Allocate memory on the device
  subroutine device_alloc(x_d, s)
    type(c_ptr), intent(inout) :: x_d
    integer(c_size_t) :: s
#ifdef HAVE_HIP
    if (hipmalloc(x_d, s) .ne. HIP_SUCCESS) then
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
#endif
  end subroutine device_free

  !> Copy data between host and device
  subroutine device_memcpy(x, x_d, n, dir)
    integer, intent(in) :: n
    class(*), intent(inout), target :: x(n)
    type(c_ptr), intent(inout) :: x_d
    integer, intent(in), value :: dir
    integer(c_size_t) :: s

    select type(x)
    type is (integer)
       s = n * 4
    type is (integer(8))
       s = n * 8
    type is (real)
       s = n * 4
    type is (double precision)
       s = n * 8
    end select
#ifdef HAVE_HIP
    if (hipMemcpy(c_loc(x), x_d, s, dir) .ne. HIP_SUCCESS) then
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

  !> Synchronize the device
  subroutine device_sync()
#ifdef HAVE_HIP
    if (hipdevicesynchronize() .ne. HIP_SUCCESS) then
       call neko_error('Error during device sync')
    end if
#endif
  end subroutine device_sync
  
end module device
