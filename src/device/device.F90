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

  !> Copy data between host and device
  interface device_memcpy
     module procedure device_memcpy_r1, device_memcpy_r2, &
          device_memcpy_r3, device_memcpy_r4
  end interface device_memcpy

  !> Map a Fortran array to a device (allocate and associate)
  interface device_map
     module procedure device_map_r1, device_map_r2, &
          device_map_r3, device_map_r4
  end interface device_map

  !> Associate a Fortran array to a (allocated) device pointer
  interface device_associate
     module procedure device_associate_r1, device_associate_r2, &
          device_associate_r3, device_associate_r4
  end interface device_associate

  !> Check if a Fortran array is assoicated with a device pointer
  interface device_associated
     module procedure device_associated_r1, device_associated_r2, &
          device_associated_r3, device_associated_r4
  end interface device_associated

  !> Return the device pointer for an associated Fortran array
  interface device_get_ptr
     module procedure device_get_ptr_r1, device_get_ptr_r2, &
          device_get_ptr_r3, device_get_ptr_r4
  end interface device_get_ptr
      
  !> Table of host to device address mappings
  type(htable_cptr_t), private :: device_addrtbl

  private :: device_memcpy_common
  
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

  !> Copy data between host and device (rank 1 arrays)
  subroutine device_memcpy_r1(x, x_d, n, dir)
    integer, intent(in) :: n
    class(*), intent(inout), target :: x(:)
    type(c_ptr), intent(inout) :: x_d
    integer, intent(in), value :: dir
    type(c_ptr) :: ptr_h
    integer(c_size_t) :: s

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

    call device_memcpy_common(ptr_h, x_d, s, dir)
    
  end subroutine device_memcpy_r1

  !> Copy data between host and device (rank 2 arrays)
  subroutine device_memcpy_r2(x, x_d, n, dir)
    integer, intent(in) :: n
    class(*), intent(inout), target :: x(:,:)
    type(c_ptr), intent(inout) :: x_d
    integer, intent(in), value :: dir
    type(c_ptr) :: ptr_h
    integer(c_size_t) :: s

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

    call device_memcpy_common(ptr_h, x_d, s, dir)
    
  end subroutine device_memcpy_r2

  !> Copy data between host and device (rank 3 arrays)
  subroutine device_memcpy_r3(x, x_d, n, dir)
    integer, intent(in) :: n
    class(*), intent(inout), target :: x(:,:,:)
    type(c_ptr), intent(inout) :: x_d
    integer, intent(in), value :: dir
    type(c_ptr) :: ptr_h
    integer(c_size_t) :: s

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

    call device_memcpy_common(ptr_h, x_d, s, dir)
    
  end subroutine device_memcpy_r3

  !> Copy data between host and device (rank 4 arrays)
  subroutine device_memcpy_r4(x, x_d, n, dir)
    integer, intent(in) :: n
    class(*), intent(inout), target :: x(:,:,:,:)
    type(c_ptr), intent(inout) :: x_d
    integer, intent(in), value :: dir
    type(c_ptr) :: ptr_h
    integer(c_size_t) :: s

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

    call device_memcpy_common(ptr_h, x_d, s, dir)
    
  end subroutine device_memcpy_r4

  !> Copy data between host and device
  subroutine device_memcpy_common(ptr_h, x_d, s, dir)
    type(c_ptr), intent(inout) :: ptr_h
    type(c_ptr), intent(inout) :: x_d
    integer(c_size_t), intent(in) :: s
    integer, intent(in), value :: dir
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

  end subroutine device_memcpy_common

  !> Associate a Fortran rank 1 array to a (allocated) device pointer
  subroutine device_associate_r1(x, x_d, n)
    integer, intent(in) :: n
    class(*), intent(inout), target :: x(:)
    type(c_ptr), intent(inout) :: x_d
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(8))       
       htbl_ptr_h%ptr = c_loc(x)
    type is (real)
       htbl_ptr_h%ptr = c_loc(x)
    type is (double precision)
       htbl_ptr_h%ptr = c_loc(x)
    class default
       call neko_error('Unknown Fortran type')
    end select

    htbl_ptr_d%ptr = x_d
    
    call device_addrtbl%set(htbl_ptr_h, htbl_ptr_d)

  end subroutine device_associate_r1

  !> Associate a Fortran rank 2 array to a (allocated) device pointer
  subroutine device_associate_r2(x, x_d, n)
    integer, intent(in) :: n
    class(*), intent(inout), target :: x(:,:)
    type(c_ptr), intent(inout) :: x_d
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(8))       
       htbl_ptr_h%ptr = c_loc(x)
    type is (real)
       htbl_ptr_h%ptr = c_loc(x)
    type is (double precision)
       htbl_ptr_h%ptr = c_loc(x)
    class default
       call neko_error('Unknown Fortran type')
    end select

    htbl_ptr_d%ptr = x_d
    
    call device_addrtbl%set(htbl_ptr_h, htbl_ptr_d)

  end subroutine device_associate_r2

  !> Associate a Fortran rank 3 array to a (allocated) device pointer
  subroutine device_associate_r3(x, x_d, n)
    integer, intent(in) :: n
    class(*), intent(inout), target :: x(:,:,:)
    type(c_ptr), intent(inout) :: x_d
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(8))       
       htbl_ptr_h%ptr = c_loc(x)
    type is (real)
       htbl_ptr_h%ptr = c_loc(x)
    type is (double precision)
       htbl_ptr_h%ptr = c_loc(x)
    class default
       call neko_error('Unknown Fortran type')
    end select

    htbl_ptr_d%ptr = x_d
    
    call device_addrtbl%set(htbl_ptr_h, htbl_ptr_d)

  end subroutine device_associate_r3

  !> Associate a Fortran rank 4 array to a (allocated) device pointer
  subroutine device_associate_r4(x, x_d, n)
    integer, intent(in) :: n
    class(*), intent(inout), target :: x(:,:,:,:)
    type(c_ptr), intent(inout) :: x_d
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(8))       
       htbl_ptr_h%ptr = c_loc(x)
    type is (real)
       htbl_ptr_h%ptr = c_loc(x)
    type is (double precision)
       htbl_ptr_h%ptr = c_loc(x)
    class default
       call neko_error('Unknown Fortran type')
    end select

    htbl_ptr_d%ptr = x_d
    
    call device_addrtbl%set(htbl_ptr_h, htbl_ptr_d)

  end subroutine device_associate_r4
  
  !> Map a Fortran rank 1 array to a device (allocate and associate)
  subroutine device_map_r1(x, x_d, n)
    integer, intent(in) :: n
    class(*), intent(inout), target :: x(:)
    type(c_ptr), intent(inout) :: x_d
    integer(c_size_t) :: s

    if (c_associated(x_d)) then
       call neko_error('Device pointer already associated')
    end if

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

    call device_alloc(x_d, s)
    call device_associate(x, x_d, n)

  end subroutine device_map_r1

  !> Map a Fortran rank 2 array to a device (allocate and associate)
  subroutine device_map_r2(x, x_d, n)
    integer, intent(in) :: n
    class(*), intent(inout), target :: x(:,:)
    type(c_ptr), intent(inout) :: x_d
    integer(c_size_t) :: s

    if (c_associated(x_d)) then
       call neko_error('Device pointer already associated')
    end if

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

    call device_alloc(x_d, s)
    call device_associate(x, x_d, n)

  end subroutine device_map_r2

  !> Map a Fortran rank 3 array to a device (allocate and associate)
  subroutine device_map_r3(x, x_d, n)
    integer, intent(in) :: n
    class(*), intent(inout), target :: x(:,:,:)
    type(c_ptr), intent(inout) :: x_d
    integer(c_size_t) :: s

    if (c_associated(x_d)) then
       call neko_error('Device pointer already associated')
    end if

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

    call device_alloc(x_d, s)
    call device_associate(x, x_d, n)

  end subroutine device_map_r3

  !> Map a Fortran rank 4 array to a device (allocate and associate)
  subroutine device_map_r4(x, x_d, n)
    integer, intent(in) :: n
    class(*), intent(inout), target :: x(:,:,:,:)
    type(c_ptr), intent(inout) :: x_d
    integer(c_size_t) :: s

    if (c_associated(x_d)) then
       call neko_error('Device pointer already associated')
    end if

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

    call device_alloc(x_d, s)
    call device_associate(x, x_d, n)

  end subroutine device_map_r4

  !> Check if a Fortran rank 1 array is assoicated with a device pointer
  function device_associated_r1(x) result(assoc)
    class(*), intent(inout), target :: x(:)
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d
    logical :: assoc

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(8))
       htbl_ptr_h%ptr = c_loc(x)
    type is (real)
       htbl_ptr_h%ptr = c_loc(x)
    type is (double precision)
       htbl_ptr_h%ptr = c_loc(x)
    class default
       call neko_error('Unknown Fortran type')
    end select

    if (device_addrtbl%get(htbl_ptr_h, htbl_ptr_d) .eq. 0) then
       assoc = .true.
    else
       assoc = .false.       
    end if
    
  end function device_associated_r1

  !> Check if a Fortran rank 2 array is assoicated with a device pointer
  function device_associated_r2(x) result(assoc)
    class(*), intent(inout), target :: x(:,:)
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d
    logical :: assoc

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(8))
       htbl_ptr_h%ptr = c_loc(x)
    type is (real)
       htbl_ptr_h%ptr = c_loc(x)
    type is (double precision)
       htbl_ptr_h%ptr = c_loc(x)
    class default
       call neko_error('Unknown Fortran type')
    end select

    if (device_addrtbl%get(htbl_ptr_h, htbl_ptr_d) .eq. 0) then
       assoc = .true.
    else
       assoc = .false.       
    end if
    
  end function device_associated_r2

  !> Check if a Fortran rank 3 array is assoicated with a device pointer
  function device_associated_r3(x) result(assoc)
    class(*), intent(inout), target :: x(:,:,:)
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d
    logical :: assoc

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(8))
       htbl_ptr_h%ptr = c_loc(x)
    type is (real)
       htbl_ptr_h%ptr = c_loc(x)
    type is (double precision)
       htbl_ptr_h%ptr = c_loc(x)
    class default
       call neko_error('Unknown Fortran type')
    end select

    if (device_addrtbl%get(htbl_ptr_h, htbl_ptr_d) .eq. 0) then
       assoc = .true.
    else
       assoc = .false.       
    end if
    
  end function device_associated_r3

  !> Check if a Fortran rank 4 array is assoicated with a device pointer
  function device_associated_r4(x) result(assoc)
    class(*), intent(inout), target :: x(:,:,:,:)
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d
    logical :: assoc

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(8))
       htbl_ptr_h%ptr = c_loc(x)
    type is (real)
       htbl_ptr_h%ptr = c_loc(x)
    type is (double precision)
       htbl_ptr_h%ptr = c_loc(x)
    class default
       call neko_error('Unknown Fortran type')
    end select

    if (device_addrtbl%get(htbl_ptr_h, htbl_ptr_d) .eq. 0) then
       assoc = .true.
    else
       assoc = .false.       
    end if
    
  end function device_associated_r4

  !> Return the device pointer for an associated Fortran rank 1 array
  function device_get_ptr_r1(x, n)
    integer, intent(in) :: n
    class(*), intent(in), target :: x(:)
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d
    type(c_ptr) :: device_get_ptr_r1

    device_get_ptr_r1 = C_NULL_PTR

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(8))
       htbl_ptr_h%ptr = c_loc(x)
    type is (real)
       htbl_ptr_h%ptr = c_loc(x)
    type is (double precision)
       htbl_ptr_h%ptr = c_loc(x)
    class default
       call neko_error('Unknown Fortran type')
    end select

    if (device_addrtbl%get(htbl_ptr_h, htbl_ptr_d) .eq. 0) then
       device_get_ptr_r1 = htbl_ptr_d%ptr
    else
       call neko_error('Array not associated with device')
    end if
  end function device_get_ptr_r1

  !> Return the device pointer for an associated Fortran rank 2 array
  function device_get_ptr_r2(x, n)
    integer, intent(in) :: n
    class(*), intent(in), target :: x(:,:)
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d
    type(c_ptr) :: device_get_ptr_r2

    device_get_ptr_r2 = C_NULL_PTR

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(8))
       htbl_ptr_h%ptr = c_loc(x)
    type is (real)
       htbl_ptr_h%ptr = c_loc(x)
    type is (double precision)
       htbl_ptr_h%ptr = c_loc(x)
    class default
       call neko_error('Unknown Fortran type')
    end select

    if (device_addrtbl%get(htbl_ptr_h, htbl_ptr_d) .eq. 0) then
       device_get_ptr_r2 = htbl_ptr_d%ptr
    else
       call neko_error('Array not associated with device')
    end if
  end function device_get_ptr_r2

  !> Return the device pointer for an associated Fortran rank 3 array
  function device_get_ptr_r3(x, n)
    integer, intent(in) :: n
    class(*), intent(in), target :: x(:,:,:)
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d
    type(c_ptr) :: device_get_ptr_r3

    device_get_ptr_r3 = C_NULL_PTR

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(8))
       htbl_ptr_h%ptr = c_loc(x)
    type is (real)
       htbl_ptr_h%ptr = c_loc(x)
    type is (double precision)
       htbl_ptr_h%ptr = c_loc(x)
    class default
       call neko_error('Unknown Fortran type')
    end select

    if (device_addrtbl%get(htbl_ptr_h, htbl_ptr_d) .eq. 0) then
       device_get_ptr_r3 = htbl_ptr_d%ptr
    else
       call neko_error('Array not associated with device')
    end if
  end function device_get_ptr_r3

  !> Return the device pointer for an associated Fortran rank 4 array
  function device_get_ptr_r4(x, n)
    integer, intent(in) :: n
    class(*), intent(in), target :: x(:,:,:,:)
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d
    type(c_ptr) :: device_get_ptr_r4

    device_get_ptr_r4 = C_NULL_PTR

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(8))
       htbl_ptr_h%ptr = c_loc(x)
    type is (real)
       htbl_ptr_h%ptr = c_loc(x)
    type is (double precision)
       htbl_ptr_h%ptr = c_loc(x)
    class default
       call neko_error('Unknown Fortran type')
    end select

    if (device_addrtbl%get(htbl_ptr_h, htbl_ptr_d) .eq. 0) then
       device_get_ptr_r4 = htbl_ptr_d%ptr
    else
       call neko_error('Array not associated with device')
    end if
  end function device_get_ptr_r4
  
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
