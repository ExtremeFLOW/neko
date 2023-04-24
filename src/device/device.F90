! Copyright (c) 2021-2022, The Neko Authors
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
!> Device abstraction, common interface for various accelerators
module device
  use num_types
  use opencl_intf
  use cuda_intf
  use hip_intf
  use htable
  use utils
  use opencl_prgm_lib
  use, intrinsic :: iso_c_binding
  implicit none

  integer, parameter :: HOST_TO_DEVICE = 1, DEVICE_TO_HOST = 2, &
       DEVICE_TO_DEVICE = 3

  !> Copy data between host and device (or device and device)
  interface device_memcpy
     module procedure device_memcpy_r1, device_memcpy_r2, &
          device_memcpy_r3, device_memcpy_r4, device_memcpy_cptr
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

  !> Deassociate a Fortran array from a device pointer
  interface device_deassociate
     module procedure device_deassociate_r1, device_deassociate_r2, &
          device_deassociate_r3, device_deassociate_r4
  end interface device_deassociate
  
  !> Return the device pointer for an associated Fortran array
  interface device_get_ptr
     module procedure device_get_ptr_r1, device_get_ptr_r2, &
          device_get_ptr_r3, device_get_ptr_r4
  end interface device_get_ptr

  !> Synchronize a device or stream
  interface device_sync
     module procedure device_sync_device, device_sync_stream
  end interface device_sync
      
  !> Table of host to device address mappings
  type(htable_cptr_t), private :: device_addrtbl

  private :: device_memcpy_common
  
contains

  subroutine device_init
#if defined(HAVE_HIP) || defined(HAVE_CUDA) || defined(HAVE_OPENCL)
    call device_addrtbl%init(64)

#if defined(HAVE_OPENCL)
    call opencl_init
#endif

#endif    
  end subroutine device_init

  subroutine device_finalize
#if defined(HAVE_HIP) || defined(HAVE_CUDA) || defined(HAVE_OPENCL)
    call device_addrtbl%free()

#if defined(HAVE_OPENCL)
    call opencl_prgm_lib_release
    call opencl_finalize
#endif

#endif
  end subroutine device_finalize

  subroutine device_name(name)
    character(len=*), intent(inout) :: name

#ifdef HAVE_HIP
    call hip_device_name(name)
#elif HAVE_CUDA
    call cuda_device_name(name)
#elif HAVE_OPENCL
    call opencl_device_name(name)
#endif
  end subroutine device_name
  
  !> Allocate memory on the device
  subroutine device_alloc(x_d, s)
    type(c_ptr), intent(inout) :: x_d
    integer(c_size_t) :: s
    integer :: ierr
#ifdef HAVE_HIP
    if (hipMalloc(x_d, s) .ne. hipSuccess) then
       call neko_error('Memory allocation on device failed')
    end if
#elif HAVE_CUDA
    if (cudamalloc(x_d, s) .ne. cudaSuccess) then
       call neko_error('Memory allocation on device failed')
    end if
#elif HAVE_OPENCL
    x_d = clCreateBuffer(glb_ctx, CL_MEM_READ_WRITE, s, C_NULL_PTR, ierr)
    if (ierr .ne. CL_SUCCESS) then
       call neko_error('Memory allocation on device failed')
    end if
#endif
  end subroutine device_alloc

  !> Deallocate memory on the device
  subroutine device_free(x_d)
    type(c_ptr), intent(inout) :: x_d
#ifdef HAVE_HIP
    if (hipfree(x_d) .ne. hipSuccess) then
       call neko_error('Memory deallocation on device failed')
    end if
#elif HAVE_CUDA
    if (cudafree(x_d) .ne. cudaSuccess) then
       call neko_error('Memory deallocation on device failed')
    end if
#elif HAVE_OPENCL
    if (clReleaseMemObject(x_d) .ne. CL_SUCCESS) then
       call neko_error('Memory deallocation on device failed')
    end if
#endif
    x_d = C_NULL_PTR
  end subroutine device_free

  !> Copy data between host and device (rank 1 arrays)
  subroutine device_memcpy_r1(x, x_d, n, dir, sync)
    integer, intent(in) :: n
    class(*), intent(inout), target :: x(:)
    type(c_ptr), intent(inout) :: x_d
    integer, intent(in), value :: dir
    logical, optional :: sync
    type(c_ptr) :: ptr_h
    integer(c_size_t) :: s
    logical :: sync_device

    if (present(sync)) then
       sync_device = sync
    else
       sync_device = .true.
    end if

    select type(x)
    type is (integer)
       s = n * 4
       ptr_h = c_loc(x)
    type is (integer(i8))       
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

    call device_memcpy_common(ptr_h, x_d, s, dir, sync_device)
    
  end subroutine device_memcpy_r1

  !> Copy data between host and device (rank 2 arrays)
  subroutine device_memcpy_r2(x, x_d, n, dir, sync)
    integer, intent(in) :: n
    class(*), intent(inout), target :: x(:,:)
    type(c_ptr), intent(inout) :: x_d
    integer, intent(in), value :: dir
    logical, optional :: sync
    type(c_ptr) :: ptr_h
    integer(c_size_t) :: s
    logical :: sync_device
    
    if (present(sync)) then
       sync_device = sync
    else
       sync_device = .true.
    end if

    select type(x)
    type is (integer)
       s = n * 4
       ptr_h = c_loc(x)
    type is (integer(i8))       
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

    call device_memcpy_common(ptr_h, x_d, s, dir, sync_device)
    
  end subroutine device_memcpy_r2

  !> Copy data between host and device (rank 3 arrays)
  subroutine device_memcpy_r3(x, x_d, n, dir, sync)
    integer, intent(in) :: n
    class(*), intent(inout), target :: x(:,:,:)
    type(c_ptr), intent(inout) :: x_d
    integer, intent(in), value :: dir
    logical, optional :: sync
    type(c_ptr) :: ptr_h
    integer(c_size_t) :: s
    logical :: sync_device

    if (present(sync)) then
       sync_device = sync
    else
       sync_device = .true.
    end if
    
    select type(x)
    type is (integer)
       s = n * 4
       ptr_h = c_loc(x)
    type is (integer(i8))       
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

    call device_memcpy_common(ptr_h, x_d, s, dir, sync_device)
    
  end subroutine device_memcpy_r3

  !> Copy data between host and device (rank 4 arrays)
  subroutine device_memcpy_r4(x, x_d, n, dir, sync)
    integer, intent(in) :: n
    class(*), intent(inout), target :: x(:,:,:,:)
    type(c_ptr), intent(inout) :: x_d
    integer, intent(in), value :: dir
    logical, optional :: sync
    type(c_ptr) :: ptr_h
    integer(c_size_t) :: s    
    logical :: sync_device

    if (present(sync)) then
       sync_device = sync
    else
       sync_device = .true.
    end if
    
    select type(x)
    type is (integer)
       s = n * 4
       ptr_h = c_loc(x)
    type is (integer(i8))       
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

    call device_memcpy_common(ptr_h, x_d, s, dir, sync_device)
    
  end subroutine device_memcpy_r4

  !> Copy data between host and device (or device and device) (c-pointers)
  !! @note For host-device copies @a dst is the host pointer and @a src is the
  !! device pointer (regardless of @a dir)
  subroutine device_memcpy_cptr(dst, src, s, dir, sync)
    type(c_ptr), intent(inout) :: dst
    type(c_ptr), intent(inout) :: src
    integer(c_size_t), intent(in) :: s
    integer, intent(in), value :: dir
    logical, optional :: sync
    logical :: sync_device

    if (present(sync)) then
       sync_device = sync
    else
       sync_device = .true.
    end if

    call device_memcpy_common(dst, src, s, dir, sync_device)
    
  end subroutine device_memcpy_cptr
  
  !> Copy data between host and device
  !! @note For device to device copies, @a ptr_h is assumed
  !! to be the dst device pointer
  subroutine device_memcpy_common(ptr_h, x_d, s, dir, sync_device)
    type(c_ptr), intent(inout) :: ptr_h
    type(c_ptr), intent(inout) :: x_d
    integer(c_size_t), intent(in) :: s
    integer, intent(in), value :: dir
    logical, intent(in) :: sync_device
#ifdef HAVE_HIP
    if (sync_device) then
       if (dir .eq. HOST_TO_DEVICE) then
          if (hipMemcpy(x_d, ptr_h, s, hipMemcpyHostToDevice) &
               .ne. hipSuccess) then
             call neko_error('Device memcpy (host-to-device) failed')
          end if
       else if (dir .eq. DEVICE_TO_HOST) then       
          if (hipMemcpy(ptr_h, x_d, s, hipMemcpyDeviceToHost) &
               .ne. hipSuccess) then
             call neko_error('Device memcpy (device-to-host) failed')
          end if
       else if (dir .eq. DEVICE_TO_DEVICE) then       
          if (hipMemcpy(ptr_h, x_d, s, hipMemcpyDeviceToDevice) &
               .ne. hipSuccess) then
             call neko_error('Device memcpy (device-to-device) failed')
          end if
       else
          call neko_error('Device memcpy failed (invalid direction')
       end if
    else
       if (dir .eq. HOST_TO_DEVICE) then
          if (hipMemcpyAsync(x_d, ptr_h, s, hipMemcpyHostToDevice) &
               .ne. hipSuccess) then
             call neko_error('Device memcpy async (host-to-device) failed')
          end if
       else if (dir .eq. DEVICE_TO_HOST) then       
          if (hipMemcpyAsync(ptr_h, x_d, s, hipMemcpyDeviceToHost) &
               .ne. hipSuccess) then
             call neko_error('Device memcpy async (device-to-host) failed')
          end if
       else if (dir .eq. DEVICE_TO_DEVICE) then       
          if (hipMemcpyAsync(ptr_h, x_d, s, hipMemcpyDeviceToDevice) &
               .ne. hipSuccess) then
             call neko_error('Device memcpy async (device-to-device) failed')
          end if
       else
          call neko_error('Device memcpy failed (invalid direction')
       end if
    end if
#elif HAVE_CUDA
    if (sync_device) then
       if (dir .eq. HOST_TO_DEVICE) then
          if (cudaMemcpy(x_d, ptr_h, s, cudaMemcpyHostToDevice) &
               .ne. cudaSuccess) then
             call neko_error('Device memcpy (host-to-device) failed')
          end if
       else if (dir .eq. DEVICE_TO_HOST) then       
          if (cudaMemcpy(ptr_h, x_d, s, cudaMemcpyDeviceToHost) &
               .ne. cudaSuccess) then
             call neko_error('Device memcpy (device-to-host) failed')
          end if
       else if (dir .eq. DEVICE_TO_DEVICE) then       
          if (cudaMemcpy(ptr_h, x_d, s, cudaMemcpyDeviceToDevice) &
               .ne. cudaSuccess) then
             call neko_error('Device memcpy (device-to-device) failed')
          end if
       else
          call neko_error('Device memcpy failed (invalid direction')
       end if
    else
       if (dir .eq. HOST_TO_DEVICE) then
          if (cudaMemcpyAsync(x_d, ptr_h, s, cudaMemcpyHostToDevice) &
               .ne. cudaSuccess) then
             call neko_error('Device memcpy async (host-to-device) failed')
          end if
       else if (dir .eq. DEVICE_TO_HOST) then       
          if (cudaMemcpyAsync(ptr_h, x_d, s, cudaMemcpyDeviceToHost) &
               .ne. cudaSuccess) then
             call neko_error('Device memcpy async (device-to-host) failed')
          end if
       else if (dir .eq. DEVICE_TO_DEVICE) then       
          if (cudaMemcpyAsync(ptr_h, x_d, s, cudaMemcpyDeviceToDevice) &
               .ne. cudaSuccess) then
             call neko_error('Device memcpy async (device-to-device) failed')
          end if
       else
          call neko_error('Device memcpy failed (invalid direction')
       end if
    end if
#elif HAVE_OPENCL
    if (sync_device) then
       if (dir .eq. HOST_TO_DEVICE) then
          if (clEnqueueWriteBuffer(glb_cmd_queue, x_d, CL_TRUE, 0_i8, s, ptr_h, &
               0, C_NULL_PTR, C_NULL_PTR) .ne. CL_SUCCESS) then
             call neko_error('Device memcpy (host-to-device) failed')
          end if
       else if (dir .eq. DEVICE_TO_HOST) then
          if (clEnqueueReadBuffer(glb_cmd_queue, x_d, CL_TRUE, 0_i8, s, ptr_h, &
               0, C_NULL_PTR, C_NULL_PTR) .ne. CL_SUCCESS) then
             call neko_error('Device memcpy (host-to-device) failed')
          end if
       else if (dir .eq. DEVICE_TO_DEVICE) then
          if (clEnqueueCopyBuffer(glb_cmd_queue, x_d, ptr_h, 0_i8, 0_i8, s, &
               0, C_NULL_PTR, C_NULL_PTR) .ne. CL_SUCCESS) then
             call neko_error('Device memcpy (device-to-device) failed')
          end if
       else
          call neko_error('Device memcpy failed (invalid direction')
       end if
    else
       if (dir .eq. HOST_TO_DEVICE) then
          if (clEnqueueWriteBuffer(glb_cmd_queue, x_d, CL_FALSE, 0_i8, s, ptr_h, &
               0, C_NULL_PTR, C_NULL_PTR) .ne. CL_SUCCESS) then
             call neko_error('Device memcpy (host-to-device) failed')
          end if
       else if (dir .eq. DEVICE_TO_HOST) then
          if (clEnqueueReadBuffer(glb_cmd_queue, x_d, CL_FALSE, 0_i8, s, ptr_h, &
               0, C_NULL_PTR, C_NULL_PTR) .ne. CL_SUCCESS) then
             call neko_error('Device memcpy (host-to-device) failed')
          end if
       else if (dir .eq. DEVICE_TO_DEVICE) then
          if (clEnqueueCopyBuffer(glb_cmd_queue, x_d, ptr_h, 0_i8, 0_i8, s, &
               0, C_NULL_PTR, C_NULL_PTR) .ne. CL_SUCCESS) then
             call neko_error('Device memcpy (device-to-device) failed')
          end if
       else
          call neko_error('Device memcpy failed (invalid direction')
       end if
    end if
#endif
  end subroutine device_memcpy_common

  !> Associate a Fortran rank 1 array to a (allocated) device pointer
  subroutine device_associate_r1(x, x_d)
    class(*), intent(inout), target :: x(:)
    type(c_ptr), intent(inout) :: x_d
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(i8))       
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
  subroutine device_associate_r2(x, x_d)
    class(*), intent(inout), target :: x(:,:)
    type(c_ptr), intent(inout) :: x_d
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(i8))       
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
  subroutine device_associate_r3(x, x_d)
    class(*), intent(inout), target :: x(:,:,:)
    type(c_ptr), intent(inout) :: x_d
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(i8))       
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
  subroutine device_associate_r4(x, x_d)
    class(*), intent(inout), target :: x(:,:,:,:)
    type(c_ptr), intent(inout) :: x_d
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(i8))       
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

  !> Deassociate a Fortran rank 1 array from a device pointer
  subroutine device_deassociate_r1(x)
    class(*), intent(inout), target :: x(:)
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(i8))       
       htbl_ptr_h%ptr = c_loc(x)
    type is (real)
       htbl_ptr_h%ptr = c_loc(x)
    type is (double precision)
       htbl_ptr_h%ptr = c_loc(x)
    class default
       call neko_error('Unknown Fortran type')
    end select

    if (device_addrtbl%get(htbl_ptr_h, htbl_ptr_d) .eq. 0) then
       call device_addrtbl%remove(htbl_ptr_h)
    end if

  end subroutine device_deassociate_r1

  !> Deassociate a Fortran rank 2 array from a device pointer
  subroutine device_deassociate_r2(x)
    class(*), intent(inout), target :: x(:,:)
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(i8))       
       htbl_ptr_h%ptr = c_loc(x)
    type is (real)
       htbl_ptr_h%ptr = c_loc(x)
    type is (double precision)
       htbl_ptr_h%ptr = c_loc(x)
    class default
       call neko_error('Unknown Fortran type')
    end select
    
    if (device_addrtbl%get(htbl_ptr_h, htbl_ptr_d) .eq. 0) then
       call device_addrtbl%remove(htbl_ptr_h)
    end if

  end subroutine device_deassociate_r2

  !> Deassociate a Fortran rank 3 array from a device pointer
  subroutine device_deassociate_r3(x)
    class(*), intent(inout), target :: x(:,:,:)
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(i8))       
       htbl_ptr_h%ptr = c_loc(x)
    type is (real)
       htbl_ptr_h%ptr = c_loc(x)
    type is (double precision)
       htbl_ptr_h%ptr = c_loc(x)
    class default
       call neko_error('Unknown Fortran type')
    end select

    if (device_addrtbl%get(htbl_ptr_h, htbl_ptr_d) .eq. 0) then
       call device_addrtbl%remove(htbl_ptr_h)
    end if

  end subroutine device_deassociate_r3

  !> Deassociate a Fortran rank 4 array from a device pointer
  subroutine device_deassociate_r4(x)
    class(*), intent(inout), target :: x(:,:,:,:)
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(i8))       
       htbl_ptr_h%ptr = c_loc(x)
    type is (real)
       htbl_ptr_h%ptr = c_loc(x)
    type is (double precision)
       htbl_ptr_h%ptr = c_loc(x)
    class default
       call neko_error('Unknown Fortran type')
    end select

    if (device_addrtbl%get(htbl_ptr_h, htbl_ptr_d) .eq. 0) then
       call device_addrtbl%remove(htbl_ptr_h)
    end if

  end subroutine device_deassociate_r4
  
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
    type is (integer(i8))
       s = n * 8
    type is (real)
       s = n * 4
    type is (double precision)
       s = n * 8
    class default
       call neko_error('Unknown Fortran type')
    end select

    call device_alloc(x_d, s)
    call device_associate(x, x_d)

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
    type is (integer(i8))
       s = n * 8
    type is (real)
       s = n * 4
    type is (double precision)
       s = n * 8
    class default
       call neko_error('Unknown Fortran type')
    end select

    call device_alloc(x_d, s)
    call device_associate(x, x_d)

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
    type is (integer(i8))
       s = n * 8
    type is (real)
       s = n * 4
    type is (double precision)
       s = n * 8
    class default
       call neko_error('Unknown Fortran type')
    end select

    call device_alloc(x_d, s)
    call device_associate(x, x_d)

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
    type is (integer(i8))
       s = n * 8
    type is (real)
       s = n * 4
    type is (double precision)
       s = n * 8
    class default
       call neko_error('Unknown Fortran type')
    end select

    call device_alloc(x_d, s)
    call device_associate(x, x_d)

  end subroutine device_map_r4

  !> Check if a Fortran rank 1 array is assoicated with a device pointer
  function device_associated_r1(x) result(assoc)
    class(*), intent(inout), target :: x(:)
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d
    logical :: assoc

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(i8))
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
    type is (integer(i8))
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
    type is (integer(i8))
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
    type is (integer(i8))
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
  function device_get_ptr_r1(x)
    class(*), intent(in), target :: x(:)
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d
    type(c_ptr) :: device_get_ptr_r1

    device_get_ptr_r1 = C_NULL_PTR

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(i8))
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
  function device_get_ptr_r2(x)
    class(*), intent(in), target :: x(:,:)
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d
    type(c_ptr) :: device_get_ptr_r2

    device_get_ptr_r2 = C_NULL_PTR

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(i8))
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
  function device_get_ptr_r3(x)
    class(*), intent(in), target :: x(:,:,:)
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d
    type(c_ptr) :: device_get_ptr_r3

    device_get_ptr_r3 = C_NULL_PTR

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(i8))
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
  function device_get_ptr_r4(x)
    class(*), intent(in), target :: x(:,:,:,:)
    type(h_cptr_t) :: htbl_ptr_h, htbl_ptr_d
    type(c_ptr) :: device_get_ptr_r4

    device_get_ptr_r4 = C_NULL_PTR

    select type(x)
    type is (integer)
       htbl_ptr_h%ptr = c_loc(x)
    type is (integer(i8))
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
  subroutine device_sync_device()
#ifdef HAVE_HIP
    if (hipDeviceSynchronize() .ne. hipSuccess) then
       call neko_error('Error during device sync')
    end if
#elif HAVE_CUDA
    if (cudaDeviceSynchronize() .ne. cudaSuccess) then
       call neko_error('Error during device sync')
    end if
#elif HAVE_OPENCL
    if (clFinish(glb_cmd_queue) .ne. CL_SUCCESS) then
       call neko_error('Error during device sync')
    end if
#endif
  end subroutine device_sync_device

  !> Synchronize a device stream
  subroutine device_sync_stream(stream)
    type(c_ptr), intent(in) :: stream
#ifdef HAVE_HIP
    if (hipStreamSynchronize(stream) .ne. hipSuccess) then
       call neko_error('Error during stream sync')
    end if
#elif HAVE_CUDA
    if (cudaStreamSynchronize(stream) .ne. cudaSuccess) then
       call neko_error('Error during stream sync')
    end if
#elif HAVE_OPENCL
    if (clFinish(stream) .ne. CL_SUCCESS) then
       call neko_error('Error during stream sync')
    end if
#endif
  end subroutine device_sync_stream

  !> Create a device stream/command queue
  subroutine device_stream_create(stream, flags)
    type(c_ptr), intent(inout) :: stream
    integer, optional :: flags
#ifdef HAVE_HIP
    if (present(flags)) then
       if (hipStreamCreateWithFlags(stream, flags) .ne. hipSuccess) then
          call neko_error('Error during stream create (w. flags)')
       end if
    else
       if (hipStreamCreate(stream) .ne. hipSuccess) then
          call neko_error('Error during stream create')
       end if
    end if
#elif HAVE_CUDA
    if (present(flags)) then
       if (cudaStreamCreateWithFlags(stream, flags) .ne. cudaSuccess) then
          call neko_error('Error during stream create (w. flags)')
       end if
    else
       if (cudaStreamCreate(stream) .ne. cudaSuccess) then
          call neko_error('Error during stream create')
       end if
    end if
#elif HAVE_OPENCL
    call neko_error('Not implemented yet')
#endif
  end subroutine device_stream_create

  !> Destroy a device stream/command queue
  subroutine device_stream_destroy(stream)
    type(c_ptr), intent(inout) :: stream
#ifdef HAVE_HIP
    if (hipStreamDestroy(stream) .ne. hipSuccess) then
       call neko_error('Error during stream destroy')
    end if
#elif HAVE_CUDA
    if (cudaStreamDestroy(stream) .ne. cudaSuccess) then
       call neko_error('Error during stream destroy')
    end if
#elif HAVE_OPENCL
    call neko_error('Not implemented yet')
#endif
  end subroutine device_stream_destroy

  !> Synchronize a device stream with an event
  subroutine device_stream_wait_event(stream, event, flags)
    type(c_ptr), intent(in) :: stream
    type(c_ptr), intent(in) :: event
    integer :: flags
#ifdef HAVE_HIP
    if (hipStreamWaitEvent(stream, event, flags) .ne. hipSuccess) then
       call neko_error('Error during stream sync')
    end if
#elif HAVE_CUDA
    if (cudaStreamWaitEvent(stream, event, flags) .ne. cudaSuccess) then
       call neko_error('Error during stream sync')
    end if
#elif HAVE_OPENCL
    call neko_error('Not implemented yet')
#endif
  end subroutine device_stream_wait_event
  
  !> Start device profiling
  subroutine device_profiler_start()
#if HAVE_CUDA
    if (cudaProfilerStart() .ne. cudaSuccess) then
       call neko_error('Error starting profiler')
    end if
#endif
  end subroutine device_profiler_start
  
  !> Stop device profiling
  subroutine device_profiler_stop()
#if HAVE_CUDA
    if (cudaProfilerStop() .ne. cudaSuccess) then
       call neko_error('Error stopping profiler')
    end if
#endif
  end subroutine device_profiler_stop

  !> Create a device event queue
  subroutine device_event_create(event, flags)
    type(c_ptr), intent(inout) :: event
    integer, optional :: flags
#ifdef HAVE_HIP
    if (present(flags)) then
       if (hipEventCreateWithFlags(event, flags) .ne. hipSuccess) then
          call neko_error('Error during event create (w. flags)')
       end if
    else
       if (hipEventCreate(event) .ne. hipSuccess) then
          call neko_error('Error during event create')
       end if
    end if
#elif HAVE_CUDA
    if (present(flags)) then
       if (cudaEventCreateWithFlags(event, flags) .ne. cudaSuccess) then
          call neko_error('Error during event create (w. flags)')
       end if
    else
       if (cudaEventCreate(event) .ne. cudaSuccess) then
          call neko_error('Error during event create')
       end if
    end if
#elif HAVE_OPENCL
    call neko_error('Not implemented yet')
#endif
  end subroutine device_event_create

  !> Destroy a device event
  subroutine device_event_destroy(event)
    type(c_ptr), intent(inout) :: event
#ifdef HAVE_HIP
    if (hipEventDestroy(event) .ne. hipSuccess) then
       call neko_error('Error during stream destroy')
    end if
#elif HAVE_CUDA
    if (cudaEventDestroy(event) .ne. cudaSuccess) then
       call neko_error('Error during stream destroy')
    end if
#elif HAVE_OPENCL
    call neko_error('Not implemented yet')
#endif
  end subroutine device_event_destroy
  
  !> Record a device event
  subroutine device_event_record(event, stream)
    type(c_ptr), intent(in) :: event
    type(c_ptr), intent(in) :: stream
#ifdef HAVE_HIP
    if (hipEventRecord(event, stream) .ne. hipSuccess) then
       call neko_error('Error recording an event')
    end if
#elif HAVE_CUDA
    if (cudaEventRecord(event, stream) .ne. cudaSuccess) then
       call neko_error('Error recording an event')
    end if
#elif HAVE_OPENCL
    call neko_error('Not implemented yet')
#endif
  end subroutine device_event_record

  !> Synchronize an event
  subroutine device_event_sync(event)
    type(c_ptr), intent(in) :: event
#ifdef HAVE_HIP
    if (hipEventSynchronize(event) .ne. hipSuccess) then
       call neko_error('Error during event sync')
    end if
#elif HAVE_CUDA
    if (cudaEventSynchronize(event) .ne. cudaSuccess) then
       call neko_error('Error during event sync')
    end if
#elif HAVE_OPENCL
    call neko_error('Not implemented yet')
#endif
  end subroutine device_event_sync

  
end module device
