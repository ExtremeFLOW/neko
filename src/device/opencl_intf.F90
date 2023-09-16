! Copyright (c) 2021-2023, The Neko Authors
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
!> Fortran OpenCL interface
module opencl_intf
  use num_types, only : i8
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_loc, c_size_t, &
                                          c_char, c_intptr_t, c_int64_t, &
                                          C_NULL_PTR, C_NULL_FUNPTR, &
                                          c_associated, c_size_t
  implicit none
  
#ifdef HAVE_OPENCL

  !> Global OpenCL command queue
  type(c_ptr), bind(c) :: glb_cmd_queue = C_NULL_PTR

  !> Aux OpenCL command queue
  type(c_ptr), bind(c) :: aux_cmd_queue = C_NULL_PTR

  !> Global OpenCL context
  type(c_ptr), bind(c) :: glb_ctx = C_NULL_PTR

  !> Global OpenCL device_id
  type(c_ptr), bind(c), target :: glb_device_id = C_NULL_PTR
  
  !> Enum Error Codes
  enum, bind(c)
     enumerator :: CL_SUCCESS = 0
     enumerator :: CL_DEVICE_NOT_FOUND = -1
     enumerator :: CL_DEVICE_NOT_AVAILABLE = -2
     enumerator :: CL_COMPILER_NOT_AVAILABLE = -3
     enumerator :: CL_MEM_OBJECT_ALLOCATION_FAILURE = -4
     enumerator :: CL_OUT_OF_RESOURCES = -5
     enumerator :: CL_OUT_OF_HOST_MEMORY = -6
     enumerator :: CL_PROFILING_INFO_NOT_AVAILABLE = -7
     enumerator :: CL_MEM_COPY_OVERLAP = -8
     enumerator :: CL_IMAGE_FORMAT_MISMATCH = -9
     enumerator :: CL_IMAGE_FORMAT_NOT_SUPPORTED = -10
     enumerator :: CL_BUILD_PROGRAM_FAILURE = -11
     enumerator :: CL_MAP_FAILURE = -12
  end enum

  !> Enum mem flags
  enum, bind(c)
     enumerator :: CL_MEM_READ_WRITE = 1
     enumerator :: CL_MEM_WRITE_ONLY = 2
     enumerator :: CL_MEM_READ_ONLY = 4
     enumerator :: CL_MEM_USE_HOST_PTR = 8
     enumerator :: CL_MEM_ALLOC_HOST_PTR = 16
     enumerator :: CL_MEM_HOST_WRITE_ONLY = 128
     enumerator :: CL_MEM_HOST_READ_ONLY = 256
     enumerator :: CL_MEM_HOST_NO_ACCESS = 512
  end enum

  !> Enum event flags
  enum, bind(c)
     enumerator :: CL_COMPLETE = 0
     enumerator :: CL_RUNNING = 1
     enumerator :: CL_SUBMITTED = 2
     enumerator :: CL_QUEUED = 3
  end enum

  !> Enum boolean
  enum, bind(c)
     enumerator :: CL_FALSE = 0
     enumerator :: CL_TRUE = 1
  end enum

  enum, bind(c)
     enumerator :: CL_CONTEXT_PLATFORM = int(Z'1084')
  end enum

  !> Enum device info
  enum, bind(c)
     enumerator :: CL_DEVICE_NAME = 4139
  end enum

  !> Device types
  integer(c_int64_t), parameter :: CL_DEVICE_TYPE_DEFAULT = 1
  integer(c_int64_t), parameter :: CL_DEVICE_TYPE_CPU = 2
  integer(c_int64_t), parameter :: CL_DEVICE_TYPE_GPU = 4
  integer(c_int64_t), parameter :: CL_DEVICE_TYPE_ACCELERATOR = 8
  integer(c_int64_t), parameter :: CL_DEVICE_TYPE_CUSTOM = 16
  integer(c_int64_t), parameter :: CL_DEVICE_TYPE_ALL = int(Z'FFFFFFFF', i8)
    
  interface
     integer (c_int) function clGetPlatformIDs(num_entries, &
          platforms, num_platforms) bind(c, name='clGetPlatformIDs')       
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: num_entries
       type(c_ptr), value :: platforms
       integer(c_int) :: num_platforms
     end function clGetPlatformIDs
  end interface

  interface
     integer (c_int) function clGetDeviceIDs(platform, device_type, &
          num_entries, devices, num_devices) bind(c, name='clGetDeviceIDs')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: platform
       integer(c_int64_t), value :: device_type
       integer(c_int), value :: num_entries
       type(c_ptr), value :: devices
       integer(c_int) :: num_devices
     end function clGetDeviceIDs
  end interface

  interface
     type (c_ptr) function clCreateContext(properties, num_devices, &
          devices, pfn_notify, user_data, ierr) bind(c, name='clCreateContext')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: properties
       integer(c_int), value :: num_devices
       type(c_ptr), value :: devices
       type(c_funptr), value :: pfn_notify
       type(c_ptr), value :: user_data
       integer(c_int) :: ierr
     end function clCreateContext
  end interface

  interface
     type(c_ptr) function clCreateCommandQueue(context, device, &
          properties, ierr)  bind(c, name='clCreateCommandQueue')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: context
       type(c_ptr), value :: device
       integer(c_int64_t), value :: properties
       integer(c_int) :: ierr
     end function clCreateCommandQueue
  end interface

  interface
     type(c_ptr) function clCreateBuffer(context, flags, size,  host_ptr, ierr) &
          bind(c, name='clCreateBuffer')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: context
       integer(c_int), value :: flags
       integer(c_size_t), value :: size
       type(c_ptr), value :: host_ptr
       integer(c_int) :: ierr
     end function clCreateBuffer
  end interface

  interface
     type(c_ptr) function clCreateUserEvent(context, ierr) &
          bind(c, name='clCreateUserEvent')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: context
       integer(c_int) :: ierr
     end function clCreateUserEvent
  end interface

  interface
     integer (c_int) function clEnqueueReadBuffer(queue, buffer, &
          blocking_read, offset, size, ptr, num_events_in_wait_list, &
          event_wait_list, event) bind(c, name='clEnqueueReadBuffer')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: queue
       type(c_ptr), value :: buffer
       integer(c_int), value :: blocking_read
       integer(c_size_t), value :: offset
       integer(c_size_t), value :: size
       type(c_ptr), value :: ptr
       integer(c_int), value :: num_events_in_wait_list
       type(c_ptr), value :: event_wait_list
       type(c_ptr), value :: event
     end function clEnqueueReadBuffer
  end interface

  interface
     integer (c_int) function clEnqueueWriteBuffer(queue, buffer, &
          blocking_write, offset, size, ptr, num_events_in_wait_list, &
          event_wait_list, event) bind(c, name='clEnqueueWriteBuffer')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: queue
       type(c_ptr), value :: buffer
       integer(c_int), value :: blocking_write
       integer(c_size_t), value :: offset
       integer(c_size_t), value :: size
       type(c_ptr), value :: ptr
       integer(c_int), value :: num_events_in_wait_list
       type(c_ptr), value :: event_wait_list
       type(c_ptr), value :: event
     end function clEnqueueWriteBuffer
  end interface

  interface
     integer (c_int) function clEnqueueCopyBuffer(queue, src_buffer, &
          dst_buffer, src_offset, dst_offset, size, num_events_in_wait_list, &
          event_wait_list, event) bind(c, name='clEnqueueCopyBuffer')
       use,intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: queue
       type(c_ptr), value :: src_buffer
       type(c_ptr), value :: dst_buffer
       integer(c_size_t), value :: src_offset
       integer(c_size_t), value :: dst_offset
       integer(c_size_t), value :: size
       integer(c_int), value :: num_events_in_wait_list
       type(c_ptr), value :: event_wait_list
       type(c_ptr), value :: event
     end function clEnqueueCopyBuffer
  end interface

  interface
     integer (c_int) function clEnqueueMarker(cmd_queue, event) &
          bind(c, name='clEnqueueMarker') 
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cmd_queue
       type(c_ptr), value :: event
     end function clEnqueueMarker
  end interface

  interface
     integer (c_int) function clEnqueueWaitForEvents(queue, &
          num_events, event_list) bind(c, name='clEnqueueWaitForEvents')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: queue
       integer(c_int), value :: num_events
       type(c_ptr), value :: event_list
     end function clEnqueueWaitForEvents
  end interface

  interface
     integer (c_int) function clWaitForEvents(num_events, event_list) &
          bind(c, name='clWaitForEvents')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: num_events
       type(c_ptr), value :: event_list
     end function clWaitForEvents
  end interface

  interface
     integer (c_int) function clSetUserEventStatus(event, status) &
          bind(c, name='clSetUserEventStatus')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: event
       integer(c_int), value :: status
     end function clSetUserEventStatus
  end interface
       
  interface
     integer (c_int) function clGetDeviceInfo(device, param_name, &
          param_value_size, param_value, param_value_size_ret) &
          bind(c, name='clGetDeviceInfo')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: device
       integer(c_int), value :: param_name
       integer(c_size_t), value :: param_value_size
       type(c_ptr), value :: param_value
       type(c_ptr), value :: param_value_size_ret
     end function clGetDeviceInfo
  end interface
       
  interface
     integer (c_int) function clReleaseContext(context) &
          bind(c, name='clReleaseContext')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: context
     end function clReleaseContext
  end interface
  
  interface
     integer (c_int) function clReleaseCommandQueue(queue) &
          bind(c, name='clReleaseCommandQueue')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: queue
     end function clReleaseCommandQueue
  end interface
  
  interface
     integer (c_int) function clReleaseDevice(device) &
          bind(c, name='clReleaseDevice')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: device
     end function clReleaseDevice
  end interface

  interface
     integer (c_int) function clReleaseProgram(prgm) &
          bind(c, name='clReleaseProgram')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: prgm
     end function clReleaseProgram
  end interface
  
  interface
     integer (c_int) function clReleaseMemObject(ptr_d) &
          bind(c, name='clReleaseMemObject')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: ptr_d
     end function clReleaseMemObject
  end interface

  interface
     integer (c_int) function clReleaseEvent(event) &
          bind(c, name='clReleaseEvent')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: event
     end function clReleaseEvent
  end interface

  interface
     integer (c_int) function clFlush(cmd_queue) &
          bind(c, name='clFlush')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cmd_queue
     end function clFlush
  end interface

    interface
     integer (c_int) function clFinish(cmd_queue) &
          bind(c, name='clFinish')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cmd_queue
     end function clFinish
  end interface
  
contains

  subroutine opencl_init
    type(c_ptr), target :: platform_id
    integer(c_int) :: num_platforms, num_devices, ierr
    integer(c_intptr_t) :: ctx_prop(3)
    integer(c_int64_t), parameter :: queue_props = 0
    integer :: i

    if (clGetPlatformIDs(1, c_loc(platform_id), &
         num_platforms) .ne. CL_SUCCESS) then
       call neko_error('Failed to get a platform id')
    end if

    if (clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, &
         c_loc(glb_device_id), num_devices) .ne. CL_SUCCESS) then
       call neko_error('Failed to get a device id')
    end if

    if (c_associated(glb_ctx)) then
       if (clReleaseContext(glb_ctx) .ne. CL_SUCCESS) then
          call neko_error('Failed to release context')
       end if
    end if

    glb_ctx = clCreateContext(C_NULL_PTR, num_devices, c_loc(glb_device_id), &
         C_NULL_FUNPTR, C_NULL_PTR, ierr)

    if (ierr .ne. CL_SUCCESS) then
       call neko_error('Failed to create an OpenCL context')
    end if

    if (c_associated(glb_cmd_queue)) then
       if (clReleaseCommandQueue(glb_cmd_queue) .ne. CL_SUCCESS) then
          call neko_error('Faield to release command queue')
       end if
    end if
    
    glb_cmd_queue = clCreateCommandQueue(glb_ctx, glb_device_id, queue_props, ierr)

    if (ierr .ne. CL_SUCCESS) then
       call neko_error('Failed to create a command queue')
    end if

  end subroutine opencl_init

  subroutine opencl_finalize

    if (c_associated(glb_ctx)) then
       if (clReleaseContext(glb_ctx) .ne. CL_SUCCESS) then
          call neko_error('Failed to release context')
       end if
       glb_ctx = C_NULL_PTR
    end if

    if (c_associated(glb_cmd_queue)) then
       if (clReleaseCommandQueue(glb_cmd_queue) .ne. CL_SUCCESS) then
          call neko_error('Faield to release command queue')
       end if
       glb_cmd_queue = C_NULL_PTR
    end if

    if (c_associated(glb_device_id)) then
       if (clReleaseDevice(glb_device_id) .ne. CL_SUCCESS) then
          call neko_error('Faield to release device')
       end if
    end if
    
  end subroutine opencl_finalize

  subroutine opencl_device_name(name)
    character(len=*), intent(inout) :: name
    character(kind=c_char, len=1024), target :: c_name
    integer(c_size_t), target :: name_len

    if (clGetDeviceInfo(glb_device_id, CL_DEVICE_NAME, int(1024, i8), &
         c_loc(c_name), c_loc(name_len)) .ne. CL_SUCCESS) then
       call neko_error('Failed to query device')
    end if
        
    name(1:name_len) = c_name(1:name_len)
    
  end subroutine opencl_device_name
  
#endif
  
end module opencl_intf
