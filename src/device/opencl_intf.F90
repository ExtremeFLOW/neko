!> Fortran OpenCL interface
module opencl_intf
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_OPENCL

  !> Global OpenCL command queue
  type(c_ptr) :: glb_cmd_queue = C_NULL_PTR

  !> Global OpenCL context
  type(c_ptr) :: glb_ctx = C_NULL_PTR
  
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

  enum, bind(c)
     enumerator :: CL_CONTEXT_PLATFORM = int(Z'1084')
  end enum

  !> Device types
  integer(c_int64_t), parameter :: CL_DEVICE_TYPE_DEFAULT = 1
  integer(c_int64_t), parameter :: CL_DEVICE_TYPE_CPU = 2
  integer(c_int64_t), parameter :: CL_DEVICE_TYPE_GPU = 4
  integer(c_int64_t), parameter :: CL_DEVICE_TYPE_ACCELERATOR = 8
  integer(c_int64_t), parameter :: CL_DEVICE_TYPE_CUSTOM = 16
  integer(c_int64_t), parameter :: CL_DEVICE_TYPE_ALL = int(Z'FFFFFFFF', 8)
    
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
     integer (c_int) function clReleaseMemObject(ptr_d) &
          bind(c, name='clReleaseMemObject')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: ptr_d
     end function clReleaseMemObject
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
    type(c_ptr), target :: platform_id, device_id
    integer(c_int) :: num_platforms, num_devices, ierr
    integer(c_intptr_t) :: ctx_prop(3)
    integer(c_int64_t), parameter :: queue_props = 0
    integer :: i

    if (clGetPlatformIDs(1, c_loc(platform_id), &
         num_platforms) .ne. CL_SUCCESS) then
       call neko_error('Faield to get a platform id')
    end if

    if (clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, &
         c_loc(device_id), num_devices) .ne. CL_SUCCESS) then
       call neko_error('Failed to get a device id')
    end if

    if (c_associated(glb_ctx)) then
       if (clReleaseContext(glb_ctx) .ne. CL_SUCCESS) then
          call neko_error('Failed to release context')
       end if
    end if

    glb_ctx = clCreateContext(C_NULL_PTR, num_devices, c_loc(device_id), &
         C_NULL_FUNPTR, C_NULL_PTR, ierr)

    if (ierr .ne. CL_SUCCESS) then
       call neko_error('Failed to create an OpenCL context')
    end if

    if (c_associated(glb_cmd_queue)) then
       if (clReleaseCommandQueue(glb_cmd_queue) .ne. CL_SUCCESS) then
          call neko_error('Faield to release command queue')
       end if
    end if
    
    glb_cmd_queue = clCreateCommandQueue(glb_ctx, device_id, queue_props, ierr)

    if (ierr .ne. CL_SUCCESS) then
       call neko_error('Failed to create a command queue')
    end if

  end subroutine opencl_init

  subroutine opencl_finalize

    if (c_associated(glb_ctx)) then
       if (clReleaseContext(glb_ctx) .ne. CL_SUCCESS) then
          call neko_error('Failed to release context')
       end if
    end if

    if (c_associated(glb_cmd_queue)) then
       if (clReleaseCommandQueue(glb_cmd_queue) .ne. CL_SUCCESS) then
          call neko_error('Faield to release command queue')
       end if
    end if
    
  end subroutine opencl_finalize
  
#endif
  
end module opencl_intf
