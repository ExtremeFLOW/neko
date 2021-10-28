!> Fortran OpenCL interface
module opencl_intf
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_OPENCL

  !> Global OpenCL command queue
  type(c_ptr) :: glb_cmd_queue = C_NULL_PTR
  
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

#endif
  
end module opencl_intf
