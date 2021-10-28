!> Fortran HIP interface
module hip_intf
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_HIP
  
  !> Enum @a hipError_t
  enum, bind(c)
     enumerator :: hipSuccess = 0
     enumerator :: hipErrorInvalidContext = 1
     enumerator :: hipErrorInvalidKernelFile = 2	
     enumerator :: hipErrorMemoryAllocation = 3
     enumerator :: hipErrorInitializationError = 4	
     enumerator :: hipErrorLaunchFailure = 5
     enumerator :: hipErrorLaunchOutOfResources = 6
     enumerator :: hipErrorInvalidDevice = 7
     enumerator :: hipErrorInvalidValue = 8	
     enumerator :: hipErrorInvalidDevicePointer = 9
     enumerator :: hipErrorInvalidMemcpyDirection = 10	
     enumerator :: hipErrorUnknown = 11
     enumerator :: hipErrorInvalidResourceHandle = 12
     enumerator :: hipErrorNotReady = 13
     enumerator :: hipErrorNoDevice = 14
     enumerator :: hipErrorPeerAccessAlreadyEnabled = 15
     enumerator :: hipErrorPeerAccessNotEnabled = 16
     enumerator :: hipErrorRuntimeMemory = 17
     enumerator :: hipErrorRuntimeOther = 18	
     enumerator :: hipErrorHostMemoryAlreadyRegistered = 19
     enumerator :: hipErrorHostMemoryNotRegistered = 20
     enumerator :: hipErrorMapBufferObjectFailed = 21	
     enumerator :: hipErrorTbd = 22
  end enum
  
  !> Enum @a hipMemcpyKind
  enum, bind(c)
     enumerator :: hipMemcpyHostToHost = 0
     enumerator :: hipMemcpyHostToDevice = 1
     enumerator :: hipMemcpyDeviceToHost = 2
     enumerator :: hipMemcpyDevicetoDevice = 3
     enumerator :: hipMemcpyDefault = 4
  end enum
  
  interface
     integer (c_int) function hipMalloc(ptr_d, s) &
          bind(c, name='hipMalloc')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: ptr_d
       integer(c_size_t), value :: s
     end function hipMalloc
  end interface

  interface
     integer (c_int) function hipFree(ptr_d) &
          bind(c, name='hipFree')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: ptr_d
     end function hipFree
  end interface
  
  interface
     integer (c_int) function hipMemcpy(ptr_dst, ptr_src, s, dir) &
          bind(c, name='hipMemcpy')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: ptr_dst, ptr_src
       integer(c_size_t), value :: s
       integer(c_int), value :: dir
     end function hipMemcpy
  end interface

  interface
     integer (c_int) function hipDeviceSynchronize() &
          bind(c, name='hipDeviceSynchronize')
       use, intrinsic :: iso_c_binding
       implicit none
     end function hipDeviceSynchronize
  end interface

#endif
 
end module hip_intf
