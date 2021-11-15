!> Fortran CUDA interface
module cuda_intf
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_CUDA

  !> Enum @a cudaError
  enum, bind(c)
     enumerator :: cudaSuccess = 0
     enumerator :: cudaErrorInvalidValue = 1
     enumerator :: cudaErrorMemoryAllocation = 2
     enumerator :: cudaErrorInitializationError = 3
  end enum

  !> Enum @a cudaMemcpyKind
  enum, bind(c)
     enumerator :: cudaMemcpyHostToHost = 0
     enumerator :: cudaMemcpyHostToDevice = 1
     enumerator :: cudaMemcpyDeviceToHost = 2
     enumerator :: cudaMemcpyDevicetoDevice = 3
     enumerator :: cudaMemcpyDefault = 4
  end enum

  interface
     integer (c_int) function cudaMalloc(ptr_d, s) &
          bind(c, name='cudaMalloc')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: ptr_d
       integer(c_size_t), value :: s
     end function cudaMalloc
  end interface

  interface
     integer (c_int) function cudaFree(ptr_d) &
          bind(c, name='cudaFree')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: ptr_d
     end function cudaFree
  end interface
  
  interface
     integer (c_int) function cudaMemcpy(ptr_dst, ptr_src, s, dir) &
          bind(c, name='cudaMemcpy')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: ptr_dst, ptr_src
       integer(c_size_t), value :: s
       integer(c_int), value :: dir
     end function cudaMemcpy
  end interface

  interface
     integer (c_int) function cudaDeviceSynchronize() &
          bind(c, name='cudaDeviceSynchronize')
       use, intrinsic :: iso_c_binding
       implicit none
     end function cudaDeviceSynchronize
  end interface

#endif
 
end module cuda_intf
