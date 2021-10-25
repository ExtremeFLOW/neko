!> Fortran CUDA interface
module cuda_intf
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_CUDA
  
  integer, parameter :: CUDA_SUCCESS = 0, CUDA_ERROR = 1

  interface
     integer (c_int) function cudamalloc(ptr_d, s) &
          bind(c, name='cudaMalloc_wrapper')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: ptr_d
       integer(c_size_t) :: s
     end function cudamalloc
  end interface

  interface
     integer (c_int) function cudafree(ptr_d) &
          bind(c, name='cudaFree_wrapper')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: ptr_d
     end function cudafree
  end interface
  
  interface
     integer (c_int) function cudamemcpy(ptr, ptr_d, s, dir) &
          bind(c, name='cudaMemcpy_wrapper')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: ptr, ptr_d
       integer(c_size_t) :: s
       integer(c_int) :: dir
     end function cudaMemcpy
  end interface

  interface
     integer (c_int) function cudadevicesynchronize() &
          bind(c, name='cudaDeviceSynchronize_wrapper')
       use, intrinsic :: iso_c_binding
       implicit none
     end function cudadevicesynchronize
  end interface

#endif
 
end module cuda_intf
