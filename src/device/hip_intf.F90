!> Fortran HIP interface
module hip_intf
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_HIP
  
  integer, parameter :: HIP_SUCCESS = 0, HIP_ERROR = 1

  interface
     integer (c_int) function hipmalloc(ptr_d, s) &
          bind(c, name='hipMalloc_wrapper')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: ptr_d
       integer(c_size_t) :: s
     end function hipmalloc
  end interface

  interface
     integer (c_int) function hipfree(ptr_d) &
          bind(c, name='hipFree_wrapper')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: ptr_d
     end function hipfree
  end interface
  
  interface
     integer (c_int) function hipmemcpy(ptr, ptr_d, s, dir) &
          bind(c, name='hipMemcpy_wrapper')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: ptr, ptr_d
       integer(c_size_t) :: s
       integer(c_int) :: dir
     end function hipMemcpy
  end interface

  interface
     integer (c_int) function hipdevicesynchronize() &
          bind(c, name='hipDeviceSynchronize_wrapper')
       use, intrinsic :: iso_c_binding
       implicit none
     end function hipdevicesynchronize
  end interface

#endif
 
end module hip_intf
