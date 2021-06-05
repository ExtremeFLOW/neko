module hip_math
  implicit none

#ifdef HAVE_HIP

  interface
     subroutine hip_copy(a_d, b_d, n) &
          bind(c, name='hip_copy')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine hip_copy
  end interface

  interface
     subroutine hip_rzero(a_d, n) &
          bind(c, name='hip_rzero')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: a_d
       integer(c_int) :: n
     end subroutine hip_rzero
  end interface

  interface
     subroutine hip_add2s1(a_d, b_d, c1, n) &
          bind(c, name='hip_add2s1')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d
       real(c_double) :: c1
       integer(c_int) :: n
     end subroutine hip_add2s1
  end interface

  interface
     subroutine hip_add2s2(a_d, b_d, c1, n) &
          bind(c, name='hip_add2s2')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d
       real(c_double) :: c1
       integer(c_int) :: n
     end subroutine hip_add2s2
  end interface

  interface
     real(c_double) function hip_glsc3(a_d, b_d, c_d, n) &
          bind(c, name='hip_glsc3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end function hip_glsc3
  end interface



#endif  
end module hip_math
