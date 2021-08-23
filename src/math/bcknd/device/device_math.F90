module device_math
  use num_types
  use utils
  use comm
  use, intrinsic :: iso_c_binding
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
     subroutine hip_invcol2(a_d, b_d, n) &
          bind(c, name='hip_invcol2')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine hip_invcol2
  end interface
  
  interface
     subroutine hip_col2(a_d, b_d, n) &
          bind(c, name='hip_col2')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end subroutine hip_col2
  end interface
  
  interface
     subroutine hip_col3(a_d, b_d, c_d, n) &
          bind(c, name='hip_col3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine hip_col3
  end interface
  
  interface
     subroutine hip_sub3(a_d, b_d, c_d, n) &
          bind(c, name='hip_sub3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine hip_sub3
  end interface

  interface
     subroutine hip_addcol3(a_d, b_d, c_d, n) &
          bind(c, name='hip_addcol3')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: a_d, b_d, c_d
       integer(c_int) :: n
     end subroutine hip_addcol3
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
  
contains

  subroutine device_copy(a_d, b_d, n)
    type(c_ptr) :: a_d, b_d
    integer :: n
#ifdef HAVE_HIP
    call hip_copy(a_d, b_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_copy

  subroutine device_rzero(a_d, n)
    type(c_ptr) :: a_d
    integer :: n
#ifdef HAVE_HIP
    call hip_rzero(a_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_rzero

  subroutine device_add2s1(a_d, b_d, c1, n)
    type(c_ptr) :: a_d, b_d
    real(kind=rp) :: c1
    integer :: n
#ifdef HAVE_HIP
    call hip_add2s1(a_d, b_d, c1, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_add2s1

  subroutine device_add2s2(a_d, b_d, c1, n)
    type(c_ptr) :: a_d, b_d
    real(kind=rp) :: c1
    integer :: n
#ifdef HAVE_HIP
    call hip_add2s2(a_d, b_d, c1, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_add2s2

  subroutine device_invcol2(a_d, b_d, n)
    type(c_ptr) :: a_d, b_d
    integer :: n
#ifdef HAVE_HIP
    call hip_invcol2(a_d, b_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_invcol2

  subroutine device_col2(a_d, b_d, n)
    type(c_ptr) :: a_d, b_d
    integer :: n
#ifdef HAVE_HIP
    call hip_col2(a_d, b_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_col2
  
  subroutine device_col3(a_d, b_d, c_d, n)
    type(c_ptr) :: a_d, b_d, c_d
    integer :: n
#ifdef HAVE_HIP
    call hip_col3(a_d, b_d, c_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_col3
  
  subroutine device_sub3(a_d, b_d, c_d, n)
    type(c_ptr) :: a_d, b_d, c_d
    integer :: n
#ifdef HAVE_HIP
    call hip_sub3(a_d, b_d, c_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_sub3

  subroutine device_addcol3(a_d, b_d, c_d, n)
    type(c_ptr) :: a_d, b_d, c_d
    integer :: n
#ifdef HAVE_HIP
    call hip_addcol3(a_d, b_d, c_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_addcol3

  function device_glsc3(a_d, b_d, c_d, n) result(res)
    type(c_ptr) :: a_d, b_d, c_d
    integer :: n, ierr
    real(kind=rp) :: res
#ifdef HAVE_HIP
    res = hip_glsc3(a_d, b_d, c_d, n)
#else
    call neko_error('No device backend configured')
#endif

    if (pe_size .gt. 1) then
       call MPI_Allreduce(MPI_IN_PLACE, res, 1, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    end if
  end function device_glsc3
  
end module device_math
