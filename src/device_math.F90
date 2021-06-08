module device_math
  use num_types
  use hip_math
  use comm
  use, intrinsic :: iso_c_binding
  implicit none

contains

  subroutine device_copy(a_d, b_d, n)
    type(c_ptr) :: a_d, b_d
    integer :: n
#ifdef HAVE_HIP
    call hip_copy(a_d, b_d, n)
#endif
  end subroutine device_copy

  subroutine device_rzero(a_d, n)
    type(c_ptr) :: a_d
    integer :: n
#ifdef HAVE_HIP
    call hip_rzero(a_d, n)
#endif
  end subroutine device_rzero

  subroutine device_add2s1(a_d, b_d, c1, n)
    type(c_ptr) :: a_d, b_d
    real(kind=rp) :: c1
    integer :: n
#ifdef HAVE_HIP
    call hip_add2s1(a_d, b_d, c1, n)    
#endif
  end subroutine device_add2s1

  subroutine device_add2s2(a_d, b_d, c1, n)
    type(c_ptr) :: a_d, b_d
    real(kind=rp) :: c1
    integer :: n
#ifdef HAVE_HIP
    call hip_add2s1(a_d, b_d, c1, n)    
#endif
  end subroutine device_add2s2

  function device_glsc3(a_d, b_d, c_d, n) result(res)
    type(c_ptr) :: a_d, b_d, c_d
    integer :: n, ierr
    real(kind=rp) :: res
#ifdef HAVE_HIP
    res = hip_glsc3(a_d, b_d, c_d, n)
#endif

    if (pe_size .gt. 1) then
       call MPI_Allreduce(MPI_IN_PLACE, res, 1, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    end if
  end function device_glsc3
  
end module device_math
