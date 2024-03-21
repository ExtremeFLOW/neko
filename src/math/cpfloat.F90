!Provides fortran wrapper for functionalities in cpfloat
!Introcudes the modules pcg_f (for random numbers) and cpfloat_f for rounding
! cpfloat is described in:
! Massimiliano Fasi and Mantas Mikaitis. CPFloat: A C Library for Simulating Low-precision Arithmetic.
! ACM Trans. Math. Softw. 49, 2, Article 18 (June 2023), 32 pages.
! https://doi.org/10.1145/3585515

module pcg_f
  use, intrinsic :: iso_fortran_env
  use, intrinsic :: iso_c_binding
  implicit none

  !same as pcg32_random_t
  type, public, bind(c) :: pcg_state_setseq_64
     integer(c_int64_t) :: state
     integer(c_int64_t) :: inc
  end type

  !same as pcg64_random_t
  type, public, bind(c) :: pcg_state_setseq_128
     integer(c_int128_t) :: state
     integer(c_int128_t) :: inc
  end type
#ifdef HAVE_CPFLOAT
  interface
    integer(c_int32_t) function pcg32_random()&
        bind(c, name="pcg32_random")
      use, intrinsic :: iso_c_binding
    end function pcg32_random
  end interface
  interface
    subroutine pcg32_srandom(initstate, initseq)& 
        bind(c, name="pcg32_srandom")
      use, intrinsic :: iso_c_binding
      integer(c_int64_t), value :: initstate
      integer(c_int64_t), value :: initseq
    end subroutine pcg32_srandom
  end interface

  interface
    integer(c_int64_t) function pcg64_random()&
        bind(c, name="pcg64_random")
      use, intrinsic :: iso_c_binding
    end function pcg64_random
  end interface
  interface
    subroutine pcg64_srandom(initstate, initseq)&
        bind(c, name="pcg64_srandom")
      use, intrinsic :: iso_c_binding
      integer(c_int128_t), value :: initstate
      integer(c_int128_t), value :: initseq
    end subroutine pcg64_srandom
  end interface
#endif
end module pcg_f


module cpfloat_f
  use, intrinsic :: iso_fortran_env
  use, intrinsic :: iso_c_binding
  use pcg_f
  use utils, only : neko_error
  implicit none

  type, public, bind(c) :: optstruct
      character(c_char) :: format(15)
      integer(c_int) :: precision
      integer(c_int) :: emax
      integer(c_int) :: subnormal
      integer(c_int) :: explim
      integer(c_int) :: round
      integer(c_int) :: flip 
      real(c_double) :: p 
      type(pcg_state_setseq_64) :: bitseed
      type(pcg_state_setseq_128) :: randseedf
      type(pcg_state_setseq_128) :: randseed
   end type optstruct

  interface
     integer(c_int) function cpfloat_intf(y,x,n,fpopts)&
          bind(c, name='cpfloat')
       use, intrinsic :: iso_c_binding
       import optstruct
       implicit none
       integer(c_size_t), value :: n
       real(kind=c_double) :: y(n), x(n)
       type(optstruct) :: fpopts
     end function cpfloat_intf
  end interface

#ifdef HAVE_CPFLOAT
  interface
     integer(c_int) function cpfloat_validate_optstruct_intf(fpopts)&
          bind(c, name='cpfloat_validate_optstruct')
       use, intrinsic :: iso_c_binding
       import optstruct
       implicit none
       type(optstruct) :: fpopts
     end function cpfloat_validate_optstruct_intf
  end interface

  interface
      type(optstruct) function cpfloat_init_optstruct_intf()&
          bind(c, name='init_optstruct')
       use, intrinsic :: iso_c_binding
       import optstruct
       implicit none
     end function cpfloat_init_optstruct_intf
  end interface
#endif
  contains

    function cpfloat_validate_optstruct(fpopts) result(ierr)
       type(optstruct) :: fpopts
       integer :: ierr
       ierr = cpfloat_validate_optstruct_intf(fpopts)
    end function cpfloat_validate_optstruct

    function cpfloat_init_optstruct() result(fpopts)
      use, intrinsic :: iso_c_binding
      implicit none
      type(optstruct) :: fpopts
#ifdef HAVE_CPFLOAT
      fpopts = cpfloat_init_optstruct_intf()
#else  
       call neko_error('Need to compile with cpfloat')
#endif
    end function cpfloat_init_optstruct

    function cpfloat(y,x,n,fpopts) result(ierr)
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_size_t) :: n
       real(c_double) :: y(n), x(n)
       type(optstruct) :: fpopts
       integer(c_int) :: ierr
#ifdef HAVE_CPFLOAT
       ierr = cpfloat_intf(y,x,n,fpopts)
#else  
       call neko_error('Need to compile with cpfloat')
#endif
       
     end function cpfloat
end module cpfloat_f


