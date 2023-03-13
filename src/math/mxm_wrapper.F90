!> Wrapper for all matrix-matrix product implementations
module mxm_wrapper
  use num_types
  use utils, only : neko_error
#ifdef HAVE_LIBXSMM
  use libxsmm
#endif
  implicit none
  private

  public :: mxm

  interface mxm_blas
     module procedure mxm_blas_sp, mxm_blas_dp, mxm_blas_qp
  end interface mxm_blas
  
  interface mxm_libxsmm
     module procedure mxm_libxsmm_sp, mxm_libxsmm_dp, mxm_libxsmm_qp
  end interface mxm_libxsmm

  private :: mxm_blas_sp, mxm_blas_dp, mxm_blas_qp
  private :: mxm_libxsmm_sp, mxm_libxsmm_dp, mxm_libxsmm_qp

contains

  !> Compute matrix-matrix product \f$ C = A \cdot B \f$
  !! for contiguously packed matrices A,B, and C.
  subroutine mxm(a,n1,b,n2,c,n3)
    integer, intent(in) :: n1, n2, n3
    real(kind=rp), intent(in) :: a(n1, n2)
    real(kind=rp), intent(in) :: b(n2, n3)
    real(kind=rp), intent(inout) :: c(n1, n3)

#ifdef HAVE_LIBXSMM
    call mxm_libxsmm(a,n1,b,n2,c,n3)
#else
    call mxm_blas(a,n1,b,n2,c,n3)
#endif

  end subroutine mxm

  subroutine mxm_blas_sp(a,n1,b,n2,c,n3)
    integer, intent(in) :: n1, n2, n3
    real(kind=sp), intent(in) :: a(n1, n2)
    real(kind=sp), intent(in) :: b(n2, n3)
    real(kind=sp), intent(inout) :: c(n1, n3)

    call sgemm('N','N',n1,n3,n2,1.0,a,n1,b,n2,0.0,c,n1)

  end subroutine mxm_blas_sp
  
  subroutine mxm_blas_dp(a,n1,b,n2,c,n3)
    integer, intent(in) :: n1, n2, n3
    real(kind=dp), intent(in) :: a(n1, n2)
    real(kind=dp), intent(in) :: b(n2, n3)
    real(kind=dp), intent(inout) :: c(n1, n3)

    call dgemm('N','N',n1,n3,n2,1d0,a,n1,b,n2,0d0,c,n1)

  end subroutine mxm_blas_dp

  subroutine mxm_blas_qp(a,n1,b,n2,c,n3)
    integer, intent(in) :: n1, n2, n3
    real(kind=qp), intent(in) :: a(n1, n2)
    real(kind=qp), intent(in) :: b(n2, n3)
    real(kind=qp), intent(inout) :: c(n1, n3)

    call neko_error('Not implemented yet!')

  end subroutine mxm_blas_qp

  subroutine mxm_libxsmm_sp(a,n1,b,n2,c,n3)
    integer, intent(in) :: n1, n2, n3
    real(kind=sp), intent(in) :: a(n1, n2)
    real(kind=sp), intent(in) :: b(n2, n3)
    real(kind=sp), intent(inout) :: c(n1, n3)
#ifdef HAVE_LIBXSMM
    type(libxsmm_smmfunction) :: xmm

    call libxsmm_dispatch(xmm, n1, n3, n2, &
         alpha=1.0, beta=0.0, prefetch=LIBXSMM_PREFETCH)
    if (libxsmm_available(xmm)) then
       call libxsmm_smmcall_abc(xmm, a, b, c)
       return
    end if
#endif       
  end subroutine mxm_libxsmm_sp
  
  subroutine mxm_libxsmm_dp(a,n1,b,n2,c,n3)
    integer, intent(in) :: n1, n2, n3
    real(kind=dp), intent(in) :: a(n1, n2)
    real(kind=dp), intent(in) :: b(n2, n3)
    real(kind=dp), intent(inout) :: c(n1, n3)
#ifdef HAVE_LIBXSMM
    type(libxsmm_dmmfunction) :: xmm

    call libxsmm_dispatch(xmm, n1, n3, n2, &
         alpha=1d0, beta=0d0, prefetch=LIBXSMM_PREFETCH)
    if (libxsmm_available(xmm)) then
       call libxsmm_dmmcall_abc(xmm, a, b, c)
       return
    end if
#endif
  end subroutine mxm_libxsmm_dp

  subroutine mxm_libxsmm_qp(a,n1,b,n2,c,n3)
    integer, intent(in) :: n1, n2, n3
    real(kind=qp), intent(in) :: a(n1, n2)
    real(kind=qp), intent(in) :: b(n2, n3)
    real(kind=qp), intent(inout) :: c(n1, n3)

    call neko_error('Not implemented yet!')

  end subroutine mxm_libxsmm_qp
  
end module mxm_wrapper
