!> Wrapper for all matrix-matrix product implementations
module mxm_wrapper
  use mxm_std
  use num_types
#ifdef HAVE_LIBXSMM
  use libxsmm, libxsmm_mmcall => libxsmm_dmmcall_abc
#endif
  implicit none

contains

  !> Compute matrix-matrix product \f$ C = A \cdot B \f$
  !! for contiguously packed matrices A,B, and C.
  subroutine mxm(a,n1,b,n2,c,n3)
    integer, intent(in) :: n1, n2, n3
    real(kind=rp), intent(in) :: a(n1, n2)
    real(kind=rp), intent(in) :: b(n2, n3)
    real(kind=rp), intent(inout) :: c(n1, n3)
#ifdef HAVE_LIBXSMM
    type(libxsmm_dmmfunction) :: xmm
    
    call libxsmm_dispatch(xmm, n1, n3, n2, &
         alpha=1d0, beta=0d0, prefetch=LIBXSMM_PREFETCH)
    if (libxsmm_available(xmm)) then
       call libxsmm_mmcall(xmm, a, b, c)
       return
    end if
#endif

    ! #ifdef BLAS_MXM
!       call dgemm('N','N',n1,n3,n2,1.0,a,n1,b,n2,0.0,c,n1)
!       goto 111
! #endif

    call mxmf2(a,n1,b,n2,c,n3)
    
      
  end subroutine mxm
  
end module mxm_wrapper
