!> Wrapper for all matrix-matrix product implementations
module mxm_wrapper
  use mxm_std
  use num_types
  implicit none

contains

  !> Compute matrix-matrix product \f$ C = A \cdot B \f$
  !! for contiguously packed matrices A,B, and C.
  subroutine mxm(a,n1,b,n2,c,n3)
    integer, intent(in) :: n1, n2, n3
    real(kind=dp), intent(inout) :: a(n1, n2)
    real(kind=dp), intent(inout) :: b(n2, n3)
    real(kind=dp), intent(inout) :: c(n1, n3)


! #ifdef XSMM
!       if ((n1*n2*n3)**(1./3) .gt. 6) then
!          call libxsmm_dgemm('N','N',n1,n3,n2,1.0,a,n1,b,n2,0.0,c,n1)
!          goto 111
!       else
!          goto 101
!       endif
! #endif

! #ifdef BLAS_MXM
!       call dgemm('N','N',n1,n3,n2,1.0,a,n1,b,n2,0.0,c,n1)
!       goto 111
! #endif

 101  call mxmf2(a,n1,b,n2,c,n3)

 111  continue
      
  end subroutine mxm
  
  !> Collect matrix-matrix product statistics
  subroutine mxm_test_all(nid,ivb)
    integer :: nid, ivb
    external mxms,mxmur2,mxmur3,mxmd,mxmfb,mxmf3,mxmu4
    external madd,mxm,mxm44,mxmf2
    integer ,parameter :: nn=24
    integer, parameter :: nt=10
    character(len=5) c(3,nt)
    real(kind=dp) :: s(nn,2,nt,3)
    real(kind=dp) :: a(nn,2,nt,3)
    integer :: k

!      call nekgsync

!      do k=1,3   ! 3 tests:  N^2 x N, NxN, NxN^2
!         call mxmtest(s(1,1, 1,k),nn,c(k, 1),mxm44 ,'mxm44',k,ivb)
!         call mxmtest(s(1,1, 2,k),nn,c(k, 2),mxms  ,' std ',k,ivb)
!         call mxmtest(s(1,1, 3,k),nn,c(k, 3),mxmur2,'mxmu2',k,ivb)
!         call mxmtest(s(1,1, 4,k),nn,c(k, 4),mxmur3,'mxmu3',k,ivb)
!         call mxmtest(s(1,1, 5,k),nn,c(k, 5),mxmd  ,'mxmd ',k,ivb)
!         call mxmtest(s(1,1, 6,k),nn,c(k, 6),mxmfb ,'mxmfb',k,ivb)
!         call mxmtest(s(1,1, 7,k),nn,c(k, 7),mxmu4 ,'mxmu4',k,ivb)
!         call mxmtest(s(1,1, 8,k),nn,c(k, 8),mxmf3 ,'mxmf3',k,ivb)
!         if (k.eq.2)  & ! Add works only for NxN case
!         call mxmtest(s(1,1, 9,k),nn,c(k, 9),madd  ,'madd ',k,ivb)
!         call mxmtest(s(1,1,10,k),nn,c(k,10),mxmf2 ,'mxmf2  ',k,ivb)
!      enddo

!      call nekgsync
!      if (nid.eq.0) call mxm_analyze(s,a,nn,c,nt,ivb)
!      call nekgsync

    end subroutine mxm_test_all

  end module mxm_wrapper
