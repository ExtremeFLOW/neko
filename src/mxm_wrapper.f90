!> Wrapper for all matrix-matrix product implementations
module mxm_wrapper
  use num_types
  implicit none

contains

  !> Compute matrix-matrix product \f$ C = A \cdot B \f$
  !! for contiguously packed matrices A,B, and C.
  subroutine mxm(a,n1,b,n2,c,n3)
    real(kind=dp), intent(inout) :: a(n1, n2)
    real(kind=dp), intent(inout) :: b(n2, n3)
    real(kind=dp), intent(inout) :: c(n1, n3)
    integer, intent(inout) :: n1, n2, n3

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
  end module mxm_wrapper
