!> Tensor operations
module tensor
  use num_types
  use mxm_wrapper
  implicit none

  interface transpose
     module procedure trsp, trsp1
  end interface transpose

contains

  !> Tensor product \f$ v =(C \otimes B \otimes A) u \f$
  subroutine tensr3(v, nv, u, nu, A, Bt, Ct, w)
    real(kind=dp), intent(inout) :: v(nv, nv, nv)
    real(kind=dp), intent(inout) :: u(nv, nv, nv)
    real(kind=dp), intent(inout) :: w(nu*nu*nv)
    real(kind=dp), intent(inout) :: A(nv, nu)
    real(kind=dp), intent(inout) :: Bt(nu, nv)
    real(kind=dp), intent(inout) :: Ct(nu, nv)
    integer :: nv
    integer :: nu
    integer :: j, k, l, nunu, nvnv, nunv
    
    nunu = nu**2
    nvnv = nv**2
    nunv = nu*nv

    !>@todo Add 2d case

    call mxm(A, nv, u, nu, v, nunu)
    k = 1
    l = 1
    do j = 1, nu
       call mxm(v(k, 1, 1), nv, Bt, nu, w(l), nv)
       k = k + nunv
       l = l + nvnv
    end do
    call mxm(w, nvnv, Ct, nu, v, nv)
  end subroutine tensr3

  !> Transpose of a rectangular tensor \f$ A = B^T \f$
  subroutine trsp(a, lda, b, ldb)
    real(kind=dp), intent(inout) :: a(lda, ldb)
    real(kind=dp), intent(in) :: b(ldb, lda)
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer :: i, j

    do j = 1, ldb
       do i = 1, lda
          a(i, j) = b(j, i)
       end do
    end do
    
  end subroutine trsp

  !> In-place transpose of a square tensor
  subroutine trsp1(a, n)
    real(kind=dp), intent(inout) :: a(n, n)
    integer, intent(in) :: n
    real(kind=dp) :: tmp
    integer :: i, j

    do j = 1, n
       do i = j + 1, n
          tmp = a(i, j)
          a(i, j) = a(j, i)
          a(j, i) = tmp
       end do
    end do
    
  end subroutine trsp1
  
end module tensor
