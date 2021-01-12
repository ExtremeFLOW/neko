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
    integer :: nv
    integer :: nu
    real(kind=dp), intent(inout) :: v(nv, nv, nv)
    real(kind=dp), intent(inout) :: u(nv, nv, nv)
    real(kind=dp), intent(inout) :: w(nu*nu*nv)
    real(kind=dp), intent(inout) :: A(nv, nu)
    real(kind=dp), intent(inout) :: Bt(nu, nv)
    real(kind=dp), intent(inout) :: Ct(nu, nv)
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
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    real(kind=dp), intent(inout) :: a(lda, ldb)
    real(kind=dp), intent(in) :: b(ldb, lda)
    integer :: i, j

    do j = 1, ldb
       do i = 1, lda
          a(i, j) = b(j, i)
       end do
    end do
    
  end subroutine trsp

  !> In-place transpose of a square tensor
  subroutine trsp1(a, n)
    integer, intent(in) :: n
    real(kind=dp), intent(inout) :: a(n, n)
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
!      computes
!              T
!     v = A u B
  subroutine tnsr2d_el(v,nv,u,nu,A,Bt)
    integer, intent(in) :: nv,nu
    real(kind=dp), intent(inout) :: v(nv*nv),u(nu*nu),A(nv,nu),Bt(nu,nv)
    real(kind=dp) :: work(0:nu*nv*nu)

    call mxm(A,nv,u,nu,work,nu)
    call mxm(work,nv,Bt,nu,v,nv)
  end subroutine tnsr2d_el


  subroutine tnsr3d_el(v,nv,u,nu,A,Bt,Ct)
    integer, intent(in) :: nv,nu
    real(kind=dp), intent(inout) :: v(nv*nv*nv),u(nu*nu*nu),A(nv,nu),Bt(nu, nv),Ct(nu,nv)
    real(kind=dp) :: work(0:(nu+nv)**3),work2(0:(nv*nv*nu))
    integer :: i, nunu, nvnu, nvnv
    nvnu = nv * nu
    nunu = nu * nu 
    nvnv = nv * nv
    call mxm(A,nv,u(1),nu,work,nunu)
    do i=0,nu-1
       call mxm(work(nvnu*i),nv,Bt,nu,work2(nv*nv*i),nv)
    enddo
    call mxm(work2,nvnv,Ct,nu,v(1),nv)
  end subroutine tnsr3d_el
 
  subroutine tnsr3d(v,nv,u,nu,A,Bt,Ct, nelv)
    integer, intent(inout) :: nv,nu, nelv
    real(kind=dp), intent(inout) :: v(nv*nv*nv,nelv),u(nu*nu*nu,nelv),A(nv,nu),Bt(nu, nv),Ct(nu,nv)
    real(kind=dp) :: work(0:(nu+nv)**3),work2(0:nu*nv*nv)
    integer :: ie, i, nunu, nvnu, nvnv
    nvnu = nv * nu
    nunu = nu * nu 
    nvnv = nv * nv
    do ie=1,nelv
       call mxm(A,nv,u(1,ie),nu,work,nunu)
       do i=0,nu-1
          call mxm(work(nvnu*i),nv,Bt,nu,work2(nv*nv*i),nv)
       enddo
       call mxm(work2,nvnv,Ct,nu,v(1,ie),nv)
    enddo
  end subroutine tnsr3d

  subroutine tnsr1_3d(v,nv,nu,A,Bt,Ct, nelv) ! v = [C (x) B (x) A] u
    integer, intent(in) :: nv,nu, nelv
    real(kind=dp), intent(inout) :: v(nv*nv*nv*nelv),A(nv,nu),Bt(nu, nv),Ct(nu,nv)
    real(kind=dp) :: work(0:(nu+nv)**3),work2(0:(nu+nv)**3)
    integer :: e,e0,ee,es, iu, iv, i, nu3, nv3
    e0=1
    es=1
    ee=nelv

    if (nv.gt.nu) then
       e0=nelv
       es=-1
       ee=1
    endif

    nu3 = nu**3
    nv3 = nv**3

    do e=e0,ee,es
       iu = 1 + (e-1)*nu3
       iv = 1 + (e-1)*nv3
       call mxm(A,nv,v(iu),nu,work,nu*nu)
       do i=0,nu-1
          call mxm(work(nv*nu*i),nv,Bt,nu,work2(nv*nv*i),nv)
       enddo
       call mxm(work2,nv*nv,Ct,nu,v(iv),nv)
    enddo
    end subroutine tnsr1_3d

end module tensor
