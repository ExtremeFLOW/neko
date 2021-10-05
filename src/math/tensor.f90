!> Tensor operations
module tensor
  use tensor_xsmm
  use tensor_cpu
  use tensor_sx
  use num_types
  use mxm_wrapper
  use neko_config
  implicit none
  private

  interface transpose
     module procedure trsp, trsp1
  end interface transpose

public tensr3, transpose, trsp, trsp1, &
     tnsr2d_el, tnsr3d_el, tnsr3d, tnsr1_3d, addtnsr

contains

  !> Tensor product \f$ v =(C \otimes B \otimes A) u \f$
  subroutine tensr3(v, nv, u, nu, A, Bt, Ct, w)
    integer :: nv
    integer :: nu
    real(kind=rp), intent(inout) :: v(nv, nv, nv)
    real(kind=rp), intent(inout) :: u(nu, nu, nu)
    real(kind=rp), intent(inout) :: w(nu*nu*nv)
    real(kind=rp), intent(inout) :: A(nv, nu)
    real(kind=rp), intent(inout) :: Bt(nu, nv)
    real(kind=rp), intent(inout) :: Ct(nu, nv)
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
    real(kind=rp), intent(inout) :: a(lda, ldb)
    real(kind=rp), intent(in) :: b(ldb, lda)
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
    real(kind=rp), intent(inout) :: a(n, n)
    real(kind=rp) :: tmp
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
  subroutine tnsr2d_el(v, nv, u, nu, A, Bt)
    integer, intent(in) :: nv, nu
    real(kind=rp), intent(inout) :: v(nv*nv), u(nu*nu)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu,nv)

    if (NEKO_BCKND_SX .eq. 1) then
       call tnsr2d_el_sx(v, nv, u, nu, A, Bt)
    else if (NEKO_BCKND_XSMM .eq. 1) then
       call tnsr2d_el_xsmm(v, nv, u, nu, A, Bt)
    else
       call tnsr2d_el_cpu(v, nv, u, nu, A, Bt)
    end if
    
  end subroutine tnsr2d_el


  subroutine tnsr3d_el(v, nv, u, nu, A, Bt, Ct)
    integer, intent(in) :: nv, nu
    real(kind=rp), intent(inout) :: v(nv*nv*nv), u(nu*nu*nu)
    real(kind=rp), intent(inout) :: A(nv,nu),Bt(nu, nv),Ct(nu,nv)

    if (NEKO_BCKND_SX .eq. 1) then
       call tnsr3d_el_sx(v, nv, u, nu, A, Bt, Ct)
    else if (NEKO_BCKND_XSMM .eq. 1) then
       call tnsr3d_el_xsmm(v, nv, u, nu, A, Bt, Ct)
    else
       call tnsr3d_el_cpu(v, nv, u, nu, A, Bt, Ct)
    end if
    
  end subroutine tnsr3d_el
 
  subroutine tnsr3d(v, nv, u, nu, A, Bt, Ct, nelv)
    integer, intent(inout) :: nv, nu, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv,nelv), u(nu*nu*nu,nelv)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu, nv), Ct(nu,nv)

    if (NEKO_BCKND_SX .eq. 1) then
       call tnsr3d_sx(v, nv, u, nu, A, Bt, Ct, nelv)
    else if (NEKO_BCKND_XSMM .eq. 1) then
       call tnsr3d_xsmm(v, nv, u, nu, A, Bt, Ct, nelv)
    else
       call tnsr3d_cpu(v, nv, u, nu, A, Bt, Ct, nelv)
    end if
    
  end subroutine tnsr3d

  subroutine tnsr1_3d(v, nv, nu, A, Bt, Ct, nelv) ! v = [C (x) B (x) A] u
    integer, intent(inout) :: nv, nu, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv*nelv)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu, nv), Ct(nu,nv)

    if (NEKO_BCKND_SX .eq. 1) then
       call tnsr1_3d_sx(v, nv, nu, A, Bt, Ct, nelv)
    else if (NEKO_BCKND_XSMM .eq. 1) then
       call tnsr1_3d_xsmm(v, nv, nu, A, Bt, Ct, nelv)
    else
       call tnsr1_3d_cpu(v, nv, nu, A, Bt, Ct, nelv)
    end if

  end subroutine tnsr1_3d

  subroutine addtnsr(s, h1, h2, h3, nx, ny, nz)

    !Map and add to S a tensor product form of the three functions H1,H2,H3.
    !This is a single element routine used for deforming geometry.
    integer, intent(in) :: nx, ny, nz
    real(kind=rp), intent(in) :: h1(nx), h2(ny), h3(nz) 
    real(kind=rp), intent(inout) ::  s(nx, ny, nz)
    real(kind=rp) :: hh
    integer :: ix, iy, iz
  
    do iz = 1,nz
       do iy = 1,ny
          hh = h2(iy)*h3(iz)
          do ix = 1,nx
             s(ix,iy,iz) = s(ix,iy,iz)+hh*h1(ix)
          end do
       end do
    end do
    
  end subroutine addtnsr
  
end module tensor
