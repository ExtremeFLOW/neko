!> Tensor operations libxsmm backend
module tensor_xsmm
  use num_types
  use mxm_wrapper
  implicit none

contains
  
  subroutine tnsr2d_el_xsmm(v, nv, u, nu, A, Bt)
    integer, intent(in) :: nv, nu
    real(kind=rp), intent(inout) :: v(nv*nv), u(nu*nu)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu,nv)
    real(kind=rp) :: work(0:nu**2*nv)

    call mxm(A, nv, u, nu, work, nu)
    call mxm(work, nv, Bt, nu, v, nv)
    
  end subroutine tnsr2d_el_xsmm

  subroutine tnsr3d_el_xsmm(v, nv, u, nu, A, Bt, Ct)
    integer, intent(in) :: nv, nu
    real(kind=rp), intent(inout) :: v(nv*nv*nv), u(nu*nu*nu)
    real(kind=rp), intent(inout) :: A(nv,nu),Bt(nu, nv),Ct(nu,nv)
    real(kind=rp) :: work(0:nu**2*nv), work2(0:nu*nv**2)
    integer :: i, nunu, nvnu, nvnv

    nvnu = nv * nu
    nunu = nu * nu 
    nvnv = nv * nv
    
    call mxm(A, nv, u(1), nu ,work, nunu)
    do i = 0,nu-1
       call mxm(work(nvnu*i), nv, Bt, nu, work2(nv*nv*i), nv)
    end do
    call mxm(work2, nvnv, Ct, nu, v(1), nv)
    
  end subroutine tnsr3d_el_xsmm
 
  subroutine tnsr3d_xsmm(v, nv, u, nu, A, Bt, Ct, nelv)
    integer, intent(inout) :: nv, nu, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv,nelv), u(nu*nu*nu,nelv)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu, nv), Ct(nu,nv)
    real(kind=rp) :: work(0:nu**2*nv), work2(0:nu*nv**2)
    integer :: ie, i, nunu, nvnu, nvnv

    nvnu = nv * nu
    nunu = nu * nu 
    nvnv = nv * nv
    
    do ie = 1,nelv
       call mxm(A, nv, u(1,ie), nu, work, nunu)
       do i = 0,nu-1
          call mxm(work(nvnu*i), nv, Bt, nu, work2(nv*nv*i), nv)
       end do
       call mxm(work2, nvnv, Ct, nu, v(1,ie), nv)
    end do
    
  end subroutine tnsr3d_xsmm

  subroutine tnsr1_3d_xsmm(v, nv, nu, A, Bt, Ct, nelv) 
    integer, intent(in) :: nv, nu, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv*nelv)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu, nv), Ct(nu,nv)
    real(kind=rp) :: work(0:nu**2*nv), work2(0:nu*nv**2)
    integer :: e, e0, ee, es, iu, iv, i, nu3, nv3

    e0 = 1
    es = 1
    ee = nelv

    if (nv.gt.nu) then
       e0 = nelv
       es = -1
       ee = 1
    endif

    nu3 = nu**3
    nv3 = nv**3

    do e = e0,ee,es
       iu = 1 + (e-1)*nu3
       iv = 1 + (e-1)*nv3
       call mxm(A, nv, v(iu), nu, work, nu*nu)
       do i = 0,nu-1
          call mxm(work(nv*nu*i), nv, Bt, nu, work2(nv*nv*i), nv)
       end do
       call mxm(work2, nv*nv, Ct, nu, v(iv), nv)
    end do
  end subroutine tnsr1_3d_xsmm

end module tensor_xsmm
