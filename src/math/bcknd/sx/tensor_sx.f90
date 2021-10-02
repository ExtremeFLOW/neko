module tensor_sx
  use num_types
  use mxm_wrapper
  implicit none

contains

  subroutine tnsr2d_el_sx(v, nv, u, nu, A, Bt)
    integer, intent(in) :: nv, nu
    real(kind=rp), intent(inout) :: v(nv*nv), u(nu*nu)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu,nv)
    real(kind=rp) :: work(0:nu**2*nv)

    call mxm(A, nv, u, nu, work, nu)
    call mxm(work, nv, Bt, nu, v, nv)
    
  end subroutine tnsr2d_el_sx
  
  subroutine tnsr3d_el_sx(v, nv, u, nu, A, Bt, Ct)
    integer, intent(in) :: nv, nu
    real(kind=rp), intent(inout) :: v(nv*nv*nv), u(nu*nu*nu)
    real(kind=rp), intent(inout) :: A(nv,nu),Bt(nu, nv),Ct(nu,nv)
    real(kind=rp) :: work(nu**2*nv), work2(nu*nv**2)
    real(kind=rp) :: tmp
    integer :: i, j, k, l, nunu, nvnu, nvnv
    integer :: ii, jj
    nvnu = nv * nu
    nunu = nu * nu 
    nvnv = nv * nv

    do j = 1, nunu
       do i = 1, nv
          ii = i + nv * (j - 1)
          tmp = 0.0_rp
          do k = 1, nu
             tmp = tmp + A(i,k) * u(k + nu * (j - 1))
          end do
          work(ii) = tmp
       end do
    end do
    
    do i = 1, nu
       do j = 1, nv
          do l = 1, nv
             ii = l + nv * (j - 1) + nvnv * (i - 1)
             tmp = 0.0_rp
             do k = 1, nu
                jj = l + nv * (k - 1) + nvnu * (i - 1)
                tmp = tmp + work(jj) * Bt(k,j)
             end do
             work2(ii) = tmp
          end do
       end do
    end do
     
    do j = 1, nv
       do i = 1, nvnv
          jj = i + nvnv * (j - 1)
          tmp = 0.0_rp
          do k = 1, nu
             ii = i + nvnv * (k - 1)
             tmp = tmp + work2(ii) * Ct(k, j)
          end do
          v(jj) = tmp
       end do
    end do
    
  end subroutine tnsr3d_el_sx

  subroutine tnsr3d_sx(v, nv, u, nu, A, Bt, Ct, nelv)
    integer, intent(in) :: nv, nu, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv,nelv), u(nu*nu*nu,nelv)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu, nv), Ct(nu,nv)
    real(kind=rp) :: work(nu**2*nv), work2(nu*nv**2), tmp
    integer :: ie, i, j, k, l, ii, jj
    integer :: nunu, nvnu, nvnv

    nvnu = nv * nu
    nunu = nu * nu 
    nvnv = nv * nv

    do ie = 1,nelv
       do j = 1, nunu
          do i = 1, nv
             ii = i + nv * (j - 1)
             tmp = 0.0_rp
             do k = 1, nu
                tmp = tmp + A(i,k) * u(k + nu * (j - 1), ie)
             end do
             work(ii) = tmp
          end do
       end do
       
       do i = 1, nu
          do j = 1, nv
             do l = 1, nv
                ii = l + nv * (j - 1) + nvnv * (i - 1)
                tmp = 0.0_rp
                do k = 1, nu
                   jj = l + nv * (k - 1) + nvnu * (i - 1)
                   tmp = tmp + work(jj) * Bt(k,j)
                end do
                work2(ii) = tmp
             end do
          end do
       end do
       
       do j = 1, nv
          do i = 1, nvnv
             jj = i + nvnv * (j - 1)
             tmp = 0.0_rp
             do k = 1, nu
                ii = i + nvnv * (k - 1)
                tmp = tmp + work2(ii) * Ct(k, j)
             end do
             v(jj, ie) = tmp
          end do
       end do
    end do
    
  end subroutine tnsr3d_sx

  subroutine tnsr1_3d_sx(v, nv, nu, A, Bt, Ct, nelv)
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
  end subroutine tnsr1_3d_sx
  
end module tensor_sx
