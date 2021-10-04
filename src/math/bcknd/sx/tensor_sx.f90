!> Tensor operations SX-Aurora backend
module tensor_sx
  use num_types
  use mxm_wrapper
  implicit none
  private

  public :: tnsr2d_el_sx, tnsr3d_el_sx, tnsr3d_sx

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

    if (nu .eq. 2 .and. nv .eq. 4) then
       call tnsr3d_nu2nv4_sx(v, u, A, Bt, Ct, nelv)
    else if (nu .eq. 4) then
       call tnsr3d_nu4_sx(v, nv, u, A, Bt, Ct, nelv)
    else
       call tnsr3d_nvnu_sx(v, nv, u, nu, A, Bt, Ct, nelv)
    end if

  end subroutine tnsr3d_sx

  subroutine tnsr3d_nvnu_sx(v, nv, u, nu, A, Bt, Ct, nelv)
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

  end subroutine tnsr3d_nvnu_sx

  subroutine tnsr3d_nu2nv4_sx(v, u, A, Bt, Ct, nelv)
    integer, parameter :: nu = 2
    integer, parameter :: nv = 4
    integer, parameter :: nunu = 4
    integer, parameter :: nvnu = 8
    integer, parameter :: nvnv = 16
    integer, intent(in) :: nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv,nelv), u(nu*nu*nu,nelv)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu, nv), Ct(nu,nv)
    real(kind=rp) :: work(nu**2*nv,nelv), work2(nu*nv**2,nelv), tmp
    integer :: ie, i, j, k, l, ii, jj
    

    do j = 1, nunu
       do i = 1, nv
          do ie = 1,nelv
             ii = i + nv * (j - 1)
             work(ii, ie) = A(i,1) * u(1 + nu * (j - 1), ie) &
                          + A(i,2) * u(2 + nu * (j - 1), ie)             
          end do
       end do
    end do
    
    do i = 1, nu
       do j = 1, nv
          do l = 1, nv
             do ie = 1,nelv
                ii = l + nv * (j - 1) + nvnv * (i - 1)
                tmp = 0.0_rp
                !NEC$ unroll_completely
                do k = 1, nu
                   jj = l + nv * (k - 1) + nvnu * (i - 1)
                   tmp = tmp + work(jj,ie) * Bt(k,j)
                end do
                work2(ii, ie) = tmp
             end do
          end do
       end do
    end do
    
    do j = 1, nv
       do i = 1, nvnv
          do ie = 1,nelv
             jj = i + nvnv * (j - 1)
             v(jj, ie) = work2(i + nvnv * (1 - 1),ie) * Ct(1, j) &
                       + work2(i + nvnv * (2 - 1),ie) * Ct(2, j) 
          end do
       end do
    end do

  end subroutine tnsr3d_nu2nv4_sx

  subroutine tnsr3d_nu4_sx(v, nv, u, A, Bt, Ct, nelv)
    integer, parameter :: nu = 4
    integer, parameter :: nunu = 16
    integer, intent(in) :: nv, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv,nelv), u(nu*nu*nu,nelv)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu, nv), Ct(nu,nv)
    real(kind=rp) :: work(nu**2*nv, nelv), work2(nu*nv**2, nelv), tmp
    integer :: ie, i, j, k, l, ii, jj
    integer :: nvnu, nvnv

    nvnu = nv * nu
    nvnv = nv * nv
    
    do j = 1, nunu
       do i = 1, nv
          do ie = 1,nelv                 
             ii = i + nv * (j - 1)
             work(ii, ie) = A(i,1) * u(1 + nu * (j - 1), ie) &
                          + A(i,2) * u(2 + nu * (j - 1), ie) &
                          + A(i,3) * u(3 + nu * (j - 1), ie) &
                          + A(i,4) * u(4 + nu * (j - 1), ie)             
          end do
       end do
    end do

    do i = 1, nu
       do j = 1, nv
          do l = 1, nv
             do ie = 1,nelv                
                ii = l + nv * (j - 1) + nvnv * (i - 1)
                tmp = 0.0_rp
                !NEC$ unroll_completely
                do k = 1, nu
                   jj = l + nv * (k - 1) + nvnu * (i - 1)
                   tmp = tmp + work(jj,ie) * Bt(k,j)
                end do
                work2(ii, ie) = tmp
             end do
          end do
       end do
    end do
    
    do j = 1, nv
       do i = 1, nvnv
          do ie = 1,nelv                 
             jj = i + nvnv * (j - 1)
             v(jj, ie) = work2(i + nvnv * (1 - 1),ie) * Ct(1, j) &
                       + work2(i + nvnv * (2 - 1),ie) * Ct(2, j) &
                       + work2(i + nvnv * (3 - 1),ie) * Ct(3, j) &
                       + work2(i + nvnv * (4 - 1),ie) * Ct(4, j) 
          end do
       end do
    end do

  end subroutine tnsr3d_nu4_sx

  subroutine tnsr1_3d_sx(v, nv, nu, A, Bt, Ct, nelv)
    integer, intent(in) :: nv, nu, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv*nelv)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu, nv), Ct(nu,nv)
    real(kind=rp) :: work(nu**2*nv), work2(nu*nv**2)
    integer :: e, e0, ee, es, iu, iv, nu3, nv3
    integer :: i, j, k, l, ii, jj, kk
    integer :: nunu, nvnu, nvnv
    real(kind=rp) :: tmp


    if (nu .eq. 4 .and. nv .eq. 2) then
       call tnsr1_3d_nu4nv2_sx(v, A, Bt, Ct, nelv)
    else
       call tnsr1_3d_nvnu_sx(v, nv, nu, A, Bt, Ct, nelv)
    end if
    
  end subroutine tnsr1_3d_sx

  subroutine tnsr1_3d_nvnu_sx(v, nv, nu, A, Bt, Ct, nelv)
    integer, intent(in) :: nv, nu, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv*nelv)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu, nv), Ct(nu,nv)
    real(kind=rp) :: work(nu**2*nv), work2(nu*nv**2)
    integer :: e, e0, ee, es, iu, iv, nu3, nv3
    integer :: i, j, k, l, ii, jj, kk
    integer :: nunu, nvnu, nvnv
    real(kind=rp) :: tmp

    nvnu = nv * nu
    nunu = nu * nu 
    nvnv = nv * nv
    
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
       iu = (e-1)*nu3
       iv = (e-1)*nv3
          
       do j = 1, nunu
          do i = 1, nv
             ii = i + nv * (j - 1)
             tmp = 0.0_rp
             do k = 1, nu
                kk = k + nu * (j - 1) + iu
                tmp = tmp + A(i,k) * v(kk)
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
             jj = i + nvnv * (j - 1) + iv
             tmp = 0.0_rp
             do k = 1, nu
                ii = i + nvnv * (k - 1)
                tmp = tmp + work2(ii) * Ct(k, j)
             end do
             v(jj) = tmp
          end do
       end do
    end do
    
  end subroutine tnsr1_3d_nvnu_sx

  subroutine tnsr1_3d_nu4nv2_sx(v, A, Bt, Ct, nelv)
    integer, parameter :: nu = 4
    integer, parameter :: nv = 2
    integer, parameter :: nunu = 16
    integer, parameter :: nvnu = 8
    integer, parameter :: nvnv = 4
    integer, parameter :: nununu = 64
    integer, parameter :: nvnvnv = 8
    integer, intent(in) :: nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv*nelv)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu, nv), Ct(nu,nv)
    real(kind=rp) :: work(nu**2*nv,nelv), work2(nu*nv**2,nelv)
    integer :: ie, e0, ee, es, iu, iv
    integer :: i, j, k, l, ii, jj, kk
    real(kind=rp) :: tmp

    do j = 1, nunu
       do i = 1, nv
          do ie = 1,nelv
             iu = (ie-1)*nununu
             ii = i + nv * (j - 1)
             work(ii, ie) = A(i,1) * v(1 + nu * (j - 1) + iu) &
                          + A(i,2) * v(2 + nu * (j - 1) + iu) &
                          + A(i,3) * v(3 + nu * (j - 1) + iu) &
                          + A(i,4) * v(4 + nu * (j - 1) + iu)
          end do
       end do
    end do
    
    do i = 1, nu
       do j = 1, nv
          do l = 1, nv
             do ie = 1,nelv
                ii = l + nv * (j - 1) + nvnv * (i - 1)
                tmp = 0.0_rp
                !NEC$ unroll_completely                
                do k = 1, nu
                   jj = l + nv * (k - 1) + nvnu * (i - 1)
                   tmp = tmp + work(jj,ie) * Bt(k,j)
                end do
                work2(ii,ie) = tmp
             end do
          end do
       end do
    end do
          
    do j = 1, nv
       do i = 1, nvnv
          do ie = 1,nelv
             iv = (ie-1)*nvnvnv
             jj = i + nvnv * (j - 1) + iv
             v(jj) = work2(i + nvnv * (1 - 1),ie) * Ct(1, j) &
                   + work2(i + nvnv * (2 - 1),ie) * Ct(2, j) &
                   + work2(i + nvnv * (3 - 1),ie) * Ct(3, j) &
                   + work2(i + nvnv * (4 - 1),ie) * Ct(4, j)                
          end do
       end do
    end do
    
  end subroutine tnsr1_3d_nu4nv2_sx
  
end module tensor_sx
