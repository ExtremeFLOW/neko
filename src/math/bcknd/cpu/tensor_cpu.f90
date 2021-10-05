module tensor_cpu
  use num_types
  use mxm_wrapper
  implicit none
  private

  public :: tnsr2d_el_cpu, tnsr3d_el_cpu, tnsr3d_cpu, tnsr1_3d_cpu
 
contains

  subroutine tnsr2d_el_cpu(v, nv, u, nu, A, Bt)
    integer, intent(in) :: nv, nu
    real(kind=rp), intent(inout) :: v(nv*nv), u(nu*nu)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu,nv)
    real(kind=rp) :: work(0:nu**2*nv)

    call mxm(A, nv, u, nu, work, nu)
    call mxm(work, nv, Bt, nu, v, nv)
    
  end subroutine tnsr2d_el_cpu
  
  subroutine tnsr3d_el_cpu(v, nv, u, nu, A, Bt, Ct)
    integer, intent(in) :: nv, nu
    real(kind=rp), intent(inout) :: v(nv*nv*nv), u(nu*nu*nu)
    real(kind=rp), intent(inout) :: A(nv,nu),Bt(nu, nv),Ct(nu,nv)

    if (nv .eq. nu) then
       select case (nv)
       case (14)
          call tnsr3d_el_n14_cpu(v, u, A, Bt, Ct)
       case (13)
          call tnsr3d_el_n13_cpu(v, u, A, Bt, Ct)
       case (12)
          call tnsr3d_el_n12_cpu(v, u, A, Bt, Ct)
       case (11)
          call tnsr3d_el_n11_cpu(v, u, A, Bt, Ct)
       case (10)
          call tnsr3d_el_n10_cpu(v, u, A, Bt, Ct)
       case (9)
          call tnsr3d_el_n9_cpu(v, u, A, Bt, Ct)
       case (8)
          call tnsr3d_el_n8_cpu(v, u, A, Bt, Ct)
       case (7)
          call tnsr3d_el_n7_cpu(v, u, A, Bt, Ct)
       case (6)
          call tnsr3d_el_n6_cpu(v, u, A, Bt, Ct)
       case (5)
          call tnsr3d_el_n5_cpu(v, u, A, Bt, Ct)
       case (4)
          call tnsr3d_el_n4_cpu(v, u, A, Bt, Ct)
       case (3)
          call tnsr3d_el_n3_cpu(v, u, A, Bt, Ct)
       case (2)
          call tnsr3d_el_n2_cpu(v, u, A, Bt, Ct)
       end select
    else
       call tnsr3d_el_nvnu_cpu(v, nv, u, nu, A, Bt, Ct)
    end if
    
  end subroutine tnsr3d_el_cpu

  subroutine tnsr3d_el_nvnu_cpu(v, nv, u, nu, A, Bt, Ct)
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
    
  end subroutine tnsr3d_el_nvnu_cpu

  subroutine tnsr3d_el_n14_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 14
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, k, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
                   + A(i,2) * u(2 + n * (j - 1)) &
                   + A(i,3) * u(3 + n * (j - 1)) &
                   + A(i,4) * u(4 + n * (j - 1)) &
                   + A(i,5) * u(5 + n * (j - 1)) &
                   + A(i,6) * u(6 + n * (j - 1)) &
                   + A(i,7) * u(7 + n * (j - 1)) &
                   + A(i,8) * u(8 + n * (j - 1)) &
                   + A(i,9) * u(9 + n * (j - 1)) &
                   + A(i,10) * u(10 + n * (j - 1)) &
                   + A(i,11) * u(11 + n * (j - 1)) &
                   + A(i,12) * u(12 + n * (j - 1)) &
                   + A(i,13) * u(13 + n * (j - 1)) &
                   + A(i,14) * u(14 + n * (j - 1))                    
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                       + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                       + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                       + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) &
                       + work(l + n * (5 - 1) + nn * (i - 1)) * Bt(5,j) &
                       + work(l + n * (6 - 1) + nn * (i - 1)) * Bt(6,j) &
                       + work(l + n * (7 - 1) + nn * (i - 1)) * Bt(7,j) &
                       + work(l + n * (8 - 1) + nn * (i - 1)) * Bt(8,j) &
                       + work(l + n * (9 - 1) + nn * (i - 1)) * Bt(9,j) &
                       + work(l + n * (10 - 1) + nn * (i - 1)) * Bt(10,j) &
                       + work(l + n * (11 - 1) + nn * (i - 1)) * Bt(11,j) &
                       + work(l + n * (12 - 1) + nn * (i - 1)) * Bt(12,j) &
                       + work(l + n * (13 - 1) + nn * (i - 1)) * Bt(13,j) &
                       + work(l + n * (14 - 1) + nn * (i - 1)) * Bt(14,j) 
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
                + work2(i + nn * (2 - 1)) * Ct(2, j) &
                + work2(i + nn * (3 - 1)) * Ct(3, j) &
                + work2(i + nn * (4 - 1)) * Ct(4, j) &
                + work2(i + nn * (5 - 1)) * Ct(5, j) &
                + work2(i + nn * (6 - 1)) * Ct(6, j) &
                + work2(i + nn * (7 - 1)) * Ct(7, j) &
                + work2(i + nn * (8 - 1)) * Ct(8, j) &
                + work2(i + nn * (9 - 1)) * Ct(9, j) &
                + work2(i + nn * (10 - 1)) * Ct(10, j) &
                + work2(i + nn * (11 - 1)) * Ct(11, j) &
                + work2(i + nn * (12 - 1)) * Ct(12, j) &
                + work2(i + nn * (13 - 1)) * Ct(13, j) &
                + work2(i + nn * (14 - 1)) * Ct(14, j) 
       end do
    end do
    
  end subroutine tnsr3d_el_n14_cpu

  subroutine tnsr3d_el_n13_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 13
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, k, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
                   + A(i,2) * u(2 + n * (j - 1)) &
                   + A(i,3) * u(3 + n * (j - 1)) &
                   + A(i,4) * u(4 + n * (j - 1)) &
                   + A(i,5) * u(5 + n * (j - 1)) &
                   + A(i,6) * u(6 + n * (j - 1)) &
                   + A(i,7) * u(7 + n * (j - 1)) &
                   + A(i,8) * u(8 + n * (j - 1)) &
                   + A(i,9) * u(9 + n * (j - 1)) &
                   + A(i,10) * u(10 + n * (j - 1)) &
                   + A(i,11) * u(11 + n * (j - 1)) &
                   + A(i,12) * u(12 + n * (j - 1)) &
                   + A(i,13) * u(13 + n * (j - 1)) 
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                       + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                       + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                       + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) &
                       + work(l + n * (5 - 1) + nn * (i - 1)) * Bt(5,j) &
                       + work(l + n * (6 - 1) + nn * (i - 1)) * Bt(6,j) &
                       + work(l + n * (7 - 1) + nn * (i - 1)) * Bt(7,j) &
                       + work(l + n * (8 - 1) + nn * (i - 1)) * Bt(8,j) &
                       + work(l + n * (9 - 1) + nn * (i - 1)) * Bt(9,j) &
                       + work(l + n * (10 - 1) + nn * (i - 1)) * Bt(10,j) &
                       + work(l + n * (11 - 1) + nn * (i - 1)) * Bt(11,j) &
                       + work(l + n * (12 - 1) + nn * (i - 1)) * Bt(12,j) &
                       + work(l + n * (13 - 1) + nn * (i - 1)) * Bt(13,j) 
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
                + work2(i + nn * (2 - 1)) * Ct(2, j) &
                + work2(i + nn * (3 - 1)) * Ct(3, j) &
                + work2(i + nn * (4 - 1)) * Ct(4, j) &
                + work2(i + nn * (5 - 1)) * Ct(5, j) &
                + work2(i + nn * (6 - 1)) * Ct(6, j) &
                + work2(i + nn * (7 - 1)) * Ct(7, j) &
                + work2(i + nn * (8 - 1)) * Ct(8, j) &
                + work2(i + nn * (9 - 1)) * Ct(9, j) &
                + work2(i + nn * (10 - 1)) * Ct(10, j) &
                + work2(i + nn * (11 - 1)) * Ct(11, j) &
                + work2(i + nn * (12 - 1)) * Ct(12, j) &
                + work2(i + nn * (13 - 1)) * Ct(13, j) 
       end do
    end do
    
  end subroutine tnsr3d_el_n13_cpu

  subroutine tnsr3d_el_n12_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 12
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, k, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
                   + A(i,2) * u(2 + n * (j - 1)) &
                   + A(i,3) * u(3 + n * (j - 1)) &
                   + A(i,4) * u(4 + n * (j - 1)) &
                   + A(i,5) * u(5 + n * (j - 1)) &
                   + A(i,6) * u(6 + n * (j - 1)) &
                   + A(i,7) * u(7 + n * (j - 1)) &
                   + A(i,8) * u(8 + n * (j - 1)) &
                   + A(i,9) * u(9 + n * (j - 1)) &
                   + A(i,10) * u(10 + n * (j - 1)) &
                   + A(i,11) * u(11 + n * (j - 1)) &
                   + A(i,12) * u(12 + n * (j - 1)) 
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                       + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                       + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                       + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) &
                       + work(l + n * (5 - 1) + nn * (i - 1)) * Bt(5,j) &
                       + work(l + n * (6 - 1) + nn * (i - 1)) * Bt(6,j) &
                       + work(l + n * (7 - 1) + nn * (i - 1)) * Bt(7,j) &
                       + work(l + n * (8 - 1) + nn * (i - 1)) * Bt(8,j) &
                       + work(l + n * (9 - 1) + nn * (i - 1)) * Bt(9,j) &
                       + work(l + n * (10 - 1) + nn * (i - 1)) * Bt(10,j) &
                       + work(l + n * (11 - 1) + nn * (i - 1)) * Bt(11,j) &
                       + work(l + n * (12 - 1) + nn * (i - 1)) * Bt(12,j) 
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
                + work2(i + nn * (2 - 1)) * Ct(2, j) &
                + work2(i + nn * (3 - 1)) * Ct(3, j) &
                + work2(i + nn * (4 - 1)) * Ct(4, j) &
                + work2(i + nn * (5 - 1)) * Ct(5, j) &
                + work2(i + nn * (6 - 1)) * Ct(6, j) &
                + work2(i + nn * (7 - 1)) * Ct(7, j) &
                + work2(i + nn * (8 - 1)) * Ct(8, j) &
                + work2(i + nn * (9 - 1)) * Ct(9, j) &
                + work2(i + nn * (10 - 1)) * Ct(10, j) &
                + work2(i + nn * (11 - 1)) * Ct(11, j) &
                + work2(i + nn * (12 - 1)) * Ct(12, j) 
       end do
    end do
    
  end subroutine tnsr3d_el_n12_cpu

  subroutine tnsr3d_el_n11_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 11
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, k, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
                   + A(i,2) * u(2 + n * (j - 1)) &
                   + A(i,3) * u(3 + n * (j - 1)) &
                   + A(i,4) * u(4 + n * (j - 1)) &
                   + A(i,5) * u(5 + n * (j - 1)) &
                   + A(i,6) * u(6 + n * (j - 1)) &
                   + A(i,7) * u(7 + n * (j - 1)) &
                   + A(i,8) * u(8 + n * (j - 1)) &
                   + A(i,9) * u(9 + n * (j - 1)) &
                   + A(i,10) * u(10 + n * (j - 1)) &
                   + A(i,11) * u(11 + n * (j - 1)) 
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                       + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                       + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                       + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) &
                       + work(l + n * (5 - 1) + nn * (i - 1)) * Bt(5,j) &
                       + work(l + n * (6 - 1) + nn * (i - 1)) * Bt(6,j) &
                       + work(l + n * (7 - 1) + nn * (i - 1)) * Bt(7,j) &
                       + work(l + n * (8 - 1) + nn * (i - 1)) * Bt(8,j) &
                       + work(l + n * (9 - 1) + nn * (i - 1)) * Bt(9,j) &
                       + work(l + n * (10 - 1) + nn * (i - 1)) * Bt(10,j) &
                       + work(l + n * (11 - 1) + nn * (i - 1)) * Bt(11,j) 
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
                + work2(i + nn * (2 - 1)) * Ct(2, j) &
                + work2(i + nn * (3 - 1)) * Ct(3, j) &
                + work2(i + nn * (4 - 1)) * Ct(4, j) &
                + work2(i + nn * (5 - 1)) * Ct(5, j) &
                + work2(i + nn * (6 - 1)) * Ct(6, j) &
                + work2(i + nn * (7 - 1)) * Ct(7, j) &
                + work2(i + nn * (8 - 1)) * Ct(8, j) &
                + work2(i + nn * (9 - 1)) * Ct(9, j) &
                + work2(i + nn * (10 - 1)) * Ct(10, j) &
                + work2(i + nn * (11 - 1)) * Ct(11, j) 
       end do
    end do
    
  end subroutine tnsr3d_el_n11_cpu

  subroutine tnsr3d_el_n10_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 10
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, k, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
                   + A(i,2) * u(2 + n * (j - 1)) &
                   + A(i,3) * u(3 + n * (j - 1)) &
                   + A(i,4) * u(4 + n * (j - 1)) &
                   + A(i,5) * u(5 + n * (j - 1)) &
                   + A(i,6) * u(6 + n * (j - 1)) &
                   + A(i,7) * u(7 + n * (j - 1)) &
                   + A(i,8) * u(8 + n * (j - 1)) &
                   + A(i,9) * u(9 + n * (j - 1)) &
                   + A(i,10) * u(10 + n * (j - 1)) 
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                       + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                       + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                       + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) &
                       + work(l + n * (5 - 1) + nn * (i - 1)) * Bt(5,j) &
                       + work(l + n * (6 - 1) + nn * (i - 1)) * Bt(6,j) &
                       + work(l + n * (7 - 1) + nn * (i - 1)) * Bt(7,j) &
                       + work(l + n * (8 - 1) + nn * (i - 1)) * Bt(8,j) &
                       + work(l + n * (9 - 1) + nn * (i - 1)) * Bt(9,j) &
                       + work(l + n * (10 - 1) + nn * (i - 1)) * Bt(10,j) 
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
                + work2(i + nn * (2 - 1)) * Ct(2, j) &
                + work2(i + nn * (3 - 1)) * Ct(3, j) &
                + work2(i + nn * (4 - 1)) * Ct(4, j) &
                + work2(i + nn * (5 - 1)) * Ct(5, j) &
                + work2(i + nn * (6 - 1)) * Ct(6, j) &
                + work2(i + nn * (7 - 1)) * Ct(7, j) &
                + work2(i + nn * (8 - 1)) * Ct(8, j) &
                + work2(i + nn * (9 - 1)) * Ct(9, j) &
                + work2(i + nn * (10 - 1)) * Ct(10, j) 
       end do
    end do
    
  end subroutine tnsr3d_el_n10_cpu

  subroutine tnsr3d_el_n9_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 9
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, k, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
                   + A(i,2) * u(2 + n * (j - 1)) &
                   + A(i,3) * u(3 + n * (j - 1)) &
                   + A(i,4) * u(4 + n * (j - 1)) &
                   + A(i,5) * u(5 + n * (j - 1)) &
                   + A(i,6) * u(6 + n * (j - 1)) &
                   + A(i,7) * u(7 + n * (j - 1)) &
                   + A(i,8) * u(8 + n * (j - 1)) &
                   + A(i,9) * u(9 + n * (j - 1)) 
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                       + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                       + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                       + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) &
                       + work(l + n * (5 - 1) + nn * (i - 1)) * Bt(5,j) &
                       + work(l + n * (6 - 1) + nn * (i - 1)) * Bt(6,j) &
                       + work(l + n * (7 - 1) + nn * (i - 1)) * Bt(7,j) &
                       + work(l + n * (8 - 1) + nn * (i - 1)) * Bt(8,j) &
                       + work(l + n * (9 - 1) + nn * (i - 1)) * Bt(9,j) 
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
                + work2(i + nn * (2 - 1)) * Ct(2, j) &
                + work2(i + nn * (3 - 1)) * Ct(3, j) &
                + work2(i + nn * (4 - 1)) * Ct(4, j) &
                + work2(i + nn * (5 - 1)) * Ct(5, j) &
                + work2(i + nn * (6 - 1)) * Ct(6, j) &
                + work2(i + nn * (7 - 1)) * Ct(7, j) &
                + work2(i + nn * (8 - 1)) * Ct(8, j) &
                + work2(i + nn * (9 - 1)) * Ct(9, j) 
       end do
    end do
    
  end subroutine tnsr3d_el_n9_cpu

  subroutine tnsr3d_el_n8_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 8
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, k, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
                   + A(i,2) * u(2 + n * (j - 1)) &
                   + A(i,3) * u(3 + n * (j - 1)) &
                   + A(i,4) * u(4 + n * (j - 1)) &
                   + A(i,5) * u(5 + n * (j - 1)) &
                   + A(i,6) * u(6 + n * (j - 1)) &
                   + A(i,7) * u(7 + n * (j - 1)) &
                   + A(i,8) * u(8 + n * (j - 1)) 
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                       + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                       + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                       + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) &
                       + work(l + n * (5 - 1) + nn * (i - 1)) * Bt(5,j) &
                       + work(l + n * (6 - 1) + nn * (i - 1)) * Bt(6,j) &
                       + work(l + n * (7 - 1) + nn * (i - 1)) * Bt(7,j) &
                       + work(l + n * (8 - 1) + nn * (i - 1)) * Bt(8,j) 
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
                + work2(i + nn * (2 - 1)) * Ct(2, j) &
                + work2(i + nn * (3 - 1)) * Ct(3, j) &
                + work2(i + nn * (4 - 1)) * Ct(4, j) &
                + work2(i + nn * (5 - 1)) * Ct(5, j) &
                + work2(i + nn * (6 - 1)) * Ct(6, j) &
                + work2(i + nn * (7 - 1)) * Ct(7, j) &
                + work2(i + nn * (8 - 1)) * Ct(8, j) 
       end do
    end do
    
  end subroutine tnsr3d_el_n8_cpu

  subroutine tnsr3d_el_n7_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 7
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, k, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
                   + A(i,2) * u(2 + n * (j - 1)) &
                   + A(i,3) * u(3 + n * (j - 1)) &
                   + A(i,4) * u(4 + n * (j - 1)) &
                   + A(i,5) * u(5 + n * (j - 1)) &
                   + A(i,6) * u(6 + n * (j - 1)) &
                   + A(i,7) * u(7 + n * (j - 1)) 
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                       + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                       + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                       + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) &
                       + work(l + n * (5 - 1) + nn * (i - 1)) * Bt(5,j) &
                       + work(l + n * (6 - 1) + nn * (i - 1)) * Bt(6,j) &
                       + work(l + n * (7 - 1) + nn * (i - 1)) * Bt(7,j) 
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
                + work2(i + nn * (2 - 1)) * Ct(2, j) &
                + work2(i + nn * (3 - 1)) * Ct(3, j) &
                + work2(i + nn * (4 - 1)) * Ct(4, j) &
                + work2(i + nn * (5 - 1)) * Ct(5, j) &
                + work2(i + nn * (6 - 1)) * Ct(6, j) &
                + work2(i + nn * (7 - 1)) * Ct(7, j) 
       end do
    end do
    
  end subroutine tnsr3d_el_n7_cpu

  subroutine tnsr3d_el_n6_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 6
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, k, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
                   + A(i,2) * u(2 + n * (j - 1)) &
                   + A(i,3) * u(3 + n * (j - 1)) &
                   + A(i,4) * u(4 + n * (j - 1)) &
                   + A(i,5) * u(5 + n * (j - 1)) &
                   + A(i,6) * u(6 + n * (j - 1)) 
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                       + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                       + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                       + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) &
                       + work(l + n * (5 - 1) + nn * (i - 1)) * Bt(5,j) &
                       + work(l + n * (6 - 1) + nn * (i - 1)) * Bt(6,j) 
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
                + work2(i + nn * (2 - 1)) * Ct(2, j) &
                + work2(i + nn * (3 - 1)) * Ct(3, j) &
                + work2(i + nn * (4 - 1)) * Ct(4, j) &
                + work2(i + nn * (5 - 1)) * Ct(5, j) &
                + work2(i + nn * (6 - 1)) * Ct(6, j) 
       end do
    end do
    
  end subroutine tnsr3d_el_n6_cpu

  subroutine tnsr3d_el_n5_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 5
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, k, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
                   + A(i,2) * u(2 + n * (j - 1)) &
                   + A(i,3) * u(3 + n * (j - 1)) &
                   + A(i,4) * u(4 + n * (j - 1)) &
                   + A(i,5) * u(5 + n * (j - 1)) 
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                       + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                       + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                       + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) &
                       + work(l + n * (5 - 1) + nn * (i - 1)) * Bt(5,j)
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
                + work2(i + nn * (2 - 1)) * Ct(2, j) &
                + work2(i + nn * (3 - 1)) * Ct(3, j) &
                + work2(i + nn * (4 - 1)) * Ct(4, j) &
                + work2(i + nn * (5 - 1)) * Ct(5, j) 
       end do
    end do
    
  end subroutine tnsr3d_el_n5_cpu

  subroutine tnsr3d_el_n4_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 4
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, k, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
                   + A(i,2) * u(2 + n * (j - 1)) &
                   + A(i,3) * u(3 + n * (j - 1)) &
                   + A(i,4) * u(4 + n * (j - 1)) 
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                       + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                       + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) &
                       + work(l + n * (4 - 1) + nn * (i - 1)) * Bt(4,j) 
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
                + work2(i + nn * (2 - 1)) * Ct(2, j) &
                + work2(i + nn * (3 - 1)) * Ct(3, j) &
                + work2(i + nn * (4 - 1)) * Ct(4, j) 
       end do
    end do
    
  end subroutine tnsr3d_el_n4_cpu

  subroutine tnsr3d_el_n3_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 3
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, k, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
                   + A(i,2) * u(2 + n * (j - 1)) &
                   + A(i,3) * u(3 + n * (j - 1)) 
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                       + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) &
                       + work(l + n * (3 - 1) + nn * (i - 1)) * Bt(3,j) 
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
                + work2(i + nn * (2 - 1)) * Ct(2, j) &
                + work2(i + nn * (3 - 1)) * Ct(3, j) 
       end do
    end do
    
  end subroutine tnsr3d_el_n3_cpu

  subroutine tnsr3d_el_n2_cpu(v, u, A, Bt, Ct)
    integer, parameter :: n = 2
    integer, parameter :: nn = n**2
    real(kind=rp), intent(inout) :: v(n*n*n), u(n*n*n)
    real(kind=rp), intent(inout) :: A(n,n), Bt(n,n), Ct(n,n)
    real(kind=rp) :: work(n**3), work2(n**3)
    integer :: i, j, k, l
    integer :: ii, jj

    do j = 1, nn
       do i = 1, n
          ii = i + n * (j - 1)
          work(ii) = A(i,1) * u(1 + n * (j - 1)) &
                   + A(i,2) * u(2 + n * (j - 1)) 
       end do
    end do
    
    do i = 1, n
       do j = 1, n
          do l = 1, n
             ii = l + n * (j - 1) + nn * (i - 1)
             work2(ii) = work(l + n * (1 - 1) + nn * (i - 1)) * Bt(1,j) &
                       + work(l + n * (2 - 1) + nn * (i - 1)) * Bt(2,j) 
          end do
       end do
    end do
     
    do j = 1, n
       do i = 1, nn
          jj = i + nn * (j - 1)
          v(jj) = work2(i + nn * (1 - 1)) * Ct(1, j) &
                + work2(i + nn * (2 - 1)) * Ct(2, j) 
       end do
    end do
    
  end subroutine tnsr3d_el_n2_cpu
  
  subroutine tnsr3d_cpu(v, nv, u, nu, A, Bt, Ct, nelv)
    integer, intent(in) :: nv, nu, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv,nelv), u(nu*nu*nu,nelv)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu, nv), Ct(nu,nv)

    if (nu .eq. 2 .and. nv .eq. 4) then
       call tnsr3d_nu2nv4_cpu(v, u, A, Bt, Ct, nelv)
    else if (nu .eq. 4) then
       call tnsr3d_nu4_cpu(v, nv, u, A, Bt, Ct, nelv)
    else
       call tnsr3d_nvnu_cpu(v, nv, u, nu, A, Bt, Ct, nelv)
    end if

  end subroutine tnsr3d_cpu
  
  subroutine tnsr3d_nvnu_cpu(v, nv, u, nu, A, Bt, Ct, nelv)
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
    
  end subroutine tnsr3d_nvnu_cpu

  subroutine tnsr3d_nu2nv4_cpu(v, u, A, Bt, Ct, nelv)
    integer, parameter :: nu = 2
    integer, parameter :: nv = 4
    integer, parameter :: nunu = 4
    integer, parameter :: nvnu = 8
    integer, parameter :: nvnv = 16
    integer, intent(in) :: nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv,nelv), u(nu*nu*nu,nelv)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu, nv), Ct(nu,nv)
    real(kind=rp) :: work(nu**2*nv), work2(nu*nv**2), tmp
    integer :: ie, i, j, k, l, ii, jj
    
    do ie = 1,nelv
       do j = 1, nunu
          do i = 1, nv
             ii = i + nv * (j - 1)
             work(ii) = A(i,1) * u(1 + nu * (j - 1), ie) &
                      + A(i,2) * u(2 + nu * (j - 1), ie)             
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
             v(jj, ie) = work2(i + nvnv * (1 - 1)) * Ct(1, j) &
                       + work2(i + nvnv * (2 - 1)) * Ct(2, j) 
          end do
       end do
    end do
    
  end subroutine tnsr3d_nu2nv4_cpu

  subroutine tnsr3d_nu4_cpu(v, nv, u, A, Bt, Ct, nelv)
    integer, parameter :: nu = 4
    integer, parameter :: nunu = 16
    integer, intent(in) :: nv, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv,nelv), u(nu*nu*nu,nelv)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu, nv), Ct(nu,nv)
    real(kind=rp) :: work(nu**2*nv), work2(nu*nv**2), tmp
    integer :: ie, i, j, k, l, ii, jj
    integer :: nvnu, nvnv

    nvnu = nv * nu
    nvnv = nv * nv
    
    do ie = 1,nelv
       do j = 1, nunu
          do i = 1, nv
             ii = i + nv * (j - 1)
             work(ii) = A(i,1) * u(1 + nu * (j - 1), ie) &
                      + A(i,2) * u(2 + nu * (j - 1), ie) &
                      + A(i,3) * u(3 + nu * (j - 1), ie) &
                      + A(i,4) * u(4 + nu * (j - 1), ie)             
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
             v(jj, ie) = work2(i + nvnv * (1 - 1)) * Ct(1, j) &
                       + work2(i + nvnv * (2 - 1)) * Ct(2, j) &
                       + work2(i + nvnv * (3 - 1)) * Ct(3, j) &
                       + work2(i + nvnv * (4 - 1)) * Ct(4, j) 
          end do
       end do
    end do

  end subroutine tnsr3d_nu4_cpu

  subroutine tnsr1_3d_cpu(v, nv, nu, A, Bt, Ct, nelv)
    integer, intent(in) :: nv, nu, nelv
    real(kind=rp), intent(inout) :: v(nv*nv*nv*nelv)
    real(kind=rp), intent(inout) :: A(nv,nu), Bt(nu, nv), Ct(nu,nv)

    if (nu .eq. 4 .and. nv .eq. 2) then
       call tnsr1_3d_nu4nv2_cpu(v, A, Bt, Ct, nelv)
    else
       call tnsr1_3d_nvnu_cpu(v, nv, nu, A, Bt, Ct, nelv)
    end if
    
  end subroutine tnsr1_3d_cpu

    subroutine tnsr1_3d_nvnu_cpu(v, nv, nu, A, Bt, Ct, nelv)
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
    
  end subroutine tnsr1_3d_nvnu_cpu

  subroutine tnsr1_3d_nu4nv2_cpu(v, A, Bt, Ct, nelv)
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
    real(kind=rp) :: work(nu**2*nv), work2(nu*nv**2)
    integer :: e, e0, ee, es, iu, iv
    integer :: i, j, k, l, ii, jj, kk
    real(kind=rp) :: tmp

    do e = 1,nelv
       iu = (e-1)*nununu
       iv = (e-1)*nvnvnv
          
       do j = 1, nunu
          do i = 1, nv
             ii = i + nv * (j - 1)
             work(ii) = A(i,1) * v(1 + nu * (j - 1) + iu) &
                      + A(i,2) * v(2 + nu * (j - 1) + iu) &
                      + A(i,3) * v(3 + nu * (j - 1) + iu) &
                      + A(i,4) * v(4 + nu * (j - 1) + iu) 
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
             v(jj) = work2(i + nvnv * (1 - 1)) * Ct(1, j) &
                   + work2(i + nvnv * (2 - 1)) * Ct(2, j) &
                   + work2(i + nvnv * (3 - 1)) * Ct(3, j) &
                   + work2(i + nvnv * (4 - 1)) * Ct(4, j)
                  
          end do
       end do
    end do
    
  end subroutine tnsr1_3d_nu4nv2_cpu
  
end module tensor_cpu
