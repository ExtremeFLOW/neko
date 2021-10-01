!> DT*X kernels
module cpu_cdtp
  use num_types
  use math
  implicit none

contains

  subroutine cpu_cdtp_lx12(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 12
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k
    
    do e = 1, nel
       
       do i = 1, lx*lx*lx
          wx(i,1,1) = ( B(i,1,1,e) * x(i,1,1,e) ) / jac(i,1,1,e)
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dr(i,1,1,e)
       end do
       
       do j = 1, lx * lx
          do i = 1, lx
             dtx(i,j,1,e) = dxt(i,1) * ta1(1,j,1) &
                          + dxt(i,2) * ta1(2,j,1) &
                          + dxt(i,3) * ta1(3,j,1) &
                          + dxt(i,4) * ta1(4,j,1) &
                          + dxt(i,5) * ta1(5,j,1) &
                          + dxt(i,6) * ta1(6,j,1) &
                          + dxt(i,7) * ta1(7,j,1) &
                          + dxt(i,8) * ta1(8,j,1) &
                          + dxt(i,9) * ta1(9,j,1) &
                          + dxt(i,10) * ta1(10,j,1) &
                          + dxt(i,11) * ta1(11,j,1) &
                          + dxt(i,12) * ta1(12,j,1) 
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * ds(i,1,1,e)
       end do
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                dtx(i,j,k,e) = dtx(i,j,k,e) &
                             + dyt(j,1) * ta1(i,1,k) &
                             + dyt(j,2) * ta1(i,2,k) &
                             + dyt(j,3) * ta1(i,3,k) &
                             + dyt(j,4) * ta1(i,4,k) &
                             + dyt(j,5) * ta1(i,5,k) &
                             + dyt(j,6) * ta1(i,6,k) &
                             + dyt(j,7) * ta1(i,7,k) &
                             + dyt(j,8) * ta1(i,8,k) &
                             + dyt(j,9) * ta1(i,9,k) &
                             + dyt(j,10) * ta1(i,10,k) &
                             + dyt(j,11) * ta1(i,11,k) &
                             + dyt(j,12) * ta1(i,12,k)
             end do
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dt(i,1,1,e)
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             dtx(i,1,k,e) = dtx(i,1,k,e) &
                          + dzt(k,1) * ta1(i,1,1) &
                          + dzt(k,2) * ta1(i,1,2) &
                          + dzt(k,3) * ta1(i,1,3) &
                          + dzt(k,4) * ta1(i,1,4) &
                          + dzt(k,5) * ta1(i,1,5) &
                          + dzt(k,6) * ta1(i,1,6) &
                          + dzt(k,7) * ta1(i,1,7) &
                          + dzt(k,8) * ta1(i,1,8) &
                          + dzt(k,9) * ta1(i,1,9) &
                          + dzt(k,10) * ta1(i,1,10) &
                          + dzt(k,11) * ta1(i,1,11) &
                          + dzt(k,12) * ta1(i,1,12)
          end do
       end do
       
    end do
  end subroutine cpu_cdtp_lx12
  
  subroutine cpu_cdtp_lx11(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 11
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k
    
    do e = 1, nel
       
       do i = 1, lx*lx*lx
          wx(i,1,1) = ( B(i,1,1,e) * x(i,1,1,e) ) / jac(i,1,1,e)
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dr(i,1,1,e)
       end do
       
       do j = 1, lx * lx
          do i = 1, lx
             dtx(i,j,1,e) = dxt(i,1) * ta1(1,j,1) &
                          + dxt(i,2) * ta1(2,j,1) &
                          + dxt(i,3) * ta1(3,j,1) &
                          + dxt(i,4) * ta1(4,j,1) &
                          + dxt(i,5) * ta1(5,j,1) &
                          + dxt(i,6) * ta1(6,j,1) &
                          + dxt(i,7) * ta1(7,j,1) &
                          + dxt(i,8) * ta1(8,j,1) &
                          + dxt(i,9) * ta1(9,j,1) &
                          + dxt(i,10) * ta1(10,j,1) &
                          + dxt(i,11) * ta1(11,j,1) 
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * ds(i,1,1,e)
       end do
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                dtx(i,j,k,e) = dtx(i,j,k,e) &
                             + dyt(j,1) * ta1(i,1,k) &
                             + dyt(j,2) * ta1(i,2,k) &
                             + dyt(j,3) * ta1(i,3,k) &
                             + dyt(j,4) * ta1(i,4,k) &
                             + dyt(j,5) * ta1(i,5,k) &
                             + dyt(j,6) * ta1(i,6,k) &
                             + dyt(j,7) * ta1(i,7,k) &
                             + dyt(j,8) * ta1(i,8,k) &
                             + dyt(j,9) * ta1(i,9,k) &
                             + dyt(j,10) * ta1(i,10,k) &
                             + dyt(j,11) * ta1(i,11,k) 
             end do
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dt(i,1,1,e)
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             dtx(i,1,k,e) = dtx(i,1,k,e) &
                          + dzt(k,1) * ta1(i,1,1) &
                          + dzt(k,2) * ta1(i,1,2) &
                          + dzt(k,3) * ta1(i,1,3) &
                          + dzt(k,4) * ta1(i,1,4) &
                          + dzt(k,5) * ta1(i,1,5) &
                          + dzt(k,6) * ta1(i,1,6) &
                          + dzt(k,7) * ta1(i,1,7) &
                          + dzt(k,8) * ta1(i,1,8) &
                          + dzt(k,9) * ta1(i,1,9) &
                          + dzt(k,10) * ta1(i,1,10) &
                          + dzt(k,11) * ta1(i,1,11) 
          end do
       end do
       
    end do
  end subroutine cpu_cdtp_lx11
  
  subroutine cpu_cdtp_lx10(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 10
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k
    
    do e = 1, nel

       do i = 1, lx*lx*lx
          wx(i,1,1) = ( B(i,1,1,e) * x(i,1,1,e) ) / jac(i,1,1,e)
       end do

       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dr(i,1,1,e)
       end do

       do j = 1, lx * lx
          do i = 1, lx
             dtx(i,j,1,e) = dxt(i,1) * ta1(1,j,1) &
                          + dxt(i,2) * ta1(2,j,1) &
                          + dxt(i,3) * ta1(3,j,1) &
                          + dxt(i,4) * ta1(4,j,1) &
                          + dxt(i,5) * ta1(5,j,1) &
                          + dxt(i,6) * ta1(6,j,1) &
                          + dxt(i,7) * ta1(7,j,1) &
                          + dxt(i,8) * ta1(8,j,1) &
                          + dxt(i,9) * ta1(9,j,1) &
                          + dxt(i,10) * ta1(10,j,1) 
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * ds(i,1,1,e)
       end do
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                dtx(i,j,k,e) = dtx(i,j,k,e) &
                             + dyt(j,1) * ta1(i,1,k) &
                             + dyt(j,2) * ta1(i,2,k) &
                             + dyt(j,3) * ta1(i,3,k) &
                             + dyt(j,4) * ta1(i,4,k) &
                             + dyt(j,5) * ta1(i,5,k) &
                             + dyt(j,6) * ta1(i,6,k) &
                             + dyt(j,7) * ta1(i,7,k) &
                             + dyt(j,8) * ta1(i,8,k) &
                             + dyt(j,9) * ta1(i,9,k) &
                             + dyt(j,10) * ta1(i,10,k) 

             end do
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dt(i,1,1,e)
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             dtx(i,1,k,e) = dtx(i,1,k,e) &
                          + dzt(k,1) * ta1(i,1,1) &
                          + dzt(k,2) * ta1(i,1,2) &
                          + dzt(k,3) * ta1(i,1,3) &
                          + dzt(k,4) * ta1(i,1,4) &
                          + dzt(k,5) * ta1(i,1,5) &
                          + dzt(k,6) * ta1(i,1,6) &
                          + dzt(k,7) * ta1(i,1,7) &
                          + dzt(k,8) * ta1(i,1,8) &
                          + dzt(k,9) * ta1(i,1,9) &
                          + dzt(k,10) * ta1(i,1,10) 
          end do
       end do
       
    end do
  end subroutine cpu_cdtp_lx10

  subroutine cpu_cdtp_lx9(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 9
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k
    
    do e = 1, nel
       
       do i = 1, lx*lx*lx
          wx(i,1,1) = ( B(i,1,1,e) * x(i,1,1,e) ) / jac(i,1,1,e)
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dr(i,1,1,e)
       end do
       
       do j = 1, lx * lx
          do i = 1, lx
             dtx(i,j,1,e) = dxt(i,1) * ta1(1,j,1) &
                          + dxt(i,2) * ta1(2,j,1) &
                          + dxt(i,3) * ta1(3,j,1) &
                          + dxt(i,4) * ta1(4,j,1) &
                          + dxt(i,5) * ta1(5,j,1) &
                          + dxt(i,6) * ta1(6,j,1) &
                          + dxt(i,7) * ta1(7,j,1) &
                          + dxt(i,8) * ta1(8,j,1) &
                          + dxt(i,9) * ta1(9,j,1) 
          end do
       end do

       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * ds(i,1,1,e)
       end do
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                dtx(i,j,k,e) = dtx(i,j,k,e) &
                             + dyt(j,1) * ta1(i,1,k) &
                             + dyt(j,2) * ta1(i,2,k) &
                             + dyt(j,3) * ta1(i,3,k) &
                             + dyt(j,4) * ta1(i,4,k) &
                             + dyt(j,5) * ta1(i,5,k) &
                             + dyt(j,6) * ta1(i,6,k) &
                             + dyt(j,7) * ta1(i,7,k) &
                             + dyt(j,8) * ta1(i,8,k) &
                             + dyt(j,9) * ta1(i,9,k) 
             end do
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dt(i,1,1,e)
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             dtx(i,1,k,e) = dtx(i,1,k,e) &
                          + dzt(k,1) * ta1(i,1,1) &
                          + dzt(k,2) * ta1(i,1,2) &
                          + dzt(k,3) * ta1(i,1,3) &
                          + dzt(k,4) * ta1(i,1,4) &
                          + dzt(k,5) * ta1(i,1,5) &
                          + dzt(k,6) * ta1(i,1,6) &
                          + dzt(k,7) * ta1(i,1,7) &
                          + dzt(k,8) * ta1(i,1,8) &
                          + dzt(k,9) * ta1(i,1,9) 
          end do
       end do
       
    end do
  end subroutine cpu_cdtp_lx9

  subroutine cpu_cdtp_lx8(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 8
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k
    
    do e = 1, nel
       
       do i = 1, lx*lx*lx
          wx(i,1,1) = ( B(i,1,1,e) * x(i,1,1,e) ) / jac(i,1,1,e)
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dr(i,1,1,e)
       end do
       
       do j = 1, lx * lx
          do i = 1, lx
             dtx(i,j,1,e) = dxt(i,1) * ta1(1,j,1) &
                          + dxt(i,2) * ta1(2,j,1) &
                          + dxt(i,3) * ta1(3,j,1) &
                          + dxt(i,4) * ta1(4,j,1) &
                          + dxt(i,5) * ta1(5,j,1) &
                          + dxt(i,6) * ta1(6,j,1) &
                          + dxt(i,7) * ta1(7,j,1) &
                          + dxt(i,8) * ta1(8,j,1) 
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * ds(i,1,1,e)
       end do
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                dtx(i,j,k,e) = dtx(i,j,k,e) &
                             + dyt(j,1) * ta1(i,1,k) &
                             + dyt(j,2) * ta1(i,2,k) &
                             + dyt(j,3) * ta1(i,3,k) &
                             + dyt(j,4) * ta1(i,4,k) &
                             + dyt(j,5) * ta1(i,5,k) &
                             + dyt(j,6) * ta1(i,6,k) &
                             + dyt(j,7) * ta1(i,7,k) &
                             + dyt(j,8) * ta1(i,8,k)
             end do
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dt(i,1,1,e)
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             dtx(i,1,k,e) = dtx(i,1,k,e) &
                         + dzt(k,1) * ta1(i,1,1) &
                         + dzt(k,2) * ta1(i,1,2) &
                         + dzt(k,3) * ta1(i,1,3) &
                         + dzt(k,4) * ta1(i,1,4) &
                         + dzt(k,5) * ta1(i,1,5) &
                         + dzt(k,6) * ta1(i,1,6) &
                         + dzt(k,7) * ta1(i,1,7) &
                         + dzt(k,8) * ta1(i,1,8)
          end do
       end do
       
    end do
  end subroutine cpu_cdtp_lx8

  subroutine cpu_cdtp_lx7(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 7
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k
    
    do e = 1, nel
       
       do i = 1, lx*lx*lx
          wx(i,1,1) = ( B(i,1,1,e) * x(i,1,1,e) ) / jac(i,1,1,e)
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dr(i,1,1,e)
       end do

       do j = 1, lx * lx
          do i = 1, lx
             dtx(i,j,1,e) = dxt(i,1) * ta1(1,j,1) &
                          + dxt(i,2) * ta1(2,j,1) &
                          + dxt(i,3) * ta1(3,j,1) &
                          + dxt(i,4) * ta1(4,j,1) &
                          + dxt(i,5) * ta1(5,j,1) &
                          + dxt(i,6) * ta1(6,j,1) &
                          + dxt(i,7) * ta1(7,j,1) 
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * ds(i,1,1,e)
       end do
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                dtx(i,j,k,e) = dtx(i,j,k,e) &
                             + dyt(j,1) * ta1(i,1,k) &
                             + dyt(j,2) * ta1(i,2,k) &
                             + dyt(j,3) * ta1(i,3,k) &
                             + dyt(j,4) * ta1(i,4,k) &
                             + dyt(j,5) * ta1(i,5,k) &
                             + dyt(j,6) * ta1(i,6,k) &
                             + dyt(j,7) * ta1(i,7,k) 
             end do
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dt(i,1,1,e)
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             dtx(i,1,k,e) = dtx(i,1,k,e) &
                          + dzt(k,1) * ta1(i,1,1) &
                          + dzt(k,2) * ta1(i,1,2) &
                          + dzt(k,3) * ta1(i,1,3) &
                          + dzt(k,4) * ta1(i,1,4) &
                          + dzt(k,5) * ta1(i,1,5) &
                          + dzt(k,6) * ta1(i,1,6) &
                          + dzt(k,7) * ta1(i,1,7) 
          end do
       end do
       
    end do
  end subroutine cpu_cdtp_lx7
  
  subroutine cpu_cdtp_lx6(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 6
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k

    do e = 1, nel
       
       do i = 1, lx*lx*lx
          wx(i,1,1) = ( B(i,1,1,e) * x(i,1,1,e) ) / jac(i,1,1,e)
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dr(i,1,1,e)
       end do
       
       do j = 1, lx * lx
          do i = 1, lx
             dtx(i,j,1,e) = dxt(i,1) * ta1(1,j,1) &
                          + dxt(i,2) * ta1(2,j,1) &
                          + dxt(i,3) * ta1(3,j,1) &
                          + dxt(i,4) * ta1(4,j,1) &
                          + dxt(i,5) * ta1(5,j,1) &
                          + dxt(i,6) * ta1(6,j,1) 
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * ds(i,1,1,e)
       end do
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                dtx(i,j,k,e) = dtx(i,j,k,e) &
                             + dyt(j,1) * ta1(i,1,k) &
                             + dyt(j,2) * ta1(i,2,k) &
                             + dyt(j,3) * ta1(i,3,k) &
                             + dyt(j,4) * ta1(i,4,k) &
                             + dyt(j,5) * ta1(i,5,k) &
                             + dyt(j,6) * ta1(i,6,k) 
             end do
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dt(i,1,1,e)
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             dtx(i,1,k,e) = dtx(i,1,k,e) &
                          + dzt(k,1) * ta1(i,1,1) &
                          + dzt(k,2) * ta1(i,1,2) &
                          + dzt(k,3) * ta1(i,1,3) &
                          + dzt(k,4) * ta1(i,1,4) &
                          + dzt(k,5) * ta1(i,1,5) &
                          + dzt(k,6) * ta1(i,1,6)              
          end do
       end do
       
    end do
  end subroutine cpu_cdtp_lx6

  subroutine cpu_cdtp_lx5(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 5
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k
    
    do e = 1, nel
       
       do i = 1, lx*lx*lx
          wx(i,1,1) = ( B(i,1,1,e) * x(i,1,1,e) ) / jac(i,1,1,e)
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dr(i,1,1,e)
       end do
       
       do j = 1, lx * lx
          do i = 1, lx
             dtx(i,j,1,e) = dxt(i,1) * ta1(1,j,1) &
                          + dxt(i,2) * ta1(2,j,1) &
                          + dxt(i,3) * ta1(3,j,1) &
                          + dxt(i,4) * ta1(4,j,1) &
                          + dxt(i,5) * ta1(5,j,1) 
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * ds(i,1,1,e)
       end do
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                dtx(i,j,k,e) = dtx(i,j,k,e) &
                             + dyt(j,1) * ta1(i,1,k) &
                             + dyt(j,2) * ta1(i,2,k) &
                             + dyt(j,3) * ta1(i,3,k) &
                             + dyt(j,4) * ta1(i,4,k) &
                             + dyt(j,5) * ta1(i,5,k) 
             end do
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dt(i,1,1,e)
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             dtx(i,1,k,e) = dtx(i,1,k,e) &
                          + dzt(k,1) * ta1(i,1,1) &
                          + dzt(k,2) * ta1(i,1,2) &
                          + dzt(k,3) * ta1(i,1,3) &
                          + dzt(k,4) * ta1(i,1,4) &
                          + dzt(k,5) * ta1(i,1,5) 
          end do
       end do
       
    end do
  end subroutine cpu_cdtp_lx5

  subroutine cpu_cdtp_lx4(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 4
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k

    do e = 1, nel
       
       do i = 1, lx*lx*lx
          wx(i,1,1) = ( B(i,1,1,e) * x(i,1,1,e) ) / jac(i,1,1,e)
       end do

       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dr(i,1,1,e)
       end do

       do j = 1, lx * lx
          do i = 1, lx
             dtx(i,j,1,e) = dxt(i,1) * ta1(1,j,1) &
                          + dxt(i,2) * ta1(2,j,1) &
                          + dxt(i,3) * ta1(3,j,1) &
                          + dxt(i,4) * ta1(4,j,1) 
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * ds(i,1,1,e)
       end do
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                dtx(i,j,k,e) = dtx(i,j,k,e) &
                             + dyt(j,1) * ta1(i,1,k) &
                             + dyt(j,2) * ta1(i,2,k) &
                             + dyt(j,3) * ta1(i,3,k) &
                             + dyt(j,4) * ta1(i,4,k) 
             end do
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dt(i,1,1,e)
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             dtx(i,1,k,e) = dtx(i,1,k,e) &
                          + dzt(k,1) * ta1(i,1,1) &
                          + dzt(k,2) * ta1(i,1,2) &
                          + dzt(k,3) * ta1(i,1,3) &
                          + dzt(k,4) * ta1(i,1,4) 
          end do
       end do
       
    end do
  end subroutine cpu_cdtp_lx4

  subroutine cpu_cdtp_lx3(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 3
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k

    do e = 1, nel
       
       do i = 1, lx*lx*lx
          wx(i,1,1) = ( B(i,1,1,e) * x(i,1,1,e) ) / jac(i,1,1,e)
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dr(i,1,1,e)
       end do

       do j = 1, lx * lx
          do i = 1, lx
             dtx(i,j,1,e) = dxt(i,1) * ta1(1,j,1) &
                          + dxt(i,2) * ta1(2,j,1) &
                          + dxt(i,3) * ta1(3,j,1)
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * ds(i,1,1,e)
       end do
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                dtx(i,j,k,e) = dtx(i,j,k,e) &
                             + dyt(j,1) * ta1(i,1,k) &
                             + dyt(j,2) * ta1(i,2,k) &
                             + dyt(j,3) * ta1(i,3,k) 
             end do
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dt(i,1,1,e)
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             dtx(i,1,k,e) = dtx(i,1,k,e) &
                          + dzt(k,1) * ta1(i,1,1) &
                          + dzt(k,2) * ta1(i,1,2) &
                          + dzt(k,3) * ta1(i,1,3) 
          end do
       end do
       
    end do
  end subroutine cpu_cdtp_lx3
  
  subroutine cpu_cdtp_lx2(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 2
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in)  :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k

    do e = 1, nel

       do i = 1, lx*lx*lx
          wx(i,1,1) = ( B(i,1,1,e) * x(i,1,1,e) ) / jac(i,1,1,e)
       end do

       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dr(i,1,1,e)
       end do

       do j = 1, lx * lx
          do i = 1, lx
             dtx(i,j,1,e) = dxt(i,1) * ta1(1,j,1) &
                          + dxt(i,2) * ta1(2,j,1) 
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * ds(i,1,1,e)
       end do
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                dtx(i,j,k,e) = dtx(i,j,k,e) &
                             + dyt(j,1) * ta1(i,1,k) &
                             + dyt(j,2) * ta1(i,2,k) 
             end do
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dt(i,1,1,e)
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             dtx(i,1,k,e) = dtx(i,1,k,e) &
                          + dzt(k,1) * ta1(i,1,1) &
                          + dzt(k,2) * ta1(i,1,2) 
          end do
       end do
       
    end do
  end subroutine cpu_cdtp_lx2

end module cpu_cdtp
