! Copyright (c) 2021-2023, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> DT*X kernels
module cpu_cdtp
  use num_types
  use math
  implicit none

contains
  
  subroutine cpu_cdtp_lx(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel, lx)
    integer, intent(in) :: nel, lx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in) :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l
    !$omp do
    do e = 1, nel
       
       do i = 1, lx*lx*lx
          wx(i,1,1) = ( B(i,1,1,e) * x(i,1,1,e) ) / jac(i,1,1,e)
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dr(i,1,1,e)
       end do
       
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             !DIR$ LOOP_INFO MIN_TRIPS(15)             
             do k = 1, lx
                tmp = tmp + dxt(i,k) * ta1(k,j,1)
             end do
             dtx(i,j,1,e) = tmp
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * ds(i,1,1,e)
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                !DIR$ LOOP_INFO MIN_TRIPS(15)                
                do l = 1, lx
                   tmp = tmp + dyt(j,l) * ta1(i,l,k)                   
                end do
                dtx(i,j,k,e) = dtx(i,j,k,e) + tmp
             end do
          end do
       end do
       
       do i = 1, lx*lx*lx
          ta1(i,1,1) = wx(i,1,1) * dt(i,1,1,e)
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             !DIR$ LOOP_INFO MIN_TRIPS(15)             
             do l = 1, lx
                tmp = tmp + dzt(k,l) * ta1(i,1,l)
             end do
             dtx(i,1,k,e) = dtx(i,1,k,e) + tmp
          end do
       end do
       
    end do
    !$omp end do
  end subroutine cpu_cdtp_lx
  
  subroutine cpu_cdtp_lx14(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 14
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in) :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k
    !$omp do
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
                          + dxt(i,12) * ta1(12,j,1) &
                          + dxt(i,13) * ta1(13,j,1) &
                          + dxt(i,14) * ta1(14,j,1) 
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
                             + dyt(j,12) * ta1(i,12,k) &
                             + dyt(j,13) * ta1(i,13,k) &
                             + dyt(j,14) * ta1(i,14,k) 
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
                          + dzt(k,12) * ta1(i,1,12) &
                          + dzt(k,13) * ta1(i,1,13) &
                          + dzt(k,14) * ta1(i,1,14) 
          end do
       end do
       
    end do
    !$omp end do
  end subroutine cpu_cdtp_lx14
  
  subroutine cpu_cdtp_lx13(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 13
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in) :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k
    !$omp do
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
                          + dxt(i,12) * ta1(12,j,1) &
                          + dxt(i,13) * ta1(13,j,1) 
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
                             + dyt(j,12) * ta1(i,12,k) &
                             + dyt(j,13) * ta1(i,13,k)
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
                          + dzt(k,12) * ta1(i,1,12) &
                          + dzt(k,13) * ta1(i,1,13) 
          end do
       end do
       
    end do
    !$omp end do
  end subroutine cpu_cdtp_lx13
  
  subroutine cpu_cdtp_lx12(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 12
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in) :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k
    !$omp do
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
    !$omp end do
  end subroutine cpu_cdtp_lx12
  
  subroutine cpu_cdtp_lx11(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 11
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in) :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k
    !$omp do 
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
    !$omp end do
  end subroutine cpu_cdtp_lx11
  
  subroutine cpu_cdtp_lx10(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 10
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in) :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k
    !$omp do
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
    !$omp end do
  end subroutine cpu_cdtp_lx10

  subroutine cpu_cdtp_lx9(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 9
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in) :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k
    !$omp do
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
    !$omp end do
  end subroutine cpu_cdtp_lx9

  subroutine cpu_cdtp_lx8(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 8
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in) :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k
    !$omp do
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
    !$omp end do
  end subroutine cpu_cdtp_lx8

  subroutine cpu_cdtp_lx7(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 7
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in) :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k
    !$omp do
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
    !$omp end do
  end subroutine cpu_cdtp_lx7
  
  subroutine cpu_cdtp_lx6(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 6
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in) :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k
    !$omp do
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
    !$omp end do
  end subroutine cpu_cdtp_lx6

  subroutine cpu_cdtp_lx5(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 5
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in) :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k
    !$omp do
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
    !$omp end do
  end subroutine cpu_cdtp_lx5

  subroutine cpu_cdtp_lx4(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 4
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in) :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k
    !$omp do
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
    !$omp end do
  end subroutine cpu_cdtp_lx4

  subroutine cpu_cdtp_lx3(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 3
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in) :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k
    !$omp do
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
    !$omp end do
  end subroutine cpu_cdtp_lx3
  
  subroutine cpu_cdtp_lx2(dtx, x, dr, ds, dt, dxt, dyt, dzt, B, jac, nel)
    integer, parameter :: lx = 2
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: x, dr, ds, dt, jac, B
    real(kind=rp), intent(in) :: dxt(lx,lx), dyt(lx,lx), dzt(lx,lx)   
    real(kind=rp), dimension(lx,lx,lx) :: wx, ta1
    integer :: e, i, j, k
    !$omp do
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
    !$omp end do
  end subroutine cpu_cdtp_lx2

end module cpu_cdtp
