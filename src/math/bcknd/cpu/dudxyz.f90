! Copyright (c) 2021, The Neko Authors
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
!> Derivative kernels
module cpu_dudxyz
  use num_types, only : rp
  implicit none

contains

  subroutine cpu_dudxyz_lx(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, lx)
    integer, intent(in) :: nel, lx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: drst
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l
    !$omp do
    do e = 1, nel
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             du(i,j,1,e) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                drst(i,j,k) =  tmp
             end do
          end do
       end do
        
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             drst(i,1,k) = tmp
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
       end do
        
    end do
    !$omp end do
  end subroutine cpu_dudxyz_lx
  
  subroutine cpu_dudxyz_lx14(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
    integer, parameter :: lx = 14
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: drst
    integer :: e, i, j, k
    !$omp do
    do e = 1, nel
       do j = 1, lx * lx
          do i = 1, lx
             du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e) &
                         + dx(i,6) * u(6,j,1,e) &
                         + dx(i,7) * u(7,j,1,e) &
                         + dx(i,8) * u(8,j,1,e) &
                         + dx(i,9) * u(9,j,1,e) &
                         + dx(i,10) * u(10,j,1,e) &
                         + dx(i,11) * u(11,j,1,e) &
                         + dx(i,12) * u(12,j,1,e) &
                         + dx(i,13) * u(13,j,1,e) &
                         + dx(i,14) * u(14,j,1,e) 
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e) &
                            + dy(j,6) * u(i,6,k,e) &
                            + dy(j,7) * u(i,7,k,e) &
                            + dy(j,8) * u(i,8,k,e) &
                            + dy(j,9) * u(i,9,k,e) &
                            + dy(j,10) * u(i,10,k,e) &
                            + dy(j,11) * u(i,11,k,e) &
                            + dy(j,12) * u(i,12,k,e) &
                            + dy(j,13) * u(i,13,k,e) &
                            + dy(j,14) * u(i,14,k,e) 
             end do
          end do
       end do
        
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) &
                         + dz(k,6) * u(i,1,6,e) &
                         + dz(k,7) * u(i,1,7,e) &
                         + dz(k,8) * u(i,1,8,e) &
                         + dz(k,9) * u(i,1,9,e) &
                         + dz(k,10) * u(i,1,10,e) &
                         + dz(k,11) * u(i,1,11,e) &
                         + dz(k,12) * u(i,1,12,e) &
                         + dz(k,13) * u(i,1,13,e) &
                         + dz(k,14) * u(i,1,14,e) 
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
       end do
        
    end do
    !$omp end do
  end subroutine cpu_dudxyz_lx14
  
  subroutine cpu_dudxyz_lx13(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
    integer, parameter :: lx = 13
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: drst
    integer :: e, i, j, k
    !$omp do
    do e = 1, nel
       do j = 1, lx * lx
          do i = 1, lx
             du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e) &
                         + dx(i,6) * u(6,j,1,e) &
                         + dx(i,7) * u(7,j,1,e) &
                         + dx(i,8) * u(8,j,1,e) &
                         + dx(i,9) * u(9,j,1,e) &
                         + dx(i,10) * u(10,j,1,e) &
                         + dx(i,11) * u(11,j,1,e) &
                         + dx(i,12) * u(12,j,1,e) &
                         + dx(i,13) * u(13,j,1,e) 
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e) &
                            + dy(j,6) * u(i,6,k,e) &
                            + dy(j,7) * u(i,7,k,e) &
                            + dy(j,8) * u(i,8,k,e) &
                            + dy(j,9) * u(i,9,k,e) &
                            + dy(j,10) * u(i,10,k,e) &
                            + dy(j,11) * u(i,11,k,e) &
                            + dy(j,12) * u(i,12,k,e) &
                            + dy(j,13) * u(i,13,k,e) 
             end do
          end do
       end do
        
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) &
                         + dz(k,6) * u(i,1,6,e) &
                         + dz(k,7) * u(i,1,7,e) &
                         + dz(k,8) * u(i,1,8,e) &
                         + dz(k,9) * u(i,1,9,e) &
                         + dz(k,10) * u(i,1,10,e) &
                         + dz(k,11) * u(i,1,11,e) &
                         + dz(k,12) * u(i,1,12,e) &
                         + dz(k,13) * u(i,1,13,e) 
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
       end do
        
    end do
    !$omp end do
  end subroutine cpu_dudxyz_lx13
  
  subroutine cpu_dudxyz_lx12(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
    integer, parameter :: lx = 12
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: drst
    integer :: e, i, j, k
    !$omp do
    do e = 1, nel
       do j = 1, lx * lx
          do i = 1, lx
             du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e) &
                         + dx(i,6) * u(6,j,1,e) &
                         + dx(i,7) * u(7,j,1,e) &
                         + dx(i,8) * u(8,j,1,e) &
                         + dx(i,9) * u(9,j,1,e) &
                         + dx(i,10) * u(10,j,1,e) &
                         + dx(i,11) * u(11,j,1,e) &
                         + dx(i,12) * u(12,j,1,e) 
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e) &
                            + dy(j,6) * u(i,6,k,e) &
                            + dy(j,7) * u(i,7,k,e) &
                            + dy(j,8) * u(i,8,k,e) &
                            + dy(j,9) * u(i,9,k,e) &
                            + dy(j,10) * u(i,10,k,e) &
                            + dy(j,11) * u(i,11,k,e) &
                            + dy(j,12) * u(i,12,k,e)
             end do
          end do
       end do
        
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) &
                         + dz(k,6) * u(i,1,6,e) &
                         + dz(k,7) * u(i,1,7,e) &
                         + dz(k,8) * u(i,1,8,e) &
                         + dz(k,9) * u(i,1,9,e) &
                         + dz(k,10) * u(i,1,10,e) &
                         + dz(k,11) * u(i,1,11,e) &
                         + dz(k,12) * u(i,1,12,e)
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
       end do
        
    end do
    !$omp end do
  end subroutine cpu_dudxyz_lx12

  subroutine cpu_dudxyz_lx11(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
    integer, parameter :: lx = 11
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: drst
    integer :: e, i, j, k
    !$omp do
    do e = 1, nel
       do j = 1, lx * lx
          do i = 1, lx
             du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e) &
                         + dx(i,6) * u(6,j,1,e) &
                         + dx(i,7) * u(7,j,1,e) &
                         + dx(i,8) * u(8,j,1,e) &
                         + dx(i,9) * u(9,j,1,e) &
                         + dx(i,10) * u(10,j,1,e) &
                         + dx(i,11) * u(11,j,1,e) 
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
       end do
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e) &
                            + dy(j,6) * u(i,6,k,e) &
                            + dy(j,7) * u(i,7,k,e) &
                            + dy(j,8) * u(i,8,k,e) &
                            + dy(j,9) * u(i,9,k,e) &
                            + dy(j,10) * u(i,10,k,e) &
                            + dy(j,11) * u(i,11,k,e)
             end do
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
       end do

       do k = 1, lx
          do i = 1, lx*lx
             drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) &
                         + dz(k,6) * u(i,1,6,e) &
                         + dz(k,7) * u(i,1,7,e) &
                         + dz(k,8) * u(i,1,8,e) &
                         + dz(k,9) * u(i,1,9,e) &
                         + dz(k,10) * u(i,1,10,e) &
                         + dz(k,11) * u(i,1,11,e) 
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
       end do
        
    end do
    !$omp end do
  end subroutine cpu_dudxyz_lx11

  subroutine cpu_dudxyz_lx10(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
    integer, parameter :: lx = 10
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: drst
    integer :: e, i, j, k
    !$omp do
    do e = 1, nel
       do j = 1, lx * lx
          do i = 1, lx
             du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e) &
                         + dx(i,6) * u(6,j,1,e) &
                         + dx(i,7) * u(7,j,1,e) &
                         + dx(i,8) * u(8,j,1,e) &
                         + dx(i,9) * u(9,j,1,e) &
                         + dx(i,10) * u(10,j,1,e) 
          end do
       end do
        
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
       end do
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e) &
                            + dy(j,6) * u(i,6,k,e) &
                            + dy(j,7) * u(i,7,k,e) &
                            + dy(j,8) * u(i,8,k,e) &
                            + dy(j,9) * u(i,9,k,e) &
                            + dy(j,10) * u(i,10,k,e) 
             end do
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) &
                         + dz(k,6) * u(i,1,6,e) &
                         + dz(k,7) * u(i,1,7,e) &
                         + dz(k,8) * u(i,1,8,e) &
                         + dz(k,9) * u(i,1,9,e) &
                         + dz(k,10) * u(i,1,10,e)
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
       end do

    end do
    !$omp end do
  end subroutine cpu_dudxyz_lx10

  subroutine cpu_dudxyz_lx9(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
    integer, parameter :: lx = 9
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: drst
    integer :: e, i, j, k
    !$omp do
    do e = 1, nel
       do j = 1, lx * lx
          do i = 1, lx
             du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e) &
                         + dx(i,6) * u(6,j,1,e) &
                         + dx(i,7) * u(7,j,1,e) &
                         + dx(i,8) * u(8,j,1,e) &
                         + dx(i,9) * u(9,j,1,e) 
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e) &
                            + dy(j,6) * u(i,6,k,e) &
                            + dy(j,7) * u(i,7,k,e) &
                            + dy(j,8) * u(i,8,k,e) &
                            + dy(j,9) * u(i,9,k,e) 
             end do
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) &
                         + dz(k,6) * u(i,1,6,e) &
                         + dz(k,7) * u(i,1,7,e) &
                         + dz(k,8) * u(i,1,8,e) &
                         + dz(k,9) * u(i,1,9,e) 
          end do
       end do
        
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
       end do

    end do
    !$omp end do
  end subroutine cpu_dudxyz_lx9

  subroutine cpu_dudxyz_lx8(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
    integer, parameter :: lx = 8
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: drst
    integer :: e, i, j, k
    !$omp do
    do e = 1, nel
       do j = 1, lx * lx
          do i = 1, lx
             du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e) &
                         + dx(i,6) * u(6,j,1,e) &
                         + dx(i,7) * u(7,j,1,e) &
                         + dx(i,8) * u(8,j,1,e) 
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e) &
                            + dy(j,6) * u(i,6,k,e) &
                            + dy(j,7) * u(i,7,k,e) &
                            + dy(j,8) * u(i,8,k,e)
             end do
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
       end do

       do k = 1, lx
          do i = 1, lx*lx
             drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) &
                         + dz(k,6) * u(i,1,6,e) &
                         + dz(k,7) * u(i,1,7,e) &
                         + dz(k,8) * u(i,1,8,e) 
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
       end do
       
    end do
    !$omp end do
  end subroutine cpu_dudxyz_lx8

  subroutine cpu_dudxyz_lx7(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
    integer, parameter :: lx = 7
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: drst
    integer :: e, i, j, k
    !$omp do
    do e = 1, nel
       do j = 1, lx * lx
          do i = 1, lx
             du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e) &
                         + dx(i,6) * u(6,j,1,e) &
                         + dx(i,7) * u(7,j,1,e) 
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e) &
                            + dy(j,6) * u(i,6,k,e) &
                            + dy(j,7) * u(i,7,k,e)
             end do
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
       end do

       do k = 1, lx
          do i = 1, lx*lx
             drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) &
                         + dz(k,6) * u(i,1,6,e) &
                         + dz(k,7) * u(i,1,7,e)
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
       end do
       
    end do
    !$omp end do
  end subroutine cpu_dudxyz_lx7

  subroutine cpu_dudxyz_lx6(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
    integer, parameter :: lx = 6
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: drst
    integer :: e, i, j, k
    !$omp do
    do e = 1, nel
       do j = 1, lx * lx
          do i = 1, lx
             du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e) &
                         + dx(i,6) * u(6,j,1,e)
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e) &
                            + dy(j,6) * u(i,6,k,e)
             end do
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
       end do

       do k = 1, lx
          do i = 1, lx*lx
             drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e) &
                         + dz(k,6) * u(i,1,6,e)
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
       end do
       
    end do
    !$omp end do
  end subroutine cpu_dudxyz_lx6

  subroutine cpu_dudxyz_lx5(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
    integer, parameter :: lx = 5
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: drst
    integer :: e, i, j, k
    !$omp do
    do e = 1, nel
       do j = 1, lx * lx
          do i = 1, lx
             du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e) &
                         + dx(i,5) * u(5,j,1,e)
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e) &
                            + dy(j,5) * u(i,5,k,e)
             end do
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e) &
                         + dz(k,5) * u(i,1,5,e)
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
       end do

    end do
    !$omp end do
  end subroutine cpu_dudxyz_lx5

  subroutine cpu_dudxyz_lx4(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
    integer, parameter :: lx = 4
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: drst
    integer :: e, i, j, k
    !$omp do
    do e = 1, nel
       do j = 1, lx * lx
          do i = 1, lx
             du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e) &
                         + dx(i,4) * u(4,j,1,e)
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
       end do
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e) &
                            + dy(j,4) * u(i,4,k,e)
             end do
          end do
       end do
        
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
       end do

       do k = 1, lx
          do i = 1, lx*lx
             drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e) &
                         + dz(k,4) * u(i,1,4,e)
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
       end do

    end do
    !$omp end do
  end subroutine cpu_dudxyz_lx4

  subroutine cpu_dudxyz_lx3(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
    integer, parameter :: lx = 3
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: drst
    integer :: e, i, j, k
    !$omp do
    do e = 1, nel
       do j = 1, lx * lx
          do i = 1, lx
             du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) &
                         + dx(i,3) * u(3,j,1,e)
          end do
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e) &
                            + dy(j,3) * u(i,3,k,e)
             end do
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
       end do

       do k = 1, lx
          do i = 1, lx*lx
             drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) &
                         + dz(k,3) * u(i,1,3,e)
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
       end do

    end do
    !$omp end do
  end subroutine cpu_dudxyz_lx3

  subroutine cpu_dudxyz_lx2(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel)
    integer, parameter :: lx = 2
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: drst
    integer :: e, i, j, k
    !$omp do
    do e = 1, nel
       do j = 1, lx * lx
          do i = 1, lx
             du(i,j,1,e) = dx(i,1) * u(1,j,1,e) &
                         + dx(i,2) * u(2,j,1,e) 
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * dr(i,1,1,e)
       end do
       
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                drst(i,j,k) = dy(j,1) * u(i,1,k,e) &
                            + dy(j,2) * u(i,2,k,e)
             end do
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * ds(i,1,1,e)
       end do
       
       do k = 1, lx
          do i = 1, lx*lx
             drst(i,1,k) = dz(k,1) * u(i,1,1,e) &
                         + dz(k,2) * u(i,1,2,e) 
          end do
       end do
       
       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) + drst(i,1,1) * dt(i,1,1,e)
       end do

       do i = 1, lx * lx * lx
          du(i,1,1,e) = du(i,1,1,e) * jacinv(i,1,1,e)
       end do

    end do
    !$omp end do
  end subroutine cpu_dudxyz_lx2

end module cpu_dudxyz
