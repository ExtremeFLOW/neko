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
!> Derivative kernels for SX-Aurora
module sx_dudxyz
  use num_types, only : rp
  use math
  implicit none
  private

  public :: sx_dudxyz_lx, sx_dudxyz_lx14, sx_dudxyz_lx13, sx_dudxyz_lx12, &
       sx_dudxyz_lx11, sx_dudxyz_lx10, sx_dudxyz_lx9, sx_dudxyz_lx8, &
       sx_dudxyz_lx7, sx_dudxyz_lx6, sx_dudxyz_lx5, sx_dudxyz_lx4, &
       sx_dudxyz_lx3, sx_dudxyz_lx2

contains

  subroutine sx_dudxyz_lx(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd, lx)
    integer, intent(in) :: nel, nd, lx
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k
    integer :: i, j, jj, kk 
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx
  
  subroutine sx_dudxyz_lx14(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 14
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k
    integer :: i, j, jj, kk 
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx14
  
  subroutine sx_dudxyz_lx13(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 13
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k
    integer :: i, j, jj, kk 
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx13

  subroutine sx_dudxyz_lx12(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 12
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k
    integer :: i, j, jj, kk 
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx12

  subroutine sx_dudxyz_lx11(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 11
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k
    integer :: i, j, jj, kk 
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx11

  subroutine sx_dudxyz_lx10(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 10
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k
    integer :: i, j, jj, kk 
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx10

  subroutine sx_dudxyz_lx9(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 9
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k
    integer :: i, j, jj, kk 
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx9

  subroutine sx_dudxyz_lx8(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 8
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k
    integer :: i, j, jj, kk 
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx8

  subroutine sx_dudxyz_lx7(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 7
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k
    integer :: i, j, jj, kk 
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx7

  subroutine sx_dudxyz_lx6(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 6
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k
    integer :: i, j, jj, kk 
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx6

  subroutine sx_dudxyz_lx5(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 5
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k
    integer :: i, j, jj, kk 
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx5

  subroutine sx_dudxyz_lx4(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 4
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k
    integer :: i, j, jj, kk 
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx4

  subroutine sx_dudxyz_lx3(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 3
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k
    integer :: i, j, jj, kk 
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx3

  subroutine sx_dudxyz_lx2(du, u, dr, ds, dt, dx, dy, dz, jacinv, nel, nd)
    integer, parameter :: lx = 2
    integer, intent(in) :: nel, nd
    real(kind=rp), dimension(lx,lx,lx,nel), intent(inout) ::  du
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) ::  u, dr, ds, dt
    real(kind=rp), dimension(lx,lx,lx,nel), intent(in) :: jacinv
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx,nel) :: drst
    integer :: e, k
    integer :: i, j, jj, kk 
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx*lx*nel
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i,kk)*u(kk,jj,1,1)
          end do
          du(i,jj,1,1) = wr
       end do
    end do

    call col2 (du, dr, nd)

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nel     
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j,kk)*u(i,kk,k,e)
                end do
                drst(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    call addcol3(du, drst, ds, nd)

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nel
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k,kk)*u(i,j,kk,e)
                end do
                drst(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    call addcol3(du, drst, dt, nd)
    call col2 (du, jacinv, nd)
  end subroutine sx_dudxyz_lx2

end module sx_dudxyz
