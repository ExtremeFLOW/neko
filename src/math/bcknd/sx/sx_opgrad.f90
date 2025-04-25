! Copyright (c) 2021, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!  * Redistributions of source code must retain the above copyright
!   notice, this list of conditions and the following disclaimer.
!
!  * Redistributions in binary form must reproduce the above
!   copyright notice, this list of conditions and the following
!   disclaimer in the documentation and/or other materials provided
!   with the distribution.
!
!  * Neither the name of the authors nor the names of its
!   contributors may be used to endorse or promote products derived
!   from this software without specific prior written permission.
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
!> Gradient kernels for SX-Aurora
submodule (opr_sx) sx_opgrad
  implicit none

contains

  module subroutine opr_sx_opgrad(ux, uy, uz, u, coef)
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: ux(coef%Xh%lxyz, coef%msh%nelv)
    real(kind=rp), intent(inout) :: uy(coef%Xh%lxyz, coef%msh%nelv)
    real(kind=rp), intent(inout) :: uz(coef%Xh%lxyz, coef%msh%nelv)
    real(kind=rp), intent(in) :: u(coef%Xh%lxyz, coef%msh%nelv)

    associate(Xh => coef%Xh, msh => coef%msh)
      select case (Xh%lx)
        case (18)
         call sx_opgrad_lx18(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
        case (17)
         call sx_opgrad_lx17(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
        case (16)
         call sx_opgrad_lx16(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
        case (15)
         call sx_opgrad_lx15(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
        case (14)
         call sx_opgrad_lx14(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
        case (13)
         call sx_opgrad_lx13(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
        case (12)
         call sx_opgrad_lx12(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
        case (11)
         call sx_opgrad_lx11(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
        case (10)
         call sx_opgrad_lx10(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
        case (9)
         call sx_opgrad_lx9(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
        case (8)
         call sx_opgrad_lx8(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
        case (7)
         call sx_opgrad_lx7(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
        case (6)
         call sx_opgrad_lx6(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
        case (5)
         call sx_opgrad_lx5(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
        case (4)
         call sx_opgrad_lx4(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
        case (3)
         call sx_opgrad_lx3(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
        case (2)
         call sx_opgrad_lx2(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv)
        case default
         call sx_opgrad_lx(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              coef%drdx, coef%dsdx, coef%dtdx, &
              coef%drdy, coef%dsdy, coef%dtdy, &
              coef%drdz, coef%dsdz, coef%dtdz, &
              Xh%w3, msh%nelv, Xh%lx)
      end select
    end associate

  end subroutine opr_sx_opgrad

  subroutine sx_opgrad_lx(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n, lx)
    integer, intent(in) :: n, lx
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          ur(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk, k,e)
                end do
                us(i, j, k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i, j, kk, e)
                end do
                ut(i, j, k, e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )

          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx

  subroutine sx_opgrad_lx18(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 18
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          ur(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk, k,e)
                end do
                us(i, j, k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i, j, kk,e)
                end do
                ut(i, j, k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )

          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx18

  subroutine sx_opgrad_lx17(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 17
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj, 1, 1)
          end do
          ur(i, jj, 1, 1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk, k,e)
                end do
                us(i, j, k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i, j, kk, e)
                end do
                ut(i, j, k, e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )

          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx17

  subroutine sx_opgrad_lx16(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 16
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          ur(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk, k,e)
                end do
                us(i, j, k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i, j, kk,e)
                end do
                ut(i, j, k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )

          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx16

  subroutine sx_opgrad_lx15(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 15
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj, 1, 1)
          end do
          ur(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk, k, e)
                end do
                us(i, j, k, e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i, j, kk, e)
                end do
                ut(i, j, k, e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )

          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx15

  subroutine sx_opgrad_lx14(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 14
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          ur(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk, k, e)
                end do
                us(i, j, k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i, j, kk,e)
                end do
                ut(i, j, k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )

          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx14

  subroutine sx_opgrad_lx13(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 13
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          ur(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk, k,e)
                end do
                us(i, j, k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i, j, kk,e)
                end do
                ut(i, j, k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )

          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx13

  subroutine sx_opgrad_lx12(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 12
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          ur(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk, k,e)
                end do
                us(i, j, k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i, j, kk,e)
                end do
                ut(i, j, k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )

          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx12

  subroutine sx_opgrad_lx11(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 11
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          ur(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk, k,e)
                end do
                us(i, j, k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i, j, kk,e)
                end do
                ut(i, j, k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )

          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx11

  subroutine sx_opgrad_lx10(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 10
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          ur(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk, k,e)
                end do
                us(i, j, k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i, j, kk,e)
                end do
                ut(i, j, k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )

          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx10

  subroutine sx_opgrad_lx9(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 9
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          ur(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk, k,e)
                end do
                us(i, j, k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i, j, kk,e)
                end do
                ut(i, j, k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )

          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx9

  subroutine sx_opgrad_lx8(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 8
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          ur(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk, k,e)
                end do
                us(i, j, k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i, j, kk,e)
                end do
                ut(i, j, k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )

          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx8

  subroutine sx_opgrad_lx7(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 7
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          ur(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk, k,e)
                end do
                us(i, j, k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i, j, kk,e)
                end do
                ut(i, j, k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )

          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx7

  subroutine sx_opgrad_lx6(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 6
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          ur(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk, k,e)
                end do
                us(i, j, k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i, j, kk,e)
                end do
                ut(i, j, k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )

          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx6

  subroutine sx_opgrad_lx5(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 5
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          ur(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk, k,e)
                end do
                us(i, j, k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i, j, kk,e)
                end do
                ut(i, j, k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )

          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx5

  subroutine sx_opgrad_lx4(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 4
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          ur(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk, k,e)
                end do
                us(i, j, k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i, j, kk,e)
                end do
                ut(i, j, k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )

          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx4

  subroutine sx_opgrad_lx3(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 3
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          ur(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk, k,e)
                end do
                us(i, j, k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i, j, kk,e)
                end do
                ut(i, j, k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )

          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx3

  subroutine sx_opgrad_lx2(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 2
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          ur(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk, k,e)
                end do
                us(i, j, k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i, j, kk,e)
                end do
                ut(i, j, k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          ux(i,1,1,e) = w3(i,1,1) * &
               ( drdx(i,1,1,e) * ur(i,1,1,e) &
               + dsdx(i,1,1,e) * us(i,1,1,e) &
               + dtdx(i,1,1,e) * ut(i,1,1,e))

          uy(i,1,1,e) = w3(i,1,1) * &
               ( dsdy(i,1,1,e) * us(i,1,1,e) &
               + drdy(i,1,1,e) * ur(i,1,1,e) &
               + dtdy(i,1,1,e) * ut(i,1,1,e) )

          uz(i,1,1,e) = w3(i,1,1) * &
               ( dtdz(i,1,1,e) * ut(i,1,1,e) &
               + drdz(i,1,1,e) * ur(i,1,1,e) &
               + dsdz(i,1,1,e) * us(i,1,1,e))
       end do
    end do
  end subroutine sx_opgrad_lx2

end submodule sx_opgrad
