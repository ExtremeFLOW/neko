! Copyright (c) 2021-2024, The Neko Authors
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
!> conv1 SX-Aurora kernels
submodule (opr_sx) sx_conv1
  implicit none

contains

  module subroutine opr_sx_conv1(du, u, vx, vy, vz, Xh, coef, nelv)
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    integer, intent(in) :: nelv
    real(kind=rp), intent(inout) :: du(Xh%lxyz, nelv)
    real(kind=rp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, nelv)
    real(kind=rp), intent(inout) :: vx(Xh%lx, Xh%ly, Xh%lz, nelv)
    real(kind=rp), intent(inout) :: vy(Xh%lx, Xh%ly, Xh%lz, nelv)
    real(kind=rp), intent(inout) :: vz(Xh%lx, Xh%ly, Xh%lz, nelv)

    select case (Xh%lx)
      case (14)
       call sx_conv1_lx14(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
      case (13)
       call sx_conv1_lx13(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
      case (12)
       call sx_conv1_lx12(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
      case (11)
       call sx_conv1_lx11(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
      case (10)
       call sx_conv1_lx10(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
      case (9)
       call sx_conv1_lx9(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
      case (8)
       call sx_conv1_lx8(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
      case (7)
       call sx_conv1_lx7(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
      case (6)
       call sx_conv1_lx6(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
      case (5)
       call sx_conv1_lx5(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
      case (4)
       call sx_conv1_lx4(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
      case (3)
       call sx_conv1_lx3(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
      case (2)
       call sx_conv1_lx2(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv)
      case default
       call sx_conv1_lx(du, u, vx, vy, vz, Xh%dx, Xh%dy, Xh%dz, &
            coef%drdx, coef%dsdx, coef%dtdx, &
            coef%drdy, coef%dsdy, coef%dtdy, &
            coef%drdz, coef%dsdz, coef%dtdz, &
            coef%jacinv, nelv, Xh%lx)
    end select

  end subroutine opr_sx_conv1

  subroutine sx_conv1_lx(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv, lx)
    integer, intent(in) :: nelv, lx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudr
    real(kind=rp), dimension(lx, lx, lx, nelv) :: duds
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudt
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * nelv
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          dudr(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i,j, kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do

  end subroutine sx_conv1_lx

  subroutine sx_conv1_lx14(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 14
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudr
    real(kind=rp), dimension(lx, lx, lx, nelv) :: duds
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudt
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * nelv
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          dudr(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i,j, kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do

  end subroutine sx_conv1_lx14

  subroutine sx_conv1_lx13(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 13
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudr
    real(kind=rp), dimension(lx, lx, lx, nelv) :: duds
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudt
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * nelv
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          dudr(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i,j, kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do

  end subroutine sx_conv1_lx13

  subroutine sx_conv1_lx12(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 12
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudr
    real(kind=rp), dimension(lx, lx, lx, nelv) :: duds
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudt
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * nelv
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          dudr(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i,j, kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do

  end subroutine sx_conv1_lx12

  subroutine sx_conv1_lx11(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 11
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudr
    real(kind=rp), dimension(lx, lx, lx, nelv) :: duds
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudt
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * nelv
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          dudr(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i,j, kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do

  end subroutine sx_conv1_lx11

  subroutine sx_conv1_lx10(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 10
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudr
    real(kind=rp), dimension(lx, lx, lx, nelv) :: duds
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudt
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * nelv
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          dudr(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i,j, kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do

  end subroutine sx_conv1_lx10

  subroutine sx_conv1_lx9(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 9
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudr
    real(kind=rp), dimension(lx, lx, lx, nelv) :: duds
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudt
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * nelv
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          dudr(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i,j, kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do

  end subroutine sx_conv1_lx9

  subroutine sx_conv1_lx8(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 8
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudr
    real(kind=rp), dimension(lx, lx, lx, nelv) :: duds
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudt
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * nelv
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          dudr(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i,j, kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do

  end subroutine sx_conv1_lx8

  subroutine sx_conv1_lx7(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 7
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudr
    real(kind=rp), dimension(lx, lx, lx, nelv) :: duds
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudt
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * nelv
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          dudr(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i,j, kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do

  end subroutine sx_conv1_lx7

  subroutine sx_conv1_lx6(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 6
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudr
    real(kind=rp), dimension(lx, lx, lx, nelv) :: duds
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudt
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * nelv
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          dudr(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i,j, kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do

  end subroutine sx_conv1_lx6

  subroutine sx_conv1_lx5(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 5
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudr
    real(kind=rp), dimension(lx, lx, lx, nelv) :: duds
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudt
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * nelv
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          dudr(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i,j, kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do

  end subroutine sx_conv1_lx5

  subroutine sx_conv1_lx4(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 4
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudr
    real(kind=rp), dimension(lx, lx, lx, nelv) :: duds
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudt
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * nelv
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          dudr(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i,j, kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do

  end subroutine sx_conv1_lx4

  subroutine sx_conv1_lx3(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 3
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudr
    real(kind=rp), dimension(lx, lx, lx, nelv) :: duds
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudt
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * nelv
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          dudr(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i,j, kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do

  end subroutine sx_conv1_lx3

  subroutine sx_conv1_lx2(du, u, vx, vy, vz, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, &
       jacinv, nelv)
    integer, parameter :: lx = 2
    integer, intent(in) :: nelv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: u, vx, vy, vz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudr
    real(kind=rp), dimension(lx, lx, lx, nelv) :: duds
    real(kind=rp), dimension(lx, lx, lx, nelv) :: dudt
    real(kind=rp) :: wr, ws, wt
    integer :: e, i, j, k, jj, kk

    do i = 1, lx
       do jj = 1, lx * lx * nelv
          wr = 0d0
          do kk = 1, lx
             wr = wr + dx(i, kk) * u(kk, jj,1,1)
          end do
          dudr(i, jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, nelv
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + dy(j, kk) * u(i, kk,k,e)
                end do
                duds(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, nelv
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dz(k, kk) * u(i,j, kk,e)
                end do
                dudt(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, nelv * lx * lx * lx
       du(i,1,1,1) = jacinv(i,1,1,1) * &
            (vx(i,1,1,1) * ( &
            drdx(i,1,1,1)*dudr(i,1,1,1) &
            + dsdx(i,1,1,1)*duds(i,1,1,1) &
            + dtdx(i,1,1,1)*dudt(i,1,1,1)) &
            + vy(i,1,1,1)*( &
            drdy(i,1,1,1)*dudr(i,1,1,1) &
            + dsdy(i,1,1,1)*duds(i,1,1,1) &
            + dtdy(i,1,1,1)*dudt(i,1,1,1)) &
            + vz(i,1,1,1)*( &
            drdz(i,1,1,1)*dudr(i,1,1,1) &
            + dsdz(i,1,1,1)*duds(i,1,1,1) &
            + dtdz(i,1,1,1)*dudt(i,1,1,1)))
    end do

  end subroutine sx_conv1_lx2

end submodule sx_conv1
