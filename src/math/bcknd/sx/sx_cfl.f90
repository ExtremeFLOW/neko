! Copyright (c) 2022-2024, The Neko Authors
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
!> CFL SX-Aurora kernels
submodule (opr_sx) sx_cfl
  implicit none

contains

  module function opr_sx_cfl(dt, u, v, w, Xh, coef, nelv) result(cfl)
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: nelv
    real(kind=rp), intent(in) :: dt
    real(kind=rp), dimension(Xh%lx, Xh%ly, Xh%lz, nelv) ::  u, v, w
    real(kind=rp) :: cfl

    select case (Xh%lx)
    case (14)
       cfl = sx_cfl_lx14(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (13)
       cfl = sx_cfl_lx13(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (12)
       cfl = sx_cfl_lx12(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (11)
       cfl = sx_cfl_lx11(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (10)
       cfl = sx_cfl_lx10(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (9)
       cfl = sx_cfl_lx9(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (8)
       cfl = sx_cfl_lx8(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (7)
       cfl = sx_cfl_lx7(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (6)
       cfl = sx_cfl_lx6(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (5)
       cfl = sx_cfl_lx5(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (4)
       cfl = sx_cfl_lx4(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (3)
       cfl = sx_cfl_lx3(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case (2)
       cfl = sx_cfl_lx2(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv)
    case default
       cfl = sx_cfl_lx(dt, u, v, w, &
           coef%drdx, coef%dsdx, coef%dtdx, &
           coef%drdy, coef%dsdy, coef%dtdy, &
           coef%drdz, coef%dsdz, coef%dtdz, &
           Xh%dr_inv, Xh%ds_inv, Xh%dt_inv, &
           coef%jacinv, nelv, Xh%lx)
    end select

  end function opr_sx_cfl
  
  function sx_cfl_lx(dt, u, v, w, drdx, dsdx, dtdx, drdy, dsdy, dtdy, &
       drdz, dsdz, dtdz, dr_inv, ds_inv, dt_inv, &
       jacinv, nelv, lx) result(cfl)
    integer, intent(in) :: nelv, lx
    real(kind=rp), intent(in) :: dt
    real(kind=rp), dimension(lx, lx, lx, nelv) ::  u, v, w
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx), intent(in) :: dr_inv, ds_inv, dt_inv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp) :: cflr, cfls, cflt, cflm
    real(kind=rp) :: ur, us, ut
    real(kind=rp) :: cfl
    integer :: i, j, k, e
    cfl = 0d0

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             do e = 1, nelv
                ur = ( u(i,j,k,e)*drdx(i,j,k,e) &
                   +   v(i,j,k,e)*drdy(i,j,k,e) &
                   +   w(i,j,k,e)*drdz(i,j,k,e) ) * jacinv(i,j,k,e)
                us = ( u(i,j,k,e)*dsdx(i,j,k,e) &
                   +   v(i,j,k,e)*dsdy(i,j,k,e) &
                   +   w(i,j,k,e)*dsdz(i,j,k,e) ) * jacinv(i,j,k,e)
                ut = ( u(i,j,k,e)*dtdx(i,j,k,e) &
                   +   v(i,j,k,e)*dtdy(i,j,k,e) &
                   +   w(i,j,k,e)*dtdz(i,j,k,e) ) * jacinv(i,j,k,e)

                cflr = abs(dt*ur*dr_inv(i))
                cfls = abs(dt*us*ds_inv(j))
                cflt = abs(dt*ut*dt_inv(k))

                cflm = cflr + cfls + cflt
                cfl  = max(cfl, cflm)
             end do
          end do
       end do
    end do

  end function sx_cfl_lx

  function sx_cfl_lx14(dt, u, v, w, drdx, dsdx, dtdx, drdy, dsdy, dtdy, &
       drdz, dsdz, dtdz, dr_inv, ds_inv, dt_inv, &
       jacinv, nelv) result(cfl)
    integer, parameter :: lx = 14
    integer, intent(in) :: nelv
    real(kind=rp), intent(in) :: dt
    real(kind=rp), dimension(lx, lx, lx, nelv) ::  u, v, w
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx), intent(in) :: dr_inv, ds_inv, dt_inv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp) :: cflr, cfls, cflt, cflm
    real(kind=rp) :: ur, us, ut
    real(kind=rp) :: cfl
    integer :: i, j, k, e
    cfl = 0d0

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             do e = 1, nelv
                ur = ( u(i,j,k,e)*drdx(i,j,k,e) &
                   +   v(i,j,k,e)*drdy(i,j,k,e) &
                   +   w(i,j,k,e)*drdz(i,j,k,e) ) * jacinv(i,j,k,e)
                us = ( u(i,j,k,e)*dsdx(i,j,k,e) &
                   +   v(i,j,k,e)*dsdy(i,j,k,e) &
                   +   w(i,j,k,e)*dsdz(i,j,k,e) ) * jacinv(i,j,k,e)
                ut = ( u(i,j,k,e)*dtdx(i,j,k,e) &
                   +   v(i,j,k,e)*dtdy(i,j,k,e) &
                   +   w(i,j,k,e)*dtdz(i,j,k,e) ) * jacinv(i,j,k,e)

                cflr = abs(dt*ur*dr_inv(i))
                cfls = abs(dt*us*ds_inv(j))
                cflt = abs(dt*ut*dt_inv(k))

                cflm = cflr + cfls + cflt
                cfl  = max(cfl, cflm)
             end do
          end do
       end do
    end do

  end function sx_cfl_lx14

  function sx_cfl_lx13(dt, u, v, w, drdx, dsdx, dtdx, drdy, dsdy, dtdy, &
       drdz, dsdz, dtdz, dr_inv, ds_inv, dt_inv, &
       jacinv, nelv) result(cfl)
    integer, parameter :: lx = 13
    integer, intent(in) :: nelv
    real(kind=rp), intent(in) :: dt
    real(kind=rp), dimension(lx, lx, lx, nelv) ::  u, v, w
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx), intent(in) :: dr_inv, ds_inv, dt_inv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp) :: cflr, cfls, cflt, cflm
    real(kind=rp) :: ur, us, ut
    real(kind=rp) :: cfl
    integer :: i, j, k, e
    cfl = 0d0

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             do e = 1, nelv
                ur = ( u(i,j,k,e)*drdx(i,j,k,e) &
                   +   v(i,j,k,e)*drdy(i,j,k,e) &
                   +   w(i,j,k,e)*drdz(i,j,k,e) ) * jacinv(i,j,k,e)
                us = ( u(i,j,k,e)*dsdx(i,j,k,e) &
                   +   v(i,j,k,e)*dsdy(i,j,k,e) &
                   +   w(i,j,k,e)*dsdz(i,j,k,e) ) * jacinv(i,j,k,e)
                ut = ( u(i,j,k,e)*dtdx(i,j,k,e) &
                   +   v(i,j,k,e)*dtdy(i,j,k,e) &
                   +   w(i,j,k,e)*dtdz(i,j,k,e) ) * jacinv(i,j,k,e)

                cflr = abs(dt*ur*dr_inv(i))
                cfls = abs(dt*us*ds_inv(j))
                cflt = abs(dt*ut*dt_inv(k))

                cflm = cflr + cfls + cflt
                cfl  = max(cfl, cflm)
             end do
          end do
       end do
    end do

  end function sx_cfl_lx13

  function sx_cfl_lx12(dt, u, v, w, drdx, dsdx, dtdx, drdy, dsdy, dtdy, &
       drdz, dsdz, dtdz, dr_inv, ds_inv, dt_inv, &
       jacinv, nelv) result(cfl)
    integer, parameter :: lx = 12
    integer, intent(in) :: nelv
    real(kind=rp), intent(in) :: dt
    real(kind=rp), dimension(lx, lx, lx, nelv) ::  u, v, w
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx), intent(in) :: dr_inv, ds_inv, dt_inv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp) :: cflr, cfls, cflt, cflm
    real(kind=rp) :: ur, us, ut
    real(kind=rp) :: cfl
    integer :: i, j, k, e
    cfl = 0d0

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             do e = 1, nelv
                ur = ( u(i,j,k,e)*drdx(i,j,k,e) &
                   +   v(i,j,k,e)*drdy(i,j,k,e) &
                   +   w(i,j,k,e)*drdz(i,j,k,e) ) * jacinv(i,j,k,e)
                us = ( u(i,j,k,e)*dsdx(i,j,k,e) &
                   +   v(i,j,k,e)*dsdy(i,j,k,e) &
                   +   w(i,j,k,e)*dsdz(i,j,k,e) ) * jacinv(i,j,k,e)
                ut = ( u(i,j,k,e)*dtdx(i,j,k,e) &
                   +   v(i,j,k,e)*dtdy(i,j,k,e) &
                   +   w(i,j,k,e)*dtdz(i,j,k,e) ) * jacinv(i,j,k,e)

                cflr = abs(dt*ur*dr_inv(i))
                cfls = abs(dt*us*ds_inv(j))
                cflt = abs(dt*ut*dt_inv(k))

                cflm = cflr + cfls + cflt
                cfl  = max(cfl, cflm)
             end do
          end do
       end do
    end do

  end function sx_cfl_lx12

  function sx_cfl_lx11(dt, u, v, w, drdx, dsdx, dtdx, drdy, dsdy, dtdy, &
       drdz, dsdz, dtdz, dr_inv, ds_inv, dt_inv, &
       jacinv, nelv) result(cfl)
    integer, parameter :: lx = 11
    integer, intent(in) :: nelv
    real(kind=rp), intent(in) :: dt
    real(kind=rp), dimension(lx, lx, lx, nelv) ::  u, v, w
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx), intent(in) :: dr_inv, ds_inv, dt_inv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp) :: cflr, cfls, cflt, cflm
    real(kind=rp) :: ur, us, ut
    real(kind=rp) :: cfl
    integer :: i, j, k, e
    cfl = 0d0

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             do e = 1, nelv
                ur = ( u(i,j,k,e)*drdx(i,j,k,e) &
                   +   v(i,j,k,e)*drdy(i,j,k,e) &
                   +   w(i,j,k,e)*drdz(i,j,k,e) ) * jacinv(i,j,k,e)
                us = ( u(i,j,k,e)*dsdx(i,j,k,e) &
                   +   v(i,j,k,e)*dsdy(i,j,k,e) &
                   +   w(i,j,k,e)*dsdz(i,j,k,e) ) * jacinv(i,j,k,e)
                ut = ( u(i,j,k,e)*dtdx(i,j,k,e) &
                   +   v(i,j,k,e)*dtdy(i,j,k,e) &
                   +   w(i,j,k,e)*dtdz(i,j,k,e) ) * jacinv(i,j,k,e)

                cflr = abs(dt*ur*dr_inv(i))
                cfls = abs(dt*us*ds_inv(j))
                cflt = abs(dt*ut*dt_inv(k))

                cflm = cflr + cfls + cflt
                cfl  = max(cfl, cflm)
             end do
          end do
       end do
    end do

  end function sx_cfl_lx11

  function sx_cfl_lx10(dt, u, v, w, drdx, dsdx, dtdx, drdy, dsdy, dtdy, &
       drdz, dsdz, dtdz, dr_inv, ds_inv, dt_inv, &
       jacinv, nelv) result(cfl)
    integer, parameter :: lx = 10
    integer, intent(in) :: nelv
    real(kind=rp), intent(in) :: dt
    real(kind=rp), dimension(lx, lx, lx, nelv) ::  u, v, w
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx), intent(in) :: dr_inv, ds_inv, dt_inv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp) :: cflr, cfls, cflt, cflm
    real(kind=rp) :: ur, us, ut
    real(kind=rp) :: cfl
    integer :: i, j, k, e
    cfl = 0d0

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             do e = 1, nelv
                ur = ( u(i,j,k,e)*drdx(i,j,k,e) &
                   +   v(i,j,k,e)*drdy(i,j,k,e) &
                   +   w(i,j,k,e)*drdz(i,j,k,e) ) * jacinv(i,j,k,e)
                us = ( u(i,j,k,e)*dsdx(i,j,k,e) &
                   +   v(i,j,k,e)*dsdy(i,j,k,e) &
                   +   w(i,j,k,e)*dsdz(i,j,k,e) ) * jacinv(i,j,k,e)
                ut = ( u(i,j,k,e)*dtdx(i,j,k,e) &
                   +   v(i,j,k,e)*dtdy(i,j,k,e) &
                   +   w(i,j,k,e)*dtdz(i,j,k,e) ) * jacinv(i,j,k,e)

                cflr = abs(dt*ur*dr_inv(i))
                cfls = abs(dt*us*ds_inv(j))
                cflt = abs(dt*ut*dt_inv(k))

                cflm = cflr + cfls + cflt
                cfl  = max(cfl, cflm)
             end do
          end do
       end do
    end do

  end function sx_cfl_lx10

  function sx_cfl_lx9(dt, u, v, w, drdx, dsdx, dtdx, drdy, dsdy, dtdy, &
       drdz, dsdz, dtdz, dr_inv, ds_inv, dt_inv, &
       jacinv, nelv) result(cfl)
    integer, parameter :: lx = 9
    integer, intent(in) :: nelv
    real(kind=rp), intent(in) :: dt
    real(kind=rp), dimension(lx, lx, lx, nelv) ::  u, v, w
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx), intent(in) :: dr_inv, ds_inv, dt_inv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp) :: cflr, cfls, cflt, cflm
    real(kind=rp) :: ur, us, ut
    real(kind=rp) :: cfl
    integer :: i, j, k, e
    cfl = 0d0

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             do e = 1, nelv
                ur = ( u(i,j,k,e)*drdx(i,j,k,e) &
                   +   v(i,j,k,e)*drdy(i,j,k,e) &
                   +   w(i,j,k,e)*drdz(i,j,k,e) ) * jacinv(i,j,k,e)
                us = ( u(i,j,k,e)*dsdx(i,j,k,e) &
                   +   v(i,j,k,e)*dsdy(i,j,k,e) &
                   +   w(i,j,k,e)*dsdz(i,j,k,e) ) * jacinv(i,j,k,e)
                ut = ( u(i,j,k,e)*dtdx(i,j,k,e) &
                   +   v(i,j,k,e)*dtdy(i,j,k,e) &
                   +   w(i,j,k,e)*dtdz(i,j,k,e) ) * jacinv(i,j,k,e)

                cflr = abs(dt*ur*dr_inv(i))
                cfls = abs(dt*us*ds_inv(j))
                cflt = abs(dt*ut*dt_inv(k))

                cflm = cflr + cfls + cflt
                cfl  = max(cfl, cflm)
             end do
          end do
       end do
    end do

  end function sx_cfl_lx9

  function sx_cfl_lx8(dt, u, v, w, drdx, dsdx, dtdx, drdy, dsdy, dtdy, &
       drdz, dsdz, dtdz, dr_inv, ds_inv, dt_inv, &
       jacinv, nelv) result(cfl)
    integer, parameter :: lx = 8
    integer, intent(in) :: nelv
    real(kind=rp), intent(in) :: dt
    real(kind=rp), dimension(lx, lx, lx, nelv) ::  u, v, w
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx), intent(in) :: dr_inv, ds_inv, dt_inv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp) :: cflr, cfls, cflt, cflm
    real(kind=rp) :: ur, us, ut
    real(kind=rp) :: cfl
    integer :: i, j, k, e
    cfl = 0d0

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             do e = 1, nelv
                ur = ( u(i,j,k,e)*drdx(i,j,k,e) &
                   +   v(i,j,k,e)*drdy(i,j,k,e) &
                   +   w(i,j,k,e)*drdz(i,j,k,e) ) * jacinv(i,j,k,e)
                us = ( u(i,j,k,e)*dsdx(i,j,k,e) &
                   +   v(i,j,k,e)*dsdy(i,j,k,e) &
                   +   w(i,j,k,e)*dsdz(i,j,k,e) ) * jacinv(i,j,k,e)
                ut = ( u(i,j,k,e)*dtdx(i,j,k,e) &
                   +   v(i,j,k,e)*dtdy(i,j,k,e) &
                   +   w(i,j,k,e)*dtdz(i,j,k,e) ) * jacinv(i,j,k,e)

                cflr = abs(dt*ur*dr_inv(i))
                cfls = abs(dt*us*ds_inv(j))
                cflt = abs(dt*ut*dt_inv(k))

                cflm = cflr + cfls + cflt
                cfl  = max(cfl, cflm)
             end do
          end do
       end do
    end do

  end function sx_cfl_lx8

  function sx_cfl_lx7(dt, u, v, w, drdx, dsdx, dtdx, drdy, dsdy, dtdy, &
       drdz, dsdz, dtdz, dr_inv, ds_inv, dt_inv, &
       jacinv, nelv) result(cfl)
    integer, parameter :: lx = 7
    integer, intent(in) :: nelv
    real(kind=rp), intent(in) :: dt
    real(kind=rp), dimension(lx, lx, lx, nelv) ::  u, v, w
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx), intent(in) :: dr_inv, ds_inv, dt_inv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp) :: cflr, cfls, cflt, cflm
    real(kind=rp) :: ur, us, ut
    real(kind=rp) :: cfl
    integer :: i, j, k, e
    cfl = 0d0

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             do e = 1, nelv
                ur = ( u(i,j,k,e)*drdx(i,j,k,e) &
                   +   v(i,j,k,e)*drdy(i,j,k,e) &
                   +   w(i,j,k,e)*drdz(i,j,k,e) ) * jacinv(i,j,k,e)
                us = ( u(i,j,k,e)*dsdx(i,j,k,e) &
                   +   v(i,j,k,e)*dsdy(i,j,k,e) &
                   +   w(i,j,k,e)*dsdz(i,j,k,e) ) * jacinv(i,j,k,e)
                ut = ( u(i,j,k,e)*dtdx(i,j,k,e) &
                   +   v(i,j,k,e)*dtdy(i,j,k,e) &
                   +   w(i,j,k,e)*dtdz(i,j,k,e) ) * jacinv(i,j,k,e)

                cflr = abs(dt*ur*dr_inv(i))
                cfls = abs(dt*us*ds_inv(j))
                cflt = abs(dt*ut*dt_inv(k))

                cflm = cflr + cfls + cflt
                cfl  = max(cfl, cflm)
             end do
          end do
       end do
    end do

  end function sx_cfl_lx7

  function sx_cfl_lx6(dt, u, v, w, drdx, dsdx, dtdx, drdy, dsdy, dtdy, &
       drdz, dsdz, dtdz, dr_inv, ds_inv, dt_inv, &
       jacinv, nelv) result(cfl)
    integer, parameter :: lx = 6
    integer, intent(in) :: nelv
    real(kind=rp), intent(in) :: dt
    real(kind=rp), dimension(lx, lx, lx, nelv) ::  u, v, w
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx), intent(in) :: dr_inv, ds_inv, dt_inv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp) :: cflr, cfls, cflt, cflm
    real(kind=rp) :: ur, us, ut
    real(kind=rp) :: cfl
    integer :: i, j, k, e
    cfl = 0d0

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             do e = 1, nelv
                ur = ( u(i,j,k,e)*drdx(i,j,k,e) &
                   +   v(i,j,k,e)*drdy(i,j,k,e) &
                   +   w(i,j,k,e)*drdz(i,j,k,e) ) * jacinv(i,j,k,e)
                us = ( u(i,j,k,e)*dsdx(i,j,k,e) &
                   +   v(i,j,k,e)*dsdy(i,j,k,e) &
                   +   w(i,j,k,e)*dsdz(i,j,k,e) ) * jacinv(i,j,k,e)
                ut = ( u(i,j,k,e)*dtdx(i,j,k,e) &
                   +   v(i,j,k,e)*dtdy(i,j,k,e) &
                   +   w(i,j,k,e)*dtdz(i,j,k,e) ) * jacinv(i,j,k,e)

                cflr = abs(dt*ur*dr_inv(i))
                cfls = abs(dt*us*ds_inv(j))
                cflt = abs(dt*ut*dt_inv(k))

                cflm = cflr + cfls + cflt
                cfl  = max(cfl, cflm)
             end do
          end do
       end do
    end do

  end function sx_cfl_lx6

  function sx_cfl_lx5(dt, u, v, w, drdx, dsdx, dtdx, drdy, dsdy, dtdy, &
       drdz, dsdz, dtdz, dr_inv, ds_inv, dt_inv, &
       jacinv, nelv) result(cfl)
    integer, parameter :: lx = 5
    integer, intent(in) :: nelv
    real(kind=rp), intent(in) :: dt
    real(kind=rp), dimension(lx, lx, lx, nelv) ::  u, v, w
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx), intent(in) :: dr_inv, ds_inv, dt_inv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp) :: cflr, cfls, cflt, cflm
    real(kind=rp) :: ur, us, ut
    real(kind=rp) :: cfl
    integer :: i, j, k, e
    cfl = 0d0

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             do e = 1, nelv
                ur = ( u(i,j,k,e)*drdx(i,j,k,e) &
                   +   v(i,j,k,e)*drdy(i,j,k,e) &
                   +   w(i,j,k,e)*drdz(i,j,k,e) ) * jacinv(i,j,k,e)
                us = ( u(i,j,k,e)*dsdx(i,j,k,e) &
                   +   v(i,j,k,e)*dsdy(i,j,k,e) &
                   +   w(i,j,k,e)*dsdz(i,j,k,e) ) * jacinv(i,j,k,e)
                ut = ( u(i,j,k,e)*dtdx(i,j,k,e) &
                   +   v(i,j,k,e)*dtdy(i,j,k,e) &
                   +   w(i,j,k,e)*dtdz(i,j,k,e) ) * jacinv(i,j,k,e)

                cflr = abs(dt*ur*dr_inv(i))
                cfls = abs(dt*us*ds_inv(j))
                cflt = abs(dt*ut*dt_inv(k))

                cflm = cflr + cfls + cflt
                cfl  = max(cfl, cflm)
             end do
          end do
       end do
    end do

  end function sx_cfl_lx5

  function sx_cfl_lx4(dt, u, v, w, drdx, dsdx, dtdx, drdy, dsdy, dtdy, &
       drdz, dsdz, dtdz, dr_inv, ds_inv, dt_inv, &
       jacinv, nelv) result(cfl)
    integer, parameter :: lx = 4
    integer, intent(in) :: nelv
    real(kind=rp), intent(in) :: dt
    real(kind=rp), dimension(lx, lx, lx, nelv) ::  u, v, w
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx), intent(in) :: dr_inv, ds_inv, dt_inv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp) :: cflr, cfls, cflt, cflm
    real(kind=rp) :: ur, us, ut
    real(kind=rp) :: cfl
    integer :: i, j, k, e
    cfl = 0d0

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             do e = 1, nelv
                ur = ( u(i,j,k,e)*drdx(i,j,k,e) &
                   +   v(i,j,k,e)*drdy(i,j,k,e) &
                   +   w(i,j,k,e)*drdz(i,j,k,e) ) * jacinv(i,j,k,e)
                us = ( u(i,j,k,e)*dsdx(i,j,k,e) &
                   +   v(i,j,k,e)*dsdy(i,j,k,e) &
                   +   w(i,j,k,e)*dsdz(i,j,k,e) ) * jacinv(i,j,k,e)
                ut = ( u(i,j,k,e)*dtdx(i,j,k,e) &
                   +   v(i,j,k,e)*dtdy(i,j,k,e) &
                   +   w(i,j,k,e)*dtdz(i,j,k,e) ) * jacinv(i,j,k,e)

                cflr = abs(dt*ur*dr_inv(i))
                cfls = abs(dt*us*ds_inv(j))
                cflt = abs(dt*ut*dt_inv(k))

                cflm = cflr + cfls + cflt
                cfl  = max(cfl, cflm)
             end do
          end do
       end do
    end do

  end function sx_cfl_lx4

  function sx_cfl_lx3(dt, u, v, w, drdx, dsdx, dtdx, drdy, dsdy, dtdy, &
       drdz, dsdz, dtdz, dr_inv, ds_inv, dt_inv, &
       jacinv, nelv) result(cfl)
    integer, parameter :: lx = 3
    integer, intent(in) :: nelv
    real(kind=rp), intent(in) :: dt
    real(kind=rp), dimension(lx, lx, lx, nelv) ::  u, v, w
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx), intent(in) :: dr_inv, ds_inv, dt_inv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp) :: cflr, cfls, cflt, cflm
    real(kind=rp) :: ur, us, ut
    real(kind=rp) :: cfl
    integer :: i, j, k, e
    cfl = 0d0

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             do e = 1, nelv
                ur = ( u(i,j,k,e)*drdx(i,j,k,e) &
                   +   v(i,j,k,e)*drdy(i,j,k,e) &
                   +   w(i,j,k,e)*drdz(i,j,k,e) ) * jacinv(i,j,k,e)
                us = ( u(i,j,k,e)*dsdx(i,j,k,e) &
                   +   v(i,j,k,e)*dsdy(i,j,k,e) &
                   +   w(i,j,k,e)*dsdz(i,j,k,e) ) * jacinv(i,j,k,e)
                ut = ( u(i,j,k,e)*dtdx(i,j,k,e) &
                   +   v(i,j,k,e)*dtdy(i,j,k,e) &
                   +   w(i,j,k,e)*dtdz(i,j,k,e) ) * jacinv(i,j,k,e)

                cflr = abs(dt*ur*dr_inv(i))
                cfls = abs(dt*us*ds_inv(j))
                cflt = abs(dt*ut*dt_inv(k))

                cflm = cflr + cfls + cflt
                cfl  = max(cfl, cflm)
             end do
          end do
       end do
    end do

  end function sx_cfl_lx3

  function sx_cfl_lx2(dt, u, v, w, drdx, dsdx, dtdx, drdy, dsdy, dtdy, &
       drdz, dsdz, dtdz, dr_inv, ds_inv, dt_inv, &
       jacinv, nelv) result(cfl)
    integer, parameter :: lx = 2
    integer, intent(in) :: nelv
    real(kind=rp), intent(in) :: dt
    real(kind=rp), dimension(lx, lx, lx, nelv) ::  u, v, w
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx), intent(in) :: dr_inv, ds_inv, dt_inv
    real(kind=rp), dimension(lx, lx, lx, nelv), intent(in) :: jacinv
    real(kind=rp) :: cflr, cfls, cflt, cflm
    real(kind=rp) :: ur, us, ut
    real(kind=rp) :: cfl
    integer :: i, j, k, e
    cfl = 0d0

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             do e = 1, nelv
                ur = ( u(i,j,k,e)*drdx(i,j,k,e) &
                   +   v(i,j,k,e)*drdy(i,j,k,e) &
                   +   w(i,j,k,e)*drdz(i,j,k,e) ) * jacinv(i,j,k,e)
                us = ( u(i,j,k,e)*dsdx(i,j,k,e) &
                   +   v(i,j,k,e)*dsdy(i,j,k,e) &
                   +   w(i,j,k,e)*dsdz(i,j,k,e) ) * jacinv(i,j,k,e)
                ut = ( u(i,j,k,e)*dtdx(i,j,k,e) &
                   +   v(i,j,k,e)*dtdy(i,j,k,e) &
                   +   w(i,j,k,e)*dtdz(i,j,k,e) ) * jacinv(i,j,k,e)

                cflr = abs(dt*ur*dr_inv(i))
                cfls = abs(dt*us*ds_inv(j))
                cflt = abs(dt*ut*dt_inv(k))

                cflm = cflr + cfls + cflt
                cfl  = max(cfl, cflm)
             end do
          end do
       end do
    end do

  end function sx_cfl_lx2

end submodule sx_cfl
