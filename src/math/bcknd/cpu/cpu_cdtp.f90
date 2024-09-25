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
!> DT*X kernels
submodule (opr_cpu) cpu_cdtp
  implicit none

contains

  module subroutine opr_cpu_cdtp(dtx, x, dr, ds, dt, coef, e_start, e_end)
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: e_start, e_end
    real(kind=rp), intent(inout) :: dtx(coef%Xh%lxyz, e_end - e_start + 1)
    real(kind=rp), intent(inout) :: x(coef%Xh%lxyz, e_end - e_start + 1)
    real(kind=rp), intent(in) :: dr(coef%Xh%lxyz, e_end - e_start + 1)
    real(kind=rp), intent(in) :: ds(coef%Xh%lxyz, e_end - e_start + 1)
    real(kind=rp), intent(in) :: dt(coef%Xh%lxyz, e_end - e_start + 1)
    integer :: e_len
    e_len = e_end - e_start + 1

    if (e_len .eq. 1) then
       call opr_cpu_cdtp_single(dtx, x, dr, ds, dt, coef, e_start)
    else
       call opr_cpu_cdtp_many(dtx, x, dr, ds, dt, coef, e_start, e_len)
    end if

  end subroutine opr_cpu_cdtp

  subroutine opr_cpu_cdtp_many(dtx, x, dr, ds, dt, coef, e_start, e_len)
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: e_start, e_len
    real(kind=rp), intent(inout) :: dtx(coef%Xh%lxyz, e_len)
    real(kind=rp), intent(inout) :: x(coef%Xh%lxyz, e_len)
    real(kind=rp), intent(in) :: dr(coef%Xh%lxyz, e_len)
    real(kind=rp), intent(in) :: ds(coef%Xh%lxyz, e_len)
    real(kind=rp), intent(in) :: dt(coef%Xh%lxyz, e_len)

    associate(Xh => coef%Xh)
      select case (Xh%lx)
        case (14)
         call cpu_cdtp_lx14(dtx, x, &
              dr(1, e_start), ds(1, e_start), dt(1, e_start), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3, e_len)
        case (13)
         call cpu_cdtp_lx13(dtx, x, &
              dr(1, e_start), ds(1, e_start), dt(1, e_start), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3, e_len)
        case (12)
         call cpu_cdtp_lx12(dtx, x, &
              dr(1, e_start), ds(1, e_start), dt(1, e_start), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3, e_len)
        case (11)
         call cpu_cdtp_lx11(dtx, x, &
              dr(1, e_start), ds(1, e_start), dt(1, e_start), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3, e_len)
        case (10)
         call cpu_cdtp_lx10(dtx, x, &
              dr(1, e_start), ds(1, e_start), dt(1, e_start), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3, e_len)
        case (9)
         call cpu_cdtp_lx9(dtx, x, &
              dr(1, e_start), ds(1, e_start), dt(1, e_start), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3, e_len)
        case (8)
         call cpu_cdtp_lx8(dtx, x, &
              dr(1, e_start), ds(1, e_start), dt(1, e_start), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3, e_len)
        case (7)
         call cpu_cdtp_lx7(dtx, x, &
              dr(1, e_start), ds(1, e_start), dt(1, e_start), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3, e_len)
        case (6)
         call cpu_cdtp_lx6(dtx, x, &
              dr(1, e_start), ds(1, e_start), dt(1, e_start), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3, e_len)
        case (5)
         call cpu_cdtp_lx5(dtx, x, &
              dr(1, e_start), ds(1, e_start), dt(1, e_start), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3, e_len)
        case (4)
         call cpu_cdtp_lx4(dtx, x, &
              dr(1, e_start), ds(1, e_start), dt(1, e_start), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3, e_len)
        case (3)
         call cpu_cdtp_lx3(dtx, x, &
              dr(1, e_start), ds(1, e_start), dt(1, e_start), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3, e_len)
        case (2)
         call cpu_cdtp_lx2(dtx, x, &
              dr(1, e_start), ds(1, e_start), dt(1, e_start), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3, e_len)
        case default
         call cpu_cdtp_lx(dtx, x, &
              dr(1, e_start), ds(1, e_start), dt(1, e_start), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3, e_len, Xh%lx)
      end select
    end associate

  end subroutine opr_cpu_cdtp_many

  subroutine opr_cpu_cdtp_single(dtx, x, dr, ds, dt, coef, e)
    integer, parameter :: e_len = 1
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: e
    real(kind=rp), intent(inout) :: dtx(coef%Xh%lxyz, e_len)
    real(kind=rp), intent(inout) :: x(coef%Xh%lxyz, e_len)
    real(kind=rp), intent(in) :: dr(coef%Xh%lxyz, e_len)
    real(kind=rp), intent(in) :: ds(coef%Xh%lxyz, e_len)
    real(kind=rp), intent(in) :: dt(coef%Xh%lxyz, e_len)

    associate(Xh => coef%Xh)
      select case (Xh%lx)
        case (14)
         call cpu_cdtp_lx14_single(dtx, x, dr(1,e), ds(1,e), dt(1,e), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3)
        case (13)
         call cpu_cdtp_lx13_single(dtx, x, dr(1,e), ds(1,e), dt(1,e), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3)
        case (12)
         call cpu_cdtp_lx12_single(dtx, x, dr(1,e), ds(1,e), dt(1,e), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3)
        case (11)
         call cpu_cdtp_lx11_single(dtx, x, dr(1,e), ds(1,e), dt(1,e), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3)
        case (10)
         call cpu_cdtp_lx10_single(dtx, x, dr(1,e), ds(1,e), dt(1,e), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3)
        case (9)
         call cpu_cdtp_lx9_single(dtx, x, dr(1,e), ds(1,e), dt(1,e), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3)
        case (8)
         call cpu_cdtp_lx8_single(dtx, x, dr(1,e), ds(1,e), dt(1,e), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3)
        case (7)
         call cpu_cdtp_lx7_single(dtx, x, dr(1,e), ds(1,e), dt(1,e), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3)
        case (6)
         call cpu_cdtp_lx6_single(dtx, x, dr(1,e), ds(1,e), dt(1,e), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3)
        case (5)
         call cpu_cdtp_lx5_single(dtx, x, dr(1,e), ds(1,e), dt(1,e), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3)
        case (4)
         call cpu_cdtp_lx4_single(dtx, x, dr(1,e), ds(1,e), dt(1,e), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3)
        case (3)
         call cpu_cdtp_lx3_single(dtx, x, dr(1,e), ds(1,e), dt(1,e), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3)
        case (2)
         call cpu_cdtp_lx2_single(dtx, x, dr(1,e), ds(1,e), dt(1,e), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3)
        case default
         call cpu_cdtp_lx_single(dtx, x, dr(1,e), ds(1,e), dt(1,e), &
              Xh%dxt, Xh%dyt, Xh%dzt, Xh%w3, Xh%lx)
      end select
    end associate

  end subroutine opr_cpu_cdtp_single

  subroutine cpu_cdtp_lx(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3, nel, lx)
    integer, intent(in) :: nel, lx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l

    do e = 1, nel

       do i = 1, lx*lx*lx
          wx(i,1,1) = x(i,1,1,e) * w3(i,1,1)
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
  end subroutine cpu_cdtp_lx

  subroutine cpu_cdtp_lx14(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3, nel)
    integer, parameter :: lx = 14
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: e, i, j, k

    do e = 1, nel

       do i = 1, lx*lx*lx
          wx(i,1,1) = x(i,1,1,e) * w3(i,1,1)
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
  end subroutine cpu_cdtp_lx14

  subroutine cpu_cdtp_lx13(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3, nel)
    integer, parameter :: lx = 13
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: e, i, j, k

    do e = 1, nel

       do i = 1, lx*lx*lx
          wx(i,1,1) = x(i,1,1,e) * w3(i,1,1)
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
  end subroutine cpu_cdtp_lx13

  subroutine cpu_cdtp_lx12(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3, nel)
    integer, parameter :: lx = 12
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: e, i, j, k

    do e = 1, nel

       do i = 1, lx*lx*lx
          wx(i,1,1) = x(i,1,1,e) * w3(i,1,1)
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

  subroutine cpu_cdtp_lx11(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3, nel)
    integer, parameter :: lx = 11
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: e, i, j, k

    do e = 1, nel

       do i = 1, lx*lx*lx
          wx(i,1,1) = x(i,1,1,e) * w3(i,1,1)
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

  subroutine cpu_cdtp_lx10(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3, nel)
    integer, parameter :: lx = 10
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: e, i, j, k

    do e = 1, nel

       do i = 1, lx*lx*lx
          wx(i,1,1) = x(i,1,1,e) * w3(i,1,1)
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

  subroutine cpu_cdtp_lx9(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3, nel)
    integer, parameter :: lx = 9
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: e, i, j, k

    do e = 1, nel

       do i = 1, lx*lx*lx
          wx(i,1,1) = x(i,1,1,e) * w3(i,1,1)
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

  subroutine cpu_cdtp_lx8(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3, nel)
    integer, parameter :: lx = 8
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: e, i, j, k

    do e = 1, nel

       do i = 1, lx*lx*lx
          wx(i,1,1) = x(i,1,1,e) * w3(i,1,1)
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

  subroutine cpu_cdtp_lx7(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3, nel)
    integer, parameter :: lx = 7
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: e, i, j, k

    do e = 1, nel

       do i = 1, lx*lx*lx
          wx(i,1,1) = x(i,1,1,e) * w3(i,1,1)
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

  subroutine cpu_cdtp_lx6(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3, nel)
    integer, parameter :: lx = 6
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: e, i, j, k

    do e = 1, nel

       do i = 1, lx*lx*lx
          wx(i,1,1) = x(i,1,1,e) * w3(i,1,1)
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

  subroutine cpu_cdtp_lx5(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3, nel)
    integer, parameter :: lx = 5
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: e, i, j, k

    do e = 1, nel

       do i = 1, lx*lx*lx
          wx(i,1,1) = x(i,1,1,e) * w3(i,1,1)
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

  subroutine cpu_cdtp_lx4(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3, nel)
    integer, parameter :: lx = 4
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: e, i, j, k

    do e = 1, nel

       do i = 1, lx*lx*lx
          wx(i,1,1) = x(i,1,1,e) * w3(i,1,1)
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

  subroutine cpu_cdtp_lx3(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3, nel)
    integer, parameter :: lx = 3
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: e, i, j, k

    do e = 1, nel

       do i = 1, lx*lx*lx
          wx(i,1,1) = x(i,1,1,e) * w3(i,1,1)
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

  subroutine cpu_cdtp_lx2(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3, nel)
    integer, parameter :: lx = 2
    integer, intent(in) :: nel
    real(kind=rp), dimension(lx, lx, lx, nel), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx, nel), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: e, i, j, k

    do e = 1, nel

       do i = 1, lx*lx*lx
          wx(i,1,1) = x(i,1,1,e) * w3(i,1,1)
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

  subroutine cpu_cdtp_lx_single(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3, lx)
    integer, intent(in) :: lx
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    real(kind=rp) :: tmp
    integer :: i, j, k, l

    do i = 1, lx*lx*lx
       wx(i,1,1) = x(i,1,1) * w3(i,1,1)
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * dr(i,1,1)
    end do

    do j = 1, lx * lx
       do i = 1, lx
          tmp = 0.0_rp
          !DIR$ LOOP_INFO MIN_TRIPS(15)
          do k = 1, lx
             tmp = tmp + dxt(i,k) * ta1(k,j,1)
          end do
          dtx(i,j,1) = tmp
       end do
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * ds(i,1,1)
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             tmp = 0.0_rp
             !DIR$ LOOP_INFO MIN_TRIPS(15)
             do l = 1, lx
                tmp = tmp + dyt(j,l) * ta1(i,l,k)
             end do
             dtx(i,j,k) = dtx(i,j,k) + tmp
          end do
       end do
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * dt(i,1,1)
    end do

    do k = 1, lx
       do i = 1, lx*lx
          tmp = 0.0_rp
          !DIR$ LOOP_INFO MIN_TRIPS(15)
          do l = 1, lx
             tmp = tmp + dzt(k,l) * ta1(i,1,l)
          end do
          dtx(i,1,k) = dtx(i,1,k) + tmp
       end do
    end do

  end subroutine cpu_cdtp_lx_single

  subroutine cpu_cdtp_lx14_single(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3)
    integer, parameter :: lx = 14
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: i, j, k

    do i = 1, lx*lx*lx
       wx(i,1,1) = x(i,1,1) * w3(i,1,1)
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * dr(i,1,1)
    end do

    do j = 1, lx * lx
       do i = 1, lx
          dtx(i,j,1) = dxt(i,1) * ta1(1,j,1) &
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
       ta1(i,1,1) = wx(i,1,1) * ds(i,1,1)
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             dtx(i,j,k) = dtx(i,j,k) &
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
       ta1(i,1,1) = wx(i,1,1) * dt(i,1,1)
    end do

    do k = 1, lx
       do i = 1, lx*lx
          dtx(i,1,k) = dtx(i,1,k) &
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

  end subroutine cpu_cdtp_lx14_single

  subroutine cpu_cdtp_lx13_single(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3)
    integer, parameter :: lx = 13
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: i, j, k

    do i = 1, lx*lx*lx
       wx(i,1,1) = x(i,1,1) * w3(i,1,1)
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * dr(i,1,1)
    end do

    do j = 1, lx * lx
       do i = 1, lx
          dtx(i,j,1) = dxt(i,1) * ta1(1,j,1) &
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
       ta1(i,1,1) = wx(i,1,1) * ds(i,1,1)
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             dtx(i,j,k) = dtx(i,j,k) &
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
       ta1(i,1,1) = wx(i,1,1) * dt(i,1,1)
    end do

    do k = 1, lx
       do i = 1, lx*lx
          dtx(i,1,k) = dtx(i,1,k) &
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

  end subroutine cpu_cdtp_lx13_single

  subroutine cpu_cdtp_lx12_single(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3)
    integer, parameter :: lx = 12
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: i, j, k

    do i = 1, lx*lx*lx
       wx(i,1,1) = x(i,1,1) * w3(i,1,1)
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * dr(i,1,1)
    end do

    do j = 1, lx * lx
       do i = 1, lx
          dtx(i,j,1) = dxt(i,1) * ta1(1,j,1) &
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
       ta1(i,1,1) = wx(i,1,1) * ds(i,1,1)
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             dtx(i,j,k) = dtx(i,j,k) &
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
       ta1(i,1,1) = wx(i,1,1) * dt(i,1,1)
    end do

    do k = 1, lx
       do i = 1, lx*lx
          dtx(i,1,k) = dtx(i,1,k) &
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

  end subroutine cpu_cdtp_lx12_single

  subroutine cpu_cdtp_lx11_single(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3)
    integer, parameter :: lx = 11
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: i, j, k

    do i = 1, lx*lx*lx
       wx(i,1,1) = x(i,1,1) * w3(i,1,1)
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * dr(i,1,1)
    end do

    do j = 1, lx * lx
       do i = 1, lx
          dtx(i,j,1) = dxt(i,1) * ta1(1,j,1) &
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
       ta1(i,1,1) = wx(i,1,1) * ds(i,1,1)
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             dtx(i,j,k) = dtx(i,j,k) &
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
       ta1(i,1,1) = wx(i,1,1) * dt(i,1,1)
    end do

    do k = 1, lx
       do i = 1, lx*lx
          dtx(i,1,k) = dtx(i,1,k) &
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

  end subroutine cpu_cdtp_lx11_single

  subroutine cpu_cdtp_lx10_single(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3)
    integer, parameter :: lx = 10
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: i, j, k

    do i = 1, lx*lx*lx
       wx(i,1,1) = x(i,1,1) * w3(i,1,1)
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * dr(i,1,1)
    end do

    do j = 1, lx * lx
       do i = 1, lx
          dtx(i,j,1) = dxt(i,1) * ta1(1,j,1) &
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
       ta1(i,1,1) = wx(i,1,1) * ds(i,1,1)
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             dtx(i,j,k) = dtx(i,j,k) &
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
       ta1(i,1,1) = wx(i,1,1) * dt(i,1,1)
    end do

    do k = 1, lx
       do i = 1, lx*lx
          dtx(i,1,k) = dtx(i,1,k) &
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

  end subroutine cpu_cdtp_lx10_single

  subroutine cpu_cdtp_lx9_single(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3)
    integer, parameter :: lx = 9
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: i, j, k

    do i = 1, lx*lx*lx
       wx(i,1,1) = x(i,1,1) * w3(i,1,1)
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * dr(i,1,1)
    end do

    do j = 1, lx * lx
       do i = 1, lx
          dtx(i,j,1) = dxt(i,1) * ta1(1,j,1) &
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
       ta1(i,1,1) = wx(i,1,1) * ds(i,1,1)
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             dtx(i,j,k) = dtx(i,j,k) &
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
       ta1(i,1,1) = wx(i,1,1) * dt(i,1,1)
    end do

    do k = 1, lx
       do i = 1, lx*lx
          dtx(i,1,k) = dtx(i,1,k) &
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

  end subroutine cpu_cdtp_lx9_single

  subroutine cpu_cdtp_lx8_single(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3)
    integer, parameter :: lx = 8
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: i, j, k

    do i = 1, lx*lx*lx
       wx(i,1,1) = x(i,1,1) * w3(i,1,1)
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * dr(i,1,1)
    end do

    do j = 1, lx * lx
       do i = 1, lx
          dtx(i,j,1) = dxt(i,1) * ta1(1,j,1) &
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
       ta1(i,1,1) = wx(i,1,1) * ds(i,1,1)
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             dtx(i,j,k) = dtx(i,j,k) &
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
       ta1(i,1,1) = wx(i,1,1) * dt(i,1,1)
    end do

    do k = 1, lx
       do i = 1, lx*lx
          dtx(i,1,k) = dtx(i,1,k) &
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

  end subroutine cpu_cdtp_lx8_single

  subroutine cpu_cdtp_lx7_single(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3)
    integer, parameter :: lx = 7
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: i, j, k

    do i = 1, lx*lx*lx
       wx(i,1,1) = x(i,1,1) * w3(i,1,1)
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * dr(i,1,1)
    end do

    do j = 1, lx * lx
       do i = 1, lx
          dtx(i,j,1) = dxt(i,1) * ta1(1,j,1) &
               + dxt(i,2) * ta1(2,j,1) &
               + dxt(i,3) * ta1(3,j,1) &
               + dxt(i,4) * ta1(4,j,1) &
               + dxt(i,5) * ta1(5,j,1) &
               + dxt(i,6) * ta1(6,j,1) &
               + dxt(i,7) * ta1(7,j,1)
       end do
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * ds(i,1,1)
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             dtx(i,j,k) = dtx(i,j,k) &
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
       ta1(i,1,1) = wx(i,1,1) * dt(i,1,1)
    end do

    do k = 1, lx
       do i = 1, lx*lx
          dtx(i,1,k) = dtx(i,1,k) &
               + dzt(k,1) * ta1(i,1,1) &
               + dzt(k,2) * ta1(i,1,2) &
               + dzt(k,3) * ta1(i,1,3) &
               + dzt(k,4) * ta1(i,1,4) &
               + dzt(k,5) * ta1(i,1,5) &
               + dzt(k,6) * ta1(i,1,6) &
               + dzt(k,7) * ta1(i,1,7)
       end do
    end do

  end subroutine cpu_cdtp_lx7_single

  subroutine cpu_cdtp_lx6_single(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3)
    integer, parameter :: lx = 6
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: i, j, k

    do i = 1, lx*lx*lx
       wx(i,1,1) = x(i,1,1) * w3(i,1,1)
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * dr(i,1,1)
    end do

    do j = 1, lx * lx
       do i = 1, lx
          dtx(i,j,1) = dxt(i,1) * ta1(1,j,1) &
               + dxt(i,2) * ta1(2,j,1) &
               + dxt(i,3) * ta1(3,j,1) &
               + dxt(i,4) * ta1(4,j,1) &
               + dxt(i,5) * ta1(5,j,1) &
               + dxt(i,6) * ta1(6,j,1)
       end do
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * ds(i,1,1)
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             dtx(i,j,k) = dtx(i,j,k) &
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
       ta1(i,1,1) = wx(i,1,1) * dt(i,1,1)
    end do

    do k = 1, lx
       do i = 1, lx*lx
          dtx(i,1,k) = dtx(i,1,k) &
               + dzt(k,1) * ta1(i,1,1) &
               + dzt(k,2) * ta1(i,1,2) &
               + dzt(k,3) * ta1(i,1,3) &
               + dzt(k,4) * ta1(i,1,4) &
               + dzt(k,5) * ta1(i,1,5) &
               + dzt(k,6) * ta1(i,1,6)
       end do
    end do

  end subroutine cpu_cdtp_lx6_single

  subroutine cpu_cdtp_lx5_single(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3)
    integer, parameter :: lx = 5
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: i, j, k

    do i = 1, lx*lx*lx
       wx(i,1,1) = x(i,1,1) * w3(i,1,1)
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * dr(i,1,1)
    end do

    do j = 1, lx * lx
       do i = 1, lx
          dtx(i,j,1) = dxt(i,1) * ta1(1,j,1) &
               + dxt(i,2) * ta1(2,j,1) &
               + dxt(i,3) * ta1(3,j,1) &
               + dxt(i,4) * ta1(4,j,1) &
               + dxt(i,5) * ta1(5,j,1)
       end do
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * ds(i,1,1)
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             dtx(i,j,k) = dtx(i,j,k) &
                  + dyt(j,1) * ta1(i,1,k) &
                  + dyt(j,2) * ta1(i,2,k) &
                  + dyt(j,3) * ta1(i,3,k) &
                  + dyt(j,4) * ta1(i,4,k) &
                  + dyt(j,5) * ta1(i,5,k)
          end do
       end do
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * dt(i,1,1)
    end do

    do k = 1, lx
       do i = 1, lx*lx
          dtx(i,1,k) = dtx(i,1,k) &
               + dzt(k,1) * ta1(i,1,1) &
               + dzt(k,2) * ta1(i,1,2) &
               + dzt(k,3) * ta1(i,1,3) &
               + dzt(k,4) * ta1(i,1,4) &
               + dzt(k,5) * ta1(i,1,5)
       end do
    end do

  end subroutine cpu_cdtp_lx5_single

  subroutine cpu_cdtp_lx4_single(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3)
    integer, parameter :: lx = 4
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: i, j, k

    do i = 1, lx*lx*lx
       wx(i,1,1) = x(i,1,1) * w3(i,1,1)
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * dr(i,1,1)
    end do

    do j = 1, lx * lx
       do i = 1, lx
          dtx(i,j,1) = dxt(i,1) * ta1(1,j,1) &
               + dxt(i,2) * ta1(2,j,1) &
               + dxt(i,3) * ta1(3,j,1) &
               + dxt(i,4) * ta1(4,j,1)
       end do
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * ds(i,1,1)
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             dtx(i,j,k) = dtx(i,j,k) &
                  + dyt(j,1) * ta1(i,1,k) &
                  + dyt(j,2) * ta1(i,2,k) &
                  + dyt(j,3) * ta1(i,3,k) &
                  + dyt(j,4) * ta1(i,4,k)
          end do
       end do
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * dt(i,1,1)
    end do

    do k = 1, lx
       do i = 1, lx*lx
          dtx(i,1,k) = dtx(i,1,k) &
               + dzt(k,1) * ta1(i,1,1) &
               + dzt(k,2) * ta1(i,1,2) &
               + dzt(k,3) * ta1(i,1,3) &
               + dzt(k,4) * ta1(i,1,4)
       end do
    end do

  end subroutine cpu_cdtp_lx4_single

  subroutine cpu_cdtp_lx3_single(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3)
    integer, parameter :: lx = 3
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: i, j, k

    do i = 1, lx*lx*lx
       wx(i,1,1) = x(i,1,1) * w3(i,1,1)
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * dr(i,1,1)
    end do

    do j = 1, lx * lx
       do i = 1, lx
          dtx(i,j,1) = dxt(i,1) * ta1(1,j,1) &
               + dxt(i,2) * ta1(2,j,1) &
               + dxt(i,3) * ta1(3,j,1)
       end do
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * ds(i,1,1)
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             dtx(i,j,k) = dtx(i,j,k) &
                  + dyt(j,1) * ta1(i,1,k) &
                  + dyt(j,2) * ta1(i,2,k) &
                  + dyt(j,3) * ta1(i,3,k)
          end do
       end do
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * dt(i,1,1)
    end do

    do k = 1, lx
       do i = 1, lx*lx
          dtx(i,1,k) = dtx(i,1,k) &
               + dzt(k,1) * ta1(i,1,1) &
               + dzt(k,2) * ta1(i,1,2) &
               + dzt(k,3) * ta1(i,1,3)
       end do
    end do

  end subroutine cpu_cdtp_lx3_single

  subroutine cpu_cdtp_lx2_single(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3)
    integer, parameter :: lx = 2
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: dtx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: x, dr, ds, dt
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp), intent(in), dimension(lx, lx) :: dxt, dyt, dzt
    real(kind=rp), dimension(lx, lx, lx) :: wx, ta1
    integer :: i, j, k

    do i = 1, lx*lx*lx
       wx(i,1,1) = x(i,1,1) * w3(i,1,1)
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * dr(i,1,1)
    end do

    do j = 1, lx * lx
       do i = 1, lx
          dtx(i,j,1) = dxt(i,1) * ta1(1,j,1) &
               + dxt(i,2) * ta1(2,j,1)
       end do
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * ds(i,1,1)
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             dtx(i,j,k) = dtx(i,j,k) &
                  + dyt(j,1) * ta1(i,1,k) &
                  + dyt(j,2) * ta1(i,2,k)
          end do
       end do
    end do

    do i = 1, lx*lx*lx
       ta1(i,1,1) = wx(i,1,1) * dt(i,1,1)
    end do

    do k = 1, lx
       do i = 1, lx*lx
          dtx(i,1,k) = dtx(i,1,k) &
               + dzt(k,1) * ta1(i,1,1) &
               + dzt(k,2) * ta1(i,1,2)
       end do
    end do

  end subroutine cpu_cdtp_lx2_single
end submodule cpu_cdtp
