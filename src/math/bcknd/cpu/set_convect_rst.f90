! Copyright (c) 2024, The Neko Authors
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
!> Gradient kernels
module cpu_set_convect_rst
  use num_types, only : rp
  implicit none

contains

  subroutine cpu_set_convect_rst_lx(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n, lx)
    integer, intent(in) :: n, lx
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    integer :: e, i

    do e = 1, n
       do i = 1, lx * lx * lx
          cr(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * drdx(i,1,1,e) &
                        + cy(i,1,1,e) * drdy(i,1,1,e) &
                        + cz(i,1,1,e) * drdz(i,1,1,e) )
          cs(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dsdx(i,1,1,e) &
                        + cy(i,1,1,e) * dsdy(i,1,1,e) &
                        + cz(i,1,1,e) * dsdz(i,1,1,e))
          ct(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dtdx(i,1,1,e) &
                        + cy(i,1,1,e) * dtdy(i,1,1,e) &
                        + cz(i,1,1,e) * dtdz(i,1,1,e))
       end do
    end do

  end subroutine cpu_set_convect_rst_lx

  subroutine cpu_set_convect_rst_lx18(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 18
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    integer :: e, i

    do e = 1, n
       do i = 1, lx * lx * lx
          cr(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * drdx(i,1,1,e) &
                        + cy(i,1,1,e) * drdy(i,1,1,e) &
                        + cz(i,1,1,e) * drdz(i,1,1,e) )
          cs(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dsdx(i,1,1,e) &
                        + cy(i,1,1,e) * dsdy(i,1,1,e) &
                        + cz(i,1,1,e) * dsdz(i,1,1,e))
          ct(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dtdx(i,1,1,e) &
                        + cy(i,1,1,e) * dtdy(i,1,1,e) &
                        + cz(i,1,1,e) * dtdz(i,1,1,e))
       end do
    end do

  end subroutine cpu_set_convect_rst_lx18

  subroutine cpu_set_convect_rst_lx17(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 17
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    integer :: e, i

    do e = 1, n
       do i = 1, lx * lx * lx
          cr(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * drdx(i,1,1,e) &
                        + cy(i,1,1,e) * drdy(i,1,1,e) &
                        + cz(i,1,1,e) * drdz(i,1,1,e) )
          cs(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dsdx(i,1,1,e) &
                        + cy(i,1,1,e) * dsdy(i,1,1,e) &
                        + cz(i,1,1,e) * dsdz(i,1,1,e))
          ct(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dtdx(i,1,1,e) &
                        + cy(i,1,1,e) * dtdy(i,1,1,e) &
                        + cz(i,1,1,e) * dtdz(i,1,1,e))
       end do
    end do

  end subroutine cpu_set_convect_rst_lx17

  subroutine cpu_set_convect_rst_lx16(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 16
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    integer :: e, i

    do e = 1, n
       do i = 1, lx * lx * lx
          cr(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * drdx(i,1,1,e) &
                        + cy(i,1,1,e) * drdy(i,1,1,e) &
                        + cz(i,1,1,e) * drdz(i,1,1,e) )
          cs(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dsdx(i,1,1,e) &
                        + cy(i,1,1,e) * dsdy(i,1,1,e) &
                        + cz(i,1,1,e) * dsdz(i,1,1,e))
          ct(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dtdx(i,1,1,e) &
                        + cy(i,1,1,e) * dtdy(i,1,1,e) &
                        + cz(i,1,1,e) * dtdz(i,1,1,e))
       end do
    end do

  end subroutine cpu_set_convect_rst_lx16

  subroutine cpu_set_convect_rst_lx15(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 15
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    integer :: e, i

    do e = 1, n
       do i = 1, lx * lx * lx
          cr(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * drdx(i,1,1,e) &
                        + cy(i,1,1,e) * drdy(i,1,1,e) &
                        + cz(i,1,1,e) * drdz(i,1,1,e) )
          cs(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dsdx(i,1,1,e) &
                        + cy(i,1,1,e) * dsdy(i,1,1,e) &
                        + cz(i,1,1,e) * dsdz(i,1,1,e))
          ct(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dtdx(i,1,1,e) &
                        + cy(i,1,1,e) * dtdy(i,1,1,e) &
                        + cz(i,1,1,e) * dtdz(i,1,1,e))
       end do
    end do

  end subroutine cpu_set_convect_rst_lx15

  subroutine cpu_set_convect_rst_lx14(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 14
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    integer :: e, i

    do e = 1, n
       do i = 1, lx * lx * lx
          cr(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * drdx(i,1,1,e) &
                        + cy(i,1,1,e) * drdy(i,1,1,e) &
                        + cz(i,1,1,e) * drdz(i,1,1,e) )
          cs(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dsdx(i,1,1,e) &
                        + cy(i,1,1,e) * dsdy(i,1,1,e) &
                        + cz(i,1,1,e) * dsdz(i,1,1,e))
          ct(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dtdx(i,1,1,e) &
                        + cy(i,1,1,e) * dtdy(i,1,1,e) &
                        + cz(i,1,1,e) * dtdz(i,1,1,e))
       end do
    end do

  end subroutine cpu_set_convect_rst_lx14

  subroutine cpu_set_convect_rst_lx13(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 13
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    integer :: e, i

    do e = 1, n
       do i = 1, lx * lx * lx
          cr(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * drdx(i,1,1,e) &
                        + cy(i,1,1,e) * drdy(i,1,1,e) &
                        + cz(i,1,1,e) * drdz(i,1,1,e) )
          cs(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dsdx(i,1,1,e) &
                        + cy(i,1,1,e) * dsdy(i,1,1,e) &
                        + cz(i,1,1,e) * dsdz(i,1,1,e))
          ct(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dtdx(i,1,1,e) &
                        + cy(i,1,1,e) * dtdy(i,1,1,e) &
                        + cz(i,1,1,e) * dtdz(i,1,1,e))
       end do
    end do

  end subroutine cpu_set_convect_rst_lx13

  subroutine cpu_set_convect_rst_lx12(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 12
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    integer :: e, i

    do e = 1, n
       do i = 1, lx * lx * lx
          cr(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * drdx(i,1,1,e) &
                        + cy(i,1,1,e) * drdy(i,1,1,e) &
                        + cz(i,1,1,e) * drdz(i,1,1,e) )
          cs(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dsdx(i,1,1,e) &
                        + cy(i,1,1,e) * dsdy(i,1,1,e) &
                        + cz(i,1,1,e) * dsdz(i,1,1,e))
          ct(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dtdx(i,1,1,e) &
                        + cy(i,1,1,e) * dtdy(i,1,1,e) &
                        + cz(i,1,1,e) * dtdz(i,1,1,e))
       end do
    end do

  end subroutine cpu_set_convect_rst_lx12

  subroutine cpu_set_convect_rst_lx11(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 11
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    integer :: e, i

    do e = 1, n
       do i = 1, lx * lx * lx
          cr(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * drdx(i,1,1,e) &
                        + cy(i,1,1,e) * drdy(i,1,1,e) &
                        + cz(i,1,1,e) * drdz(i,1,1,e) )
          cs(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dsdx(i,1,1,e) &
                        + cy(i,1,1,e) * dsdy(i,1,1,e) &
                        + cz(i,1,1,e) * dsdz(i,1,1,e))
          ct(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dtdx(i,1,1,e) &
                        + cy(i,1,1,e) * dtdy(i,1,1,e) &
                        + cz(i,1,1,e) * dtdz(i,1,1,e))
       end do
    end do

  end subroutine cpu_set_convect_rst_lx11

  subroutine cpu_set_convect_rst_lx10(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 10
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    integer :: e, i

    do e = 1, n
       do i = 1, lx * lx * lx
          cr(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * drdx(i,1,1,e) &
                        + cy(i,1,1,e) * drdy(i,1,1,e) &
                        + cz(i,1,1,e) * drdz(i,1,1,e) )
          cs(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dsdx(i,1,1,e) &
                        + cy(i,1,1,e) * dsdy(i,1,1,e) &
                        + cz(i,1,1,e) * dsdz(i,1,1,e))
          ct(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dtdx(i,1,1,e) &
                        + cy(i,1,1,e) * dtdy(i,1,1,e) &
                        + cz(i,1,1,e) * dtdz(i,1,1,e))
       end do
    end do

  end subroutine cpu_set_convect_rst_lx10

  subroutine cpu_set_convect_rst_lx9(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 9
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    integer :: e, i

    do e = 1, n
       do i = 1, lx * lx * lx
          cr(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * drdx(i,1,1,e) &
                        + cy(i,1,1,e) * drdy(i,1,1,e) &
                        + cz(i,1,1,e) * drdz(i,1,1,e) )
          cs(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dsdx(i,1,1,e) &
                        + cy(i,1,1,e) * dsdy(i,1,1,e) &
                        + cz(i,1,1,e) * dsdz(i,1,1,e))
          ct(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dtdx(i,1,1,e) &
                        + cy(i,1,1,e) * dtdy(i,1,1,e) &
                        + cz(i,1,1,e) * dtdz(i,1,1,e))
       end do
    end do

  end subroutine cpu_set_convect_rst_lx9

  subroutine cpu_set_convect_rst_lx8(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 8
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    integer :: e, i

    do e = 1, n
       do i = 1, lx * lx * lx
          cr(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * drdx(i,1,1,e) &
                        + cy(i,1,1,e) * drdy(i,1,1,e) &
                        + cz(i,1,1,e) * drdz(i,1,1,e) )
          cs(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dsdx(i,1,1,e) &
                        + cy(i,1,1,e) * dsdy(i,1,1,e) &
                        + cz(i,1,1,e) * dsdz(i,1,1,e))
          ct(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dtdx(i,1,1,e) &
                        + cy(i,1,1,e) * dtdy(i,1,1,e) &
                        + cz(i,1,1,e) * dtdz(i,1,1,e))
       end do
    end do

  end subroutine cpu_set_convect_rst_lx8

  subroutine cpu_set_convect_rst_lx7(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 7
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    integer :: e, i

    do e = 1, n
       do i = 1, lx * lx * lx
          cr(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * drdx(i,1,1,e) &
                        + cy(i,1,1,e) * drdy(i,1,1,e) &
                        + cz(i,1,1,e) * drdz(i,1,1,e) )
          cs(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dsdx(i,1,1,e) &
                        + cy(i,1,1,e) * dsdy(i,1,1,e) &
                        + cz(i,1,1,e) * dsdz(i,1,1,e))
          ct(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dtdx(i,1,1,e) &
                        + cy(i,1,1,e) * dtdy(i,1,1,e) &
                        + cz(i,1,1,e) * dtdz(i,1,1,e))
       end do
    end do

  end subroutine cpu_set_convect_rst_lx7

  subroutine cpu_set_convect_rst_lx6(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 6
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    integer :: e, i

    do e = 1, n
       do i = 1, lx * lx * lx
          cr(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * drdx(i,1,1,e) &
                        + cy(i,1,1,e) * drdy(i,1,1,e) &
                        + cz(i,1,1,e) * drdz(i,1,1,e) )
          cs(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dsdx(i,1,1,e) &
                        + cy(i,1,1,e) * dsdy(i,1,1,e) &
                        + cz(i,1,1,e) * dsdz(i,1,1,e))
          ct(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dtdx(i,1,1,e) &
                        + cy(i,1,1,e) * dtdy(i,1,1,e) &
                        + cz(i,1,1,e) * dtdz(i,1,1,e))
       end do
    end do

  end subroutine cpu_set_convect_rst_lx6

  subroutine cpu_set_convect_rst_lx5(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 5
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    integer :: e, i

    do e = 1, n
       do i = 1, lx * lx * lx
          cr(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * drdx(i,1,1,e) &
                        + cy(i,1,1,e) * drdy(i,1,1,e) &
                        + cz(i,1,1,e) * drdz(i,1,1,e) )
          cs(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dsdx(i,1,1,e) &
                        + cy(i,1,1,e) * dsdy(i,1,1,e) &
                        + cz(i,1,1,e) * dsdz(i,1,1,e))
          ct(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dtdx(i,1,1,e) &
                        + cy(i,1,1,e) * dtdy(i,1,1,e) &
                        + cz(i,1,1,e) * dtdz(i,1,1,e))
       end do
    end do

  end subroutine cpu_set_convect_rst_lx5

  subroutine cpu_set_convect_rst_lx4(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 4
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    integer :: e, i

    do e = 1, n
       do i = 1, lx * lx * lx
          cr(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * drdx(i,1,1,e) &
                        + cy(i,1,1,e) * drdy(i,1,1,e) &
                        + cz(i,1,1,e) * drdz(i,1,1,e) )
          cs(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dsdx(i,1,1,e) &
                        + cy(i,1,1,e) * dsdy(i,1,1,e) &
                        + cz(i,1,1,e) * dsdz(i,1,1,e))
          ct(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dtdx(i,1,1,e) &
                        + cy(i,1,1,e) * dtdy(i,1,1,e) &
                        + cz(i,1,1,e) * dtdz(i,1,1,e))
       end do
    end do

  end subroutine cpu_set_convect_rst_lx4

  subroutine cpu_set_convect_rst_lx3(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 3
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    integer :: e, i

    do e = 1, n
       do i = 1, lx * lx * lx
          cr(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * drdx(i,1,1,e) &
                        + cy(i,1,1,e) * drdy(i,1,1,e) &
                        + cz(i,1,1,e) * drdz(i,1,1,e) )
          cs(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dsdx(i,1,1,e) &
                        + cy(i,1,1,e) * dsdy(i,1,1,e) &
                        + cz(i,1,1,e) * dsdz(i,1,1,e))
          ct(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dtdx(i,1,1,e) &
                        + cy(i,1,1,e) * dtdy(i,1,1,e) &
                        + cz(i,1,1,e) * dtdz(i,1,1,e))
       end do
    end do

  end subroutine cpu_set_convect_rst_lx3

  subroutine cpu_set_convect_rst_lx2(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 2
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    integer :: e, i

    do e = 1, n
       do i = 1, lx * lx * lx
          cr(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * drdx(i,1,1,e) &
                        + cy(i,1,1,e) * drdy(i,1,1,e) &
                        + cz(i,1,1,e) * drdz(i,1,1,e) )
          cs(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dsdx(i,1,1,e) &
                        + cy(i,1,1,e) * dsdy(i,1,1,e) &
                        + cz(i,1,1,e) * dsdz(i,1,1,e))
          ct(i,1,1,e) = w3(i,1,1) &
                      * ( cx(i,1,1,e) * dtdx(i,1,1,e) &
                        + cy(i,1,1,e) * dtdy(i,1,1,e) &
                        + cz(i,1,1,e) * dtdz(i,1,1,e))
       end do
    end do

  end subroutine cpu_set_convect_rst_lx2

end module cpu_set_convect_rst
