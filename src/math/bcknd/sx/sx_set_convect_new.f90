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
module sx_set_convect_new
  use num_types, only : rp
  implicit none

contains

  subroutine sx_set_convect_new_lx(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n, lx)
    integer, intent(in) :: n, lx
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    integer :: e, i

    do i = 1, lx * lx * lx
       do e = 1, n
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

  end subroutine sx_set_convect_new_lx

  subroutine sx_set_convect_new_lx18(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 18
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    integer :: e, i

    do i = 1, lx * lx * lx
       do e = 1, n
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

  end subroutine sx_set_convect_new_lx18

  subroutine sx_set_convect_new_lx17(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 17
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    integer :: e, i

    do i = 1, lx * lx * lx
       do e = 1, n
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

  end subroutine sx_set_convect_new_lx17

  subroutine sx_set_convect_new_lx16(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 16
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    integer :: e, i

    do i = 1, lx * lx * lx
       do e = 1, n
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

  end subroutine sx_set_convect_new_lx16

  subroutine sx_set_convect_new_lx15(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 15
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    integer :: e, i

    do i = 1, lx * lx * lx
       do e = 1, n
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

  end subroutine sx_set_convect_new_lx15

  subroutine sx_set_convect_new_lx14(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 14
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    integer :: e, i

    do i = 1, lx * lx * lx
       do e = 1, n
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

  end subroutine sx_set_convect_new_lx14

  subroutine sx_set_convect_new_lx13(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 13
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    integer :: e, i

    do i = 1, lx * lx * lx
       do e = 1, n
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

  end subroutine sx_set_convect_new_lx13

  subroutine sx_set_convect_new_lx12(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 12
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    integer :: e, i

    do i = 1, lx * lx * lx
       do e = 1, n
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

  end subroutine sx_set_convect_new_lx12

  subroutine sx_set_convect_new_lx11(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 11
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    integer :: e, i

    do i = 1, lx * lx * lx
       do e = 1, n
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

  end subroutine sx_set_convect_new_lx11

  subroutine sx_set_convect_new_lx10(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 10
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    integer :: e, i

    do i = 1, lx * lx * lx
       do e = 1, n
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

  end subroutine sx_set_convect_new_lx10

  subroutine sx_set_convect_new_lx9(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 9
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    integer :: e, i

    do i = 1, lx * lx * lx
       do e = 1, n
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

  end subroutine sx_set_convect_new_lx9

  subroutine sx_set_convect_new_lx8(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 8
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    integer :: e, i

    do i = 1, lx * lx * lx
       do e = 1, n
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

  end subroutine sx_set_convect_new_lx8 

  subroutine sx_set_convect_new_lx7(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 7
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    integer :: e, i

    do i = 1, lx * lx * lx
       do e = 1, n
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

  end subroutine sx_set_convect_new_lx7

  subroutine sx_set_convect_new_lx6(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 6
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    integer :: e, i

    do i = 1, lx * lx * lx
       do e = 1, n
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

end subroutine sx_set_convect_new_lx6

subroutine sx_set_convect_new_lx5(cr, cs, ct, cx, cy, cz, &
     drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
  integer, parameter :: lx = 5
  integer, intent(in) :: n
  real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: cr, cs, ct
  real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: cx, cy, cz
  real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
  real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
  real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz
  real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
  integer :: e, i

  do i = 1, lx * lx * lx
     do e = 1, n
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

  end subroutine sx_set_convect_new_lx5
  
  subroutine sx_set_convect_new_lx4(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 4
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    integer :: e, i

    do i = 1, lx * lx * lx
       do e = 1, n
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

  end subroutine sx_set_convect_new_lx4
  
  subroutine sx_set_convect_new_lx3(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 3
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    integer :: e, i

    do i = 1, lx * lx * lx
       do e = 1, n
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

  end subroutine sx_set_convect_new_lx3 

  subroutine sx_set_convect_new_lx2(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n)
    integer, parameter :: lx = 2
    integer, intent(in) :: n
    real(kind=rp), dimension(lx,lx,lx,n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx,lx,lx,n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx), intent(in) :: w3(lx,lx,lx)
    integer :: e, i

    do i = 1, lx * lx * lx
       do e = 1, n
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

  end subroutine sx_set_convect_new_lx2

end module sx_set_convect_new
