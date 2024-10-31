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
submodule (opr_sx) sx_set_convect_rst
  use num_types, only : rp
  implicit none

contains

  module subroutine opr_sx_set_convect_rst(cr, cs, ct, cx, cy, cz, Xh, coef)
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    real(kind=rp), dimension(Xh%lxyz, coef%msh%nelv), &
                   intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(Xh%lxyz, coef%msh%nelv), &
                   intent(in) :: cx, cy, cz
    associate(drdx => coef%drdx, drdy => coef%drdy, drdz => coef%drdz, &
      dsdx => coef%dsdx, dsdy => coef%dsdy, dsdz => coef%dsdz, &
      dtdx => coef%dtdx, dtdy => coef%dtdy, dtdz => coef%dtdz, &
      nelv => coef%msh%nelv, lx => Xh%lx, w3 => Xh%w3)

      select case (lx)
      case (18)
         call sx_set_convect_rst_lx18(cr, cs, ct, cx, cy, cz, drdx, &
              dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case (17)
         call sx_set_convect_rst_lx17(cr, cs, ct, cx, cy, cz, drdx, &
              dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case (16)
         call sx_set_convect_rst_lx16(cr, cs, ct, cx, cy, cz, drdx, &
              dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case (15)
         call sx_set_convect_rst_lx15(cr, cs, ct, cx, cy, cz, drdx, &
              dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case (14)
         call sx_set_convect_rst_lx14(cr, cs, ct, cx, cy, cz, drdx, &
              dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case (13)
         call sx_set_convect_rst_lx13(cr, cs, ct, cx, cy, cz, drdx, &
              dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case (12)
         call sx_set_convect_rst_lx12(cr, cs, ct, cx, cy, cz, drdx, &
              dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case (11)
         call sx_set_convect_rst_lx11(cr, cs, ct, cx, cy, cz, drdx, &
              dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case (10)
         call sx_set_convect_rst_lx10(cr, cs, ct, cx, cy, cz, drdx, &
              dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case (9)
         call sx_set_convect_rst_lx9(cr, cs, ct, cx, cy, cz, drdx, &
              dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case (8)
         call sx_set_convect_rst_lx8(cr, cs, ct, cx, cy, cz, drdx, &
              dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case (7)
         call sx_set_convect_rst_lx7(cr, cs, ct, cx, cy, cz, drdx, &
              dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case (6)
         call sx_set_convect_rst_lx6(cr, cs, ct, cx, cy, cz, drdx, &
              dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case (5)
         call sx_set_convect_rst_lx5(cr, cs, ct, cx, cy, cz, drdx, &
              dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case (4)
         call sx_set_convect_rst_lx4(cr, cs, ct, cx, cy, cz, drdx, &
              dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case (3)
         call sx_set_convect_rst_lx3(cr, cs, ct, cx, cy, cz, drdx, &
              dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case (2)
         call sx_set_convect_rst_lx2(cr, cs, ct, cx, cy, cz, drdx, &
              dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv)
      case default
         call sx_set_convect_rst_lx(cr, cs, ct, cx, cy, cz, drdx, &
              dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, nelv, lx)
      end select
    end associate
  end subroutine opr_sx_set_convect_rst

  subroutine sx_set_convect_rst_lx(cr, cs, ct, cx, cy, cz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n, lx)
    integer, intent(in) :: n, lx
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: cr, cs, ct
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cx, cy, cz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
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

  end subroutine sx_set_convect_rst_lx

  subroutine sx_set_convect_rst_lx18(cr, cs, ct, cx, cy, cz, &
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

  end subroutine sx_set_convect_rst_lx18

  subroutine sx_set_convect_rst_lx17(cr, cs, ct, cx, cy, cz, &
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

  end subroutine sx_set_convect_rst_lx17

  subroutine sx_set_convect_rst_lx16(cr, cs, ct, cx, cy, cz, &
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

  end subroutine sx_set_convect_rst_lx16

  subroutine sx_set_convect_rst_lx15(cr, cs, ct, cx, cy, cz, &
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

  end subroutine sx_set_convect_rst_lx15

  subroutine sx_set_convect_rst_lx14(cr, cs, ct, cx, cy, cz, &
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

  end subroutine sx_set_convect_rst_lx14

  subroutine sx_set_convect_rst_lx13(cr, cs, ct, cx, cy, cz, &
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

  end subroutine sx_set_convect_rst_lx13

  subroutine sx_set_convect_rst_lx12(cr, cs, ct, cx, cy, cz, &
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

  end subroutine sx_set_convect_rst_lx12

  subroutine sx_set_convect_rst_lx11(cr, cs, ct, cx, cy, cz, &
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

  end subroutine sx_set_convect_rst_lx11

  subroutine sx_set_convect_rst_lx10(cr, cs, ct, cx, cy, cz, &
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

  end subroutine sx_set_convect_rst_lx10

  subroutine sx_set_convect_rst_lx9(cr, cs, ct, cx, cy, cz, &
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

  end subroutine sx_set_convect_rst_lx9

  subroutine sx_set_convect_rst_lx8(cr, cs, ct, cx, cy, cz, &
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

  end subroutine sx_set_convect_rst_lx8

  subroutine sx_set_convect_rst_lx7(cr, cs, ct, cx, cy, cz, &
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

  end subroutine sx_set_convect_rst_lx7

  subroutine sx_set_convect_rst_lx6(cr, cs, ct, cx, cy, cz, &
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

end subroutine sx_set_convect_rst_lx6

subroutine sx_set_convect_rst_lx5(cr, cs, ct, cx, cy, cz, &
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

  end subroutine sx_set_convect_rst_lx5

  subroutine sx_set_convect_rst_lx4(cr, cs, ct, cx, cy, cz, &
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

  end subroutine sx_set_convect_rst_lx4

  subroutine sx_set_convect_rst_lx3(cr, cs, ct, cx, cy, cz, &
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

  end subroutine sx_set_convect_rst_lx3

  subroutine sx_set_convect_rst_lx2(cr, cs, ct, cx, cy, cz, &
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

  end subroutine sx_set_convect_rst_lx2

end submodule sx_set_convect_rst
