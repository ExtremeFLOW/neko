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
!> Gradient kernels
submodule (opr_cpu) cpu_opgrad
  implicit none

contains

  module subroutine opr_cpu_opgrad(ux, uy, uz, u, coef, e_start, e_end)
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: e_start, e_end
    real(kind=rp), dimension(coef%Xh%lxyz, e_end - e_start + 1), &
         intent(inout) :: ux
    real(kind=rp), dimension(coef%Xh%lxyz, e_end - e_start + 1), &
         intent(inout) :: uy
    real(kind=rp), dimension(coef%Xh%lxyz, e_end - e_start + 1), &
         intent(inout) :: uz
    real(kind=rp), dimension(coef%Xh%lxyz, e_end - e_start + 1), &
         intent(in) :: u
    integer :: e_len

    e_len = e_end-e_start+1

    if (e_len .eq. 1) then
       call opr_cpu_opgrad_single(ux, uy, uz, u, coef, e_start)
    else
       call opr_cpu_opgrad_many(ux, uy, uz, u, coef, e_start, e_len)
    end if
  end subroutine opr_cpu_opgrad

  subroutine opr_cpu_opgrad_many(ux, uy, uz, u, coef, e_start, e_len)
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: e_start, e_len
    real(kind=rp), dimension(coef%Xh%lxyz, e_len), intent(inout) :: ux
    real(kind=rp), dimension(coef%Xh%lxyz, e_len), intent(inout) :: uy
    real(kind=rp), dimension(coef%Xh%lxyz, e_len), intent(inout) :: uz
    real(kind=rp), dimension(coef%Xh%lxyz, e_len), intent(in) :: u
    
    associate(Xh => coef%Xh, msh => coef%msh, &
         drdx => coef%drdx, drdy => coef%drdy, drdz => coef%drdz, &
         dsdx => coef%dsdx, dsdy => coef%dsdy, dsdz => coef%dsdz, &
         dtdx => coef%dtdx, dtdy => coef%dtdy, dtdz => coef%dtdz)

      select case (Xh%lx)
      case (18)
         call cpu_opgrad_lx18(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)
      case (17)
         call cpu_opgrad_lx17(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)
      case (16)
         call cpu_opgrad_lx16(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)
      case (15)
         call cpu_opgrad_lx15(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)
      case (14)
         call cpu_opgrad_lx14(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)
      case (13)
         call cpu_opgrad_lx13(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)
      case (12)
         call cpu_opgrad_lx12(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)
      case (11)
         call cpu_opgrad_lx11(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)
      case (10)
         call cpu_opgrad_lx10(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)

      case (9)
         call cpu_opgrad_lx9(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)
      case (8)
         call cpu_opgrad_lx8(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)
      case (7)
         call cpu_opgrad_lx7(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)
      case (6)
         call cpu_opgrad_lx6(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)
      case (5)
         call cpu_opgrad_lx5(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)
      case (4)
         call cpu_opgrad_lx4(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)
      case (3)
         call cpu_opgrad_lx3(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)
      case (2)
         call cpu_opgrad_lx2(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len)
      case default
         call cpu_opgrad_lx(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e_start), dsdx(1,1,1,e_start), dtdx(1,1,1,e_start), &
              drdy(1,1,1,e_start), dsdy(1,1,1,e_start), dtdy(1,1,1,e_start), &
              drdz(1,1,1,e_start), dsdz(1,1,1,e_start), dtdz(1,1,1,e_start), &
              Xh%w3, e_len, Xh%lx)
      end select
    end associate

  end subroutine opr_cpu_opgrad_many
  
  subroutine cpu_opgrad_lx(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, n, lx)
    integer, intent(in) :: n, lx
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             ur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                us(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             ut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                       + dsdx(i,1,1,e) * us(i,1,1) &
                       + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine cpu_opgrad_lx

  subroutine cpu_opgrad_lx18(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: e, i, j, k

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
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
                        + dx(i,14) * u(14,j,1,e) &
                        + dx(i,15) * u(15,j,1,e) &
                        + dx(i,16) * u(16,j,1,e) &
                        + dx(i,17) * u(17,j,1,e) &
                        + dx(i,18) * u(18,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
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
                          + dy(j,14) * u(i,14,k,e) &
                          + dy(j,15) * u(i,15,k,e) &
                          + dy(j,16) * u(i,16,k,e) &
                          + dy(j,17) * u(i,17,k,e) &
                          + dy(j,18) * u(i,18,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
                       + dz(k,14) * u(i,1,14,e) &
                       + dz(k,15) * u(i,1,15,e) &
                       + dz(k,16) * u(i,1,16,e) &
                       + dz(k,17) * u(i,1,17,e) &
                       + dz(k,18) * u(i,1,18,e)
          end do
       end do

       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                       + dsdx(i,1,1,e) * us(i,1,1) &
                       + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine cpu_opgrad_lx18

  subroutine cpu_opgrad_lx17(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: e, i, j, k

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
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
                        + dx(i,14) * u(14,j,1,e) &
                        + dx(i,15) * u(15,j,1,e) &
                        + dx(i,16) * u(16,j,1,e) &
                        + dx(i,17) * u(17,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
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
                          + dy(j,14) * u(i,14,k,e) &
                          + dy(j,15) * u(i,15,k,e) &
                          + dy(j,16) * u(i,16,k,e) &
                          + dy(j,17) * u(i,17,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
                       + dz(k,14) * u(i,1,14,e) &
                       + dz(k,15) * u(i,1,15,e) &
                       + dz(k,16) * u(i,1,16,e) &
                       + dz(k,17) * u(i,1,17,e)
          end do
       end do

       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                       + dsdx(i,1,1,e) * us(i,1,1) &
                       + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine cpu_opgrad_lx17

  subroutine cpu_opgrad_lx16(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: e, i, j, k

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
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
                        + dx(i,14) * u(14,j,1,e) &
                        + dx(i,15) * u(15,j,1,e) &
                        + dx(i,16) * u(16,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
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
                          + dy(j,14) * u(i,14,k,e) &
                          + dy(j,15) * u(i,15,k,e) &
                          + dy(j,16) * u(i,16,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
                       + dz(k,14) * u(i,1,14,e) &
                       + dz(k,15) * u(i,1,15,e) &
                       + dz(k,16) * u(i,1,16,e)
          end do
       end do

       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                       + dsdx(i,1,1,e) * us(i,1,1) &
                       + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine cpu_opgrad_lx16

  subroutine cpu_opgrad_lx15(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: e, i, j, k

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
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
                        + dx(i,14) * u(14,j,1,e) &
                        + dx(i,15) * u(15,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
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
                          + dy(j,14) * u(i,14,k,e) &
                          + dy(j,15) * u(i,15,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
                       + dz(k,14) * u(i,1,14,e) &
                       + dz(k,15) * u(i,1,15,e)
          end do
       end do

       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                       + dsdx(i,1,1,e) * us(i,1,1) &
                       + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine cpu_opgrad_lx15

  subroutine cpu_opgrad_lx14(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: e, i, j, k

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
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

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
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

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                       + dsdx(i,1,1,e) * us(i,1,1) &
                       + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine cpu_opgrad_lx14

  subroutine cpu_opgrad_lx13(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: e, i, j, k

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
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

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
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

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                       + dsdx(i,1,1,e) * us(i,1,1) &
                       + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine cpu_opgrad_lx13

  subroutine cpu_opgrad_lx12(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: e, i, j, k

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
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

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
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

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                       + dsdx(i,1,1,e) * us(i,1,1) &
                       + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine cpu_opgrad_lx12

  subroutine cpu_opgrad_lx11(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: e, i, j, k

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
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

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
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

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine cpu_opgrad_lx11

  subroutine cpu_opgrad_lx10(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: e, i, j, k

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
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

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
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

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                         + dsdx(i,1,1,e) * us(i,1,1) &
                         + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine cpu_opgrad_lx10

  subroutine cpu_opgrad_lx9(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: e, i, j, k

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
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

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
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

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine cpu_opgrad_lx9

  subroutine cpu_opgrad_lx8(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: e, i, j, k

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
                        + dx(i,2) * u(2,j,1,e) &
                        + dx(i,3) * u(3,j,1,e) &
                        + dx(i,4) * u(4,j,1,e) &
                        + dx(i,5) * u(5,j,1,e) &
                        + dx(i,6) * u(6,j,1,e) &
                        + dx(i,7) * u(7,j,1,e) &
                        + dx(i,8) * u(8,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
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

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
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
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine cpu_opgrad_lx8

  subroutine cpu_opgrad_lx7(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: e, i, j, k

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
                        + dx(i,2) * u(2,j,1,e) &
                        + dx(i,3) * u(3,j,1,e) &
                        + dx(i,4) * u(4,j,1,e) &
                        + dx(i,5) * u(5,j,1,e) &
                        + dx(i,6) * u(6,j,1,e) &
                        + dx(i,7) * u(7,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
                          + dy(j,2) * u(i,2,k,e) &
                          + dy(j,3) * u(i,3,k,e) &
                          + dy(j,4) * u(i,4,k,e) &
                          + dy(j,5) * u(i,5,k,e) &
                          + dy(j,6) * u(i,6,k,e) &
                          + dy(j,7) * u(i,7,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
                       + dz(k,2) * u(i,1,2,e) &
                       + dz(k,3) * u(i,1,3,e) &
                       + dz(k,4) * u(i,1,4,e) &
                       + dz(k,5) * u(i,1,5,e) &
                       + dz(k,6) * u(i,1,6,e) &
                       + dz(k,7) * u(i,1,7,e)
          end do
       end do

       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine cpu_opgrad_lx7

  subroutine cpu_opgrad_lx6(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: e, i, j, k

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
                        + dx(i,2) * u(2,j,1,e) &
                        + dx(i,3) * u(3,j,1,e) &
                        + dx(i,4) * u(4,j,1,e) &
                        + dx(i,5) * u(5,j,1,e) &
                        + dx(i,6) * u(6,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
                          + dy(j,2) * u(i,2,k,e) &
                          + dy(j,3) * u(i,3,k,e) &
                          + dy(j,4) * u(i,4,k,e) &
                          + dy(j,5) * u(i,5,k,e) &
                          + dy(j,6) * u(i,6,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
                       + dz(k,2) * u(i,1,2,e) &
                       + dz(k,3) * u(i,1,3,e) &
                       + dz(k,4) * u(i,1,4,e) &
                       + dz(k,5) * u(i,1,5,e) &
                       + dz(k,6) * u(i,1,6,e)
          end do
       end do

       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine cpu_opgrad_lx6

  subroutine cpu_opgrad_lx5(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: e, i, j, k

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
                        + dx(i,2) * u(2,j,1,e) &
                        + dx(i,3) * u(3,j,1,e) &
                        + dx(i,4) * u(4,j,1,e) &
                        + dx(i,5) * u(5,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
                          + dy(j,2) * u(i,2,k,e) &
                          + dy(j,3) * u(i,3,k,e) &
                          + dy(j,4) * u(i,4,k,e) &
                          + dy(j,5) * u(i,5,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
                       + dz(k,2) * u(i,1,2,e) &
                       + dz(k,3) * u(i,1,3,e) &
                       + dz(k,4) * u(i,1,4,e) &
                       + dz(k,5) * u(i,1,5,e)
          end do
       end do

       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine cpu_opgrad_lx5

  subroutine cpu_opgrad_lx4(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: e, i, j, k

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
                        + dx(i,2) * u(2,j,1,e) &
                        + dx(i,3) * u(3,j,1,e) &
                        + dx(i,4) * u(4,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
                          + dy(j,2) * u(i,2,k,e) &
                          + dy(j,3) * u(i,3,k,e) &
                          + dy(j,4) * u(i,4,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
                       + dz(k,2) * u(i,1,2,e) &
                       + dz(k,3) * u(i,1,3,e) &
                       + dz(k,4) * u(i,1,4,e)
          end do
       end do

       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine cpu_opgrad_lx4

  subroutine cpu_opgrad_lx3(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: e, i, j, k

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
                        + dx(i,2) * u(2,j,1,e) &
                        + dx(i,3) * u(3,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
                          + dy(j,2) * u(i,2,k,e) &
                          + dy(j,3) * u(i,3,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
                       + dz(k,2) * u(i,1,2,e) &
                       + dz(k,3) * u(i,1,3,e)
          end do
       end do

       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine cpu_opgrad_lx3

  subroutine cpu_opgrad_lx2(ux, uy, uz, u, dx, dy, dz, &
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
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: e, i, j, k

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             ur(i,j,1) = dx(i,1) * u(1,j,1,e) &
                        + dx(i,2) * u(2,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                us(i,j,k) = dy(j,1) * u(i,1,k,e) &
                          + dy(j,2) * u(i,2,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             ut(i,1,k) = dz(k,1) * u(i,1,1,e) &
                       + dz(k,2) * u(i,1,2,e)
          end do
       end do

       do i = 1, lx * lx * lx
          ux(i,1,1,e) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          uy(i,1,1,e) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          uz(i,1,1,e) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )
       end do
    end do
  end subroutine cpu_opgrad_lx2

  subroutine opr_cpu_opgrad_single(ux, uy, uz, u, coef, e)
    integer, parameter :: e_len = 1
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: e
    real(kind=rp), dimension(coef%Xh%lxyz, e_len), intent(inout) :: ux
    real(kind=rp), dimension(coef%Xh%lxyz, e_len), intent(inout) :: uy
    real(kind=rp), dimension(coef%Xh%lxyz, e_len), intent(inout) :: uz
    real(kind=rp), dimension(coef%Xh%lxyz, e_len), intent(in) :: u
    
    associate(Xh => coef%Xh, msh => coef%msh, &
         drdx => coef%drdx, drdy => coef%drdy, drdz => coef%drdz, &
         dsdx => coef%dsdx, dsdy => coef%dsdy, dsdz => coef%dsdz, &
         dtdx => coef%dtdx, dtdy => coef%dtdy, dtdz => coef%dtdz)

      select case (Xh%lx)
      case (18)
         call cpu_opgrad_lx18_single(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e), &
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e), &
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e), &
              Xh%w3)
      case (17)
         call cpu_opgrad_lx17_single(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e), &
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e), &
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e), &
              Xh%w3)
      case (16)
         call cpu_opgrad_lx16_single(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e), &
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e), &
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e), &
              Xh%w3)
      case (15)
         call cpu_opgrad_lx15_single(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e), &
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e), &
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e), &
              Xh%w3)
      case (14)
         call cpu_opgrad_lx14_single(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e), &
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e), &
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e), &
              Xh%w3)
      case (13)
         call cpu_opgrad_lx13_single(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e), &
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e), &
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e), &
              Xh%w3)
      case (12)
         call cpu_opgrad_lx12_single(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e), &
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e), &
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e), &
              Xh%w3)
      case (11)
         call cpu_opgrad_lx11_single(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e), &
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e), &
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e), &
              Xh%w3)
      case (10)
         call cpu_opgrad_lx10_single(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e), &
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e), &
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e), &
              Xh%w3)

      case (9)
         call cpu_opgrad_lx9_single(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e), &
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e), &
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e), &
              Xh%w3)
      case (8)
         call cpu_opgrad_lx8_single(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e), &
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e), &
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e), &
              Xh%w3)
      case (7)
         call cpu_opgrad_lx7_single(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e), &
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e), &
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e), &
              Xh%w3)
      case (6)
         call cpu_opgrad_lx6_single(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e), &
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e), &
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e), &
              Xh%w3)
      case (5)
         call cpu_opgrad_lx5_single(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e), &
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e), &
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e), &
              Xh%w3)
      case (4)
         call cpu_opgrad_lx4_single(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e), &
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e), &
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e), &
              Xh%w3)
      case (3)
         call cpu_opgrad_lx3_single(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e), &
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e), &
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e), &
              Xh%w3)
      case (2)
         call cpu_opgrad_lx2_single(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e), &
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e), &
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e), &
              Xh%w3)
      case default
         call cpu_opgrad_lx_single(ux, uy, uz, u, &
              Xh%dx, Xh%dy, Xh%dz, &
              drdx(1,1,1,e), dsdx(1,1,1,e), dtdx(1,1,1,e), &
              drdy(1,1,1,e), dsdy(1,1,1,e), dtdy(1,1,1,e), &
              drdz(1,1,1,e), dsdz(1,1,1,e), dtdz(1,1,1,e), &
              Xh%w3, Xh%lx)
      end select
    end associate

  end subroutine opr_cpu_opgrad_single
  
  subroutine cpu_opgrad_lx_single(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, lx)
    integer, intent(in) :: lx
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: tmp
    integer :: i, j, k, l

    do j = 1, lx * lx
       do i = 1, lx
          tmp = 0.0_rp
          do k = 1, lx
             tmp = tmp + dx(i,k) * u(k,j,1)
          end do
          ur(i,j,1) = tmp
       end do
    end do
    
    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dy(j,l) * u(i,l,k)
             end do
             us(i,j,k) = tmp
          end do
       end do
    end do
    
    do k = 1, lx
       do i = 1, lx*lx
          tmp = 0.0_rp
          do l = 1, lx
             tmp = tmp + dz(k,l) * u(i,1,l)
          end do
          ut(i,1,k) = tmp
       end do
    end do
    
    do i = 1, lx * lx * lx
       ux(i,1,1) = w3(i,1,1) &
                   * ( drdx(i,1,1) * ur(i,1,1) &
                     + dsdx(i,1,1) * us(i,1,1) &
                     + dtdx(i,1,1) * ut(i,1,1) )
       uy(i,1,1) = w3(i,1,1) &
                   * ( dsdy(i,1,1) * us(i,1,1) &
                     + drdy(i,1,1) * ur(i,1,1) &
                     + dtdy(i,1,1) * ut(i,1,1) )
       uz(i,1,1) = w3(i,1,1) &
                   * ( dtdz(i,1,1) * ut(i,1,1) &
                     + drdz(i,1,1) * ur(i,1,1) &
                     + dsdz(i,1,1) * us(i,1,1) )
    end do

  end subroutine cpu_opgrad_lx_single

  subroutine cpu_opgrad_lx18_single(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3)
    integer, parameter :: lx = 18
      real(kind=rp), dimension(lx, lx, lx), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          ur(i,j,1) = dx(i,1) * u(1,j,1) &
                    + dx(i,2) * u(2,j,1) &
                    + dx(i,3) * u(3,j,1) &
                    + dx(i,4) * u(4,j,1) &
                    + dx(i,5) * u(5,j,1) &
                    + dx(i,6) * u(6,j,1) &
                    + dx(i,7) * u(7,j,1) &
                    + dx(i,8) * u(8,j,1) &
                    + dx(i,9) * u(9,j,1) &
                    + dx(i,10) * u(10,j,1) &
                    + dx(i,11) * u(11,j,1) &
                    + dx(i,12) * u(12,j,1) &
                    + dx(i,13) * u(13,j,1) &
                    + dx(i,14) * u(14,j,1) &
                    + dx(i,15) * u(15,j,1) &
                    + dx(i,16) * u(16,j,1) &
                    + dx(i,17) * u(17,j,1) &
                    + dx(i,18) * u(18,j,1)
       end do
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             us(i,j,k) = dy(j,1) * u(i,1,k) &
                       + dy(j,2) * u(i,2,k) &
                       + dy(j,3) * u(i,3,k) &
                       + dy(j,4) * u(i,4,k) &
                       + dy(j,5) * u(i,5,k) &
                       + dy(j,6) * u(i,6,k) &
                       + dy(j,7) * u(i,7,k) &
                       + dy(j,8) * u(i,8,k) &
                       + dy(j,9) * u(i,9,k) &
                       + dy(j,10) * u(i,10,k) &
                       + dy(j,11) * u(i,11,k) &
                       + dy(j,12) * u(i,12,k) &
                       + dy(j,13) * u(i,13,k) &
                       + dy(j,14) * u(i,14,k) &
                       + dy(j,15) * u(i,15,k) &
                       + dy(j,16) * u(i,16,k) &
                       + dy(j,17) * u(i,17,k) &
                       + dy(j,18) * u(i,18,k)
          end do
       end do
    end do

    do k = 1, lx
       do i = 1, lx*lx
          ut(i,1,k) = dz(k,1) * u(i,1,1) &
                    + dz(k,2) * u(i,1,2) &
                    + dz(k,3) * u(i,1,3) &
                    + dz(k,4) * u(i,1,4) &
                    + dz(k,5) * u(i,1,5) &
                    + dz(k,6) * u(i,1,6) &
                    + dz(k,7) * u(i,1,7) &
                    + dz(k,8) * u(i,1,8) &
                    + dz(k,9) * u(i,1,9) &
                    + dz(k,10) * u(i,1,10) &
                    + dz(k,11) * u(i,1,11) &
                    + dz(k,12) * u(i,1,12) &
                    + dz(k,13) * u(i,1,13) &
                    + dz(k,14) * u(i,1,14) &
                    + dz(k,15) * u(i,1,15) &
                    + dz(k,16) * u(i,1,16) &
                    + dz(k,17) * u(i,1,17) &
                    + dz(k,18) * u(i,1,18)
       end do
    end do
    
    do i = 1, lx * lx * lx
       ux(i,1,1) = w3(i,1,1) &
                   * ( drdx(i,1,1) * ur(i,1,1) &
                     + dsdx(i,1,1) * us(i,1,1) &
                     + dtdx(i,1,1) * ut(i,1,1) )
       uy(i,1,1) = w3(i,1,1) &
                   * ( dsdy(i,1,1) * us(i,1,1) &
                     + drdy(i,1,1) * ur(i,1,1) &
                     + dtdy(i,1,1) * ut(i,1,1) )
       uz(i,1,1) = w3(i,1,1) &
                   * ( dtdz(i,1,1) * ut(i,1,1) &
                     + drdz(i,1,1) * ur(i,1,1) &
                     + dsdz(i,1,1) * us(i,1,1) )
    end do

  end subroutine cpu_opgrad_lx18_single

  subroutine cpu_opgrad_lx17_single(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3)
    integer, parameter :: lx = 17
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          ur(i,j,1) = dx(i,1) * u(1,j,1) &
                    + dx(i,2) * u(2,j,1) &
                    + dx(i,3) * u(3,j,1) &
                    + dx(i,4) * u(4,j,1) &
                    + dx(i,5) * u(5,j,1) &
                    + dx(i,6) * u(6,j,1) &
                    + dx(i,7) * u(7,j,1) &
                    + dx(i,8) * u(8,j,1) &
                    + dx(i,9) * u(9,j,1) &
                    + dx(i,10) * u(10,j,1) &
                    + dx(i,11) * u(11,j,1) &
                    + dx(i,12) * u(12,j,1) &
                    + dx(i,13) * u(13,j,1) &
                    + dx(i,14) * u(14,j,1) &
                    + dx(i,15) * u(15,j,1) &
                    + dx(i,16) * u(16,j,1) &
                    + dx(i,17) * u(17,j,1)
       end do
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             us(i,j,k) = dy(j,1) * u(i,1,k) &
                       + dy(j,2) * u(i,2,k) &
                       + dy(j,3) * u(i,3,k) &
                       + dy(j,4) * u(i,4,k) &
                       + dy(j,5) * u(i,5,k) &
                       + dy(j,6) * u(i,6,k) &
                       + dy(j,7) * u(i,7,k) &
                       + dy(j,8) * u(i,8,k) &
                       + dy(j,9) * u(i,9,k) &
                       + dy(j,10) * u(i,10,k) &
                       + dy(j,11) * u(i,11,k) &
                       + dy(j,12) * u(i,12,k) &
                       + dy(j,13) * u(i,13,k) &
                       + dy(j,14) * u(i,14,k) &
                       + dy(j,15) * u(i,15,k) &
                       + dy(j,16) * u(i,16,k) &
                       + dy(j,17) * u(i,17,k)
          end do
       end do
    end do
    
    do k = 1, lx
       do i = 1, lx*lx
          ut(i,1,k) = dz(k,1) * u(i,1,1) &
                    + dz(k,2) * u(i,1,2) &
                    + dz(k,3) * u(i,1,3) &
                    + dz(k,4) * u(i,1,4) &
                    + dz(k,5) * u(i,1,5) &
                    + dz(k,6) * u(i,1,6) &
                    + dz(k,7) * u(i,1,7) &
                    + dz(k,8) * u(i,1,8) &
                    + dz(k,9) * u(i,1,9) &
                    + dz(k,10) * u(i,1,10) &
                    + dz(k,11) * u(i,1,11) &
                    + dz(k,12) * u(i,1,12) &
                    + dz(k,13) * u(i,1,13) &
                    + dz(k,14) * u(i,1,14) &
                    + dz(k,15) * u(i,1,15) &
                    + dz(k,16) * u(i,1,16) &
                    + dz(k,17) * u(i,1,17)
       end do
    end do
    
    do i = 1, lx * lx * lx
       ux(i,1,1) = w3(i,1,1) &
                   * ( drdx(i,1,1) * ur(i,1,1) &
                     + dsdx(i,1,1) * us(i,1,1) &
                     + dtdx(i,1,1) * ut(i,1,1) )
       uy(i,1,1) = w3(i,1,1) &
                   * ( dsdy(i,1,1) * us(i,1,1) &
                     + drdy(i,1,1) * ur(i,1,1) &
                     + dtdy(i,1,1) * ut(i,1,1) )
       uz(i,1,1) = w3(i,1,1) &
                   * ( dtdz(i,1,1) * ut(i,1,1) &
                     + drdz(i,1,1) * ur(i,1,1) &
                     + dsdz(i,1,1) * us(i,1,1) )
    end do
  end subroutine cpu_opgrad_lx17_single

  subroutine cpu_opgrad_lx16_single(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3)
    integer, parameter :: lx = 16
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer ::  i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          ur(i,j,1) = dx(i,1) * u(1,j,1) &
                    + dx(i,2) * u(2,j,1) &
                    + dx(i,3) * u(3,j,1) &
                    + dx(i,4) * u(4,j,1) &
                    + dx(i,5) * u(5,j,1) &
                    + dx(i,6) * u(6,j,1) &
                    + dx(i,7) * u(7,j,1) &
                    + dx(i,8) * u(8,j,1) &
                    + dx(i,9) * u(9,j,1) &
                    + dx(i,10) * u(10,j,1) &
                    + dx(i,11) * u(11,j,1) &
                    + dx(i,12) * u(12,j,1) &
                    + dx(i,13) * u(13,j,1) &
                    + dx(i,14) * u(14,j,1) &
                    + dx(i,15) * u(15,j,1) &
                    + dx(i,16) * u(16,j,1)
       end do
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             us(i,j,k) = dy(j,1) * u(i,1,k) &
                       + dy(j,2) * u(i,2,k) &
                       + dy(j,3) * u(i,3,k) &
                       + dy(j,4) * u(i,4,k) &
                       + dy(j,5) * u(i,5,k) &
                       + dy(j,6) * u(i,6,k) &
                       + dy(j,7) * u(i,7,k) &
                       + dy(j,8) * u(i,8,k) &
                       + dy(j,9) * u(i,9,k) &
                       + dy(j,10) * u(i,10,k) &
                       + dy(j,11) * u(i,11,k) &
                       + dy(j,12) * u(i,12,k) &
                       + dy(j,13) * u(i,13,k) &
                       + dy(j,14) * u(i,14,k) &
                       + dy(j,15) * u(i,15,k) &
                       + dy(j,16) * u(i,16,k)
          end do
       end do
    end do

    do k = 1, lx
       do i = 1, lx*lx
          ut(i,1,k) = dz(k,1) * u(i,1,1) &
                    + dz(k,2) * u(i,1,2) &
                    + dz(k,3) * u(i,1,3) &
                    + dz(k,4) * u(i,1,4) &
                    + dz(k,5) * u(i,1,5) &
                    + dz(k,6) * u(i,1,6) &
                    + dz(k,7) * u(i,1,7) &
                    + dz(k,8) * u(i,1,8) &
                    + dz(k,9) * u(i,1,9) &
                    + dz(k,10) * u(i,1,10) &
                    + dz(k,11) * u(i,1,11) &
                    + dz(k,12) * u(i,1,12) &
                    + dz(k,13) * u(i,1,13) &
                    + dz(k,14) * u(i,1,14) &
                    + dz(k,15) * u(i,1,15) &
                    + dz(k,16) * u(i,1,16)
       end do
    end do
    
    do i = 1, lx * lx * lx
       ux(i,1,1) = w3(i,1,1) &
                   * ( drdx(i,1,1) * ur(i,1,1) &
                     + dsdx(i,1,1) * us(i,1,1) &
                     + dtdx(i,1,1) * ut(i,1,1) )
       uy(i,1,1) = w3(i,1,1) &
                   * ( dsdy(i,1,1) * us(i,1,1) &
                     + drdy(i,1,1) * ur(i,1,1) &
                     + dtdy(i,1,1) * ut(i,1,1) )
       uz(i,1,1) = w3(i,1,1) &
                   * ( dtdz(i,1,1) * ut(i,1,1) &
                     + drdz(i,1,1) * ur(i,1,1) &
                     + dsdz(i,1,1) * us(i,1,1) )
    end do
  end subroutine cpu_opgrad_lx16_single

  subroutine cpu_opgrad_lx15_single(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3)
    integer, parameter :: lx = 15
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          ur(i,j,1) = dx(i,1) * u(1,j,1) &
                    + dx(i,2) * u(2,j,1) &
                    + dx(i,3) * u(3,j,1) &
                    + dx(i,4) * u(4,j,1) &
                    + dx(i,5) * u(5,j,1) &
                    + dx(i,6) * u(6,j,1) &
                    + dx(i,7) * u(7,j,1) &
                    + dx(i,8) * u(8,j,1) &
                    + dx(i,9) * u(9,j,1) &
                    + dx(i,10) * u(10,j,1) &
                    + dx(i,11) * u(11,j,1) &
                    + dx(i,12) * u(12,j,1) &
                    + dx(i,13) * u(13,j,1) &
                    + dx(i,14) * u(14,j,1) &
                    + dx(i,15) * u(15,j,1)
       end do
    end do
    
    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             us(i,j,k) = dy(j,1) * u(i,1,k) &
                       + dy(j,2) * u(i,2,k) &
                       + dy(j,3) * u(i,3,k) &
                       + dy(j,4) * u(i,4,k) &
                       + dy(j,5) * u(i,5,k) &
                       + dy(j,6) * u(i,6,k) &
                       + dy(j,7) * u(i,7,k) &
                       + dy(j,8) * u(i,8,k) &
                       + dy(j,9) * u(i,9,k) &
                       + dy(j,10) * u(i,10,k) &
                       + dy(j,11) * u(i,11,k) &
                       + dy(j,12) * u(i,12,k) &
                       + dy(j,13) * u(i,13,k) &
                       + dy(j,14) * u(i,14,k) &
                       + dy(j,15) * u(i,15,k)
          end do
       end do
    end do
    
    do k = 1, lx
       do i = 1, lx*lx
          ut(i,1,k) = dz(k,1) * u(i,1,1) &
                    + dz(k,2) * u(i,1,2) &
                    + dz(k,3) * u(i,1,3) &
                    + dz(k,4) * u(i,1,4) &
                    + dz(k,5) * u(i,1,5) &
                    + dz(k,6) * u(i,1,6) &
                    + dz(k,7) * u(i,1,7) &
                    + dz(k,8) * u(i,1,8) &
                    + dz(k,9) * u(i,1,9) &
                    + dz(k,10) * u(i,1,10) &
                    + dz(k,11) * u(i,1,11) &
                    + dz(k,12) * u(i,1,12) &
                    + dz(k,13) * u(i,1,13) &
                    + dz(k,14) * u(i,1,14) &
                    + dz(k,15) * u(i,1,15)
       end do
    end do
    
    do i = 1, lx * lx * lx
       ux(i,1,1) = w3(i,1,1) &
                   * ( drdx(i,1,1) * ur(i,1,1) &
                     + dsdx(i,1,1) * us(i,1,1) &
                     + dtdx(i,1,1) * ut(i,1,1) )
       uy(i,1,1) = w3(i,1,1) &
                   * ( dsdy(i,1,1) * us(i,1,1) &
                     + drdy(i,1,1) * ur(i,1,1) &
                     + dtdy(i,1,1) * ut(i,1,1) )
       uz(i,1,1) = w3(i,1,1) &
                   * ( dtdz(i,1,1) * ut(i,1,1) &
                     + drdz(i,1,1) * ur(i,1,1) &
                     + dsdz(i,1,1) * us(i,1,1) )
    end do
  end subroutine cpu_opgrad_lx15_single

  subroutine cpu_opgrad_lx14_single(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3)
    integer, parameter :: lx = 14
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          ur(i,j,1) = dx(i,1) * u(1,j,1) &
                    + dx(i,2) * u(2,j,1) &
                    + dx(i,3) * u(3,j,1) &
                    + dx(i,4) * u(4,j,1) &
                    + dx(i,5) * u(5,j,1) &
                    + dx(i,6) * u(6,j,1) &
                    + dx(i,7) * u(7,j,1) &
                    + dx(i,8) * u(8,j,1) &
                    + dx(i,9) * u(9,j,1) &
                    + dx(i,10) * u(10,j,1) &
                    + dx(i,11) * u(11,j,1) &
                    + dx(i,12) * u(12,j,1) &
                    + dx(i,13) * u(13,j,1) &
                    + dx(i,14) * u(14,j,1)
       end do
    end do
    
    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             us(i,j,k) = dy(j,1) * u(i,1,k) &
                       + dy(j,2) * u(i,2,k) &
                       + dy(j,3) * u(i,3,k) &
                       + dy(j,4) * u(i,4,k) &
                       + dy(j,5) * u(i,5,k) &
                       + dy(j,6) * u(i,6,k) &
                       + dy(j,7) * u(i,7,k) &
                       + dy(j,8) * u(i,8,k) &
                       + dy(j,9) * u(i,9,k) &
                       + dy(j,10) * u(i,10,k) &
                       + dy(j,11) * u(i,11,k) &
                       + dy(j,12) * u(i,12,k) &
                       + dy(j,13) * u(i,13,k) &
                       + dy(j,14) * u(i,14,k)
          end do
       end do
    end do
    
    do k = 1, lx
       do i = 1, lx*lx
          ut(i,1,k) = dz(k,1) * u(i,1,1) &
                    + dz(k,2) * u(i,1,2) &
                    + dz(k,3) * u(i,1,3) &
                    + dz(k,4) * u(i,1,4) &
                    + dz(k,5) * u(i,1,5) &
                    + dz(k,6) * u(i,1,6) &
                    + dz(k,7) * u(i,1,7) &
                    + dz(k,8) * u(i,1,8) &
                    + dz(k,9) * u(i,1,9) &
                    + dz(k,10) * u(i,1,10) &
                    + dz(k,11) * u(i,1,11) &
                    + dz(k,12) * u(i,1,12) &
                    + dz(k,13) * u(i,1,13) &
                    + dz(k,14) * u(i,1,14)
       end do
    end do

    do i = 1, lx * lx * lx
       ux(i,1,1) = w3(i,1,1) &
                   * ( drdx(i,1,1) * ur(i,1,1) &
                     + dsdx(i,1,1) * us(i,1,1) &
                     + dtdx(i,1,1) * ut(i,1,1) )
       uy(i,1,1) = w3(i,1,1) &
                   * ( dsdy(i,1,1) * us(i,1,1) &
                     + drdy(i,1,1) * ur(i,1,1) &
                     + dtdy(i,1,1) * ut(i,1,1) )
       uz(i,1,1) = w3(i,1,1) &
                   * ( dtdz(i,1,1) * ut(i,1,1) &
                     + drdz(i,1,1) * ur(i,1,1) &
                     + dsdz(i,1,1) * us(i,1,1) )
    end do
  end subroutine cpu_opgrad_lx14_single

  subroutine cpu_opgrad_lx13_single(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3)
    integer, parameter :: lx = 13
      real(kind=rp), dimension(lx, lx, lx), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          ur(i,j,1) = dx(i,1) * u(1,j,1) &
                    + dx(i,2) * u(2,j,1) &
                    + dx(i,3) * u(3,j,1) &
                    + dx(i,4) * u(4,j,1) &
                    + dx(i,5) * u(5,j,1) &
                    + dx(i,6) * u(6,j,1) &
                    + dx(i,7) * u(7,j,1) &
                    + dx(i,8) * u(8,j,1) &
                    + dx(i,9) * u(9,j,1) &
                    + dx(i,10) * u(10,j,1) &
                    + dx(i,11) * u(11,j,1) &
                    + dx(i,12) * u(12,j,1) &
                    + dx(i,13) * u(13,j,1)
       end do
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             us(i,j,k) = dy(j,1) * u(i,1,k) &
                       + dy(j,2) * u(i,2,k) &
                       + dy(j,3) * u(i,3,k) &
                       + dy(j,4) * u(i,4,k) &
                       + dy(j,5) * u(i,5,k) &
                       + dy(j,6) * u(i,6,k) &
                       + dy(j,7) * u(i,7,k) &
                       + dy(j,8) * u(i,8,k) &
                       + dy(j,9) * u(i,9,k) &
                       + dy(j,10) * u(i,10,k) &
                       + dy(j,11) * u(i,11,k) &
                       + dy(j,12) * u(i,12,k) &
                       + dy(j,13) * u(i,13,k)
          end do
       end do
    end do

    do k = 1, lx
       do i = 1, lx*lx
          ut(i,1,k) = dz(k,1) * u(i,1,1) &
                    + dz(k,2) * u(i,1,2) &
                    + dz(k,3) * u(i,1,3) &
                    + dz(k,4) * u(i,1,4) &
                    + dz(k,5) * u(i,1,5) &
                    + dz(k,6) * u(i,1,6) &
                    + dz(k,7) * u(i,1,7) &
                    + dz(k,8) * u(i,1,8) &
                    + dz(k,9) * u(i,1,9) &
                    + dz(k,10) * u(i,1,10) &
                    + dz(k,11) * u(i,1,11) &
                    + dz(k,12) * u(i,1,12) &
                    + dz(k,13) * u(i,1,13)
       end do
    end do

    do i = 1, lx * lx * lx
       ux(i,1,1) = w3(i,1,1) &
                   * ( drdx(i,1,1) * ur(i,1,1) &
                     + dsdx(i,1,1) * us(i,1,1) &
                     + dtdx(i,1,1) * ut(i,1,1) )
       uy(i,1,1) = w3(i,1,1) &
                   * ( dsdy(i,1,1) * us(i,1,1) &
                     + drdy(i,1,1) * ur(i,1,1) &
                     + dtdy(i,1,1) * ut(i,1,1) )
       uz(i,1,1) = w3(i,1,1) &
                   * ( dtdz(i,1,1) * ut(i,1,1) &
                     + drdz(i,1,1) * ur(i,1,1) &
                     + dsdz(i,1,1) * us(i,1,1) )
    end do
  end subroutine cpu_opgrad_lx13_single

  subroutine cpu_opgrad_lx12_single(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3)
    integer, parameter :: lx = 12
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          ur(i,j,1) = dx(i,1) * u(1,j,1) &
                    + dx(i,2) * u(2,j,1) &
                    + dx(i,3) * u(3,j,1) &
                    + dx(i,4) * u(4,j,1) &
                    + dx(i,5) * u(5,j,1) &
                    + dx(i,6) * u(6,j,1) &
                    + dx(i,7) * u(7,j,1) &
                    + dx(i,8) * u(8,j,1) &
                    + dx(i,9) * u(9,j,1) &
                    + dx(i,10) * u(10,j,1) &
                    + dx(i,11) * u(11,j,1) &
                    + dx(i,12) * u(12,j,1)
       end do
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             us(i,j,k) = dy(j,1) * u(i,1,k) &
                       + dy(j,2) * u(i,2,k) &
                       + dy(j,3) * u(i,3,k) &
                       + dy(j,4) * u(i,4,k) &
                       + dy(j,5) * u(i,5,k) &
                       + dy(j,6) * u(i,6,k) &
                       + dy(j,7) * u(i,7,k) &
                       + dy(j,8) * u(i,8,k) &
                       + dy(j,9) * u(i,9,k) &
                       + dy(j,10) * u(i,10,k) &
                       + dy(j,11) * u(i,11,k) &
                       + dy(j,12) * u(i,12,k)
          end do
       end do
    end do

    do k = 1, lx
       do i = 1, lx*lx
          ut(i,1,k) = dz(k,1) * u(i,1,1) &
                    + dz(k,2) * u(i,1,2) &
                    + dz(k,3) * u(i,1,3) &
                    + dz(k,4) * u(i,1,4) &
                    + dz(k,5) * u(i,1,5) &
                    + dz(k,6) * u(i,1,6) &
                    + dz(k,7) * u(i,1,7) &
                    + dz(k,8) * u(i,1,8) &
                    + dz(k,9) * u(i,1,9) &
                    + dz(k,10) * u(i,1,10) &
                    + dz(k,11) * u(i,1,11) &
                    + dz(k,12) * u(i,1,12)
       end do
    end do
       
    do i = 1, lx * lx * lx
       ux(i,1,1) = w3(i,1,1) &
                   * ( drdx(i,1,1) * ur(i,1,1) &
                     + dsdx(i,1,1) * us(i,1,1) &
                     + dtdx(i,1,1) * ut(i,1,1) )
       uy(i,1,1) = w3(i,1,1) &
                   * ( dsdy(i,1,1) * us(i,1,1) &
                     + drdy(i,1,1) * ur(i,1,1) &
                     + dtdy(i,1,1) * ut(i,1,1) )
       uz(i,1,1) = w3(i,1,1) &
                   * ( dtdz(i,1,1) * ut(i,1,1) &
                     + drdz(i,1,1) * ur(i,1,1) &
                     + dsdz(i,1,1) * us(i,1,1) )
    end do
  end subroutine cpu_opgrad_lx12_single

  subroutine cpu_opgrad_lx11_single(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3)
    integer, parameter :: lx = 11
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          ur(i,j,1) = dx(i,1) * u(1,j,1) &
                    + dx(i,2) * u(2,j,1) &
                    + dx(i,3) * u(3,j,1) &
                    + dx(i,4) * u(4,j,1) &
                    + dx(i,5) * u(5,j,1) &
                    + dx(i,6) * u(6,j,1) &
                    + dx(i,7) * u(7,j,1) &
                    + dx(i,8) * u(8,j,1) &
                    + dx(i,9) * u(9,j,1) &
                    + dx(i,10) * u(10,j,1) &
                    + dx(i,11) * u(11,j,1)
       end do
    end do
    
    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             us(i,j,k) = dy(j,1) * u(i,1,k) &
                       + dy(j,2) * u(i,2,k) &
                       + dy(j,3) * u(i,3,k) &
                       + dy(j,4) * u(i,4,k) &
                       + dy(j,5) * u(i,5,k) &
                       + dy(j,6) * u(i,6,k) &
                       + dy(j,7) * u(i,7,k) &
                       + dy(j,8) * u(i,8,k) &
                       + dy(j,9) * u(i,9,k) &
                       + dy(j,10) * u(i,10,k) &
                       + dy(j,11) * u(i,11,k)
          end do
       end do
    end do
    
    do k = 1, lx
       do i = 1, lx*lx
          ut(i,1,k) = dz(k,1) * u(i,1,1) &
                    + dz(k,2) * u(i,1,2) &
                    + dz(k,3) * u(i,1,3) &
                    + dz(k,4) * u(i,1,4) &
                    + dz(k,5) * u(i,1,5) &
                    + dz(k,6) * u(i,1,6) &
                    + dz(k,7) * u(i,1,7) &
                    + dz(k,8) * u(i,1,8) &
                    + dz(k,9) * u(i,1,9) &
                    + dz(k,10) * u(i,1,10) &
                    + dz(k,11) * u(i,1,11)
       end do
    end do
    
    do i = 1, lx * lx * lx
       ux(i,1,1) = w3(i,1,1) &
                   * ( drdx(i,1,1) * ur(i,1,1) &
                     + dsdx(i,1,1) * us(i,1,1) &
                     + dtdx(i,1,1) * ut(i,1,1) )
       uy(i,1,1) = w3(i,1,1) &
                   * ( dsdy(i,1,1) * us(i,1,1) &
                     + drdy(i,1,1) * ur(i,1,1) &
                     + dtdy(i,1,1) * ut(i,1,1) )
       uz(i,1,1) = w3(i,1,1) &
                   * ( dtdz(i,1,1) * ut(i,1,1) &
                     + drdz(i,1,1) * ur(i,1,1) &
                     + dsdz(i,1,1) * us(i,1,1) )
    end do
  end subroutine cpu_opgrad_lx11_single

  subroutine cpu_opgrad_lx10_single(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3)
    integer, parameter :: lx = 10
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          ur(i,j,1) = dx(i,1) * u(1,j,1) &
                    + dx(i,2) * u(2,j,1) &
                    + dx(i,3) * u(3,j,1) &
                    + dx(i,4) * u(4,j,1) &
                    + dx(i,5) * u(5,j,1) &
                    + dx(i,6) * u(6,j,1) &
                    + dx(i,7) * u(7,j,1) &
                    + dx(i,8) * u(8,j,1) &
                    + dx(i,9) * u(9,j,1) &
                    + dx(i,10) * u(10,j,1)
       end do
    end do
    
    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             us(i,j,k) = dy(j,1) * u(i,1,k) &
                       + dy(j,2) * u(i,2,k) &
                       + dy(j,3) * u(i,3,k) &
                       + dy(j,4) * u(i,4,k) &
                       + dy(j,5) * u(i,5,k) &
                       + dy(j,6) * u(i,6,k) &
                       + dy(j,7) * u(i,7,k) &
                       + dy(j,8) * u(i,8,k) &
                       + dy(j,9) * u(i,9,k) &
                       + dy(j,10) * u(i,10,k)
          end do
       end do
    end do
       
    do k = 1, lx
       do i = 1, lx*lx
          ut(i,1,k) = dz(k,1) * u(i,1,1) &
                    + dz(k,2) * u(i,1,2) &
                    + dz(k,3) * u(i,1,3) &
                    + dz(k,4) * u(i,1,4) &
                    + dz(k,5) * u(i,1,5) &
                    + dz(k,6) * u(i,1,6) &
                    + dz(k,7) * u(i,1,7) &
                    + dz(k,8) * u(i,1,8) &
                    + dz(k,9) * u(i,1,9) &
                    + dz(k,10) * u(i,1,10)
       end do
    end do
    
    do i = 1, lx * lx * lx
       ux(i,1,1) = w3(i,1,1) &
                   * ( drdx(i,1,1) * ur(i,1,1) &
                     + dsdx(i,1,1) * us(i,1,1) &
                     + dtdx(i,1,1) * ut(i,1,1) )
       uy(i,1,1) = w3(i,1,1) &
                   * ( dsdy(i,1,1) * us(i,1,1) &
                     + drdy(i,1,1) * ur(i,1,1) &
                     + dtdy(i,1,1) * ut(i,1,1) )
       uz(i,1,1) = w3(i,1,1) &
                   * ( dtdz(i,1,1) * ut(i,1,1) &
                     + drdz(i,1,1) * ur(i,1,1) &
                     + dsdz(i,1,1) * us(i,1,1) )
    end do
  end subroutine cpu_opgrad_lx10_single

  subroutine cpu_opgrad_lx9_single(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3)
    integer, parameter :: lx = 9
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          ur(i,j,1) = dx(i,1) * u(1,j,1) &
                    + dx(i,2) * u(2,j,1) &
                    + dx(i,3) * u(3,j,1) &
                    + dx(i,4) * u(4,j,1) &
                    + dx(i,5) * u(5,j,1) &
                    + dx(i,6) * u(6,j,1) &
                    + dx(i,7) * u(7,j,1) &
                    + dx(i,8) * u(8,j,1) &
                    + dx(i,9) * u(9,j,1)
       end do
    end do
    
    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             us(i,j,k) = dy(j,1) * u(i,1,k) &
                       + dy(j,2) * u(i,2,k) &
                       + dy(j,3) * u(i,3,k) &
                       + dy(j,4) * u(i,4,k) &
                       + dy(j,5) * u(i,5,k) &
                       + dy(j,6) * u(i,6,k) &
                       + dy(j,7) * u(i,7,k) &
                       + dy(j,8) * u(i,8,k) &
                       + dy(j,9) * u(i,9,k)
          end do
       end do
    end do
    
    do k = 1, lx
       do i = 1, lx*lx
          ut(i,1,k) = dz(k,1) * u(i,1,1) &
                    + dz(k,2) * u(i,1,2) &
                    + dz(k,3) * u(i,1,3) &
                    + dz(k,4) * u(i,1,4) &
                    + dz(k,5) * u(i,1,5) &
                    + dz(k,6) * u(i,1,6) &
                    + dz(k,7) * u(i,1,7) &
                    + dz(k,8) * u(i,1,8) &
                    + dz(k,9) * u(i,1,9)
       end do
    end do
       
    do i = 1, lx * lx * lx
       ux(i,1,1) = w3(i,1,1) &
                   * ( drdx(i,1,1) * ur(i,1,1) &
                     + dsdx(i,1,1) * us(i,1,1) &
                     + dtdx(i,1,1) * ut(i,1,1) )
       uy(i,1,1) = w3(i,1,1) &
                   * ( dsdy(i,1,1) * us(i,1,1) &
                     + drdy(i,1,1) * ur(i,1,1) &
                     + dtdy(i,1,1) * ut(i,1,1) )
       uz(i,1,1) = w3(i,1,1) &
                   * ( dtdz(i,1,1) * ut(i,1,1) &
                     + drdz(i,1,1) * ur(i,1,1) &
                     + dsdz(i,1,1) * us(i,1,1) )
    end do
  end subroutine cpu_opgrad_lx9_single

  subroutine cpu_opgrad_lx8_single(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3)
    integer, parameter :: lx = 8
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          ur(i,j,1) = dx(i,1) * u(1,j,1) &
                    + dx(i,2) * u(2,j,1) &
                    + dx(i,3) * u(3,j,1) &
                    + dx(i,4) * u(4,j,1) &
                    + dx(i,5) * u(5,j,1) &
                    + dx(i,6) * u(6,j,1) &
                    + dx(i,7) * u(7,j,1) &
                    + dx(i,8) * u(8,j,1)
       end do
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             us(i,j,k) = dy(j,1) * u(i,1,k) &
                       + dy(j,2) * u(i,2,k) &
                       + dy(j,3) * u(i,3,k) &
                       + dy(j,4) * u(i,4,k) &
                       + dy(j,5) * u(i,5,k) &
                       + dy(j,6) * u(i,6,k) &
                       + dy(j,7) * u(i,7,k) &
                       + dy(j,8) * u(i,8,k)
          end do
       end do
    end do
    
    do k = 1, lx
       do i = 1, lx*lx
          ut(i,1,k) = dz(k,1) * u(i,1,1) &
                    + dz(k,2) * u(i,1,2) &
                    + dz(k,3) * u(i,1,3) &
                    + dz(k,4) * u(i,1,4) &
                    + dz(k,5) * u(i,1,5) &
                    + dz(k,6) * u(i,1,6) &
                    + dz(k,7) * u(i,1,7) &
                    + dz(k,8) * u(i,1,8)
       end do
    end do
    
    do i = 1, lx * lx * lx
       ux(i,1,1) = w3(i,1,1) &
                   * ( drdx(i,1,1) * ur(i,1,1) &
                     + dsdx(i,1,1) * us(i,1,1) &
                     + dtdx(i,1,1) * ut(i,1,1) )
       uy(i,1,1) = w3(i,1,1) &
                   * ( dsdy(i,1,1) * us(i,1,1) &
                     + drdy(i,1,1) * ur(i,1,1) &
                     + dtdy(i,1,1) * ut(i,1,1) )
       uz(i,1,1) = w3(i,1,1) &
                   * ( dtdz(i,1,1) * ut(i,1,1) &
                     + drdz(i,1,1) * ur(i,1,1) &
                     + dsdz(i,1,1) * us(i,1,1) )
    end do
  end subroutine cpu_opgrad_lx8_single

  subroutine cpu_opgrad_lx7_single(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3)
    integer, parameter :: lx = 7
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          ur(i,j,1) = dx(i,1) * u(1,j,1) &
                    + dx(i,2) * u(2,j,1) &
                    + dx(i,3) * u(3,j,1) &
                    + dx(i,4) * u(4,j,1) &
                    + dx(i,5) * u(5,j,1) &
                    + dx(i,6) * u(6,j,1) &
                    + dx(i,7) * u(7,j,1)
       end do
    end do
    
    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             us(i,j,k) = dy(j,1) * u(i,1,k) &
                       + dy(j,2) * u(i,2,k) &
                       + dy(j,3) * u(i,3,k) &
                       + dy(j,4) * u(i,4,k) &
                       + dy(j,5) * u(i,5,k) &
                       + dy(j,6) * u(i,6,k) &
                       + dy(j,7) * u(i,7,k)
          end do
       end do
    end do

    do k = 1, lx
       do i = 1, lx*lx
          ut(i,1,k) = dz(k,1) * u(i,1,1) &
                    + dz(k,2) * u(i,1,2) &
                    + dz(k,3) * u(i,1,3) &
                    + dz(k,4) * u(i,1,4) &
                    + dz(k,5) * u(i,1,5) &
                    + dz(k,6) * u(i,1,6) &
                    + dz(k,7) * u(i,1,7)
       end do
    end do

    do i = 1, lx * lx * lx
       ux(i,1,1) = w3(i,1,1) &
                   * ( drdx(i,1,1) * ur(i,1,1) &
                     + dsdx(i,1,1) * us(i,1,1) &
                     + dtdx(i,1,1) * ut(i,1,1) )
       uy(i,1,1) = w3(i,1,1) &
                   * ( dsdy(i,1,1) * us(i,1,1) &
                     + drdy(i,1,1) * ur(i,1,1) &
                     + dtdy(i,1,1) * ut(i,1,1) )
       uz(i,1,1) = w3(i,1,1) &
                   * ( dtdz(i,1,1) * ut(i,1,1) &
                     + drdz(i,1,1) * ur(i,1,1) &
                     + dsdz(i,1,1) * us(i,1,1) )
    end do
  end subroutine cpu_opgrad_lx7_single

  subroutine cpu_opgrad_lx6_single(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3)
    integer, parameter :: lx = 6
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          ur(i,j,1) = dx(i,1) * u(1,j,1) &
                    + dx(i,2) * u(2,j,1) &
                    + dx(i,3) * u(3,j,1) &
                    + dx(i,4) * u(4,j,1) &
                    + dx(i,5) * u(5,j,1) &
                    + dx(i,6) * u(6,j,1)
       end do
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             us(i,j,k) = dy(j,1) * u(i,1,k) &
                       + dy(j,2) * u(i,2,k) &
                       + dy(j,3) * u(i,3,k) &
                       + dy(j,4) * u(i,4,k) &
                       + dy(j,5) * u(i,5,k) &
                       + dy(j,6) * u(i,6,k)
          end do
       end do
    end do

    do k = 1, lx
       do i = 1, lx*lx
          ut(i,1,k) = dz(k,1) * u(i,1,1) &
                    + dz(k,2) * u(i,1,2) &
                    + dz(k,3) * u(i,1,3) &
                    + dz(k,4) * u(i,1,4) &
                    + dz(k,5) * u(i,1,5) &
                    + dz(k,6) * u(i,1,6)
       end do
    end do
    
    do i = 1, lx * lx * lx
       ux(i,1,1) = w3(i,1,1) &
                   * ( drdx(i,1,1) * ur(i,1,1) &
                     + dsdx(i,1,1) * us(i,1,1) &
                     + dtdx(i,1,1) * ut(i,1,1) )
       uy(i,1,1) = w3(i,1,1) &
                   * ( dsdy(i,1,1) * us(i,1,1) &
                     + drdy(i,1,1) * ur(i,1,1) &
                     + dtdy(i,1,1) * ut(i,1,1) )
       uz(i,1,1) = w3(i,1,1) &
                   * ( dtdz(i,1,1) * ut(i,1,1) &
                     + drdz(i,1,1) * ur(i,1,1) &
                     + dsdz(i,1,1) * us(i,1,1) )
    end do
  end subroutine cpu_opgrad_lx6_single

  subroutine cpu_opgrad_lx5_single(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3)
    integer, parameter :: lx = 5
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          ur(i,j,1) = dx(i,1) * u(1,j,1) &
                    + dx(i,2) * u(2,j,1) &
                    + dx(i,3) * u(3,j,1) &
                    + dx(i,4) * u(4,j,1) &
                    + dx(i,5) * u(5,j,1)
       end do
    end do

    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             us(i,j,k) = dy(j,1) * u(i,1,k) &
                       + dy(j,2) * u(i,2,k) &
                       + dy(j,3) * u(i,3,k) &
                       + dy(j,4) * u(i,4,k) &
                       + dy(j,5) * u(i,5,k)
          end do
       end do
    end do

    do k = 1, lx
       do i = 1, lx*lx
          ut(i,1,k) = dz(k,1) * u(i,1,1) &
                    + dz(k,2) * u(i,1,2) &
                    + dz(k,3) * u(i,1,3) &
                    + dz(k,4) * u(i,1,4) &
                    + dz(k,5) * u(i,1,5)
       end do
    end do
    
    do i = 1, lx * lx * lx
       ux(i,1,1) = w3(i,1,1) &
                   * ( drdx(i,1,1) * ur(i,1,1) &
                     + dsdx(i,1,1) * us(i,1,1) &
                     + dtdx(i,1,1) * ut(i,1,1) )
       uy(i,1,1) = w3(i,1,1) &
                   * ( dsdy(i,1,1) * us(i,1,1) &
                     + drdy(i,1,1) * ur(i,1,1) &
                     + dtdy(i,1,1) * ut(i,1,1) )
       uz(i,1,1) = w3(i,1,1) &
                   * ( dtdz(i,1,1) * ut(i,1,1) &
                     + drdz(i,1,1) * ur(i,1,1) &
                     + dsdz(i,1,1) * us(i,1,1) )
    end do
  end subroutine cpu_opgrad_lx5_single

  subroutine cpu_opgrad_lx4_single(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3)
    integer, parameter :: lx = 4
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          ur(i,j,1) = dx(i,1) * u(1,j,1) &
                    + dx(i,2) * u(2,j,1) &
                    + dx(i,3) * u(3,j,1) &
                    + dx(i,4) * u(4,j,1)
       end do
    end do
       
    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             us(i,j,k) = dy(j,1) * u(i,1,k) &
                       + dy(j,2) * u(i,2,k) &
                       + dy(j,3) * u(i,3,k) &
                       + dy(j,4) * u(i,4,k)
          end do
       end do
    end do
    
    do k = 1, lx
       do i = 1, lx*lx
          ut(i,1,k) = dz(k,1) * u(i,1,1) &
                    + dz(k,2) * u(i,1,2) &
                    + dz(k,3) * u(i,1,3) &
                    + dz(k,4) * u(i,1,4)
       end do
    end do

    do i = 1, lx * lx * lx
       ux(i,1,1) = w3(i,1,1) &
                   * ( drdx(i,1,1) * ur(i,1,1) &
                     + dsdx(i,1,1) * us(i,1,1) &
                     + dtdx(i,1,1) * ut(i,1,1) )
       uy(i,1,1) = w3(i,1,1) &
                   * ( dsdy(i,1,1) * us(i,1,1) &
                     + drdy(i,1,1) * ur(i,1,1) &
                     + dtdy(i,1,1) * ut(i,1,1) )
       uz(i,1,1) = w3(i,1,1) &
                   * ( dtdz(i,1,1) * ut(i,1,1) &
                     + drdz(i,1,1) * ur(i,1,1) &
                     + dsdz(i,1,1) * us(i,1,1) )
    end do
  end subroutine cpu_opgrad_lx4_single

  subroutine cpu_opgrad_lx3_single(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3)
    integer, parameter :: lx = 3
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          ur(i,j,1) = dx(i,1) * u(1,j,1) &
                    + dx(i,2) * u(2,j,1) &
                    + dx(i,3) * u(3,j,1)
       end do
    end do
    
    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             us(i,j,k) = dy(j,1) * u(i,1,k) &
                       + dy(j,2) * u(i,2,k) &
                       + dy(j,3) * u(i,3,k)
          end do
       end do
    end do
    
    do k = 1, lx
       do i = 1, lx*lx
          ut(i,1,k) = dz(k,1) * u(i,1,1) &
                    + dz(k,2) * u(i,1,2) &
                    + dz(k,3) * u(i,1,3)
       end do
    end do
       
    do i = 1, lx * lx * lx
       ux(i,1,1) = w3(i,1,1) &
                   * ( drdx(i,1,1) * ur(i,1,1) &
                     + dsdx(i,1,1) * us(i,1,1) &
                     + dtdx(i,1,1) * ut(i,1,1) )
       uy(i,1,1) = w3(i,1,1) &
                   * ( dsdy(i,1,1) * us(i,1,1) &
                     + drdy(i,1,1) * ur(i,1,1) &
                     + dtdy(i,1,1) * ut(i,1,1) )
       uz(i,1,1) = w3(i,1,1) &
                   * ( dtdz(i,1,1) * ut(i,1,1) &
                     + drdz(i,1,1) * ur(i,1,1) &
                     + dsdz(i,1,1) * us(i,1,1) )
    end do
  end subroutine cpu_opgrad_lx3_single

  subroutine cpu_opgrad_lx2_single(ux, uy, uz, u, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3)
    integer, parameter :: lx = 2
    real(kind=rp), dimension(lx, lx, lx), intent(inout) :: ux, uy, uz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: u
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    integer :: i, j, k

    do j = 1, lx * lx
       do i = 1, lx
          ur(i,j,1) = dx(i,1) * u(1,j,1) &
                    + dx(i,2) * u(2,j,1)
       end do
    end do
    
    do k = 1, lx
       do j = 1, lx
          do i = 1, lx
             us(i,j,k) = dy(j,1) * u(i,1,k) &
                       + dy(j,2) * u(i,2,k)
          end do
       end do
    end do
    
    do k = 1, lx
       do i = 1, lx*lx
          ut(i,1,k) = dz(k,1) * u(i,1,1) &
                    + dz(k,2) * u(i,1,2)
       end do
    end do
    
    do i = 1, lx * lx * lx
       ux(i,1,1) = w3(i,1,1) &
                   * ( drdx(i,1,1) * ur(i,1,1) &
                     + dsdx(i,1,1) * us(i,1,1) &
                     + dtdx(i,1,1) * ut(i,1,1) )
       uy(i,1,1) = w3(i,1,1) &
                   * ( dsdy(i,1,1) * us(i,1,1) &
                     + drdy(i,1,1) * ur(i,1,1) &
                     + dtdy(i,1,1) * ut(i,1,1) )
       uz(i,1,1) = w3(i,1,1) &
                   * ( dtdz(i,1,1) * ut(i,1,1) &
                     + drdz(i,1,1) * ur(i,1,1) &
                     + dsdz(i,1,1) * us(i,1,1) )
    end do

  end subroutine cpu_opgrad_lx2_single
  
end submodule cpu_opgrad
