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
!> Lambda2 kernels for SX-Aurora
submodule (opr_sx) sx_lambda2
  use math, only : pi
  implicit none

contains

  module subroutine opr_sx_lambda2(lambda2, u, v, w, coef)
    type(coef_t), intent(in) :: coef
    type(field_t), intent(inout) :: lambda2
    type(field_t), intent(in) :: u, v, w

    associate(Xh => coef%Xh, msh => coef%msh)
      select case(Xh%lx)
      case (18)
         call sx_lambda2_lx18(lambda2%x, u%x, v%x, w%x, &
                              Xh%dx, Xh%dy, Xh%dz, &
                              coef%drdx, coef%dsdx, coef%dtdx, &
                              coef%drdy, coef%dsdy, coef%dtdy, &
                              coef%drdz, coef%dsdz, coef%dtdz, &
                              Xh%w3, coef%B, msh%nelv)
      case (17)
         call sx_lambda2_lx17(lambda2%x, u%x, v%x, w%x, &
                              Xh%dx, Xh%dy, Xh%dz, &
                              coef%drdx, coef%dsdx, coef%dtdx, &
                              coef%drdy, coef%dsdy, coef%dtdy, &
                              coef%drdz, coef%dsdz, coef%dtdz, &
                              Xh%w3, coef%B, msh%nelv)
      case (16)
         call sx_lambda2_lx16(lambda2%x, u%x, v%x, w%x, &
                              Xh%dx, Xh%dy, Xh%dz, &
                              coef%drdx, coef%dsdx, coef%dtdx, &
                              coef%drdy, coef%dsdy, coef%dtdy, &
                              coef%drdz, coef%dsdz, coef%dtdz, &
                              Xh%w3, coef%B, msh%nelv)
      case (15)
         call sx_lambda2_lx15(lambda2%x, u%x, v%x, w%x, &
                              Xh%dx, Xh%dy, Xh%dz, &
                              coef%drdx, coef%dsdx, coef%dtdx, &
                              coef%drdy, coef%dsdy, coef%dtdy, &
                              coef%drdz, coef%dsdz, coef%dtdz, &
                              Xh%w3, coef%B, msh%nelv)
      case (14)
         call sx_lambda2_lx14(lambda2%x, u%x, v%x, w%x, &
                              Xh%dx, Xh%dy, Xh%dz, &
                              coef%drdx, coef%dsdx, coef%dtdx, &
                              coef%drdy, coef%dsdy, coef%dtdy, &
                              coef%drdz, coef%dsdz, coef%dtdz, &
                              Xh%w3, coef%B, msh%nelv)
      case (13)
         call sx_lambda2_lx13(lambda2%x, u%x, v%x, w%x, &
                              Xh%dx, Xh%dy, Xh%dz, &
                              coef%drdx, coef%dsdx, coef%dtdx, &
                              coef%drdy, coef%dsdy, coef%dtdy, &
                              coef%drdz, coef%dsdz, coef%dtdz, &
                              Xh%w3, coef%B, msh%nelv)
      case (12)
         call sx_lambda2_lx12(lambda2%x, u%x, v%x, w%x, &
                              Xh%dx, Xh%dy, Xh%dz, &
                              coef%drdx, coef%dsdx, coef%dtdx, &
                              coef%drdy, coef%dsdy, coef%dtdy, &
                              coef%drdz, coef%dsdz, coef%dtdz, &
                              Xh%w3, coef%B, msh%nelv)
      case (11)
         call sx_lambda2_lx11(lambda2%x, u%x, v%x, w%x, &
                              Xh%dx, Xh%dy, Xh%dz, &
                              coef%drdx, coef%dsdx, coef%dtdx, &
                              coef%drdy, coef%dsdy, coef%dtdy, &
                              coef%drdz, coef%dsdz, coef%dtdz, &
                              Xh%w3, coef%B, msh%nelv)
      case (10)
         call sx_lambda2_lx10(lambda2%x, u%x, v%x, w%x, &
                              Xh%dx, Xh%dy, Xh%dz, &
                              coef%drdx, coef%dsdx, coef%dtdx, &
                              coef%drdy, coef%dsdy, coef%dtdy, &
                              coef%drdz, coef%dsdz, coef%dtdz, &
                              Xh%w3, coef%B, msh%nelv)
      case (9)
         call sx_lambda2_lx9(lambda2%x, u%x, v%x, w%x, &
                             Xh%dx, Xh%dy, Xh%dz, &
                             coef%drdx, coef%dsdx, coef%dtdx, &
                             coef%drdy, coef%dsdy, coef%dtdy, &
                             coef%drdz, coef%dsdz, coef%dtdz, &
                             Xh%w3, coef%B, msh%nelv)
      case (8)
         call sx_lambda2_lx8(lambda2%x, u%x, v%x, w%x, &
                             Xh%dx, Xh%dy, Xh%dz, &
                             coef%drdx, coef%dsdx, coef%dtdx, &
                             coef%drdy, coef%dsdy, coef%dtdy, &
                             coef%drdz, coef%dsdz, coef%dtdz, &
                             Xh%w3, coef%B, msh%nelv)
      case (7)
         call sx_lambda2_lx7(lambda2%x, u%x, v%x, w%x, &
                             Xh%dx, Xh%dy, Xh%dz, &
                             coef%drdx, coef%dsdx, coef%dtdx, &
                             coef%drdy, coef%dsdy, coef%dtdy, &
                             coef%drdz, coef%dsdz, coef%dtdz, &
                             Xh%w3, coef%B, msh%nelv)
      case (6)
         call sx_lambda2_lx6(lambda2%x, u%x, v%x, w%x, &
                             Xh%dx, Xh%dy, Xh%dz, &
                             coef%drdx, coef%dsdx, coef%dtdx, &
                             coef%drdy, coef%dsdy, coef%dtdy, &
                             coef%drdz, coef%dsdz, coef%dtdz, &
                             Xh%w3, coef%B, msh%nelv)
      case (5)
         call sx_lambda2_lx5(lambda2%x, u%x, v%x, w%x, &
                             Xh%dx, Xh%dy, Xh%dz, &
                             coef%drdx, coef%dsdx, coef%dtdx, &
                             coef%drdy, coef%dsdy, coef%dtdy, &
                             coef%drdz, coef%dsdz, coef%dtdz, &
                             Xh%w3, coef%B, msh%nelv)
      case (4)
         call sx_lambda2_lx4(lambda2%x, u%x, v%x, w%x, &
                             Xh%dx, Xh%dy, Xh%dz, &
                             coef%drdx, coef%dsdx, coef%dtdx, &
                             coef%drdy, coef%dsdy, coef%dtdy, &
                             coef%drdz, coef%dsdz, coef%dtdz, &
                             Xh%w3, coef%B, msh%nelv)
      case (3)
         call sx_lambda2_lx3(lambda2%x, u%x, v%x, w%x, &
                             Xh%dx, Xh%dy, Xh%dz, &
                             coef%drdx, coef%dsdx, coef%dtdx, &
                             coef%drdy, coef%dsdy, coef%dtdy, &
                             coef%drdz, coef%dsdz, coef%dtdz, &
                             Xh%w3, coef%B, msh%nelv)
      case (2)
         call sx_lambda2_lx2(lambda2%x, u%x, v%x, w%x, &
                             Xh%dx, Xh%dy, Xh%dz, &
                             coef%drdx, coef%dsdx, coef%dtdx, &
                             coef%drdy, coef%dsdy, coef%dtdy, &
                             coef%drdz, coef%dsdz, coef%dtdz, &
                             Xh%w3, coef%B, msh%nelv)
      case default
         call sx_lambda2_lx(lambda2%x, u%x, v%x, w%x, &
                            Xh%dx, Xh%dy, Xh%dz, &
                            coef%drdx, coef%dsdx, coef%dtdx, &
                            coef%drdy, coef%dsdy, coef%dtdy, &
                            coef%drdz, coef%dsdz, coef%dtdz, &
                            Xh%w3, coef%B, msh%nelv, Xh%lx)
      end select
    end associate
    
  end subroutine opr_sx_lambda2  
  
  subroutine sx_lambda2_lx(lambda2, u, v, w, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, cB, n, lx)
    integer, intent(in) :: n, lx
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: lambda2
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u    
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: v
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: w
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cB
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: grad(lx*lx*lx,3,3)
    real(kind=rp) :: eigen(3), B, C, D, q, r, theta, l2
    real(kind=rp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=rp) :: a11, a22, a33, a12, a13, a23
    real(kind=rp) :: msk1, msk2, msk3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
    real(kind=rp) :: tmp1, tmp2, tmp3
    integer :: e, i, j, k, l
    
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do k = 1, lx
                tmp1 = tmp1 + dx(i,k) * u(k,j,1,e)
                tmp2 = tmp2 + dx(i,k) * v(k,j,1,e)
                tmp3 = tmp3 + dx(i,k) * w(k,j,1,e)
             end do
             ur(i,j,1) = tmp1
             vr(i,j,1) = tmp2
             wr(i,j,1) = tmp3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp1 = 0.0_rp
                tmp2 = 0.0_rp
                tmp3 = 0.0_rp
                do l = 1, lx
                   tmp1 = tmp1 + dy(j,l) * u(i,l,k,e)
                   tmp2 = tmp2 + dy(j,l) * v(i,l,k,e)
                   tmp3 = tmp3 + dy(j,l) * w(i,l,k,e)
                end do
                us(i,j,k) = tmp1
                vs(i,j,k) = tmp2
                ws(i,j,k) = tmp3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do l = 1, lx
                tmp1 = tmp1 + dz(k,l) * u(i,1,l,e)
                tmp2 = tmp2 + dz(k,l) * v(i,1,l,e)
                tmp3 = tmp3 + dz(k,l) * w(i,1,l,e)
             end do
             ut(i,1,k) = tmp1
             vt(i,1,k) = tmp2
             wt(i,1,k) = tmp3
          end do
       end do

       do i = 1, lx * lx * lx
          grad(1,1,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          grad(1,1,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          grad(1,1,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )

          grad(1,2,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * vr(i,1,1) &
                        + dsdx(i,1,1,e) * vs(i,1,1) &
                        + dtdx(i,1,1,e) * vt(i,1,1) )
          grad(1,2,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * vs(i,1,1) &
                        + drdy(i,1,1,e) * vr(i,1,1) &
                        + dtdy(i,1,1,e) * vt(i,1,1) )
          grad(1,2,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * vt(i,1,1) &
                        + drdz(i,1,1,e) * vr(i,1,1) &
                        + dsdz(i,1,1,e) * vs(i,1,1) )

          grad(1,3,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * wr(i,1,1) &
                        + dsdx(i,1,1,e) * ws(i,1,1) &
                        + dtdx(i,1,1,e) * wt(i,1,1) )
          grad(1,3,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * ws(i,1,1) &
                        + drdy(i,1,1,e) * wr(i,1,1) &
                        + dtdy(i,1,1,e) * wt(i,1,1) )
          grad(1,3,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * wt(i,1,1) &
                        + drdz(i,1,1,e) * wr(i,1,1) &
                        + dsdz(i,1,1,e) * ws(i,1,1) )          
       end do
       

       do i = 1, lx * lx * lx
          s11 = grad(i,1,1)
          s22 = grad(i,2,2)
          s33 = grad(i,3,3)

          
          s12 = 0.5*(grad(i,1,2) + grad(i,2,1))
          s13 = 0.5*(grad(i,1,3) + grad(i,3,1))
          s23 = 0.5*(grad(i,2,3) + grad(i,3,2))
          
          o12 = 0.5*(grad(i,1,2) - grad(i,2,1))
          o13 = 0.5*(grad(i,1,3) - grad(i,3,1))
          o23 = 0.5*(grad(i,2,3) - grad(i,3,2))

          a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
          a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23
          a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23

          a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
          a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
          a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


          B = -(a11 + a22 + a33)
          C = -(a12*a12 + a13*a13 + a23*a23 &
              - a11 * a22 - a11 * a33 - a22 * a33)
          D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 &
              - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33)


          q = (3.0 * C - B*B) / 9.0
          r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0
          theta = acos( r / sqrt(-q*q*q) )
          
          eigen(1) = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0
          eigen(2) = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0
          eigen(3) = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0

          msk1 = merge(1.0_rp, 0.0_rp, eigen(2) .le. eigen(1) &
                       .and. eigen(1) .le. eigen(3) .or.  eigen(3) &
                       .le. eigen(1) .and. eigen(1) .le. eigen(2) )
          msk2 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(2) &
                       .and. eigen(2) .le. eigen(3) .or. eigen(3) &
                       .le. eigen(2) .and. eigen(2) .le. eigen(1))
          msk3 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(3) &
                       .and. eigen(3) .le. eigen(2) .or. eigen(2) &
                       .le. eigen(3) .and. eigen(3) .le. eigen(1))

          l2 = msk1 * eigen(1) + msk2 * eigen(2) + msk3 * eigen(3)

          lambda2(i,1,1,e) = l2/(cB(i,1,1,e)**2)
       end do
    end do
  end subroutine sx_lambda2_lx

  subroutine sx_lambda2_lx18(lambda2, u, v, w, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, cB, n)
    integer, parameter :: lx = 18
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: lambda2
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u    
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: v
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: w
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cB
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: grad(lx*lx*lx,3,3)
    real(kind=rp) :: eigen(3), B, C, D, q, r, theta, l2
    real(kind=rp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=rp) :: a11, a22, a33, a12, a13, a23
    real(kind=rp) :: msk1, msk2, msk3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
    real(kind=rp) :: tmp1, tmp2, tmp3
    integer :: e, i, j, k, l
    
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do k = 1, lx
                tmp1 = tmp1 + dx(i,k) * u(k,j,1,e)
                tmp2 = tmp2 + dx(i,k) * v(k,j,1,e)
                tmp3 = tmp3 + dx(i,k) * w(k,j,1,e)
             end do
             ur(i,j,1) = tmp1
             vr(i,j,1) = tmp2
             wr(i,j,1) = tmp3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp1 = 0.0_rp
                tmp2 = 0.0_rp
                tmp3 = 0.0_rp
                do l = 1, lx
                   tmp1 = tmp1 + dy(j,l) * u(i,l,k,e)
                   tmp2 = tmp2 + dy(j,l) * v(i,l,k,e)
                   tmp3 = tmp3 + dy(j,l) * w(i,l,k,e)
                end do
                us(i,j,k) = tmp1
                vs(i,j,k) = tmp2
                ws(i,j,k) = tmp3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do l = 1, lx
                tmp1 = tmp1 + dz(k,l) * u(i,1,l,e)
                tmp2 = tmp2 + dz(k,l) * v(i,1,l,e)
                tmp3 = tmp3 + dz(k,l) * w(i,1,l,e)
             end do
             ut(i,1,k) = tmp1
             vt(i,1,k) = tmp2
             wt(i,1,k) = tmp3
          end do
       end do

       do i = 1, lx * lx * lx
          grad(1,1,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          grad(1,1,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          grad(1,1,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )

          grad(1,2,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * vr(i,1,1) &
                        + dsdx(i,1,1,e) * vs(i,1,1) &
                        + dtdx(i,1,1,e) * vt(i,1,1) )
          grad(1,2,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * vs(i,1,1) &
                        + drdy(i,1,1,e) * vr(i,1,1) &
                        + dtdy(i,1,1,e) * vt(i,1,1) )
          grad(1,2,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * vt(i,1,1) &
                        + drdz(i,1,1,e) * vr(i,1,1) &
                        + dsdz(i,1,1,e) * vs(i,1,1) )

          grad(1,3,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * wr(i,1,1) &
                        + dsdx(i,1,1,e) * ws(i,1,1) &
                        + dtdx(i,1,1,e) * wt(i,1,1) )
          grad(1,3,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * ws(i,1,1) &
                        + drdy(i,1,1,e) * wr(i,1,1) &
                        + dtdy(i,1,1,e) * wt(i,1,1) )
          grad(1,3,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * wt(i,1,1) &
                        + drdz(i,1,1,e) * wr(i,1,1) &
                        + dsdz(i,1,1,e) * ws(i,1,1) )          
       end do
       

       do i = 1, lx * lx * lx
          s11 = grad(i,1,1)
          s22 = grad(i,2,2)
          s33 = grad(i,3,3)

          
          s12 = 0.5*(grad(i,1,2) + grad(i,2,1))
          s13 = 0.5*(grad(i,1,3) + grad(i,3,1))
          s23 = 0.5*(grad(i,2,3) + grad(i,3,2))
          
          o12 = 0.5*(grad(i,1,2) - grad(i,2,1))
          o13 = 0.5*(grad(i,1,3) - grad(i,3,1))
          o23 = 0.5*(grad(i,2,3) - grad(i,3,2))

          a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
          a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23
          a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23

          a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
          a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
          a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


          B = -(a11 + a22 + a33)
          C = -(a12*a12 + a13*a13 + a23*a23 &
              - a11 * a22 - a11 * a33 - a22 * a33)
          D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 &
              - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33)


          q = (3.0 * C - B*B) / 9.0
          r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0
          theta = acos( r / sqrt(-q*q*q) )
          
          eigen(1) = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0
          eigen(2) = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0
          eigen(3) = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0

          msk1 = merge(1.0_rp, 0.0_rp, eigen(2) .le. eigen(1) &
                       .and. eigen(1) .le. eigen(3) .or.  eigen(3) &
                       .le. eigen(1) .and. eigen(1) .le. eigen(2) )
          msk2 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(2) &
                       .and. eigen(2) .le. eigen(3) .or. eigen(3) &
                       .le. eigen(2) .and. eigen(2) .le. eigen(1))
          msk3 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(3) &
                       .and. eigen(3) .le. eigen(2) .or. eigen(2) &
                       .le. eigen(3) .and. eigen(3) .le. eigen(1))

          l2 = msk1 * eigen(1) + msk2 * eigen(2) + msk3 * eigen(3)

          lambda2(i,1,1,e) = l2/(cB(i,1,1,e)**2)
       end do
    end do
  end subroutine sx_lambda2_lx18

  subroutine sx_lambda2_lx17(lambda2, u, v, w, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, cB, n)
    integer, parameter :: lx = 17
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: lambda2
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u    
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: v
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: w
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cB
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: grad(lx*lx*lx,3,3)
    real(kind=rp) :: eigen(3), B, C, D, q, r, theta, l2
    real(kind=rp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=rp) :: a11, a22, a33, a12, a13, a23
    real(kind=rp) :: msk1, msk2, msk3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
    real(kind=rp) :: tmp1, tmp2, tmp3
    integer :: e, i, j, k, l
    
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do k = 1, lx
                tmp1 = tmp1 + dx(i,k) * u(k,j,1,e)
                tmp2 = tmp2 + dx(i,k) * v(k,j,1,e)
                tmp3 = tmp3 + dx(i,k) * w(k,j,1,e)
             end do
             ur(i,j,1) = tmp1
             vr(i,j,1) = tmp2
             wr(i,j,1) = tmp3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp1 = 0.0_rp
                tmp2 = 0.0_rp
                tmp3 = 0.0_rp
                do l = 1, lx
                   tmp1 = tmp1 + dy(j,l) * u(i,l,k,e)
                   tmp2 = tmp2 + dy(j,l) * v(i,l,k,e)
                   tmp3 = tmp3 + dy(j,l) * w(i,l,k,e)
                end do
                us(i,j,k) = tmp1
                vs(i,j,k) = tmp2
                ws(i,j,k) = tmp3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do l = 1, lx
                tmp1 = tmp1 + dz(k,l) * u(i,1,l,e)
                tmp2 = tmp2 + dz(k,l) * v(i,1,l,e)
                tmp3 = tmp3 + dz(k,l) * w(i,1,l,e)
             end do
             ut(i,1,k) = tmp1
             vt(i,1,k) = tmp2
             wt(i,1,k) = tmp3
          end do
       end do

       do i = 1, lx * lx * lx
          grad(1,1,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          grad(1,1,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          grad(1,1,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )

          grad(1,2,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * vr(i,1,1) &
                        + dsdx(i,1,1,e) * vs(i,1,1) &
                        + dtdx(i,1,1,e) * vt(i,1,1) )
          grad(1,2,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * vs(i,1,1) &
                        + drdy(i,1,1,e) * vr(i,1,1) &
                        + dtdy(i,1,1,e) * vt(i,1,1) )
          grad(1,2,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * vt(i,1,1) &
                        + drdz(i,1,1,e) * vr(i,1,1) &
                        + dsdz(i,1,1,e) * vs(i,1,1) )

          grad(1,3,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * wr(i,1,1) &
                        + dsdx(i,1,1,e) * ws(i,1,1) &
                        + dtdx(i,1,1,e) * wt(i,1,1) )
          grad(1,3,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * ws(i,1,1) &
                        + drdy(i,1,1,e) * wr(i,1,1) &
                        + dtdy(i,1,1,e) * wt(i,1,1) )
          grad(1,3,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * wt(i,1,1) &
                        + drdz(i,1,1,e) * wr(i,1,1) &
                        + dsdz(i,1,1,e) * ws(i,1,1) )          
       end do
       

       do i = 1, lx * lx * lx
          s11 = grad(i,1,1)
          s22 = grad(i,2,2)
          s33 = grad(i,3,3)

          
          s12 = 0.5*(grad(i,1,2) + grad(i,2,1))
          s13 = 0.5*(grad(i,1,3) + grad(i,3,1))
          s23 = 0.5*(grad(i,2,3) + grad(i,3,2))
          
          o12 = 0.5*(grad(i,1,2) - grad(i,2,1))
          o13 = 0.5*(grad(i,1,3) - grad(i,3,1))
          o23 = 0.5*(grad(i,2,3) - grad(i,3,2))

          a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
          a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23
          a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23

          a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
          a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
          a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


          B = -(a11 + a22 + a33)
          C = -(a12*a12 + a13*a13 + a23*a23 &
              - a11 * a22 - a11 * a33 - a22 * a33)
          D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 &
              - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33)


          q = (3.0 * C - B*B) / 9.0
          r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0
          theta = acos( r / sqrt(-q*q*q) )
          
          eigen(1) = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0
          eigen(2) = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0
          eigen(3) = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0

          msk1 = merge(1.0_rp, 0.0_rp, eigen(2) .le. eigen(1) &
                       .and. eigen(1) .le. eigen(3) .or.  eigen(3) &
                       .le. eigen(1) .and. eigen(1) .le. eigen(2) )
          msk2 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(2) &
                       .and. eigen(2) .le. eigen(3) .or. eigen(3) &
                       .le. eigen(2) .and. eigen(2) .le. eigen(1))
          msk3 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(3) &
                       .and. eigen(3) .le. eigen(2) .or. eigen(2) &
                       .le. eigen(3) .and. eigen(3) .le. eigen(1))

          l2 = msk1 * eigen(1) + msk2 * eigen(2) + msk3 * eigen(3)

          lambda2(i,1,1,e) = l2/(cB(i,1,1,e)**2)
       end do
    end do
  end subroutine sx_lambda2_lx17

  subroutine sx_lambda2_lx16(lambda2, u, v, w, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, cB, n)
    integer, parameter :: lx = 16
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: lambda2
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u    
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: v
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: w
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cB
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: grad(lx*lx*lx,3,3)
    real(kind=rp) :: eigen(3), B, C, D, q, r, theta, l2
    real(kind=rp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=rp) :: a11, a22, a33, a12, a13, a23
    real(kind=rp) :: msk1, msk2, msk3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
    real(kind=rp) :: tmp1, tmp2, tmp3
    integer :: e, i, j, k, l
    
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do k = 1, lx
                tmp1 = tmp1 + dx(i,k) * u(k,j,1,e)
                tmp2 = tmp2 + dx(i,k) * v(k,j,1,e)
                tmp3 = tmp3 + dx(i,k) * w(k,j,1,e)
             end do
             ur(i,j,1) = tmp1
             vr(i,j,1) = tmp2
             wr(i,j,1) = tmp3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp1 = 0.0_rp
                tmp2 = 0.0_rp
                tmp3 = 0.0_rp
                do l = 1, lx
                   tmp1 = tmp1 + dy(j,l) * u(i,l,k,e)
                   tmp2 = tmp2 + dy(j,l) * v(i,l,k,e)
                   tmp3 = tmp3 + dy(j,l) * w(i,l,k,e)
                end do
                us(i,j,k) = tmp1
                vs(i,j,k) = tmp2
                ws(i,j,k) = tmp3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do l = 1, lx
                tmp1 = tmp1 + dz(k,l) * u(i,1,l,e)
                tmp2 = tmp2 + dz(k,l) * v(i,1,l,e)
                tmp3 = tmp3 + dz(k,l) * w(i,1,l,e)
             end do
             ut(i,1,k) = tmp1
             vt(i,1,k) = tmp2
             wt(i,1,k) = tmp3
          end do
       end do

       do i = 1, lx * lx * lx
          grad(1,1,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          grad(1,1,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          grad(1,1,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )

          grad(1,2,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * vr(i,1,1) &
                        + dsdx(i,1,1,e) * vs(i,1,1) &
                        + dtdx(i,1,1,e) * vt(i,1,1) )
          grad(1,2,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * vs(i,1,1) &
                        + drdy(i,1,1,e) * vr(i,1,1) &
                        + dtdy(i,1,1,e) * vt(i,1,1) )
          grad(1,2,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * vt(i,1,1) &
                        + drdz(i,1,1,e) * vr(i,1,1) &
                        + dsdz(i,1,1,e) * vs(i,1,1) )

          grad(1,3,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * wr(i,1,1) &
                        + dsdx(i,1,1,e) * ws(i,1,1) &
                        + dtdx(i,1,1,e) * wt(i,1,1) )
          grad(1,3,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * ws(i,1,1) &
                        + drdy(i,1,1,e) * wr(i,1,1) &
                        + dtdy(i,1,1,e) * wt(i,1,1) )
          grad(1,3,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * wt(i,1,1) &
                        + drdz(i,1,1,e) * wr(i,1,1) &
                        + dsdz(i,1,1,e) * ws(i,1,1) )          
       end do
       

       do i = 1, lx * lx * lx
          s11 = grad(i,1,1)
          s22 = grad(i,2,2)
          s33 = grad(i,3,3)

          
          s12 = 0.5*(grad(i,1,2) + grad(i,2,1))
          s13 = 0.5*(grad(i,1,3) + grad(i,3,1))
          s23 = 0.5*(grad(i,2,3) + grad(i,3,2))
          
          o12 = 0.5*(grad(i,1,2) - grad(i,2,1))
          o13 = 0.5*(grad(i,1,3) - grad(i,3,1))
          o23 = 0.5*(grad(i,2,3) - grad(i,3,2))

          a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
          a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23
          a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23

          a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
          a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
          a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


          B = -(a11 + a22 + a33)
          C = -(a12*a12 + a13*a13 + a23*a23 &
              - a11 * a22 - a11 * a33 - a22 * a33)
          D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 &
              - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33)


          q = (3.0 * C - B*B) / 9.0
          r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0
          theta = acos( r / sqrt(-q*q*q) )
          
          eigen(1) = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0
          eigen(2) = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0
          eigen(3) = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0

          msk1 = merge(1.0_rp, 0.0_rp, eigen(2) .le. eigen(1) &
                       .and. eigen(1) .le. eigen(3) .or.  eigen(3) &
                       .le. eigen(1) .and. eigen(1) .le. eigen(2) )
          msk2 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(2) &
                       .and. eigen(2) .le. eigen(3) .or. eigen(3) &
                       .le. eigen(2) .and. eigen(2) .le. eigen(1))
          msk3 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(3) &
                       .and. eigen(3) .le. eigen(2) .or. eigen(2) &
                       .le. eigen(3) .and. eigen(3) .le. eigen(1))

          l2 = msk1 * eigen(1) + msk2 * eigen(2) + msk3 * eigen(3)

          lambda2(i,1,1,e) = l2/(cB(i,1,1,e)**2)
       end do
    end do
  end subroutine sx_lambda2_lx16

  subroutine sx_lambda2_lx15(lambda2, u, v, w, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, cB, n)
    integer, parameter :: lx = 15
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: lambda2
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u    
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: v
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: w
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cB
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: grad(lx*lx*lx,3,3)
    real(kind=rp) :: eigen(3), B, C, D, q, r, theta, l2
    real(kind=rp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=rp) :: a11, a22, a33, a12, a13, a23
    real(kind=rp) :: msk1, msk2, msk3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
    real(kind=rp) :: tmp1, tmp2, tmp3
    integer :: e, i, j, k, l
    
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do k = 1, lx
                tmp1 = tmp1 + dx(i,k) * u(k,j,1,e)
                tmp2 = tmp2 + dx(i,k) * v(k,j,1,e)
                tmp3 = tmp3 + dx(i,k) * w(k,j,1,e)
             end do
             ur(i,j,1) = tmp1
             vr(i,j,1) = tmp2
             wr(i,j,1) = tmp3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp1 = 0.0_rp
                tmp2 = 0.0_rp
                tmp3 = 0.0_rp
                do l = 1, lx
                   tmp1 = tmp1 + dy(j,l) * u(i,l,k,e)
                   tmp2 = tmp2 + dy(j,l) * v(i,l,k,e)
                   tmp3 = tmp3 + dy(j,l) * w(i,l,k,e)
                end do
                us(i,j,k) = tmp1
                vs(i,j,k) = tmp2
                ws(i,j,k) = tmp3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do l = 1, lx
                tmp1 = tmp1 + dz(k,l) * u(i,1,l,e)
                tmp2 = tmp2 + dz(k,l) * v(i,1,l,e)
                tmp3 = tmp3 + dz(k,l) * w(i,1,l,e)
             end do
             ut(i,1,k) = tmp1
             vt(i,1,k) = tmp2
             wt(i,1,k) = tmp3
          end do
       end do

       do i = 1, lx * lx * lx
          grad(1,1,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          grad(1,1,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          grad(1,1,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )

          grad(1,2,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * vr(i,1,1) &
                        + dsdx(i,1,1,e) * vs(i,1,1) &
                        + dtdx(i,1,1,e) * vt(i,1,1) )
          grad(1,2,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * vs(i,1,1) &
                        + drdy(i,1,1,e) * vr(i,1,1) &
                        + dtdy(i,1,1,e) * vt(i,1,1) )
          grad(1,2,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * vt(i,1,1) &
                        + drdz(i,1,1,e) * vr(i,1,1) &
                        + dsdz(i,1,1,e) * vs(i,1,1) )

          grad(1,3,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * wr(i,1,1) &
                        + dsdx(i,1,1,e) * ws(i,1,1) &
                        + dtdx(i,1,1,e) * wt(i,1,1) )
          grad(1,3,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * ws(i,1,1) &
                        + drdy(i,1,1,e) * wr(i,1,1) &
                        + dtdy(i,1,1,e) * wt(i,1,1) )
          grad(1,3,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * wt(i,1,1) &
                        + drdz(i,1,1,e) * wr(i,1,1) &
                        + dsdz(i,1,1,e) * ws(i,1,1) )          
       end do
       

       do i = 1, lx * lx * lx
          s11 = grad(i,1,1)
          s22 = grad(i,2,2)
          s33 = grad(i,3,3)

          
          s12 = 0.5*(grad(i,1,2) + grad(i,2,1))
          s13 = 0.5*(grad(i,1,3) + grad(i,3,1))
          s23 = 0.5*(grad(i,2,3) + grad(i,3,2))
          
          o12 = 0.5*(grad(i,1,2) - grad(i,2,1))
          o13 = 0.5*(grad(i,1,3) - grad(i,3,1))
          o23 = 0.5*(grad(i,2,3) - grad(i,3,2))

          a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
          a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23
          a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23

          a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
          a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
          a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


          B = -(a11 + a22 + a33)
          C = -(a12*a12 + a13*a13 + a23*a23 &
              - a11 * a22 - a11 * a33 - a22 * a33)
          D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 &
              - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33)


          q = (3.0 * C - B*B) / 9.0
          r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0
          theta = acos( r / sqrt(-q*q*q) )
          
          eigen(1) = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0
          eigen(2) = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0
          eigen(3) = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0

          msk1 = merge(1.0_rp, 0.0_rp, eigen(2) .le. eigen(1) &
                       .and. eigen(1) .le. eigen(3) .or.  eigen(3) &
                       .le. eigen(1) .and. eigen(1) .le. eigen(2) )
          msk2 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(2) &
                       .and. eigen(2) .le. eigen(3) .or. eigen(3) &
                       .le. eigen(2) .and. eigen(2) .le. eigen(1))
          msk3 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(3) &
                       .and. eigen(3) .le. eigen(2) .or. eigen(2) &
                       .le. eigen(3) .and. eigen(3) .le. eigen(1))

          l2 = msk1 * eigen(1) + msk2 * eigen(2) + msk3 * eigen(3)

          lambda2(i,1,1,e) = l2/(cB(i,1,1,e)**2)
       end do
    end do
  end subroutine sx_lambda2_lx15

  subroutine sx_lambda2_lx14(lambda2, u, v, w, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, cB, n)
    integer, parameter :: lx = 14
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: lambda2
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u    
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: v
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: w
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cB
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: grad(lx*lx*lx,3,3)
    real(kind=rp) :: eigen(3), B, C, D, q, r, theta, l2
    real(kind=rp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=rp) :: a11, a22, a33, a12, a13, a23
    real(kind=rp) :: msk1, msk2, msk3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
    real(kind=rp) :: tmp1, tmp2, tmp3
    integer :: e, i, j, k, l
    
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do k = 1, lx
                tmp1 = tmp1 + dx(i,k) * u(k,j,1,e)
                tmp2 = tmp2 + dx(i,k) * v(k,j,1,e)
                tmp3 = tmp3 + dx(i,k) * w(k,j,1,e)
             end do
             ur(i,j,1) = tmp1
             vr(i,j,1) = tmp2
             wr(i,j,1) = tmp3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp1 = 0.0_rp
                tmp2 = 0.0_rp
                tmp3 = 0.0_rp
                do l = 1, lx
                   tmp1 = tmp1 + dy(j,l) * u(i,l,k,e)
                   tmp2 = tmp2 + dy(j,l) * v(i,l,k,e)
                   tmp3 = tmp3 + dy(j,l) * w(i,l,k,e)
                end do
                us(i,j,k) = tmp1
                vs(i,j,k) = tmp2
                ws(i,j,k) = tmp3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do l = 1, lx
                tmp1 = tmp1 + dz(k,l) * u(i,1,l,e)
                tmp2 = tmp2 + dz(k,l) * v(i,1,l,e)
                tmp3 = tmp3 + dz(k,l) * w(i,1,l,e)
             end do
             ut(i,1,k) = tmp1
             vt(i,1,k) = tmp2
             wt(i,1,k) = tmp3
          end do
       end do

       do i = 1, lx * lx * lx
          grad(1,1,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          grad(1,1,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          grad(1,1,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )

          grad(1,2,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * vr(i,1,1) &
                        + dsdx(i,1,1,e) * vs(i,1,1) &
                        + dtdx(i,1,1,e) * vt(i,1,1) )
          grad(1,2,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * vs(i,1,1) &
                        + drdy(i,1,1,e) * vr(i,1,1) &
                        + dtdy(i,1,1,e) * vt(i,1,1) )
          grad(1,2,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * vt(i,1,1) &
                        + drdz(i,1,1,e) * vr(i,1,1) &
                        + dsdz(i,1,1,e) * vs(i,1,1) )

          grad(1,3,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * wr(i,1,1) &
                        + dsdx(i,1,1,e) * ws(i,1,1) &
                        + dtdx(i,1,1,e) * wt(i,1,1) )
          grad(1,3,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * ws(i,1,1) &
                        + drdy(i,1,1,e) * wr(i,1,1) &
                        + dtdy(i,1,1,e) * wt(i,1,1) )
          grad(1,3,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * wt(i,1,1) &
                        + drdz(i,1,1,e) * wr(i,1,1) &
                        + dsdz(i,1,1,e) * ws(i,1,1) )          
       end do
       

       do i = 1, lx * lx * lx
          s11 = grad(i,1,1)
          s22 = grad(i,2,2)
          s33 = grad(i,3,3)

          
          s12 = 0.5*(grad(i,1,2) + grad(i,2,1))
          s13 = 0.5*(grad(i,1,3) + grad(i,3,1))
          s23 = 0.5*(grad(i,2,3) + grad(i,3,2))
          
          o12 = 0.5*(grad(i,1,2) - grad(i,2,1))
          o13 = 0.5*(grad(i,1,3) - grad(i,3,1))
          o23 = 0.5*(grad(i,2,3) - grad(i,3,2))

          a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
          a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23
          a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23

          a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
          a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
          a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


          B = -(a11 + a22 + a33)
          C = -(a12*a12 + a13*a13 + a23*a23 &
              - a11 * a22 - a11 * a33 - a22 * a33)
          D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 &
              - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33)


          q = (3.0 * C - B*B) / 9.0
          r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0
          theta = acos( r / sqrt(-q*q*q) )
          
          eigen(1) = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0
          eigen(2) = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0
          eigen(3) = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0

          msk1 = merge(1.0_rp, 0.0_rp, eigen(2) .le. eigen(1) &
                       .and. eigen(1) .le. eigen(3) .or.  eigen(3) &
                       .le. eigen(1) .and. eigen(1) .le. eigen(2) )
          msk2 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(2) &
                       .and. eigen(2) .le. eigen(3) .or. eigen(3) &
                       .le. eigen(2) .and. eigen(2) .le. eigen(1))
          msk3 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(3) &
                       .and. eigen(3) .le. eigen(2) .or. eigen(2) &
                       .le. eigen(3) .and. eigen(3) .le. eigen(1))

          l2 = msk1 * eigen(1) + msk2 * eigen(2) + msk3 * eigen(3)

          lambda2(i,1,1,e) = l2/(cB(i,1,1,e)**2)
       end do
    end do
  end subroutine sx_lambda2_lx14

  subroutine sx_lambda2_lx13(lambda2, u, v, w, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, cB, n)
    integer, parameter :: lx = 13
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: lambda2
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u    
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: v
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: w
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cB
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: grad(lx*lx*lx,3,3)
    real(kind=rp) :: eigen(3), B, C, D, q, r, theta, l2
    real(kind=rp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=rp) :: a11, a22, a33, a12, a13, a23
    real(kind=rp) :: msk1, msk2, msk3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
    real(kind=rp) :: tmp1, tmp2, tmp3
    integer :: e, i, j, k, l
    
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do k = 1, lx
                tmp1 = tmp1 + dx(i,k) * u(k,j,1,e)
                tmp2 = tmp2 + dx(i,k) * v(k,j,1,e)
                tmp3 = tmp3 + dx(i,k) * w(k,j,1,e)
             end do
             ur(i,j,1) = tmp1
             vr(i,j,1) = tmp2
             wr(i,j,1) = tmp3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp1 = 0.0_rp
                tmp2 = 0.0_rp
                tmp3 = 0.0_rp
                do l = 1, lx
                   tmp1 = tmp1 + dy(j,l) * u(i,l,k,e)
                   tmp2 = tmp2 + dy(j,l) * v(i,l,k,e)
                   tmp3 = tmp3 + dy(j,l) * w(i,l,k,e)
                end do
                us(i,j,k) = tmp1
                vs(i,j,k) = tmp2
                ws(i,j,k) = tmp3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do l = 1, lx
                tmp1 = tmp1 + dz(k,l) * u(i,1,l,e)
                tmp2 = tmp2 + dz(k,l) * v(i,1,l,e)
                tmp3 = tmp3 + dz(k,l) * w(i,1,l,e)
             end do
             ut(i,1,k) = tmp1
             vt(i,1,k) = tmp2
             wt(i,1,k) = tmp3
          end do
       end do

       do i = 1, lx * lx * lx
          grad(1,1,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          grad(1,1,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          grad(1,1,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )

          grad(1,2,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * vr(i,1,1) &
                        + dsdx(i,1,1,e) * vs(i,1,1) &
                        + dtdx(i,1,1,e) * vt(i,1,1) )
          grad(1,2,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * vs(i,1,1) &
                        + drdy(i,1,1,e) * vr(i,1,1) &
                        + dtdy(i,1,1,e) * vt(i,1,1) )
          grad(1,2,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * vt(i,1,1) &
                        + drdz(i,1,1,e) * vr(i,1,1) &
                        + dsdz(i,1,1,e) * vs(i,1,1) )

          grad(1,3,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * wr(i,1,1) &
                        + dsdx(i,1,1,e) * ws(i,1,1) &
                        + dtdx(i,1,1,e) * wt(i,1,1) )
          grad(1,3,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * ws(i,1,1) &
                        + drdy(i,1,1,e) * wr(i,1,1) &
                        + dtdy(i,1,1,e) * wt(i,1,1) )
          grad(1,3,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * wt(i,1,1) &
                        + drdz(i,1,1,e) * wr(i,1,1) &
                        + dsdz(i,1,1,e) * ws(i,1,1) )          
       end do
       

       do i = 1, lx * lx * lx
          s11 = grad(i,1,1)
          s22 = grad(i,2,2)
          s33 = grad(i,3,3)

          
          s12 = 0.5*(grad(i,1,2) + grad(i,2,1))
          s13 = 0.5*(grad(i,1,3) + grad(i,3,1))
          s23 = 0.5*(grad(i,2,3) + grad(i,3,2))
          
          o12 = 0.5*(grad(i,1,2) - grad(i,2,1))
          o13 = 0.5*(grad(i,1,3) - grad(i,3,1))
          o23 = 0.5*(grad(i,2,3) - grad(i,3,2))

          a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
          a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23
          a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23

          a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
          a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
          a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


          B = -(a11 + a22 + a33)
          C = -(a12*a12 + a13*a13 + a23*a23 &
              - a11 * a22 - a11 * a33 - a22 * a33)
          D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 &
              - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33)


          q = (3.0 * C - B*B) / 9.0
          r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0
          theta = acos( r / sqrt(-q*q*q) )
          
          eigen(1) = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0
          eigen(2) = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0
          eigen(3) = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0

          msk1 = merge(1.0_rp, 0.0_rp, eigen(2) .le. eigen(1) &
                       .and. eigen(1) .le. eigen(3) .or.  eigen(3) &
                       .le. eigen(1) .and. eigen(1) .le. eigen(2) )
          msk2 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(2) &
                       .and. eigen(2) .le. eigen(3) .or. eigen(3) &
                       .le. eigen(2) .and. eigen(2) .le. eigen(1))
          msk3 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(3) &
                       .and. eigen(3) .le. eigen(2) .or. eigen(2) &
                       .le. eigen(3) .and. eigen(3) .le. eigen(1))

          l2 = msk1 * eigen(1) + msk2 * eigen(2) + msk3 * eigen(3)

          lambda2(i,1,1,e) = l2/(cB(i,1,1,e)**2)
       end do
    end do
  end subroutine sx_lambda2_lx13

  subroutine sx_lambda2_lx12(lambda2, u, v, w, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, cB, n)
    integer, parameter :: lx = 12
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: lambda2
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u    
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: v
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: w
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cB
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: grad(lx*lx*lx,3,3)
    real(kind=rp) :: eigen(3), B, C, D, q, r, theta, l2
    real(kind=rp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=rp) :: a11, a22, a33, a12, a13, a23
    real(kind=rp) :: msk1, msk2, msk3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
    real(kind=rp) :: tmp1, tmp2, tmp3
    integer :: e, i, j, k, l
    
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do k = 1, lx
                tmp1 = tmp1 + dx(i,k) * u(k,j,1,e)
                tmp2 = tmp2 + dx(i,k) * v(k,j,1,e)
                tmp3 = tmp3 + dx(i,k) * w(k,j,1,e)
             end do
             ur(i,j,1) = tmp1
             vr(i,j,1) = tmp2
             wr(i,j,1) = tmp3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp1 = 0.0_rp
                tmp2 = 0.0_rp
                tmp3 = 0.0_rp
                do l = 1, lx
                   tmp1 = tmp1 + dy(j,l) * u(i,l,k,e)
                   tmp2 = tmp2 + dy(j,l) * v(i,l,k,e)
                   tmp3 = tmp3 + dy(j,l) * w(i,l,k,e)
                end do
                us(i,j,k) = tmp1
                vs(i,j,k) = tmp2
                ws(i,j,k) = tmp3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do l = 1, lx
                tmp1 = tmp1 + dz(k,l) * u(i,1,l,e)
                tmp2 = tmp2 + dz(k,l) * v(i,1,l,e)
                tmp3 = tmp3 + dz(k,l) * w(i,1,l,e)
             end do
             ut(i,1,k) = tmp1
             vt(i,1,k) = tmp2
             wt(i,1,k) = tmp3
          end do
       end do

       do i = 1, lx * lx * lx
          grad(1,1,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          grad(1,1,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          grad(1,1,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )

          grad(1,2,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * vr(i,1,1) &
                        + dsdx(i,1,1,e) * vs(i,1,1) &
                        + dtdx(i,1,1,e) * vt(i,1,1) )
          grad(1,2,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * vs(i,1,1) &
                        + drdy(i,1,1,e) * vr(i,1,1) &
                        + dtdy(i,1,1,e) * vt(i,1,1) )
          grad(1,2,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * vt(i,1,1) &
                        + drdz(i,1,1,e) * vr(i,1,1) &
                        + dsdz(i,1,1,e) * vs(i,1,1) )

          grad(1,3,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * wr(i,1,1) &
                        + dsdx(i,1,1,e) * ws(i,1,1) &
                        + dtdx(i,1,1,e) * wt(i,1,1) )
          grad(1,3,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * ws(i,1,1) &
                        + drdy(i,1,1,e) * wr(i,1,1) &
                        + dtdy(i,1,1,e) * wt(i,1,1) )
          grad(1,3,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * wt(i,1,1) &
                        + drdz(i,1,1,e) * wr(i,1,1) &
                        + dsdz(i,1,1,e) * ws(i,1,1) )          
       end do
       

       do i = 1, lx * lx * lx
          s11 = grad(i,1,1)
          s22 = grad(i,2,2)
          s33 = grad(i,3,3)

          
          s12 = 0.5*(grad(i,1,2) + grad(i,2,1))
          s13 = 0.5*(grad(i,1,3) + grad(i,3,1))
          s23 = 0.5*(grad(i,2,3) + grad(i,3,2))
          
          o12 = 0.5*(grad(i,1,2) - grad(i,2,1))
          o13 = 0.5*(grad(i,1,3) - grad(i,3,1))
          o23 = 0.5*(grad(i,2,3) - grad(i,3,2))

          a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
          a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23
          a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23

          a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
          a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
          a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


          B = -(a11 + a22 + a33)
          C = -(a12*a12 + a13*a13 + a23*a23 &
              - a11 * a22 - a11 * a33 - a22 * a33)
          D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 &
              - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33)


          q = (3.0 * C - B*B) / 9.0
          r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0
          theta = acos( r / sqrt(-q*q*q) )
          
          eigen(1) = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0
          eigen(2) = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0
          eigen(3) = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0

          msk1 = merge(1.0_rp, 0.0_rp, eigen(2) .le. eigen(1) &
                       .and. eigen(1) .le. eigen(3) .or.  eigen(3) &
                       .le. eigen(1) .and. eigen(1) .le. eigen(2) )
          msk2 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(2) &
                       .and. eigen(2) .le. eigen(3) .or. eigen(3) &
                       .le. eigen(2) .and. eigen(2) .le. eigen(1))
          msk3 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(3) &
                       .and. eigen(3) .le. eigen(2) .or. eigen(2) &
                       .le. eigen(3) .and. eigen(3) .le. eigen(1))

          l2 = msk1 * eigen(1) + msk2 * eigen(2) + msk3 * eigen(3)

          lambda2(i,1,1,e) = l2/(cB(i,1,1,e)**2)
       end do
    end do
  end subroutine sx_lambda2_lx12

  subroutine sx_lambda2_lx11(lambda2, u, v, w, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, cB, n)
    integer, parameter :: lx = 11
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: lambda2
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u    
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: v
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: w
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cB
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: grad(lx*lx*lx,3,3)
    real(kind=rp) :: eigen(3), B, C, D, q, r, theta, l2
    real(kind=rp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=rp) :: a11, a22, a33, a12, a13, a23
    real(kind=rp) :: msk1, msk2, msk3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
    real(kind=rp) :: tmp1, tmp2, tmp3
    integer :: e, i, j, k, l
    
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do k = 1, lx
                tmp1 = tmp1 + dx(i,k) * u(k,j,1,e)
                tmp2 = tmp2 + dx(i,k) * v(k,j,1,e)
                tmp3 = tmp3 + dx(i,k) * w(k,j,1,e)
             end do
             ur(i,j,1) = tmp1
             vr(i,j,1) = tmp2
             wr(i,j,1) = tmp3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp1 = 0.0_rp
                tmp2 = 0.0_rp
                tmp3 = 0.0_rp
                do l = 1, lx
                   tmp1 = tmp1 + dy(j,l) * u(i,l,k,e)
                   tmp2 = tmp2 + dy(j,l) * v(i,l,k,e)
                   tmp3 = tmp3 + dy(j,l) * w(i,l,k,e)
                end do
                us(i,j,k) = tmp1
                vs(i,j,k) = tmp2
                ws(i,j,k) = tmp3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do l = 1, lx
                tmp1 = tmp1 + dz(k,l) * u(i,1,l,e)
                tmp2 = tmp2 + dz(k,l) * v(i,1,l,e)
                tmp3 = tmp3 + dz(k,l) * w(i,1,l,e)
             end do
             ut(i,1,k) = tmp1
             vt(i,1,k) = tmp2
             wt(i,1,k) = tmp3
          end do
       end do

       do i = 1, lx * lx * lx
          grad(1,1,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          grad(1,1,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          grad(1,1,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )

          grad(1,2,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * vr(i,1,1) &
                        + dsdx(i,1,1,e) * vs(i,1,1) &
                        + dtdx(i,1,1,e) * vt(i,1,1) )
          grad(1,2,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * vs(i,1,1) &
                        + drdy(i,1,1,e) * vr(i,1,1) &
                        + dtdy(i,1,1,e) * vt(i,1,1) )
          grad(1,2,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * vt(i,1,1) &
                        + drdz(i,1,1,e) * vr(i,1,1) &
                        + dsdz(i,1,1,e) * vs(i,1,1) )

          grad(1,3,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * wr(i,1,1) &
                        + dsdx(i,1,1,e) * ws(i,1,1) &
                        + dtdx(i,1,1,e) * wt(i,1,1) )
          grad(1,3,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * ws(i,1,1) &
                        + drdy(i,1,1,e) * wr(i,1,1) &
                        + dtdy(i,1,1,e) * wt(i,1,1) )
          grad(1,3,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * wt(i,1,1) &
                        + drdz(i,1,1,e) * wr(i,1,1) &
                        + dsdz(i,1,1,e) * ws(i,1,1) )          
       end do
       

       do i = 1, lx * lx * lx
          s11 = grad(i,1,1)
          s22 = grad(i,2,2)
          s33 = grad(i,3,3)

          
          s12 = 0.5*(grad(i,1,2) + grad(i,2,1))
          s13 = 0.5*(grad(i,1,3) + grad(i,3,1))
          s23 = 0.5*(grad(i,2,3) + grad(i,3,2))
          
          o12 = 0.5*(grad(i,1,2) - grad(i,2,1))
          o13 = 0.5*(grad(i,1,3) - grad(i,3,1))
          o23 = 0.5*(grad(i,2,3) - grad(i,3,2))

          a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
          a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23
          a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23

          a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
          a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
          a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


          B = -(a11 + a22 + a33)
          C = -(a12*a12 + a13*a13 + a23*a23 &
              - a11 * a22 - a11 * a33 - a22 * a33)
          D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 &
              - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33)


          q = (3.0 * C - B*B) / 9.0
          r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0
          theta = acos( r / sqrt(-q*q*q) )
          
          eigen(1) = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0
          eigen(2) = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0
          eigen(3) = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0

          msk1 = merge(1.0_rp, 0.0_rp, eigen(2) .le. eigen(1) &
                       .and. eigen(1) .le. eigen(3) .or.  eigen(3) &
                       .le. eigen(1) .and. eigen(1) .le. eigen(2) )
          msk2 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(2) &
                       .and. eigen(2) .le. eigen(3) .or. eigen(3) &
                       .le. eigen(2) .and. eigen(2) .le. eigen(1))
          msk3 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(3) &
                       .and. eigen(3) .le. eigen(2) .or. eigen(2) &
                       .le. eigen(3) .and. eigen(3) .le. eigen(1))

          l2 = msk1 * eigen(1) + msk2 * eigen(2) + msk3 * eigen(3)

          lambda2(i,1,1,e) = l2/(cB(i,1,1,e)**2)
       end do
    end do
  end subroutine sx_lambda2_lx11

  subroutine sx_lambda2_lx10(lambda2, u, v, w, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, cB, n)
    integer, parameter :: lx = 10
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: lambda2
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u    
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: v
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: w
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cB
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: grad(lx*lx*lx,3,3)
    real(kind=rp) :: eigen(3), B, C, D, q, r, theta, l2
    real(kind=rp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=rp) :: a11, a22, a33, a12, a13, a23
    real(kind=rp) :: msk1, msk2, msk3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
    real(kind=rp) :: tmp1, tmp2, tmp3
    integer :: e, i, j, k, l
    
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do k = 1, lx
                tmp1 = tmp1 + dx(i,k) * u(k,j,1,e)
                tmp2 = tmp2 + dx(i,k) * v(k,j,1,e)
                tmp3 = tmp3 + dx(i,k) * w(k,j,1,e)
             end do
             ur(i,j,1) = tmp1
             vr(i,j,1) = tmp2
             wr(i,j,1) = tmp3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp1 = 0.0_rp
                tmp2 = 0.0_rp
                tmp3 = 0.0_rp
                do l = 1, lx
                   tmp1 = tmp1 + dy(j,l) * u(i,l,k,e)
                   tmp2 = tmp2 + dy(j,l) * v(i,l,k,e)
                   tmp3 = tmp3 + dy(j,l) * w(i,l,k,e)
                end do
                us(i,j,k) = tmp1
                vs(i,j,k) = tmp2
                ws(i,j,k) = tmp3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do l = 1, lx
                tmp1 = tmp1 + dz(k,l) * u(i,1,l,e)
                tmp2 = tmp2 + dz(k,l) * v(i,1,l,e)
                tmp3 = tmp3 + dz(k,l) * w(i,1,l,e)
             end do
             ut(i,1,k) = tmp1
             vt(i,1,k) = tmp2
             wt(i,1,k) = tmp3
          end do
       end do

       do i = 1, lx * lx * lx
          grad(1,1,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          grad(1,1,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          grad(1,1,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )

          grad(1,2,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * vr(i,1,1) &
                        + dsdx(i,1,1,e) * vs(i,1,1) &
                        + dtdx(i,1,1,e) * vt(i,1,1) )
          grad(1,2,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * vs(i,1,1) &
                        + drdy(i,1,1,e) * vr(i,1,1) &
                        + dtdy(i,1,1,e) * vt(i,1,1) )
          grad(1,2,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * vt(i,1,1) &
                        + drdz(i,1,1,e) * vr(i,1,1) &
                        + dsdz(i,1,1,e) * vs(i,1,1) )

          grad(1,3,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * wr(i,1,1) &
                        + dsdx(i,1,1,e) * ws(i,1,1) &
                        + dtdx(i,1,1,e) * wt(i,1,1) )
          grad(1,3,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * ws(i,1,1) &
                        + drdy(i,1,1,e) * wr(i,1,1) &
                        + dtdy(i,1,1,e) * wt(i,1,1) )
          grad(1,3,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * wt(i,1,1) &
                        + drdz(i,1,1,e) * wr(i,1,1) &
                        + dsdz(i,1,1,e) * ws(i,1,1) )          
       end do
       

       do i = 1, lx * lx * lx
          s11 = grad(i,1,1)
          s22 = grad(i,2,2)
          s33 = grad(i,3,3)

          
          s12 = 0.5*(grad(i,1,2) + grad(i,2,1))
          s13 = 0.5*(grad(i,1,3) + grad(i,3,1))
          s23 = 0.5*(grad(i,2,3) + grad(i,3,2))
          
          o12 = 0.5*(grad(i,1,2) - grad(i,2,1))
          o13 = 0.5*(grad(i,1,3) - grad(i,3,1))
          o23 = 0.5*(grad(i,2,3) - grad(i,3,2))

          a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
          a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23
          a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23

          a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
          a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
          a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


          B = -(a11 + a22 + a33)
          C = -(a12*a12 + a13*a13 + a23*a23 &
              - a11 * a22 - a11 * a33 - a22 * a33)
          D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 &
              - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33)


          q = (3.0 * C - B*B) / 9.0
          r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0
          theta = acos( r / sqrt(-q*q*q) )
          
          eigen(1) = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0
          eigen(2) = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0
          eigen(3) = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0

          msk1 = merge(1.0_rp, 0.0_rp, eigen(2) .le. eigen(1) &
                       .and. eigen(1) .le. eigen(3) .or.  eigen(3) &
                       .le. eigen(1) .and. eigen(1) .le. eigen(2) )
          msk2 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(2) &
                       .and. eigen(2) .le. eigen(3) .or. eigen(3) &
                       .le. eigen(2) .and. eigen(2) .le. eigen(1))
          msk3 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(3) &
                       .and. eigen(3) .le. eigen(2) .or. eigen(2) &
                       .le. eigen(3) .and. eigen(3) .le. eigen(1))

          l2 = msk1 * eigen(1) + msk2 * eigen(2) + msk3 * eigen(3)

          lambda2(i,1,1,e) = l2/(cB(i,1,1,e)**2)
       end do
    end do
  end subroutine sx_lambda2_lx10

  subroutine sx_lambda2_lx9(lambda2, u, v, w, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, cB, n)
    integer, parameter :: lx = 9
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: lambda2
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u    
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: v
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: w
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cB
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: grad(lx*lx*lx,3,3)
    real(kind=rp) :: eigen(3), B, C, D, q, r, theta, l2
    real(kind=rp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=rp) :: a11, a22, a33, a12, a13, a23
    real(kind=rp) :: msk1, msk2, msk3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
    real(kind=rp) :: tmp1, tmp2, tmp3
    integer :: e, i, j, k, l
    
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do k = 1, lx
                tmp1 = tmp1 + dx(i,k) * u(k,j,1,e)
                tmp2 = tmp2 + dx(i,k) * v(k,j,1,e)
                tmp3 = tmp3 + dx(i,k) * w(k,j,1,e)
             end do
             ur(i,j,1) = tmp1
             vr(i,j,1) = tmp2
             wr(i,j,1) = tmp3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp1 = 0.0_rp
                tmp2 = 0.0_rp
                tmp3 = 0.0_rp
                do l = 1, lx
                   tmp1 = tmp1 + dy(j,l) * u(i,l,k,e)
                   tmp2 = tmp2 + dy(j,l) * v(i,l,k,e)
                   tmp3 = tmp3 + dy(j,l) * w(i,l,k,e)
                end do
                us(i,j,k) = tmp1
                vs(i,j,k) = tmp2
                ws(i,j,k) = tmp3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do l = 1, lx
                tmp1 = tmp1 + dz(k,l) * u(i,1,l,e)
                tmp2 = tmp2 + dz(k,l) * v(i,1,l,e)
                tmp3 = tmp3 + dz(k,l) * w(i,1,l,e)
             end do
             ut(i,1,k) = tmp1
             vt(i,1,k) = tmp2
             wt(i,1,k) = tmp3
          end do
       end do

       do i = 1, lx * lx * lx
          grad(1,1,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          grad(1,1,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          grad(1,1,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )

          grad(1,2,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * vr(i,1,1) &
                        + dsdx(i,1,1,e) * vs(i,1,1) &
                        + dtdx(i,1,1,e) * vt(i,1,1) )
          grad(1,2,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * vs(i,1,1) &
                        + drdy(i,1,1,e) * vr(i,1,1) &
                        + dtdy(i,1,1,e) * vt(i,1,1) )
          grad(1,2,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * vt(i,1,1) &
                        + drdz(i,1,1,e) * vr(i,1,1) &
                        + dsdz(i,1,1,e) * vs(i,1,1) )

          grad(1,3,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * wr(i,1,1) &
                        + dsdx(i,1,1,e) * ws(i,1,1) &
                        + dtdx(i,1,1,e) * wt(i,1,1) )
          grad(1,3,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * ws(i,1,1) &
                        + drdy(i,1,1,e) * wr(i,1,1) &
                        + dtdy(i,1,1,e) * wt(i,1,1) )
          grad(1,3,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * wt(i,1,1) &
                        + drdz(i,1,1,e) * wr(i,1,1) &
                        + dsdz(i,1,1,e) * ws(i,1,1) )          
       end do
       

       do i = 1, lx * lx * lx
          s11 = grad(i,1,1)
          s22 = grad(i,2,2)
          s33 = grad(i,3,3)

          
          s12 = 0.5*(grad(i,1,2) + grad(i,2,1))
          s13 = 0.5*(grad(i,1,3) + grad(i,3,1))
          s23 = 0.5*(grad(i,2,3) + grad(i,3,2))
          
          o12 = 0.5*(grad(i,1,2) - grad(i,2,1))
          o13 = 0.5*(grad(i,1,3) - grad(i,3,1))
          o23 = 0.5*(grad(i,2,3) - grad(i,3,2))

          a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
          a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23
          a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23

          a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
          a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
          a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


          B = -(a11 + a22 + a33)
          C = -(a12*a12 + a13*a13 + a23*a23 &
              - a11 * a22 - a11 * a33 - a22 * a33)
          D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 &
              - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33)


          q = (3.0 * C - B*B) / 9.0
          r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0
          theta = acos( r / sqrt(-q*q*q) )
          
          eigen(1) = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0
          eigen(2) = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0
          eigen(3) = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0

          msk1 = merge(1.0_rp, 0.0_rp, eigen(2) .le. eigen(1) &
                       .and. eigen(1) .le. eigen(3) .or.  eigen(3) &
                       .le. eigen(1) .and. eigen(1) .le. eigen(2) )
          msk2 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(2) &
                       .and. eigen(2) .le. eigen(3) .or. eigen(3) &
                       .le. eigen(2) .and. eigen(2) .le. eigen(1))
          msk3 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(3) &
                       .and. eigen(3) .le. eigen(2) .or. eigen(2) &
                       .le. eigen(3) .and. eigen(3) .le. eigen(1))

          l2 = msk1 * eigen(1) + msk2 * eigen(2) + msk3 * eigen(3)

          lambda2(i,1,1,e) = l2/(cB(i,1,1,e)**2)
       end do
    end do
  end subroutine sx_lambda2_lx9

  subroutine sx_lambda2_lx8(lambda2, u, v, w, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, cB, n)
    integer, parameter :: lx = 8
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: lambda2
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u    
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: v
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: w
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cB
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: grad(lx*lx*lx,3,3)
    real(kind=rp) :: eigen(3), B, C, D, q, r, theta, l2
    real(kind=rp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=rp) :: a11, a22, a33, a12, a13, a23
    real(kind=rp) :: msk1, msk2, msk3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
    real(kind=rp) :: tmp1, tmp2, tmp3
    integer :: e, i, j, k, l
    
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do k = 1, lx
                tmp1 = tmp1 + dx(i,k) * u(k,j,1,e)
                tmp2 = tmp2 + dx(i,k) * v(k,j,1,e)
                tmp3 = tmp3 + dx(i,k) * w(k,j,1,e)
             end do
             ur(i,j,1) = tmp1
             vr(i,j,1) = tmp2
             wr(i,j,1) = tmp3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp1 = 0.0_rp
                tmp2 = 0.0_rp
                tmp3 = 0.0_rp
                do l = 1, lx
                   tmp1 = tmp1 + dy(j,l) * u(i,l,k,e)
                   tmp2 = tmp2 + dy(j,l) * v(i,l,k,e)
                   tmp3 = tmp3 + dy(j,l) * w(i,l,k,e)
                end do
                us(i,j,k) = tmp1
                vs(i,j,k) = tmp2
                ws(i,j,k) = tmp3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do l = 1, lx
                tmp1 = tmp1 + dz(k,l) * u(i,1,l,e)
                tmp2 = tmp2 + dz(k,l) * v(i,1,l,e)
                tmp3 = tmp3 + dz(k,l) * w(i,1,l,e)
             end do
             ut(i,1,k) = tmp1
             vt(i,1,k) = tmp2
             wt(i,1,k) = tmp3
          end do
       end do

       do i = 1, lx * lx * lx
          grad(1,1,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          grad(1,1,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          grad(1,1,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )

          grad(1,2,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * vr(i,1,1) &
                        + dsdx(i,1,1,e) * vs(i,1,1) &
                        + dtdx(i,1,1,e) * vt(i,1,1) )
          grad(1,2,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * vs(i,1,1) &
                        + drdy(i,1,1,e) * vr(i,1,1) &
                        + dtdy(i,1,1,e) * vt(i,1,1) )
          grad(1,2,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * vt(i,1,1) &
                        + drdz(i,1,1,e) * vr(i,1,1) &
                        + dsdz(i,1,1,e) * vs(i,1,1) )

          grad(1,3,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * wr(i,1,1) &
                        + dsdx(i,1,1,e) * ws(i,1,1) &
                        + dtdx(i,1,1,e) * wt(i,1,1) )
          grad(1,3,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * ws(i,1,1) &
                        + drdy(i,1,1,e) * wr(i,1,1) &
                        + dtdy(i,1,1,e) * wt(i,1,1) )
          grad(1,3,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * wt(i,1,1) &
                        + drdz(i,1,1,e) * wr(i,1,1) &
                        + dsdz(i,1,1,e) * ws(i,1,1) )          
       end do
       

       do i = 1, lx * lx * lx
          s11 = grad(i,1,1)
          s22 = grad(i,2,2)
          s33 = grad(i,3,3)

          
          s12 = 0.5*(grad(i,1,2) + grad(i,2,1))
          s13 = 0.5*(grad(i,1,3) + grad(i,3,1))
          s23 = 0.5*(grad(i,2,3) + grad(i,3,2))
          
          o12 = 0.5*(grad(i,1,2) - grad(i,2,1))
          o13 = 0.5*(grad(i,1,3) - grad(i,3,1))
          o23 = 0.5*(grad(i,2,3) - grad(i,3,2))

          a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
          a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23
          a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23

          a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
          a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
          a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


          B = -(a11 + a22 + a33)
          C = -(a12*a12 + a13*a13 + a23*a23 &
              - a11 * a22 - a11 * a33 - a22 * a33)
          D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 &
              - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33)


          q = (3.0 * C - B*B) / 9.0
          r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0
          theta = acos( r / sqrt(-q*q*q) )
          
          eigen(1) = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0
          eigen(2) = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0
          eigen(3) = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0

          msk1 = merge(1.0_rp, 0.0_rp, eigen(2) .le. eigen(1) &
                       .and. eigen(1) .le. eigen(3) .or.  eigen(3) &
                       .le. eigen(1) .and. eigen(1) .le. eigen(2) )
          msk2 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(2) &
                       .and. eigen(2) .le. eigen(3) .or. eigen(3) &
                       .le. eigen(2) .and. eigen(2) .le. eigen(1))
          msk3 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(3) &
                       .and. eigen(3) .le. eigen(2) .or. eigen(2) &
                       .le. eigen(3) .and. eigen(3) .le. eigen(1))

          l2 = msk1 * eigen(1) + msk2 * eigen(2) + msk3 * eigen(3)

          lambda2(i,1,1,e) = l2/(cB(i,1,1,e)**2)
       end do
    end do
  end subroutine sx_lambda2_lx8

  subroutine sx_lambda2_lx7(lambda2, u, v, w, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, cB, n)
    integer, parameter :: lx = 7
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: lambda2
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u    
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: v
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: w
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cB
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: grad(lx*lx*lx,3,3)
    real(kind=rp) :: eigen(3), B, C, D, q, r, theta, l2
    real(kind=rp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=rp) :: a11, a22, a33, a12, a13, a23
    real(kind=rp) :: msk1, msk2, msk3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
    real(kind=rp) :: tmp1, tmp2, tmp3
    integer :: e, i, j, k, l
    
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do k = 1, lx
                tmp1 = tmp1 + dx(i,k) * u(k,j,1,e)
                tmp2 = tmp2 + dx(i,k) * v(k,j,1,e)
                tmp3 = tmp3 + dx(i,k) * w(k,j,1,e)
             end do
             ur(i,j,1) = tmp1
             vr(i,j,1) = tmp2
             wr(i,j,1) = tmp3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp1 = 0.0_rp
                tmp2 = 0.0_rp
                tmp3 = 0.0_rp
                do l = 1, lx
                   tmp1 = tmp1 + dy(j,l) * u(i,l,k,e)
                   tmp2 = tmp2 + dy(j,l) * v(i,l,k,e)
                   tmp3 = tmp3 + dy(j,l) * w(i,l,k,e)
                end do
                us(i,j,k) = tmp1
                vs(i,j,k) = tmp2
                ws(i,j,k) = tmp3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do l = 1, lx
                tmp1 = tmp1 + dz(k,l) * u(i,1,l,e)
                tmp2 = tmp2 + dz(k,l) * v(i,1,l,e)
                tmp3 = tmp3 + dz(k,l) * w(i,1,l,e)
             end do
             ut(i,1,k) = tmp1
             vt(i,1,k) = tmp2
             wt(i,1,k) = tmp3
          end do
       end do

       do i = 1, lx * lx * lx
          grad(1,1,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          grad(1,1,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          grad(1,1,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )

          grad(1,2,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * vr(i,1,1) &
                        + dsdx(i,1,1,e) * vs(i,1,1) &
                        + dtdx(i,1,1,e) * vt(i,1,1) )
          grad(1,2,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * vs(i,1,1) &
                        + drdy(i,1,1,e) * vr(i,1,1) &
                        + dtdy(i,1,1,e) * vt(i,1,1) )
          grad(1,2,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * vt(i,1,1) &
                        + drdz(i,1,1,e) * vr(i,1,1) &
                        + dsdz(i,1,1,e) * vs(i,1,1) )

          grad(1,3,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * wr(i,1,1) &
                        + dsdx(i,1,1,e) * ws(i,1,1) &
                        + dtdx(i,1,1,e) * wt(i,1,1) )
          grad(1,3,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * ws(i,1,1) &
                        + drdy(i,1,1,e) * wr(i,1,1) &
                        + dtdy(i,1,1,e) * wt(i,1,1) )
          grad(1,3,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * wt(i,1,1) &
                        + drdz(i,1,1,e) * wr(i,1,1) &
                        + dsdz(i,1,1,e) * ws(i,1,1) )          
       end do
       

       do i = 1, lx * lx * lx
          s11 = grad(i,1,1)
          s22 = grad(i,2,2)
          s33 = grad(i,3,3)

          
          s12 = 0.5*(grad(i,1,2) + grad(i,2,1))
          s13 = 0.5*(grad(i,1,3) + grad(i,3,1))
          s23 = 0.5*(grad(i,2,3) + grad(i,3,2))
          
          o12 = 0.5*(grad(i,1,2) - grad(i,2,1))
          o13 = 0.5*(grad(i,1,3) - grad(i,3,1))
          o23 = 0.5*(grad(i,2,3) - grad(i,3,2))

          a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
          a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23
          a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23

          a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
          a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
          a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


          B = -(a11 + a22 + a33)
          C = -(a12*a12 + a13*a13 + a23*a23 &
              - a11 * a22 - a11 * a33 - a22 * a33)
          D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 &
              - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33)


          q = (3.0 * C - B*B) / 9.0
          r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0
          theta = acos( r / sqrt(-q*q*q) )
          
          eigen(1) = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0
          eigen(2) = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0
          eigen(3) = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0

          msk1 = merge(1.0_rp, 0.0_rp, eigen(2) .le. eigen(1) &
                       .and. eigen(1) .le. eigen(3) .or.  eigen(3) &
                       .le. eigen(1) .and. eigen(1) .le. eigen(2) )
          msk2 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(2) &
                       .and. eigen(2) .le. eigen(3) .or. eigen(3) &
                       .le. eigen(2) .and. eigen(2) .le. eigen(1))
          msk3 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(3) &
                       .and. eigen(3) .le. eigen(2) .or. eigen(2) &
                       .le. eigen(3) .and. eigen(3) .le. eigen(1))

          l2 = msk1 * eigen(1) + msk2 * eigen(2) + msk3 * eigen(3)

          lambda2(i,1,1,e) = l2/(cB(i,1,1,e)**2)
       end do
    end do
  end subroutine sx_lambda2_lx7

  subroutine sx_lambda2_lx6(lambda2, u, v, w, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, cB, n)
    integer, parameter :: lx = 6
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: lambda2
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u    
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: v
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: w
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cB
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: grad(lx*lx*lx,3,3)
    real(kind=rp) :: eigen(3), B, C, D, q, r, theta, l2
    real(kind=rp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=rp) :: a11, a22, a33, a12, a13, a23
    real(kind=rp) :: msk1, msk2, msk3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
    real(kind=rp) :: tmp1, tmp2, tmp3
    integer :: e, i, j, k, l
    
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do k = 1, lx
                tmp1 = tmp1 + dx(i,k) * u(k,j,1,e)
                tmp2 = tmp2 + dx(i,k) * v(k,j,1,e)
                tmp3 = tmp3 + dx(i,k) * w(k,j,1,e)
             end do
             ur(i,j,1) = tmp1
             vr(i,j,1) = tmp2
             wr(i,j,1) = tmp3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp1 = 0.0_rp
                tmp2 = 0.0_rp
                tmp3 = 0.0_rp
                do l = 1, lx
                   tmp1 = tmp1 + dy(j,l) * u(i,l,k,e)
                   tmp2 = tmp2 + dy(j,l) * v(i,l,k,e)
                   tmp3 = tmp3 + dy(j,l) * w(i,l,k,e)
                end do
                us(i,j,k) = tmp1
                vs(i,j,k) = tmp2
                ws(i,j,k) = tmp3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do l = 1, lx
                tmp1 = tmp1 + dz(k,l) * u(i,1,l,e)
                tmp2 = tmp2 + dz(k,l) * v(i,1,l,e)
                tmp3 = tmp3 + dz(k,l) * w(i,1,l,e)
             end do
             ut(i,1,k) = tmp1
             vt(i,1,k) = tmp2
             wt(i,1,k) = tmp3
          end do
       end do

       do i = 1, lx * lx * lx
          grad(1,1,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          grad(1,1,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          grad(1,1,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )

          grad(1,2,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * vr(i,1,1) &
                        + dsdx(i,1,1,e) * vs(i,1,1) &
                        + dtdx(i,1,1,e) * vt(i,1,1) )
          grad(1,2,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * vs(i,1,1) &
                        + drdy(i,1,1,e) * vr(i,1,1) &
                        + dtdy(i,1,1,e) * vt(i,1,1) )
          grad(1,2,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * vt(i,1,1) &
                        + drdz(i,1,1,e) * vr(i,1,1) &
                        + dsdz(i,1,1,e) * vs(i,1,1) )

          grad(1,3,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * wr(i,1,1) &
                        + dsdx(i,1,1,e) * ws(i,1,1) &
                        + dtdx(i,1,1,e) * wt(i,1,1) )
          grad(1,3,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * ws(i,1,1) &
                        + drdy(i,1,1,e) * wr(i,1,1) &
                        + dtdy(i,1,1,e) * wt(i,1,1) )
          grad(1,3,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * wt(i,1,1) &
                        + drdz(i,1,1,e) * wr(i,1,1) &
                        + dsdz(i,1,1,e) * ws(i,1,1) )          
       end do
       

       do i = 1, lx * lx * lx
          s11 = grad(i,1,1)
          s22 = grad(i,2,2)
          s33 = grad(i,3,3)

          
          s12 = 0.5*(grad(i,1,2) + grad(i,2,1))
          s13 = 0.5*(grad(i,1,3) + grad(i,3,1))
          s23 = 0.5*(grad(i,2,3) + grad(i,3,2))
          
          o12 = 0.5*(grad(i,1,2) - grad(i,2,1))
          o13 = 0.5*(grad(i,1,3) - grad(i,3,1))
          o23 = 0.5*(grad(i,2,3) - grad(i,3,2))

          a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
          a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23
          a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23

          a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
          a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
          a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


          B = -(a11 + a22 + a33)
          C = -(a12*a12 + a13*a13 + a23*a23 &
              - a11 * a22 - a11 * a33 - a22 * a33)
          D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 &
              - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33)


          q = (3.0 * C - B*B) / 9.0
          r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0
          theta = acos( r / sqrt(-q*q*q) )
          
          eigen(1) = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0
          eigen(2) = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0
          eigen(3) = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0

          msk1 = merge(1.0_rp, 0.0_rp, eigen(2) .le. eigen(1) &
                       .and. eigen(1) .le. eigen(3) .or.  eigen(3) &
                       .le. eigen(1) .and. eigen(1) .le. eigen(2) )
          msk2 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(2) &
                       .and. eigen(2) .le. eigen(3) .or. eigen(3) &
                       .le. eigen(2) .and. eigen(2) .le. eigen(1))
          msk3 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(3) &
                       .and. eigen(3) .le. eigen(2) .or. eigen(2) &
                       .le. eigen(3) .and. eigen(3) .le. eigen(1))

          l2 = msk1 * eigen(1) + msk2 * eigen(2) + msk3 * eigen(3)

          lambda2(i,1,1,e) = l2/(cB(i,1,1,e)**2)
       end do
    end do
  end subroutine sx_lambda2_lx6

  subroutine sx_lambda2_lx5(lambda2, u, v, w, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, cB, n)
    integer, parameter :: lx = 5
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: lambda2
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u    
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: v
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: w
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cB
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: grad(lx*lx*lx,3,3)
    real(kind=rp) :: eigen(3), B, C, D, q, r, theta, l2
    real(kind=rp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=rp) :: a11, a22, a33, a12, a13, a23
    real(kind=rp) :: msk1, msk2, msk3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
    real(kind=rp) :: tmp1, tmp2, tmp3
    integer :: e, i, j, k, l
    
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do k = 1, lx
                tmp1 = tmp1 + dx(i,k) * u(k,j,1,e)
                tmp2 = tmp2 + dx(i,k) * v(k,j,1,e)
                tmp3 = tmp3 + dx(i,k) * w(k,j,1,e)
             end do
             ur(i,j,1) = tmp1
             vr(i,j,1) = tmp2
             wr(i,j,1) = tmp3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp1 = 0.0_rp
                tmp2 = 0.0_rp
                tmp3 = 0.0_rp
                do l = 1, lx
                   tmp1 = tmp1 + dy(j,l) * u(i,l,k,e)
                   tmp2 = tmp2 + dy(j,l) * v(i,l,k,e)
                   tmp3 = tmp3 + dy(j,l) * w(i,l,k,e)
                end do
                us(i,j,k) = tmp1
                vs(i,j,k) = tmp2
                ws(i,j,k) = tmp3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do l = 1, lx
                tmp1 = tmp1 + dz(k,l) * u(i,1,l,e)
                tmp2 = tmp2 + dz(k,l) * v(i,1,l,e)
                tmp3 = tmp3 + dz(k,l) * w(i,1,l,e)
             end do
             ut(i,1,k) = tmp1
             vt(i,1,k) = tmp2
             wt(i,1,k) = tmp3
          end do
       end do

       do i = 1, lx * lx * lx
          grad(1,1,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          grad(1,1,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          grad(1,1,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )

          grad(1,2,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * vr(i,1,1) &
                        + dsdx(i,1,1,e) * vs(i,1,1) &
                        + dtdx(i,1,1,e) * vt(i,1,1) )
          grad(1,2,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * vs(i,1,1) &
                        + drdy(i,1,1,e) * vr(i,1,1) &
                        + dtdy(i,1,1,e) * vt(i,1,1) )
          grad(1,2,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * vt(i,1,1) &
                        + drdz(i,1,1,e) * vr(i,1,1) &
                        + dsdz(i,1,1,e) * vs(i,1,1) )

          grad(1,3,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * wr(i,1,1) &
                        + dsdx(i,1,1,e) * ws(i,1,1) &
                        + dtdx(i,1,1,e) * wt(i,1,1) )
          grad(1,3,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * ws(i,1,1) &
                        + drdy(i,1,1,e) * wr(i,1,1) &
                        + dtdy(i,1,1,e) * wt(i,1,1) )
          grad(1,3,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * wt(i,1,1) &
                        + drdz(i,1,1,e) * wr(i,1,1) &
                        + dsdz(i,1,1,e) * ws(i,1,1) )          
       end do
       

       do i = 1, lx * lx * lx
          s11 = grad(i,1,1)
          s22 = grad(i,2,2)
          s33 = grad(i,3,3)

          
          s12 = 0.5*(grad(i,1,2) + grad(i,2,1))
          s13 = 0.5*(grad(i,1,3) + grad(i,3,1))
          s23 = 0.5*(grad(i,2,3) + grad(i,3,2))
          
          o12 = 0.5*(grad(i,1,2) - grad(i,2,1))
          o13 = 0.5*(grad(i,1,3) - grad(i,3,1))
          o23 = 0.5*(grad(i,2,3) - grad(i,3,2))

          a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
          a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23
          a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23

          a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
          a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
          a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


          B = -(a11 + a22 + a33)
          C = -(a12*a12 + a13*a13 + a23*a23 &
              - a11 * a22 - a11 * a33 - a22 * a33)
          D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 &
              - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33)


          q = (3.0 * C - B*B) / 9.0
          r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0
          theta = acos( r / sqrt(-q*q*q) )
          
          eigen(1) = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0
          eigen(2) = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0
          eigen(3) = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0

          msk1 = merge(1.0_rp, 0.0_rp, eigen(2) .le. eigen(1) &
                       .and. eigen(1) .le. eigen(3) .or.  eigen(3) &
                       .le. eigen(1) .and. eigen(1) .le. eigen(2) )
          msk2 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(2) &
                       .and. eigen(2) .le. eigen(3) .or. eigen(3) &
                       .le. eigen(2) .and. eigen(2) .le. eigen(1))
          msk3 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(3) &
                       .and. eigen(3) .le. eigen(2) .or. eigen(2) &
                       .le. eigen(3) .and. eigen(3) .le. eigen(1))

          l2 = msk1 * eigen(1) + msk2 * eigen(2) + msk3 * eigen(3)

          lambda2(i,1,1,e) = l2/(cB(i,1,1,e)**2)
       end do
    end do
  end subroutine sx_lambda2_lx5

  subroutine sx_lambda2_lx4(lambda2, u, v, w, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, cB, n)
    integer, parameter :: lx = 4
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: lambda2
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u    
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: v
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: w
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cB
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: grad(lx*lx*lx,3,3)
    real(kind=rp) :: eigen(3), B, C, D, q, r, theta, l2
    real(kind=rp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=rp) :: a11, a22, a33, a12, a13, a23
    real(kind=rp) :: msk1, msk2, msk3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
    real(kind=rp) :: tmp1, tmp2, tmp3
    integer :: e, i, j, k, l
    
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do k = 1, lx
                tmp1 = tmp1 + dx(i,k) * u(k,j,1,e)
                tmp2 = tmp2 + dx(i,k) * v(k,j,1,e)
                tmp3 = tmp3 + dx(i,k) * w(k,j,1,e)
             end do
             ur(i,j,1) = tmp1
             vr(i,j,1) = tmp2
             wr(i,j,1) = tmp3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp1 = 0.0_rp
                tmp2 = 0.0_rp
                tmp3 = 0.0_rp
                do l = 1, lx
                   tmp1 = tmp1 + dy(j,l) * u(i,l,k,e)
                   tmp2 = tmp2 + dy(j,l) * v(i,l,k,e)
                   tmp3 = tmp3 + dy(j,l) * w(i,l,k,e)
                end do
                us(i,j,k) = tmp1
                vs(i,j,k) = tmp2
                ws(i,j,k) = tmp3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do l = 1, lx
                tmp1 = tmp1 + dz(k,l) * u(i,1,l,e)
                tmp2 = tmp2 + dz(k,l) * v(i,1,l,e)
                tmp3 = tmp3 + dz(k,l) * w(i,1,l,e)
             end do
             ut(i,1,k) = tmp1
             vt(i,1,k) = tmp2
             wt(i,1,k) = tmp3
          end do
       end do

       do i = 1, lx * lx * lx
          grad(1,1,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          grad(1,1,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          grad(1,1,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )

          grad(1,2,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * vr(i,1,1) &
                        + dsdx(i,1,1,e) * vs(i,1,1) &
                        + dtdx(i,1,1,e) * vt(i,1,1) )
          grad(1,2,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * vs(i,1,1) &
                        + drdy(i,1,1,e) * vr(i,1,1) &
                        + dtdy(i,1,1,e) * vt(i,1,1) )
          grad(1,2,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * vt(i,1,1) &
                        + drdz(i,1,1,e) * vr(i,1,1) &
                        + dsdz(i,1,1,e) * vs(i,1,1) )

          grad(1,3,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * wr(i,1,1) &
                        + dsdx(i,1,1,e) * ws(i,1,1) &
                        + dtdx(i,1,1,e) * wt(i,1,1) )
          grad(1,3,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * ws(i,1,1) &
                        + drdy(i,1,1,e) * wr(i,1,1) &
                        + dtdy(i,1,1,e) * wt(i,1,1) )
          grad(1,3,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * wt(i,1,1) &
                        + drdz(i,1,1,e) * wr(i,1,1) &
                        + dsdz(i,1,1,e) * ws(i,1,1) )          
       end do
       

       do i = 1, lx * lx * lx
          s11 = grad(i,1,1)
          s22 = grad(i,2,2)
          s33 = grad(i,3,3)

          
          s12 = 0.5*(grad(i,1,2) + grad(i,2,1))
          s13 = 0.5*(grad(i,1,3) + grad(i,3,1))
          s23 = 0.5*(grad(i,2,3) + grad(i,3,2))
          
          o12 = 0.5*(grad(i,1,2) - grad(i,2,1))
          o13 = 0.5*(grad(i,1,3) - grad(i,3,1))
          o23 = 0.5*(grad(i,2,3) - grad(i,3,2))

          a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
          a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23
          a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23

          a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
          a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
          a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


          B = -(a11 + a22 + a33)
          C = -(a12*a12 + a13*a13 + a23*a23 &
              - a11 * a22 - a11 * a33 - a22 * a33)
          D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 &
              - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33)


          q = (3.0 * C - B*B) / 9.0
          r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0
          theta = acos( r / sqrt(-q*q*q) )
          
          eigen(1) = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0
          eigen(2) = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0
          eigen(3) = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0

          msk1 = merge(1.0_rp, 0.0_rp, eigen(2) .le. eigen(1) &
                       .and. eigen(1) .le. eigen(3) .or.  eigen(3) &
                       .le. eigen(1) .and. eigen(1) .le. eigen(2) )
          msk2 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(2) &
                       .and. eigen(2) .le. eigen(3) .or. eigen(3) &
                       .le. eigen(2) .and. eigen(2) .le. eigen(1))
          msk3 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(3) &
                       .and. eigen(3) .le. eigen(2) .or. eigen(2) &
                       .le. eigen(3) .and. eigen(3) .le. eigen(1))

          l2 = msk1 * eigen(1) + msk2 * eigen(2) + msk3 * eigen(3)

          lambda2(i,1,1,e) = l2/(cB(i,1,1,e)**2)
       end do
    end do
  end subroutine sx_lambda2_lx4

  subroutine sx_lambda2_lx3(lambda2, u, v, w, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, cB, n)
    integer, parameter :: lx = 3
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: lambda2
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u    
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: v
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: w
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cB
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: grad(lx*lx*lx,3,3)
    real(kind=rp) :: eigen(3), B, C, D, q, r, theta, l2
    real(kind=rp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=rp) :: a11, a22, a33, a12, a13, a23
    real(kind=rp) :: msk1, msk2, msk3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
    real(kind=rp) :: tmp1, tmp2, tmp3
    integer :: e, i, j, k, l
    
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do k = 1, lx
                tmp1 = tmp1 + dx(i,k) * u(k,j,1,e)
                tmp2 = tmp2 + dx(i,k) * v(k,j,1,e)
                tmp3 = tmp3 + dx(i,k) * w(k,j,1,e)
             end do
             ur(i,j,1) = tmp1
             vr(i,j,1) = tmp2
             wr(i,j,1) = tmp3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp1 = 0.0_rp
                tmp2 = 0.0_rp
                tmp3 = 0.0_rp
                do l = 1, lx
                   tmp1 = tmp1 + dy(j,l) * u(i,l,k,e)
                   tmp2 = tmp2 + dy(j,l) * v(i,l,k,e)
                   tmp3 = tmp3 + dy(j,l) * w(i,l,k,e)
                end do
                us(i,j,k) = tmp1
                vs(i,j,k) = tmp2
                ws(i,j,k) = tmp3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do l = 1, lx
                tmp1 = tmp1 + dz(k,l) * u(i,1,l,e)
                tmp2 = tmp2 + dz(k,l) * v(i,1,l,e)
                tmp3 = tmp3 + dz(k,l) * w(i,1,l,e)
             end do
             ut(i,1,k) = tmp1
             vt(i,1,k) = tmp2
             wt(i,1,k) = tmp3
          end do
       end do

       do i = 1, lx * lx * lx
          grad(1,1,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          grad(1,1,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          grad(1,1,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )

          grad(1,2,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * vr(i,1,1) &
                        + dsdx(i,1,1,e) * vs(i,1,1) &
                        + dtdx(i,1,1,e) * vt(i,1,1) )
          grad(1,2,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * vs(i,1,1) &
                        + drdy(i,1,1,e) * vr(i,1,1) &
                        + dtdy(i,1,1,e) * vt(i,1,1) )
          grad(1,2,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * vt(i,1,1) &
                        + drdz(i,1,1,e) * vr(i,1,1) &
                        + dsdz(i,1,1,e) * vs(i,1,1) )

          grad(1,3,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * wr(i,1,1) &
                        + dsdx(i,1,1,e) * ws(i,1,1) &
                        + dtdx(i,1,1,e) * wt(i,1,1) )
          grad(1,3,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * ws(i,1,1) &
                        + drdy(i,1,1,e) * wr(i,1,1) &
                        + dtdy(i,1,1,e) * wt(i,1,1) )
          grad(1,3,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * wt(i,1,1) &
                        + drdz(i,1,1,e) * wr(i,1,1) &
                        + dsdz(i,1,1,e) * ws(i,1,1) )          
       end do
       

       do i = 1, lx * lx * lx
          s11 = grad(i,1,1)
          s22 = grad(i,2,2)
          s33 = grad(i,3,3)

          
          s12 = 0.5*(grad(i,1,2) + grad(i,2,1))
          s13 = 0.5*(grad(i,1,3) + grad(i,3,1))
          s23 = 0.5*(grad(i,2,3) + grad(i,3,2))
          
          o12 = 0.5*(grad(i,1,2) - grad(i,2,1))
          o13 = 0.5*(grad(i,1,3) - grad(i,3,1))
          o23 = 0.5*(grad(i,2,3) - grad(i,3,2))

          a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
          a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23
          a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23

          a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
          a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
          a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


          B = -(a11 + a22 + a33)
          C = -(a12*a12 + a13*a13 + a23*a23 &
              - a11 * a22 - a11 * a33 - a22 * a33)
          D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 &
              - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33)


          q = (3.0 * C - B*B) / 9.0
          r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0
          theta = acos( r / sqrt(-q*q*q) )
          
          eigen(1) = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0
          eigen(2) = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0
          eigen(3) = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0

          msk1 = merge(1.0_rp, 0.0_rp, eigen(2) .le. eigen(1) &
                       .and. eigen(1) .le. eigen(3) .or.  eigen(3) &
                       .le. eigen(1) .and. eigen(1) .le. eigen(2) )
          msk2 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(2) &
                       .and. eigen(2) .le. eigen(3) .or. eigen(3) &
                       .le. eigen(2) .and. eigen(2) .le. eigen(1))
          msk3 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(3) &
                       .and. eigen(3) .le. eigen(2) .or. eigen(2) &
                       .le. eigen(3) .and. eigen(3) .le. eigen(1))

          l2 = msk1 * eigen(1) + msk2 * eigen(2) + msk3 * eigen(3)

          lambda2(i,1,1,e) = l2/(cB(i,1,1,e)**2)
       end do
    end do
  end subroutine sx_lambda2_lx3

  subroutine sx_lambda2_lx2(lambda2, u, v, w, dx, dy, dz, &
       drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz, w3, cB, n)
    integer, parameter :: lx = 2
    integer, intent(in) :: n
    real(kind=rp), dimension(lx, lx, lx, n), intent(inout) :: lambda2
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: u    
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: v
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: w
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdx, dsdx, dtdx
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdy, dsdy, dtdy
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: drdz, dsdz, dtdz
    real(kind=rp), dimension(lx, lx, lx, n), intent(in) :: cB
    real(kind=rp), dimension(lx, lx, lx), intent(in) :: w3
    real(kind=rp) :: grad(lx*lx*lx,3,3)
    real(kind=rp) :: eigen(3), B, C, D, q, r, theta, l2
    real(kind=rp) :: s11, s22, s33, s12, s13, s23, o12, o13, o23
    real(kind=rp) :: a11, a22, a33, a12, a13, a23
    real(kind=rp) :: msk1, msk2, msk3
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
    real(kind=rp) :: tmp1, tmp2, tmp3
    integer :: e, i, j, k, l
    
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do k = 1, lx
                tmp1 = tmp1 + dx(i,k) * u(k,j,1,e)
                tmp2 = tmp2 + dx(i,k) * v(k,j,1,e)
                tmp3 = tmp3 + dx(i,k) * w(k,j,1,e)
             end do
             ur(i,j,1) = tmp1
             vr(i,j,1) = tmp2
             wr(i,j,1) = tmp3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp1 = 0.0_rp
                tmp2 = 0.0_rp
                tmp3 = 0.0_rp
                do l = 1, lx
                   tmp1 = tmp1 + dy(j,l) * u(i,l,k,e)
                   tmp2 = tmp2 + dy(j,l) * v(i,l,k,e)
                   tmp3 = tmp3 + dy(j,l) * w(i,l,k,e)
                end do
                us(i,j,k) = tmp1
                vs(i,j,k) = tmp2
                ws(i,j,k) = tmp3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp1 = 0.0_rp
             tmp2 = 0.0_rp
             tmp3 = 0.0_rp
             do l = 1, lx
                tmp1 = tmp1 + dz(k,l) * u(i,1,l,e)
                tmp2 = tmp2 + dz(k,l) * v(i,1,l,e)
                tmp3 = tmp3 + dz(k,l) * w(i,1,l,e)
             end do
             ut(i,1,k) = tmp1
             vt(i,1,k) = tmp2
             wt(i,1,k) = tmp3
          end do
       end do

       do i = 1, lx * lx * lx
          grad(1,1,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * ur(i,1,1) &
                        + dsdx(i,1,1,e) * us(i,1,1) &
                        + dtdx(i,1,1,e) * ut(i,1,1) )
          grad(1,1,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * us(i,1,1) &
                        + drdy(i,1,1,e) * ur(i,1,1) &
                        + dtdy(i,1,1,e) * ut(i,1,1) )
          grad(1,1,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * ut(i,1,1) &
                        + drdz(i,1,1,e) * ur(i,1,1) &
                        + dsdz(i,1,1,e) * us(i,1,1) )

          grad(1,2,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * vr(i,1,1) &
                        + dsdx(i,1,1,e) * vs(i,1,1) &
                        + dtdx(i,1,1,e) * vt(i,1,1) )
          grad(1,2,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * vs(i,1,1) &
                        + drdy(i,1,1,e) * vr(i,1,1) &
                        + dtdy(i,1,1,e) * vt(i,1,1) )
          grad(1,2,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * vt(i,1,1) &
                        + drdz(i,1,1,e) * vr(i,1,1) &
                        + dsdz(i,1,1,e) * vs(i,1,1) )

          grad(1,3,1) = w3(i,1,1) &
                      * ( drdx(i,1,1,e) * wr(i,1,1) &
                        + dsdx(i,1,1,e) * ws(i,1,1) &
                        + dtdx(i,1,1,e) * wt(i,1,1) )
          grad(1,3,2) = w3(i,1,1) &
                      * ( dsdy(i,1,1,e) * ws(i,1,1) &
                        + drdy(i,1,1,e) * wr(i,1,1) &
                        + dtdy(i,1,1,e) * wt(i,1,1) )
          grad(1,3,3) = w3(i,1,1) &
                      * ( dtdz(i,1,1,e) * wt(i,1,1) &
                        + drdz(i,1,1,e) * wr(i,1,1) &
                        + dsdz(i,1,1,e) * ws(i,1,1) )          
       end do
       

       do i = 1, lx * lx * lx
          s11 = grad(i,1,1)
          s22 = grad(i,2,2)
          s33 = grad(i,3,3)

          
          s12 = 0.5*(grad(i,1,2) + grad(i,2,1))
          s13 = 0.5*(grad(i,1,3) + grad(i,3,1))
          s23 = 0.5*(grad(i,2,3) + grad(i,3,2))
          
          o12 = 0.5*(grad(i,1,2) - grad(i,2,1))
          o13 = 0.5*(grad(i,1,3) - grad(i,3,1))
          o23 = 0.5*(grad(i,2,3) - grad(i,3,2))

          a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
          a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23
          a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23

          a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
          a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13
          a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23


          B = -(a11 + a22 + a33)
          C = -(a12*a12 + a13*a13 + a23*a23 &
              - a11 * a22 - a11 * a33 - a22 * a33)
          D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 &
              - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33)


          q = (3.0 * C - B*B) / 9.0
          r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0
          theta = acos( r / sqrt(-q*q*q) )
          
          eigen(1) = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0
          eigen(2) = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0
          eigen(3) = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0

          msk1 = merge(1.0_rp, 0.0_rp, eigen(2) .le. eigen(1) &
                       .and. eigen(1) .le. eigen(3) .or.  eigen(3) &
                       .le. eigen(1) .and. eigen(1) .le. eigen(2) )
          msk2 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(2) &
                       .and. eigen(2) .le. eigen(3) .or. eigen(3) &
                       .le. eigen(2) .and. eigen(2) .le. eigen(1))
          msk3 = merge(1.0_rp, 0.0_rp, eigen(1) .le. eigen(3) &
                       .and. eigen(3) .le. eigen(2) .or. eigen(2) &
                       .le. eigen(3) .and. eigen(3) .le. eigen(1))

          l2 = msk1 * eigen(1) + msk2 * eigen(2) + msk3 * eigen(3)

          lambda2(i,1,1,e) = l2/(cB(i,1,1,e)**2)
       end do
    end do
  end subroutine sx_lambda2_lx2
  
end submodule sx_lambda2


