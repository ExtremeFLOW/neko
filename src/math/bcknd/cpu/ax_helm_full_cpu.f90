! Copyright (c) 2023-2024, The Neko Authors
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
module ax_helm_full_cpu
  use ax_helm_full, only : ax_helm_full_t
  use num_types, only : rp
  use coefs, only : coef_t
  use space, only : space_t
  use mesh, only : mesh_t
  use math, only : addcol4
  use utils, only : neko_error
  implicit none
  private

  !> CPU matrix-vector product for a Helmholtz problem with full stress tensor.
  type, public, extends(ax_helm_full_t) :: ax_helm_full_cpu_t
   contains
     !> Compute the product.
     procedure, pass(this) :: compute_vector => ax_helm_full_compute_vector
  end type ax_helm_full_cpu_t

contains

  !> Compute \f$ Ax \f$ inside a Krylov method, taking 3
  !! components of a vector field in a coupled manner.
  !! @param au Result for the first component of the vector.
  !! @param av Result for the first component of the vector.
  !! @param aw Result for the first component of the vector.
  !! @param u The first component of the vector.
  !! @param v The second component of the vector.
  !! @param w The third component of the vector.
  !! @param coef Coefficients.
  !! @param msh Mesh.
  !! @param Xh Function space \f$ X_h \f$.
  subroutine ax_helm_full_compute_vector(this, au, av, aw, u, v, w, coef, msh,&
                                         Xh)
    class(ax_helm_full_cpu_t), intent(in) :: this
    type(mesh_t), intent(inout) :: msh
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: v(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: au(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: av(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: aw(Xh%lx, Xh%ly, Xh%lz, msh%nelv)

    select case (Xh%lx)
    case (14)
       call ax_helm_stress_lx14(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, &
            coef%drdx, coef%drdy, coef%drdz, coef%dsdx, coef%dsdy, coef%dsdz, &
            coef%dtdx, coef%dtdy, coef%dtdz, &
            coef%jacinv, Xh%w3, msh%nelv)
    case (13)
       call ax_helm_stress_lx13(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, &
            coef%drdx, coef%drdy, coef%drdz, coef%dsdx, coef%dsdy, coef%dsdz, &
            coef%dtdx, coef%dtdy, coef%dtdz, &
            coef%jacinv, Xh%w3, msh%nelv)
    case (12)
       call ax_helm_stress_lx12(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, &
            coef%drdx, coef%drdy, coef%drdz, coef%dsdx, coef%dsdy, coef%dsdz, &
            coef%dtdx, coef%dtdy, coef%dtdz, &
            coef%jacinv, Xh%w3, msh%nelv)
    case (11)
       call ax_helm_stress_lx11(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, &
            coef%drdx, coef%drdy, coef%drdz, coef%dsdx, coef%dsdy, coef%dsdz, &
            coef%dtdx, coef%dtdy, coef%dtdz, &
            coef%jacinv, Xh%w3, msh%nelv)
    case (10)
       call ax_helm_stress_lx10(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, &
            coef%drdx, coef%drdy, coef%drdz, coef%dsdx, coef%dsdy, coef%dsdz, &
            coef%dtdx, coef%dtdy, coef%dtdz, &
            coef%jacinv, Xh%w3, msh%nelv)
    case (9)
       call ax_helm_stress_lx9(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, &
            coef%drdx, coef%drdy, coef%drdz, coef%dsdx, coef%dsdy, coef%dsdz, &
            coef%dtdx, coef%dtdy, coef%dtdz, &
            coef%jacinv, Xh%w3, msh%nelv)
    case (8)
       call ax_helm_stress_lx8(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, &
            coef%drdx, coef%drdy, coef%drdz, coef%dsdx, coef%dsdy, coef%dsdz, &
            coef%dtdx, coef%dtdy, coef%dtdz, &
            coef%jacinv, Xh%w3, msh%nelv)
    case (7)
       call ax_helm_stress_lx7(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, &
            coef%drdx, coef%drdy, coef%drdz, coef%dsdx, coef%dsdy, coef%dsdz, &
            coef%dtdx, coef%dtdy, coef%dtdz, &
            coef%jacinv, Xh%w3, msh%nelv)
    case (6)
       call ax_helm_stress_lx6(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, &
            coef%drdx, coef%drdy, coef%drdz, coef%dsdx, coef%dsdy, coef%dsdz, &
            coef%dtdx, coef%dtdy, coef%dtdz, &
            coef%jacinv, Xh%w3, msh%nelv)
    case (5)
       call ax_helm_stress_lx5(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, &
            coef%drdx, coef%drdy, coef%drdz, coef%dsdx, coef%dsdy, coef%dsdz, &
            coef%dtdx, coef%dtdy, coef%dtdz, &
            coef%jacinv, Xh%w3, msh%nelv)
    case (4)
       call ax_helm_stress_lx4(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, &
            coef%drdx, coef%drdy, coef%drdz, coef%dsdx, coef%dsdy, coef%dsdz, &
            coef%dtdx, coef%dtdy, coef%dtdz, &
            coef%jacinv, Xh%w3, msh%nelv)
    case (3)
       call ax_helm_stress_lx3(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, &
            coef%drdx, coef%drdy, coef%drdz, coef%dsdx, coef%dsdy, coef%dsdz, &
            coef%dtdx, coef%dtdy, coef%dtdz, &
            coef%jacinv, Xh%w3, msh%nelv)
    case (2)
       call ax_helm_stress_lx2(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, &
            coef%drdx, coef%drdy, coef%drdz, coef%dsdx, coef%dsdy, coef%dsdz, &
            coef%dtdx, coef%dtdy, coef%dtdz, &
            coef%jacinv, Xh%w3, msh%nelv)
    case default
       call ax_helm_stress_lx(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, &
            coef%drdx, coef%drdy, coef%drdz, coef%dsdx, coef%dsdy, coef%dsdz, &
            coef%dtdx, coef%dtdy, coef%dtdz, &
            coef%jacinv, Xh%w3, msh%nelv, Xh%lx)
    end select

    if (coef%ifh2) then
       call addcol4 (au, coef%h2, coef%B, u, coef%dof%size())
       call addcol4 (av, coef%h2, coef%B, v, coef%dof%size())
       call addcol4 (aw, coef%h2, coef%B, w, coef%dof%size())
    end if

  end subroutine ax_helm_full_compute_vector

  subroutine ax_helm_stress_lx(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, h2, drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
       jacinv, weights3, n, lx)

    integer, intent(in) :: n, lx
    real(kind=rp), intent(inout) :: u(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: v(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdz(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: jacinv(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp), intent(in) :: weights3(lx, lx, lx)

    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    real(kind=rp) :: wvr(lx, lx, lx)
    real(kind=rp) :: wvs(lx, lx, lx)
    real(kind=rp) :: wvt(lx, lx, lx)
    real(kind=rp) :: wwr(lx, lx, lx)
    real(kind=rp) :: wws(lx, lx, lx)
    real(kind=rp) :: wwt(lx, lx, lx)

    integer :: e, i, j, k, l
    real(kind=rp) :: dj, t1, t2, t3
    real(kind=rp) :: s11, s12, s13, s21, s22, s23, s31, s32, s33
    real(kind=rp) :: u1, u2, u3, v1, v2, v3, w1, w2, w3

    do e = 1, n

       do j = 1, lx * lx
          do i = 1, lx
             t1 = 0.0_rp
             t2 = 0.0_rp
             t3 = 0.0_rp
             do k = 1, lx
                t1 = t1 + Dx(i,k) * u(k,j,1,e)
                t2 = t2 + Dx(i,k) * v(k,j,1,e)
                t3 = t3 + Dx(i,k) * w(k,j,1,e)
             end do
             wur(i,j,1) = t1
             wvr(i,j,1) = t2
             wwr(i,j,1) = t3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                t1 = 0.0_rp
                t2 = 0.0_rp
                t3 = 0.0_rp
                do l = 1, lx
                   t1 = t1 + Dy(j,l) * u(i,l,k,e)
                   t2 = t2 + Dy(j,l) * v(i,l,k,e)
                   t3 = t3 + Dy(j,l) * w(i,l,k,e)
                end do
                wus(i,j,k) = t1
                wvs(i,j,k) = t2
                wws(i,j,k) = t3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             t1 = 0.0_rp
             t2 = 0.0_rp
             t3 = 0.0_rp
             do l = 1, lx
                t1 = t1 + Dz(k,l) * u(i,1,l,e)
                t2 = t2 + Dz(k,l) * v(i,1,l,e)
                t3 = t3 + Dz(k,l) * w(i,1,l,e)
             end do
             wut(i,1,k) = t1
             wvt(i,1,k) = t2
             wwt(i,1,k) = t3
          end do
       end do

       do i = 1, lx*lx*lx

          u1 = wur(i,1,1) * drdx(i,1,1,e) + wus(i,1,1) * dsdx(i,1,1,e) &
                                          + wut(i,1,1) * dtdx(i,1,1,e)
          u2 = wur(i,1,1) * drdy(i,1,1,e) + wus(i,1,1) * dsdy(i,1,1,e) &
                                          + wut(i,1,1) * dtdy(i,1,1,e)
          u3 = wur(i,1,1) * drdz(i,1,1,e) + wus(i,1,1) * dsdz(i,1,1,e) &
                                          + wut(i,1,1) * dtdz(i,1,1,e)

          v1 = wvr(i,1,1) * drdx(i,1,1,e) + wvs(i,1,1) * dsdx(i,1,1,e) &
                                          + wvt(i,1,1) * dtdx(i,1,1,e)
          v2 = wvr(i,1,1) * drdy(i,1,1,e) + wvs(i,1,1) * dsdy(i,1,1,e) &
                                          + wvt(i,1,1) * dtdy(i,1,1,e)
          v3 = wvr(i,1,1) * drdz(i,1,1,e) + wvs(i,1,1) * dsdz(i,1,1,e) &
                                          + wvt(i,1,1) * dtdz(i,1,1,e)

          w1 = wwr(i,1,1) * drdx(i,1,1,e) + wws(i,1,1) * dsdx(i,1,1,e) &
                                          + wwt(i,1,1) * dtdx(i,1,1,e)
          w2 = wwr(i,1,1) * drdy(i,1,1,e) + wws(i,1,1) * dsdy(i,1,1,e) &
                                          + wwt(i,1,1) * dtdy(i,1,1,e)
          w3 = wwr(i,1,1) * drdz(i,1,1,e) + wws(i,1,1) * dsdz(i,1,1,e) &
                                          + wwt(i,1,1) * dtdz(i,1,1,e)

          dj  = h1(i,1,1,e) * weights3(i,1,1) * jacinv(i,1,1,e)
          s11 = dj*(u1 + u1)
          s12 = dj*(u2 + v1)
          s13 = dj*(u3 + w1)
          s21 = dj*(v1 + u2)
          s22 = dj*(v2 + v2)
          s23 = dj*(v3 + w2)
          s31 = dj*(w1 + u3)
          s32 = dj*(w2 + v3)
          s33 = dj*(w3 + w3)

          wur(i,1,1) = drdx(i,1,1,e)*s11 + drdy(i,1,1,e)*s12 + drdz(i,1,1,e)*s13
          wus(i,1,1) = dsdx(i,1,1,e)*s11 + dsdy(i,1,1,e)*s12 + dsdz(i,1,1,e)*s13
          wut(i,1,1) = dtdx(i,1,1,e)*s11 + dtdy(i,1,1,e)*s12 + dtdz(i,1,1,e)*s13

          wvr(i,1,1) = drdx(i,1,1,e)*s21 + drdy(i,1,1,e)*s22 + drdz(i,1,1,e)*s23
          wvs(i,1,1) = dsdx(i,1,1,e)*s21 + dsdy(i,1,1,e)*s22 + dsdz(i,1,1,e)*s23
          wvt(i,1,1) = dtdx(i,1,1,e)*s21 + dtdy(i,1,1,e)*s22 + dtdz(i,1,1,e)*s23

          wwr(i,1,1) = drdx(i,1,1,e)*s31 + drdy(i,1,1,e)*s32 + drdz(i,1,1,e)*s33
          wws(i,1,1) = dsdx(i,1,1,e)*s31 + dsdy(i,1,1,e)*s32 + dsdz(i,1,1,e)*s33
          wwt(i,1,1) = dtdx(i,1,1,e)*s31 + dtdy(i,1,1,e)*s32 + dtdz(i,1,1,e)*s33
       end do

       do j = 1, lx*lx
          do i = 1, lx
             t1 = 0.0_rp
             t2 = 0.0_rp
             t3 = 0.0_rp
             do k = 1, lx
                t1 = t1 + Dxt(i,k) * wur(k,j,1)
                t2 = t2 + Dxt(i,k) * wvr(k,j,1)
                t3 = t3 + Dxt(i,k) * wwr(k,j,1)
             end do
             au(i,j,1,e) = t1
             av(i,j,1,e) = t2
             aw(i,j,1,e) = t3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                t1 = 0.0_rp
                t2 = 0.0_rp
                t3 = 0.0_rp
                do l = 1, lx
                   t1 = t1 + Dyt(j,l) * wus(i,l,k)
                   t2 = t2 + Dyt(j,l) * wvs(i,l,k)
                   t3 = t3 + Dyt(j,l) * wws(i,l,k)
                end do
                au(i,j,k,e) = au(i,j,k,e) + t1
                av(i,j,k,e) = av(i,j,k,e) + t2
                aw(i,j,k,e) = aw(i,j,k,e) + t3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             t1 = 0.0_rp
             t2 = 0.0_rp
             t3 = 0.0_rp
             do l = 1, lx
                t1 = t1 + Dzt(k,l) * wut(i,1,l)
                t2 = t2 + Dzt(k,l) * wvt(i,1,l)
                t3 = t3 + Dzt(k,l) * wwt(i,1,l)
             end do
             au(i,1,k,e) = au(i,1,k,e) + t1
             av(i,1,k,e) = av(i,1,k,e) + t2
             aw(i,1,k,e) = aw(i,1,k,e) + t3
          end do
       end do
    end do

  end subroutine ax_helm_stress_lx

  subroutine ax_helm_stress_lx14(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, h2, drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
       jacinv, weights3, n)
    integer, parameter :: lx = 14
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: u(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: v(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdz(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: jacinv(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp), intent(in) :: weights3(lx, lx, lx)

    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    real(kind=rp) :: wvr(lx, lx, lx)
    real(kind=rp) :: wvs(lx, lx, lx)
    real(kind=rp) :: wvt(lx, lx, lx)
    real(kind=rp) :: wwr(lx, lx, lx)
    real(kind=rp) :: wws(lx, lx, lx)
    real(kind=rp) :: wwt(lx, lx, lx)

    integer :: e, i, j, k, l
    real(kind=rp) :: dj
    real(kind=rp) :: s11, s12, s13, s21, s22, s23, s31, s32, s33
    real(kind=rp) :: u1, u2, u3, v1, v2, v3, w1, w2, w3

    do e = 1, n

       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) &
                        + Dx(i,4) * u(4,j,1,e) &
                        + Dx(i,5) * u(5,j,1,e) &
                        + Dx(i,6) * u(6,j,1,e) &
                        + Dx(i,7) * u(7,j,1,e) &
                        + Dx(i,8) * u(8,j,1,e) &
                        + Dx(i,9) * u(9,j,1,e) &
                        + Dx(i,10) * u(10,j,1,e) &
                        + Dx(i,11) * u(11,j,1,e) &
                        + Dx(i,12) * u(12,j,1,e) &
                        + Dx(i,13) * u(13,j,1,e) &
                        + Dx(i,14) * u(14,j,1,e)

             wvr(i,j,1) = Dx(i,1) * v(1,j,1,e) &
                        + Dx(i,2) * v(2,j,1,e) &
                        + Dx(i,3) * v(3,j,1,e) &
                        + Dx(i,4) * v(4,j,1,e) &
                        + Dx(i,5) * v(5,j,1,e) &
                        + Dx(i,6) * v(6,j,1,e) &
                        + Dx(i,7) * v(7,j,1,e) &
                        + Dx(i,8) * v(8,j,1,e) &
                        + Dx(i,9) * v(9,j,1,e) &
                        + Dx(i,10) * v(10,j,1,e) &
                        + Dx(i,11) * v(11,j,1,e) &
                        + Dx(i,12) * v(12,j,1,e) &
                        + Dx(i,13) * v(13,j,1,e) &
                        + Dx(i,14) * v(14,j,1,e)

             wwr(i,j,1) = Dx(i,1) * w(1,j,1,e) &
                        + Dx(i,2) * w(2,j,1,e) &
                        + Dx(i,3) * w(3,j,1,e) &
                        + Dx(i,4) * w(4,j,1,e) &
                        + Dx(i,5) * w(5,j,1,e) &
                        + Dx(i,6) * w(6,j,1,e) &
                        + Dx(i,7) * w(7,j,1,e) &
                        + Dx(i,8) * w(8,j,1,e) &
                        + Dx(i,9) * w(9,j,1,e) &
                        + Dx(i,10) * w(10,j,1,e) &
                        + Dx(i,11) * w(11,j,1,e) &
                        + Dx(i,12) * w(12,j,1,e) &
                        + Dx(i,13) * w(13,j,1,e) &
                        + Dx(i,14) * w(14,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) &
                           + Dy(j,4) * u(i,4,k,e) &
                           + Dy(j,5) * u(i,5,k,e) &
                           + Dy(j,6) * u(i,6,k,e) &
                           + Dy(j,7) * u(i,7,k,e) &
                           + Dy(j,8) * u(i,8,k,e) &
                           + Dy(j,9) * u(i,9,k,e) &
                           + Dy(j,10) * u(i,10,k,e) &
                           + Dy(j,11) * u(i,11,k,e) &
                           + Dy(j,12) * u(i,12,k,e) &
                           + Dy(j,13) * u(i,13,k,e) &
                           + Dy(j,14) * u(i,14,k,e)


                wvs(i,j,k) = Dy(j,1) * v(i,1,k,e) &
                           + Dy(j,2) * v(i,2,k,e) &
                           + Dy(j,3) * v(i,3,k,e) &
                           + Dy(j,4) * v(i,4,k,e) &
                           + Dy(j,5) * v(i,5,k,e) &
                           + Dy(j,6) * v(i,6,k,e) &
                           + Dy(j,7) * v(i,7,k,e) &
                           + Dy(j,8) * v(i,8,k,e) &
                           + Dy(j,9) * v(i,9,k,e) &
                           + Dy(j,10) * v(i,10,k,e) &
                           + Dy(j,11) * v(i,11,k,e) &
                           + Dy(j,12) * v(i,12,k,e) &
                           + Dy(j,13) * v(i,13,k,e) &
                           + Dy(j,14) * v(i,14,k,e)

                wws(i,j,k) = Dy(j,1) * w(i,1,k,e) &
                           + Dy(j,2) * w(i,2,k,e) &
                           + Dy(j,3) * w(i,3,k,e) &
                           + Dy(j,4) * w(i,4,k,e) &
                           + Dy(j,5) * w(i,5,k,e) &
                           + Dy(j,6) * w(i,6,k,e) &
                           + Dy(j,7) * w(i,7,k,e) &
                           + Dy(j,8) * w(i,8,k,e) &
                           + Dy(j,9) * w(i,9,k,e) &
                           + Dy(j,10) * w(i,10,k,e) &
                           + Dy(j,11) * w(i,11,k,e) &
                           + Dy(j,12) * w(i,12,k,e) &
                           + Dy(j,13) * w(i,13,k,e) &
                           + Dy(j,14) * w(i,14,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) &
                        + Dz(k,4) * u(i,1,4,e) &
                        + Dz(k,5) * u(i,1,5,e) &
                        + Dz(k,6) * u(i,1,6,e) &
                        + Dz(k,7) * u(i,1,7,e) &
                        + Dz(k,8) * u(i,1,8,e) &
                        + Dz(k,9) * u(i,1,9,e) &
                        + Dz(k,10) * u(i,1,10,e) &
                        + Dz(k,11) * u(i,1,11,e) &
                        + Dz(k,12) * u(i,1,12,e) &
                        + Dz(k,13) * u(i,1,13,e) &
                        + Dz(k,14) * u(i,1,14,e)

             wvt(i,1,k) = Dz(k,1) * v(i,1,1,e) &
                        + Dz(k,2) * v(i,1,2,e) &
                        + Dz(k,3) * v(i,1,3,e) &
                        + Dz(k,4) * v(i,1,4,e) &
                        + Dz(k,5) * v(i,1,5,e) &
                        + Dz(k,6) * v(i,1,6,e) &
                        + Dz(k,7) * v(i,1,7,e) &
                        + Dz(k,8) * v(i,1,8,e) &
                        + Dz(k,9) * v(i,1,9,e) &
                        + Dz(k,10) * v(i,1,10,e) &
                        + Dz(k,11) * v(i,1,11,e) &
                        + Dz(k,12) * v(i,1,12,e) &
                        + Dz(k,13) * v(i,1,13,e) &
                        + Dz(k,14) * v(i,1,14,e)

             wwt(i,1,k) = Dz(k,1) * w(i,1,1,e) &
                        + Dz(k,2) * w(i,1,2,e) &
                        + Dz(k,3) * w(i,1,3,e) &
                        + Dz(k,4) * w(i,1,4,e) &
                        + Dz(k,5) * w(i,1,5,e) &
                        + Dz(k,6) * w(i,1,6,e) &
                        + Dz(k,7) * w(i,1,7,e) &
                        + Dz(k,8) * w(i,1,8,e) &
                        + Dz(k,9) * w(i,1,9,e) &
                        + Dz(k,10) * w(i,1,10,e) &
                        + Dz(k,11) * w(i,1,11,e) &
                        + Dz(k,12) * w(i,1,12,e) &
                        + Dz(k,13) * w(i,1,13,e) &
                        + Dz(k,14) * w(i,1,14,e)
          end do
       end do

       do i = 1, lx*lx*lx

          u1 = wur(i,1,1) * drdx(i,1,1,e) + wus(i,1,1) * dsdx(i,1,1,e) &
                                          + wut(i,1,1) * dtdx(i,1,1,e)
          u2 = wur(i,1,1) * drdy(i,1,1,e) + wus(i,1,1) * dsdy(i,1,1,e) &
                                          + wut(i,1,1) * dtdy(i,1,1,e)
          u3 = wur(i,1,1) * drdz(i,1,1,e) + wus(i,1,1) * dsdz(i,1,1,e) &
                                          + wut(i,1,1) * dtdz(i,1,1,e)

          v1 = wvr(i,1,1) * drdx(i,1,1,e) + wvs(i,1,1) * dsdx(i,1,1,e) &
                                          + wvt(i,1,1) * dtdx(i,1,1,e)
          v2 = wvr(i,1,1) * drdy(i,1,1,e) + wvs(i,1,1) * dsdy(i,1,1,e) &
                                          + wvt(i,1,1) * dtdy(i,1,1,e)
          v3 = wvr(i,1,1) * drdz(i,1,1,e) + wvs(i,1,1) * dsdz(i,1,1,e) &
                                          + wvt(i,1,1) * dtdz(i,1,1,e)

          w1 = wwr(i,1,1) * drdx(i,1,1,e) + wws(i,1,1) * dsdx(i,1,1,e) &
                                          + wwt(i,1,1) * dtdx(i,1,1,e)
          w2 = wwr(i,1,1) * drdy(i,1,1,e) + wws(i,1,1) * dsdy(i,1,1,e) &
                                          + wwt(i,1,1) * dtdy(i,1,1,e)
          w3 = wwr(i,1,1) * drdz(i,1,1,e) + wws(i,1,1) * dsdz(i,1,1,e) &
                                          + wwt(i,1,1) * dtdz(i,1,1,e)

          dj  = h1(i,1,1,e) * weights3(i,1,1) * jacinv(i,1,1,e)
          s11 = dj*(u1 + u1)
          s12 = dj*(u2 + v1)
          s13 = dj*(u3 + w1)
          s21 = dj*(v1 + u2)
          s22 = dj*(v2 + v2)
          s23 = dj*(v3 + w2)
          s31 = dj*(w1 + u3)
          s32 = dj*(w2 + v3)
          s33 = dj*(w3 + w3)

          wur(i,1,1) = drdx(i,1,1,e)*s11 + drdy(i,1,1,e)*s12 + drdz(i,1,1,e)*s13
          wus(i,1,1) = dsdx(i,1,1,e)*s11 + dsdy(i,1,1,e)*s12 + dsdz(i,1,1,e)*s13
          wut(i,1,1) = dtdx(i,1,1,e)*s11 + dtdy(i,1,1,e)*s12 + dtdz(i,1,1,e)*s13

          wvr(i,1,1) = drdx(i,1,1,e)*s21 + drdy(i,1,1,e)*s22 + drdz(i,1,1,e)*s23
          wvs(i,1,1) = dsdx(i,1,1,e)*s21 + dsdy(i,1,1,e)*s22 + dsdz(i,1,1,e)*s23
          wvt(i,1,1) = dtdx(i,1,1,e)*s21 + dtdy(i,1,1,e)*s22 + dtdz(i,1,1,e)*s23

          wwr(i,1,1) = drdx(i,1,1,e)*s31 + drdy(i,1,1,e)*s32 + drdz(i,1,1,e)*s33
          wws(i,1,1) = dsdx(i,1,1,e)*s31 + dsdy(i,1,1,e)*s32 + dsdz(i,1,1,e)*s33
          wwt(i,1,1) = dtdx(i,1,1,e)*s31 + dtdy(i,1,1,e)*s32 + dtdz(i,1,1,e)*s33
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * wur(1,j,1) &
                         + Dxt(i,2) * wur(2,j,1) &
                         + Dxt(i,3) * wur(3,j,1) &
                         + Dxt(i,4) * wur(4,j,1) &
                         + Dxt(i,5) * wur(5,j,1) &
                         + Dxt(i,6) * wur(6,j,1) &
                         + Dxt(i,7) * wur(7,j,1) &
                         + Dxt(i,8) * wur(8,j,1) &
                         + Dxt(i,9) * wur(9,j,1) &
                         + Dxt(i,10) * wur(10,j,1) &
                         + Dxt(i,11) * wur(11,j,1) &
                         + Dxt(i,12) * wur(12,j,1) &
                         + Dxt(i,13) * wur(13,j,1) &
                         + Dxt(i,14) * wur(14,j,1)

             av(i,j,1,e) = Dxt(i,1) * wvr(1,j,1) &
                         + Dxt(i,2) * wvr(2,j,1) &
                         + Dxt(i,3) * wvr(3,j,1) &
                         + Dxt(i,4) * wvr(4,j,1) &
                         + Dxt(i,5) * wvr(5,j,1) &
                         + Dxt(i,6) * wvr(6,j,1) &
                         + Dxt(i,7) * wvr(7,j,1) &
                         + Dxt(i,8) * wvr(8,j,1) &
                         + Dxt(i,9) * wvr(9,j,1) &
                         + Dxt(i,10) * wvr(10,j,1) &
                         + Dxt(i,11) * wvr(11,j,1) &
                         + Dxt(i,12) * wvr(12,j,1) &
                         + Dxt(i,13) * wvr(13,j,1) &
                         + Dxt(i,14) * wvr(14,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wwr(1,j,1) &
                         + Dxt(i,2) * wwr(2,j,1) &
                         + Dxt(i,3) * wwr(3,j,1) &
                         + Dxt(i,4) * wwr(4,j,1) &
                         + Dxt(i,5) * wwr(5,j,1) &
                         + Dxt(i,6) * wwr(6,j,1) &
                         + Dxt(i,7) * wwr(7,j,1) &
                         + Dxt(i,8) * wwr(8,j,1) &
                         + Dxt(i,9) * wwr(9,j,1) &
                         + Dxt(i,10) * wwr(10,j,1) &
                         + Dxt(i,11) * wwr(11,j,1) &
                         + Dxt(i,12) * wwr(12,j,1) &
                         + Dxt(i,13) * wwr(13,j,1) &
                         + Dxt(i,14) * wwr(14,j,1)

          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * wus(i,1,k) &
                            + Dyt(j,2) * wus(i,2,k) &
                            + Dyt(j,3) * wus(i,3,k) &
                            + Dyt(j,4) * wus(i,4,k) &
                            + Dyt(j,5) * wus(i,5,k) &
                            + Dyt(j,6) * wus(i,6,k) &
                            + Dyt(j,7) * wus(i,7,k) &
                            + Dyt(j,8) * wus(i,8,k) &
                            + Dyt(j,9) * wus(i,9,k) &
                            + Dyt(j,10) * wus(i,10,k) &
                            + Dyt(j,11) * wus(i,11,k) &
                            + Dyt(j,12) * wus(i,12,k) &
                            + Dyt(j,13) * wus(i,13,k) &
                            + Dyt(j,14) * wus(i,14,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * wvs(i,1,k) &
                            + Dyt(j,2) * wvs(i,2,k) &
                            + Dyt(j,3) * wvs(i,3,k) &
                            + Dyt(j,4) * wvs(i,4,k) &
                            + Dyt(j,5) * wvs(i,5,k) &
                            + Dyt(j,6) * wvs(i,6,k) &
                            + Dyt(j,7) * wvs(i,7,k) &
                            + Dyt(j,8) * wvs(i,8,k) &
                            + Dyt(j,9) * wvs(i,9,k) &
                            + Dyt(j,10) * wvs(i,10,k) &
                            + Dyt(j,11) * wvs(i,11,k) &
                            + Dyt(j,12) * wvs(i,12,k) &
                            + Dyt(j,13) * wvs(i,13,k) &
                            + Dyt(j,14) * wvs(i,14,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * wws(i,1,k) &
                            + Dyt(j,2) * wws(i,2,k) &
                            + Dyt(j,3) * wws(i,3,k) &
                            + Dyt(j,4) * wws(i,4,k) &
                            + Dyt(j,5) * wws(i,5,k) &
                            + Dyt(j,6) * wws(i,6,k) &
                            + Dyt(j,7) * wws(i,7,k) &
                            + Dyt(j,8) * wws(i,8,k) &
                            + Dyt(j,9) * wws(i,9,k) &
                            + Dyt(j,10) * wws(i,10,k) &
                            + Dyt(j,11) * wws(i,11,k) &
                            + Dyt(j,12) * wws(i,12,k) &
                            + Dyt(j,13) * wws(i,13,k) &
                            + Dyt(j,14) * wws(i,14,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * wut(i,1,1) &
                         + Dzt(k,2) * wut(i,1,2) &
                         + Dzt(k,3) * wut(i,1,3) &
                         + Dzt(k,4) * wut(i,1,4) &
                         + Dzt(k,5) * wut(i,1,5) &
                         + Dzt(k,6) * wut(i,1,6) &
                         + Dzt(k,7) * wut(i,1,7) &
                         + Dzt(k,8) * wut(i,1,8) &
                         + Dzt(k,9) * wut(i,1,9) &
                         + Dzt(k,10) * wut(i,1,10) &
                         + Dzt(k,11) * wut(i,1,11) &
                         + Dzt(k,12) * wut(i,1,12) &
                         + Dzt(k,13) * wut(i,1,13) &
                         + Dzt(k,14) * wut(i,1,14)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * wvt(i,1,1) &
                         + Dzt(k,2) * wvt(i,1,2) &
                         + Dzt(k,3) * wvt(i,1,3) &
                         + Dzt(k,4) * wvt(i,1,4) &
                         + Dzt(k,5) * wvt(i,1,5) &
                         + Dzt(k,6) * wvt(i,1,6) &
                         + Dzt(k,7) * wvt(i,1,7) &
                         + Dzt(k,8) * wvt(i,1,8) &
                         + Dzt(k,9) * wvt(i,1,9) &
                         + Dzt(k,10) * wvt(i,1,10) &
                         + Dzt(k,11) * wvt(i,1,11) &
                         + Dzt(k,12) * wvt(i,1,12) &
                         + Dzt(k,13) * wvt(i,1,13) &
                         + Dzt(k,14) * wvt(i,1,14)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wwt(i,1,1) &
                         + Dzt(k,2) * wwt(i,1,2) &
                         + Dzt(k,3) * wwt(i,1,3) &
                         + Dzt(k,4) * wwt(i,1,4) &
                         + Dzt(k,5) * wwt(i,1,5) &
                         + Dzt(k,6) * wwt(i,1,6) &
                         + Dzt(k,7) * wwt(i,1,7) &
                         + Dzt(k,8) * wwt(i,1,8) &
                         + Dzt(k,9) * wwt(i,1,9) &
                         + Dzt(k,10) * wwt(i,1,10) &
                         + Dzt(k,11) * wwt(i,1,11) &
                         + Dzt(k,12) * wwt(i,1,12) &
                         + Dzt(k,13) * wwt(i,1,13) &
                         + Dzt(k,14) * wwt(i,1,14)
          end do
       end do

    end do

  end subroutine ax_helm_stress_lx14

  subroutine ax_helm_stress_lx13(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, &
       Dzt, h1, h2, drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
       jacinv, weights3, n)
    integer, parameter :: lx = 13
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: u(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: v(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdz(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: jacinv(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp), intent(in) :: weights3(lx, lx, lx)

    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    real(kind=rp) :: wvr(lx, lx, lx)
    real(kind=rp) :: wvs(lx, lx, lx)
    real(kind=rp) :: wvt(lx, lx, lx)
    real(kind=rp) :: wwr(lx, lx, lx)
    real(kind=rp) :: wws(lx, lx, lx)
    real(kind=rp) :: wwt(lx, lx, lx)

    integer :: e, i, j, k, l
    real(kind=rp) :: dj
    real(kind=rp) :: s11, s12, s13, s21, s22, s23, s31, s32, s33
    real(kind=rp) :: u1, u2, u3, v1, v2, v3, w1, w2, w3

    do e = 1, n

       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) &
                        + Dx(i,4) * u(4,j,1,e) &
                        + Dx(i,5) * u(5,j,1,e) &
                        + Dx(i,6) * u(6,j,1,e) &
                        + Dx(i,7) * u(7,j,1,e) &
                        + Dx(i,8) * u(8,j,1,e) &
                        + Dx(i,9) * u(9,j,1,e) &
                        + Dx(i,10) * u(10,j,1,e) &
                        + Dx(i,11) * u(11,j,1,e) &
                        + Dx(i,12) * u(12,j,1,e) &
                        + Dx(i,13) * u(13,j,1,e)

             wvr(i,j,1) = Dx(i,1) * v(1,j,1,e) &
                        + Dx(i,2) * v(2,j,1,e) &
                        + Dx(i,3) * v(3,j,1,e) &
                        + Dx(i,4) * v(4,j,1,e) &
                        + Dx(i,5) * v(5,j,1,e) &
                        + Dx(i,6) * v(6,j,1,e) &
                        + Dx(i,7) * v(7,j,1,e) &
                        + Dx(i,8) * v(8,j,1,e) &
                        + Dx(i,9) * v(9,j,1,e) &
                        + Dx(i,10) * v(10,j,1,e) &
                        + Dx(i,11) * v(11,j,1,e) &
                        + Dx(i,12) * v(12,j,1,e) &
                        + Dx(i,13) * v(13,j,1,e)

             wwr(i,j,1) = Dx(i,1) * w(1,j,1,e) &
                        + Dx(i,2) * w(2,j,1,e) &
                        + Dx(i,3) * w(3,j,1,e) &
                        + Dx(i,4) * w(4,j,1,e) &
                        + Dx(i,5) * w(5,j,1,e) &
                        + Dx(i,6) * w(6,j,1,e) &
                        + Dx(i,7) * w(7,j,1,e) &
                        + Dx(i,8) * w(8,j,1,e) &
                        + Dx(i,9) * w(9,j,1,e) &
                        + Dx(i,10) * w(10,j,1,e) &
                        + Dx(i,11) * w(11,j,1,e) &
                        + Dx(i,12) * w(12,j,1,e) &
                        + Dx(i,13) * w(13,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) &
                           + Dy(j,4) * u(i,4,k,e) &
                           + Dy(j,5) * u(i,5,k,e) &
                           + Dy(j,6) * u(i,6,k,e) &
                           + Dy(j,7) * u(i,7,k,e) &
                           + Dy(j,8) * u(i,8,k,e) &
                           + Dy(j,9) * u(i,9,k,e) &
                           + Dy(j,10) * u(i,10,k,e) &
                           + Dy(j,11) * u(i,11,k,e) &
                           + Dy(j,12) * u(i,12,k,e) &
                           + Dy(j,13) * u(i,13,k,e)


                wvs(i,j,k) = Dy(j,1) * v(i,1,k,e) &
                           + Dy(j,2) * v(i,2,k,e) &
                           + Dy(j,3) * v(i,3,k,e) &
                           + Dy(j,4) * v(i,4,k,e) &
                           + Dy(j,5) * v(i,5,k,e) &
                           + Dy(j,6) * v(i,6,k,e) &
                           + Dy(j,7) * v(i,7,k,e) &
                           + Dy(j,8) * v(i,8,k,e) &
                           + Dy(j,9) * v(i,9,k,e) &
                           + Dy(j,10) * v(i,10,k,e) &
                           + Dy(j,11) * v(i,11,k,e) &
                           + Dy(j,12) * v(i,12,k,e) &
                           + Dy(j,13) * v(i,13,k,e)

                wws(i,j,k) = Dy(j,1) * w(i,1,k,e) &
                           + Dy(j,2) * w(i,2,k,e) &
                           + Dy(j,3) * w(i,3,k,e) &
                           + Dy(j,4) * w(i,4,k,e) &
                           + Dy(j,5) * w(i,5,k,e) &
                           + Dy(j,6) * w(i,6,k,e) &
                           + Dy(j,7) * w(i,7,k,e) &
                           + Dy(j,8) * w(i,8,k,e) &
                           + Dy(j,9) * w(i,9,k,e) &
                           + Dy(j,10) * w(i,10,k,e) &
                           + Dy(j,11) * w(i,11,k,e) &
                           + Dy(j,12) * w(i,12,k,e) &
                           + Dy(j,13) * w(i,13,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) &
                        + Dz(k,4) * u(i,1,4,e) &
                        + Dz(k,5) * u(i,1,5,e) &
                        + Dz(k,6) * u(i,1,6,e) &
                        + Dz(k,7) * u(i,1,7,e) &
                        + Dz(k,8) * u(i,1,8,e) &
                        + Dz(k,9) * u(i,1,9,e) &
                        + Dz(k,10) * u(i,1,10,e) &
                        + Dz(k,11) * u(i,1,11,e) &
                        + Dz(k,12) * u(i,1,12,e) &
                        + Dz(k,13) * u(i,1,13,e)

             wvt(i,1,k) = Dz(k,1) * v(i,1,1,e) &
                        + Dz(k,2) * v(i,1,2,e) &
                        + Dz(k,3) * v(i,1,3,e) &
                        + Dz(k,4) * v(i,1,4,e) &
                        + Dz(k,5) * v(i,1,5,e) &
                        + Dz(k,6) * v(i,1,6,e) &
                        + Dz(k,7) * v(i,1,7,e) &
                        + Dz(k,8) * v(i,1,8,e) &
                        + Dz(k,9) * v(i,1,9,e) &
                        + Dz(k,10) * v(i,1,10,e) &
                        + Dz(k,11) * v(i,1,11,e) &
                        + Dz(k,12) * v(i,1,12,e) &
                        + Dz(k,13) * v(i,1,13,e)

             wwt(i,1,k) = Dz(k,1) * w(i,1,1,e) &
                        + Dz(k,2) * w(i,1,2,e) &
                        + Dz(k,3) * w(i,1,3,e) &
                        + Dz(k,4) * w(i,1,4,e) &
                        + Dz(k,5) * w(i,1,5,e) &
                        + Dz(k,6) * w(i,1,6,e) &
                        + Dz(k,7) * w(i,1,7,e) &
                        + Dz(k,8) * w(i,1,8,e) &
                        + Dz(k,9) * w(i,1,9,e) &
                        + Dz(k,10) * w(i,1,10,e) &
                        + Dz(k,11) * w(i,1,11,e) &
                        + Dz(k,12) * w(i,1,12,e) &
                        + Dz(k,13) * w(i,1,13,e)
          end do
       end do

       do i = 1, lx*lx*lx

          u1 = wur(i,1,1) * drdx(i,1,1,e) + wus(i,1,1) * dsdx(i,1,1,e) &
                                          + wut(i,1,1) * dtdx(i,1,1,e)
          u2 = wur(i,1,1) * drdy(i,1,1,e) + wus(i,1,1) * dsdy(i,1,1,e) &
                                          + wut(i,1,1) * dtdy(i,1,1,e)
          u3 = wur(i,1,1) * drdz(i,1,1,e) + wus(i,1,1) * dsdz(i,1,1,e) &
                                          + wut(i,1,1) * dtdz(i,1,1,e)

          v1 = wvr(i,1,1) * drdx(i,1,1,e) + wvs(i,1,1) * dsdx(i,1,1,e) &
                                          + wvt(i,1,1) * dtdx(i,1,1,e)
          v2 = wvr(i,1,1) * drdy(i,1,1,e) + wvs(i,1,1) * dsdy(i,1,1,e) &
                                          + wvt(i,1,1) * dtdy(i,1,1,e)
          v3 = wvr(i,1,1) * drdz(i,1,1,e) + wvs(i,1,1) * dsdz(i,1,1,e) &
                                          + wvt(i,1,1) * dtdz(i,1,1,e)

          w1 = wwr(i,1,1) * drdx(i,1,1,e) + wws(i,1,1) * dsdx(i,1,1,e) &
                                          + wwt(i,1,1) * dtdx(i,1,1,e)
          w2 = wwr(i,1,1) * drdy(i,1,1,e) + wws(i,1,1) * dsdy(i,1,1,e) &
                                          + wwt(i,1,1) * dtdy(i,1,1,e)
          w3 = wwr(i,1,1) * drdz(i,1,1,e) + wws(i,1,1) * dsdz(i,1,1,e) &
                                          + wwt(i,1,1) * dtdz(i,1,1,e)

          dj  = h1(i,1,1,e) * weights3(i,1,1) * jacinv(i,1,1,e)
          s11 = dj*(u1 + u1)
          s12 = dj*(u2 + v1)
          s13 = dj*(u3 + w1)
          s21 = dj*(v1 + u2)
          s22 = dj*(v2 + v2)
          s23 = dj*(v3 + w2)
          s31 = dj*(w1 + u3)
          s32 = dj*(w2 + v3)
          s33 = dj*(w3 + w3)

          wur(i,1,1) = drdx(i,1,1,e)*s11 + drdy(i,1,1,e)*s12 + drdz(i,1,1,e)*s13
          wus(i,1,1) = dsdx(i,1,1,e)*s11 + dsdy(i,1,1,e)*s12 + dsdz(i,1,1,e)*s13
          wut(i,1,1) = dtdx(i,1,1,e)*s11 + dtdy(i,1,1,e)*s12 + dtdz(i,1,1,e)*s13

          wvr(i,1,1) = drdx(i,1,1,e)*s21 + drdy(i,1,1,e)*s22 + drdz(i,1,1,e)*s23
          wvs(i,1,1) = dsdx(i,1,1,e)*s21 + dsdy(i,1,1,e)*s22 + dsdz(i,1,1,e)*s23
          wvt(i,1,1) = dtdx(i,1,1,e)*s21 + dtdy(i,1,1,e)*s22 + dtdz(i,1,1,e)*s23

          wwr(i,1,1) = drdx(i,1,1,e)*s31 + drdy(i,1,1,e)*s32 + drdz(i,1,1,e)*s33
          wws(i,1,1) = dsdx(i,1,1,e)*s31 + dsdy(i,1,1,e)*s32 + dsdz(i,1,1,e)*s33
          wwt(i,1,1) = dtdx(i,1,1,e)*s31 + dtdy(i,1,1,e)*s32 + dtdz(i,1,1,e)*s33
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * wur(1,j,1) &
                         + Dxt(i,2) * wur(2,j,1) &
                         + Dxt(i,3) * wur(3,j,1) &
                         + Dxt(i,4) * wur(4,j,1) &
                         + Dxt(i,5) * wur(5,j,1) &
                         + Dxt(i,6) * wur(6,j,1) &
                         + Dxt(i,7) * wur(7,j,1) &
                         + Dxt(i,8) * wur(8,j,1) &
                         + Dxt(i,9) * wur(9,j,1) &
                         + Dxt(i,10) * wur(10,j,1) &
                         + Dxt(i,11) * wur(11,j,1) &
                         + Dxt(i,12) * wur(12,j,1) &
                         + Dxt(i,13) * wur(13,j,1)

             av(i,j,1,e) = Dxt(i,1) * wvr(1,j,1) &
                         + Dxt(i,2) * wvr(2,j,1) &
                         + Dxt(i,3) * wvr(3,j,1) &
                         + Dxt(i,4) * wvr(4,j,1) &
                         + Dxt(i,5) * wvr(5,j,1) &
                         + Dxt(i,6) * wvr(6,j,1) &
                         + Dxt(i,7) * wvr(7,j,1) &
                         + Dxt(i,8) * wvr(8,j,1) &
                         + Dxt(i,9) * wvr(9,j,1) &
                         + Dxt(i,10) * wvr(10,j,1) &
                         + Dxt(i,11) * wvr(11,j,1) &
                         + Dxt(i,12) * wvr(12,j,1) &
                         + Dxt(i,13) * wvr(13,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wwr(1,j,1) &
                         + Dxt(i,2) * wwr(2,j,1) &
                         + Dxt(i,3) * wwr(3,j,1) &
                         + Dxt(i,4) * wwr(4,j,1) &
                         + Dxt(i,5) * wwr(5,j,1) &
                         + Dxt(i,6) * wwr(6,j,1) &
                         + Dxt(i,7) * wwr(7,j,1) &
                         + Dxt(i,8) * wwr(8,j,1) &
                         + Dxt(i,9) * wwr(9,j,1) &
                         + Dxt(i,10) * wwr(10,j,1) &
                         + Dxt(i,11) * wwr(11,j,1) &
                         + Dxt(i,12) * wwr(12,j,1) &
                         + Dxt(i,13) * wwr(13,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * wus(i,1,k) &
                            + Dyt(j,2) * wus(i,2,k) &
                            + Dyt(j,3) * wus(i,3,k) &
                            + Dyt(j,4) * wus(i,4,k) &
                            + Dyt(j,5) * wus(i,5,k) &
                            + Dyt(j,6) * wus(i,6,k) &
                            + Dyt(j,7) * wus(i,7,k) &
                            + Dyt(j,8) * wus(i,8,k) &
                            + Dyt(j,9) * wus(i,9,k) &
                            + Dyt(j,10) * wus(i,10,k) &
                            + Dyt(j,11) * wus(i,11,k) &
                            + Dyt(j,12) * wus(i,12,k) &
                            + Dyt(j,13) * wus(i,13,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * wvs(i,1,k) &
                            + Dyt(j,2) * wvs(i,2,k) &
                            + Dyt(j,3) * wvs(i,3,k) &
                            + Dyt(j,4) * wvs(i,4,k) &
                            + Dyt(j,5) * wvs(i,5,k) &
                            + Dyt(j,6) * wvs(i,6,k) &
                            + Dyt(j,7) * wvs(i,7,k) &
                            + Dyt(j,8) * wvs(i,8,k) &
                            + Dyt(j,9) * wvs(i,9,k) &
                            + Dyt(j,10) * wvs(i,10,k) &
                            + Dyt(j,11) * wvs(i,11,k) &
                            + Dyt(j,12) * wvs(i,12,k) &
                            + Dyt(j,13) * wvs(i,13,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * wws(i,1,k) &
                            + Dyt(j,2) * wws(i,2,k) &
                            + Dyt(j,3) * wws(i,3,k) &
                            + Dyt(j,4) * wws(i,4,k) &
                            + Dyt(j,5) * wws(i,5,k) &
                            + Dyt(j,6) * wws(i,6,k) &
                            + Dyt(j,7) * wws(i,7,k) &
                            + Dyt(j,8) * wws(i,8,k) &
                            + Dyt(j,9) * wws(i,9,k) &
                            + Dyt(j,10) * wws(i,10,k) &
                            + Dyt(j,11) * wws(i,11,k) &
                            + Dyt(j,12) * wws(i,12,k) &
                            + Dyt(j,13) * wws(i,13,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * wut(i,1,1) &
                         + Dzt(k,2) * wut(i,1,2) &
                         + Dzt(k,3) * wut(i,1,3) &
                         + Dzt(k,4) * wut(i,1,4) &
                         + Dzt(k,5) * wut(i,1,5) &
                         + Dzt(k,6) * wut(i,1,6) &
                         + Dzt(k,7) * wut(i,1,7) &
                         + Dzt(k,8) * wut(i,1,8) &
                         + Dzt(k,9) * wut(i,1,9) &
                         + Dzt(k,10) * wut(i,1,10) &
                         + Dzt(k,11) * wut(i,1,11) &
                         + Dzt(k,12) * wut(i,1,12) &
                         + Dzt(k,13) * wut(i,1,13)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * wvt(i,1,1) &
                         + Dzt(k,2) * wvt(i,1,2) &
                         + Dzt(k,3) * wvt(i,1,3) &
                         + Dzt(k,4) * wvt(i,1,4) &
                         + Dzt(k,5) * wvt(i,1,5) &
                         + Dzt(k,6) * wvt(i,1,6) &
                         + Dzt(k,7) * wvt(i,1,7) &
                         + Dzt(k,8) * wvt(i,1,8) &
                         + Dzt(k,9) * wvt(i,1,9) &
                         + Dzt(k,10) * wvt(i,1,10) &
                         + Dzt(k,11) * wvt(i,1,11) &
                         + Dzt(k,12) * wvt(i,1,12) &
                         + Dzt(k,13) * wvt(i,1,13)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wwt(i,1,1) &
                         + Dzt(k,2) * wwt(i,1,2) &
                         + Dzt(k,3) * wwt(i,1,3) &
                         + Dzt(k,4) * wwt(i,1,4) &
                         + Dzt(k,5) * wwt(i,1,5) &
                         + Dzt(k,6) * wwt(i,1,6) &
                         + Dzt(k,7) * wwt(i,1,7) &
                         + Dzt(k,8) * wwt(i,1,8) &
                         + Dzt(k,9) * wwt(i,1,9) &
                         + Dzt(k,10) * wwt(i,1,10) &
                         + Dzt(k,11) * wwt(i,1,11) &
                         + Dzt(k,12) * wwt(i,1,12) &
                         + Dzt(k,13) * wwt(i,1,13)
          end do
       end do

    end do

  end subroutine ax_helm_stress_lx13

  subroutine ax_helm_stress_lx12(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, &
       Dzt, h1, h2, drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
       jacinv, weights3, n)
    integer, parameter :: lx = 12
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: u(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: v(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdz(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: jacinv(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp), intent(in) :: weights3(lx, lx, lx)

    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    real(kind=rp) :: wvr(lx, lx, lx)
    real(kind=rp) :: wvs(lx, lx, lx)
    real(kind=rp) :: wvt(lx, lx, lx)
    real(kind=rp) :: wwr(lx, lx, lx)
    real(kind=rp) :: wws(lx, lx, lx)
    real(kind=rp) :: wwt(lx, lx, lx)

    integer :: e, i, j, k, l
    real(kind=rp) :: dj
    real(kind=rp) :: s11, s12, s13, s21, s22, s23, s31, s32, s33
    real(kind=rp) :: u1, u2, u3, v1, v2, v3, w1, w2, w3

    do e = 1, n

       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) &
                        + Dx(i,4) * u(4,j,1,e) &
                        + Dx(i,5) * u(5,j,1,e) &
                        + Dx(i,6) * u(6,j,1,e) &
                        + Dx(i,7) * u(7,j,1,e) &
                        + Dx(i,8) * u(8,j,1,e) &
                        + Dx(i,9) * u(9,j,1,e) &
                        + Dx(i,10) * u(10,j,1,e) &
                        + Dx(i,11) * u(11,j,1,e) &
                        + Dx(i,12) * u(12,j,1,e)

             wvr(i,j,1) = Dx(i,1) * v(1,j,1,e) &
                        + Dx(i,2) * v(2,j,1,e) &
                        + Dx(i,3) * v(3,j,1,e) &
                        + Dx(i,4) * v(4,j,1,e) &
                        + Dx(i,5) * v(5,j,1,e) &
                        + Dx(i,6) * v(6,j,1,e) &
                        + Dx(i,7) * v(7,j,1,e) &
                        + Dx(i,8) * v(8,j,1,e) &
                        + Dx(i,9) * v(9,j,1,e) &
                        + Dx(i,10) * v(10,j,1,e) &
                        + Dx(i,11) * v(11,j,1,e) &
                        + Dx(i,12) * v(12,j,1,e)

             wwr(i,j,1) = Dx(i,1) * w(1,j,1,e) &
                        + Dx(i,2) * w(2,j,1,e) &
                        + Dx(i,3) * w(3,j,1,e) &
                        + Dx(i,4) * w(4,j,1,e) &
                        + Dx(i,5) * w(5,j,1,e) &
                        + Dx(i,6) * w(6,j,1,e) &
                        + Dx(i,7) * w(7,j,1,e) &
                        + Dx(i,8) * w(8,j,1,e) &
                        + Dx(i,9) * w(9,j,1,e) &
                        + Dx(i,10) * w(10,j,1,e) &
                        + Dx(i,11) * w(11,j,1,e) &
                        + Dx(i,12) * w(12,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) &
                           + Dy(j,4) * u(i,4,k,e) &
                           + Dy(j,5) * u(i,5,k,e) &
                           + Dy(j,6) * u(i,6,k,e) &
                           + Dy(j,7) * u(i,7,k,e) &
                           + Dy(j,8) * u(i,8,k,e) &
                           + Dy(j,9) * u(i,9,k,e) &
                           + Dy(j,10) * u(i,10,k,e) &
                           + Dy(j,11) * u(i,11,k,e) &
                           + Dy(j,12) * u(i,12,k,e)


                wvs(i,j,k) = Dy(j,1) * v(i,1,k,e) &
                           + Dy(j,2) * v(i,2,k,e) &
                           + Dy(j,3) * v(i,3,k,e) &
                           + Dy(j,4) * v(i,4,k,e) &
                           + Dy(j,5) * v(i,5,k,e) &
                           + Dy(j,6) * v(i,6,k,e) &
                           + Dy(j,7) * v(i,7,k,e) &
                           + Dy(j,8) * v(i,8,k,e) &
                           + Dy(j,9) * v(i,9,k,e) &
                           + Dy(j,10) * v(i,10,k,e) &
                           + Dy(j,11) * v(i,11,k,e) &
                           + Dy(j,12) * v(i,12,k,e)

                wws(i,j,k) = Dy(j,1) * w(i,1,k,e) &
                           + Dy(j,2) * w(i,2,k,e) &
                           + Dy(j,3) * w(i,3,k,e) &
                           + Dy(j,4) * w(i,4,k,e) &
                           + Dy(j,5) * w(i,5,k,e) &
                           + Dy(j,6) * w(i,6,k,e) &
                           + Dy(j,7) * w(i,7,k,e) &
                           + Dy(j,8) * w(i,8,k,e) &
                           + Dy(j,9) * w(i,9,k,e) &
                           + Dy(j,10) * w(i,10,k,e) &
                           + Dy(j,11) * w(i,11,k,e) &
                           + Dy(j,12) * w(i,12,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) &
                        + Dz(k,4) * u(i,1,4,e) &
                        + Dz(k,5) * u(i,1,5,e) &
                        + Dz(k,6) * u(i,1,6,e) &
                        + Dz(k,7) * u(i,1,7,e) &
                        + Dz(k,8) * u(i,1,8,e) &
                        + Dz(k,9) * u(i,1,9,e) &
                        + Dz(k,10) * u(i,1,10,e) &
                        + Dz(k,11) * u(i,1,11,e) &
                        + Dz(k,12) * u(i,1,12,e)

             wvt(i,1,k) = Dz(k,1) * v(i,1,1,e) &
                        + Dz(k,2) * v(i,1,2,e) &
                        + Dz(k,3) * v(i,1,3,e) &
                        + Dz(k,4) * v(i,1,4,e) &
                        + Dz(k,5) * v(i,1,5,e) &
                        + Dz(k,6) * v(i,1,6,e) &
                        + Dz(k,7) * v(i,1,7,e) &
                        + Dz(k,8) * v(i,1,8,e) &
                        + Dz(k,9) * v(i,1,9,e) &
                        + Dz(k,10) * v(i,1,10,e) &
                        + Dz(k,11) * v(i,1,11,e) &
                        + Dz(k,12) * v(i,1,12,e)

             wwt(i,1,k) = Dz(k,1) * w(i,1,1,e) &
                        + Dz(k,2) * w(i,1,2,e) &
                        + Dz(k,3) * w(i,1,3,e) &
                        + Dz(k,4) * w(i,1,4,e) &
                        + Dz(k,5) * w(i,1,5,e) &
                        + Dz(k,6) * w(i,1,6,e) &
                        + Dz(k,7) * w(i,1,7,e) &
                        + Dz(k,8) * w(i,1,8,e) &
                        + Dz(k,9) * w(i,1,9,e) &
                        + Dz(k,10) * w(i,1,10,e) &
                        + Dz(k,11) * w(i,1,11,e) &
                        + Dz(k,12) * w(i,1,12,e)
          end do
       end do

       do i = 1, lx*lx*lx

          u1 = wur(i,1,1) * drdx(i,1,1,e) + wus(i,1,1) * dsdx(i,1,1,e) &
                                          + wut(i,1,1) * dtdx(i,1,1,e)
          u2 = wur(i,1,1) * drdy(i,1,1,e) + wus(i,1,1) * dsdy(i,1,1,e) &
                                          + wut(i,1,1) * dtdy(i,1,1,e)
          u3 = wur(i,1,1) * drdz(i,1,1,e) + wus(i,1,1) * dsdz(i,1,1,e) &
                                          + wut(i,1,1) * dtdz(i,1,1,e)

          v1 = wvr(i,1,1) * drdx(i,1,1,e) + wvs(i,1,1) * dsdx(i,1,1,e) &
                                          + wvt(i,1,1) * dtdx(i,1,1,e)
          v2 = wvr(i,1,1) * drdy(i,1,1,e) + wvs(i,1,1) * dsdy(i,1,1,e) &
                                          + wvt(i,1,1) * dtdy(i,1,1,e)
          v3 = wvr(i,1,1) * drdz(i,1,1,e) + wvs(i,1,1) * dsdz(i,1,1,e) &
                                          + wvt(i,1,1) * dtdz(i,1,1,e)

          w1 = wwr(i,1,1) * drdx(i,1,1,e) + wws(i,1,1) * dsdx(i,1,1,e) &
                                          + wwt(i,1,1) * dtdx(i,1,1,e)
          w2 = wwr(i,1,1) * drdy(i,1,1,e) + wws(i,1,1) * dsdy(i,1,1,e) &
                                          + wwt(i,1,1) * dtdy(i,1,1,e)
          w3 = wwr(i,1,1) * drdz(i,1,1,e) + wws(i,1,1) * dsdz(i,1,1,e) &
                                          + wwt(i,1,1) * dtdz(i,1,1,e)

          dj  = h1(i,1,1,e) * weights3(i,1,1) * jacinv(i,1,1,e)
          s11 = dj*(u1 + u1)
          s12 = dj*(u2 + v1)
          s13 = dj*(u3 + w1)
          s21 = dj*(v1 + u2)
          s22 = dj*(v2 + v2)
          s23 = dj*(v3 + w2)
          s31 = dj*(w1 + u3)
          s32 = dj*(w2 + v3)
          s33 = dj*(w3 + w3)

          wur(i,1,1) = drdx(i,1,1,e)*s11 + drdy(i,1,1,e)*s12 + drdz(i,1,1,e)*s13
          wus(i,1,1) = dsdx(i,1,1,e)*s11 + dsdy(i,1,1,e)*s12 + dsdz(i,1,1,e)*s13
          wut(i,1,1) = dtdx(i,1,1,e)*s11 + dtdy(i,1,1,e)*s12 + dtdz(i,1,1,e)*s13

          wvr(i,1,1) = drdx(i,1,1,e)*s21 + drdy(i,1,1,e)*s22 + drdz(i,1,1,e)*s23
          wvs(i,1,1) = dsdx(i,1,1,e)*s21 + dsdy(i,1,1,e)*s22 + dsdz(i,1,1,e)*s23
          wvt(i,1,1) = dtdx(i,1,1,e)*s21 + dtdy(i,1,1,e)*s22 + dtdz(i,1,1,e)*s23

          wwr(i,1,1) = drdx(i,1,1,e)*s31 + drdy(i,1,1,e)*s32 + drdz(i,1,1,e)*s33
          wws(i,1,1) = dsdx(i,1,1,e)*s31 + dsdy(i,1,1,e)*s32 + dsdz(i,1,1,e)*s33
          wwt(i,1,1) = dtdx(i,1,1,e)*s31 + dtdy(i,1,1,e)*s32 + dtdz(i,1,1,e)*s33
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * wur(1,j,1) &
                         + Dxt(i,2) * wur(2,j,1) &
                         + Dxt(i,3) * wur(3,j,1) &
                         + Dxt(i,4) * wur(4,j,1) &
                         + Dxt(i,5) * wur(5,j,1) &
                         + Dxt(i,6) * wur(6,j,1) &
                         + Dxt(i,7) * wur(7,j,1) &
                         + Dxt(i,8) * wur(8,j,1) &
                         + Dxt(i,9) * wur(9,j,1) &
                         + Dxt(i,10) * wur(10,j,1) &
                         + Dxt(i,11) * wur(11,j,1) &
                         + Dxt(i,12) * wur(12,j,1)

             av(i,j,1,e) = Dxt(i,1) * wvr(1,j,1) &
                         + Dxt(i,2) * wvr(2,j,1) &
                         + Dxt(i,3) * wvr(3,j,1) &
                         + Dxt(i,4) * wvr(4,j,1) &
                         + Dxt(i,5) * wvr(5,j,1) &
                         + Dxt(i,6) * wvr(6,j,1) &
                         + Dxt(i,7) * wvr(7,j,1) &
                         + Dxt(i,8) * wvr(8,j,1) &
                         + Dxt(i,9) * wvr(9,j,1) &
                         + Dxt(i,10) * wvr(10,j,1) &
                         + Dxt(i,11) * wvr(11,j,1) &
                         + Dxt(i,12) * wvr(12,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wwr(1,j,1) &
                         + Dxt(i,2) * wwr(2,j,1) &
                         + Dxt(i,3) * wwr(3,j,1) &
                         + Dxt(i,4) * wwr(4,j,1) &
                         + Dxt(i,5) * wwr(5,j,1) &
                         + Dxt(i,6) * wwr(6,j,1) &
                         + Dxt(i,7) * wwr(7,j,1) &
                         + Dxt(i,8) * wwr(8,j,1) &
                         + Dxt(i,9) * wwr(9,j,1) &
                         + Dxt(i,10) * wwr(10,j,1) &
                         + Dxt(i,11) * wwr(11,j,1) &
                         + Dxt(i,12) * wwr(12,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * wus(i,1,k) &
                            + Dyt(j,2) * wus(i,2,k) &
                            + Dyt(j,3) * wus(i,3,k) &
                            + Dyt(j,4) * wus(i,4,k) &
                            + Dyt(j,5) * wus(i,5,k) &
                            + Dyt(j,6) * wus(i,6,k) &
                            + Dyt(j,7) * wus(i,7,k) &
                            + Dyt(j,8) * wus(i,8,k) &
                            + Dyt(j,9) * wus(i,9,k) &
                            + Dyt(j,10) * wus(i,10,k) &
                            + Dyt(j,11) * wus(i,11,k) &
                            + Dyt(j,12) * wus(i,12,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * wvs(i,1,k) &
                            + Dyt(j,2) * wvs(i,2,k) &
                            + Dyt(j,3) * wvs(i,3,k) &
                            + Dyt(j,4) * wvs(i,4,k) &
                            + Dyt(j,5) * wvs(i,5,k) &
                            + Dyt(j,6) * wvs(i,6,k) &
                            + Dyt(j,7) * wvs(i,7,k) &
                            + Dyt(j,8) * wvs(i,8,k) &
                            + Dyt(j,9) * wvs(i,9,k) &
                            + Dyt(j,10) * wvs(i,10,k) &
                            + Dyt(j,11) * wvs(i,11,k) &
                            + Dyt(j,12) * wvs(i,12,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * wws(i,1,k) &
                            + Dyt(j,2) * wws(i,2,k) &
                            + Dyt(j,3) * wws(i,3,k) &
                            + Dyt(j,4) * wws(i,4,k) &
                            + Dyt(j,5) * wws(i,5,k) &
                            + Dyt(j,6) * wws(i,6,k) &
                            + Dyt(j,7) * wws(i,7,k) &
                            + Dyt(j,8) * wws(i,8,k) &
                            + Dyt(j,9) * wws(i,9,k) &
                            + Dyt(j,10) * wws(i,10,k) &
                            + Dyt(j,11) * wws(i,11,k) &
                            + Dyt(j,12) * wws(i,12,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * wut(i,1,1) &
                         + Dzt(k,2) * wut(i,1,2) &
                         + Dzt(k,3) * wut(i,1,3) &
                         + Dzt(k,4) * wut(i,1,4) &
                         + Dzt(k,5) * wut(i,1,5) &
                         + Dzt(k,6) * wut(i,1,6) &
                         + Dzt(k,7) * wut(i,1,7) &
                         + Dzt(k,8) * wut(i,1,8) &
                         + Dzt(k,9) * wut(i,1,9) &
                         + Dzt(k,10) * wut(i,1,10) &
                         + Dzt(k,11) * wut(i,1,11) &
                         + Dzt(k,12) * wut(i,1,12)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * wvt(i,1,1) &
                         + Dzt(k,2) * wvt(i,1,2) &
                         + Dzt(k,3) * wvt(i,1,3) &
                         + Dzt(k,4) * wvt(i,1,4) &
                         + Dzt(k,5) * wvt(i,1,5) &
                         + Dzt(k,6) * wvt(i,1,6) &
                         + Dzt(k,7) * wvt(i,1,7) &
                         + Dzt(k,8) * wvt(i,1,8) &
                         + Dzt(k,9) * wvt(i,1,9) &
                         + Dzt(k,10) * wvt(i,1,10) &
                         + Dzt(k,11) * wvt(i,1,11) &
                         + Dzt(k,12) * wvt(i,1,12)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wwt(i,1,1) &
                         + Dzt(k,2) * wwt(i,1,2) &
                         + Dzt(k,3) * wwt(i,1,3) &
                         + Dzt(k,4) * wwt(i,1,4) &
                         + Dzt(k,5) * wwt(i,1,5) &
                         + Dzt(k,6) * wwt(i,1,6) &
                         + Dzt(k,7) * wwt(i,1,7) &
                         + Dzt(k,8) * wwt(i,1,8) &
                         + Dzt(k,9) * wwt(i,1,9) &
                         + Dzt(k,10) * wwt(i,1,10) &
                         + Dzt(k,11) * wwt(i,1,11) &
                         + Dzt(k,12) * wwt(i,1,12)
          end do
       end do

    end do

  end subroutine ax_helm_stress_lx12

  subroutine ax_helm_stress_lx11(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, h2, drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
       jacinv, weights3, n)
    integer, parameter :: lx = 11
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: u(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: v(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdz(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: jacinv(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp), intent(in) :: weights3(lx, lx, lx)

    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    real(kind=rp) :: wvr(lx, lx, lx)
    real(kind=rp) :: wvs(lx, lx, lx)
    real(kind=rp) :: wvt(lx, lx, lx)
    real(kind=rp) :: wwr(lx, lx, lx)
    real(kind=rp) :: wws(lx, lx, lx)
    real(kind=rp) :: wwt(lx, lx, lx)

    integer :: e, i, j, k, l
    real(kind=rp) :: dj
    real(kind=rp) :: s11, s12, s13, s21, s22, s23, s31, s32, s33
    real(kind=rp) :: u1, u2, u3, v1, v2, v3, w1, w2, w3

    do e = 1, n

       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) &
                        + Dx(i,4) * u(4,j,1,e) &
                        + Dx(i,5) * u(5,j,1,e) &
                        + Dx(i,6) * u(6,j,1,e) &
                        + Dx(i,7) * u(7,j,1,e) &
                        + Dx(i,8) * u(8,j,1,e) &
                        + Dx(i,9) * u(9,j,1,e) &
                        + Dx(i,10) * u(10,j,1,e) &
                        + Dx(i,11) * u(11,j,1,e)

             wvr(i,j,1) = Dx(i,1) * v(1,j,1,e) &
                        + Dx(i,2) * v(2,j,1,e) &
                        + Dx(i,3) * v(3,j,1,e) &
                        + Dx(i,4) * v(4,j,1,e) &
                        + Dx(i,5) * v(5,j,1,e) &
                        + Dx(i,6) * v(6,j,1,e) &
                        + Dx(i,7) * v(7,j,1,e) &
                        + Dx(i,8) * v(8,j,1,e) &
                        + Dx(i,9) * v(9,j,1,e) &
                        + Dx(i,10) * v(10,j,1,e) &
                        + Dx(i,11) * v(11,j,1,e)

             wwr(i,j,1) = Dx(i,1) * w(1,j,1,e) &
                        + Dx(i,2) * w(2,j,1,e) &
                        + Dx(i,3) * w(3,j,1,e) &
                        + Dx(i,4) * w(4,j,1,e) &
                        + Dx(i,5) * w(5,j,1,e) &
                        + Dx(i,6) * w(6,j,1,e) &
                        + Dx(i,7) * w(7,j,1,e) &
                        + Dx(i,8) * w(8,j,1,e) &
                        + Dx(i,9) * w(9,j,1,e) &
                        + Dx(i,10) * w(10,j,1,e) &
                        + Dx(i,11) * w(11,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) &
                           + Dy(j,4) * u(i,4,k,e) &
                           + Dy(j,5) * u(i,5,k,e) &
                           + Dy(j,6) * u(i,6,k,e) &
                           + Dy(j,7) * u(i,7,k,e) &
                           + Dy(j,8) * u(i,8,k,e) &
                           + Dy(j,9) * u(i,9,k,e) &
                           + Dy(j,10) * u(i,10,k,e) &
                           + Dy(j,11) * u(i,11,k,e)

                wvs(i,j,k) = Dy(j,1) * v(i,1,k,e) &
                           + Dy(j,2) * v(i,2,k,e) &
                           + Dy(j,3) * v(i,3,k,e) &
                           + Dy(j,4) * v(i,4,k,e) &
                           + Dy(j,5) * v(i,5,k,e) &
                           + Dy(j,6) * v(i,6,k,e) &
                           + Dy(j,7) * v(i,7,k,e) &
                           + Dy(j,8) * v(i,8,k,e) &
                           + Dy(j,9) * v(i,9,k,e) &
                           + Dy(j,10) * v(i,10,k,e) &
                           + Dy(j,11) * v(i,11,k,e)

                wws(i,j,k) = Dy(j,1) * w(i,1,k,e) &
                           + Dy(j,2) * w(i,2,k,e) &
                           + Dy(j,3) * w(i,3,k,e) &
                           + Dy(j,4) * w(i,4,k,e) &
                           + Dy(j,5) * w(i,5,k,e) &
                           + Dy(j,6) * w(i,6,k,e) &
                           + Dy(j,7) * w(i,7,k,e) &
                           + Dy(j,8) * w(i,8,k,e) &
                           + Dy(j,9) * w(i,9,k,e) &
                           + Dy(j,10) * w(i,10,k,e) &
                           + Dy(j,11) * w(i,11,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) &
                        + Dz(k,4) * u(i,1,4,e) &
                        + Dz(k,5) * u(i,1,5,e) &
                        + Dz(k,6) * u(i,1,6,e) &
                        + Dz(k,7) * u(i,1,7,e) &
                        + Dz(k,8) * u(i,1,8,e) &
                        + Dz(k,9) * u(i,1,9,e) &
                        + Dz(k,10) * u(i,1,10,e) &
                        + Dz(k,11) * u(i,1,11,e)

             wvt(i,1,k) = Dz(k,1) * v(i,1,1,e) &
                        + Dz(k,2) * v(i,1,2,e) &
                        + Dz(k,3) * v(i,1,3,e) &
                        + Dz(k,4) * v(i,1,4,e) &
                        + Dz(k,5) * v(i,1,5,e) &
                        + Dz(k,6) * v(i,1,6,e) &
                        + Dz(k,7) * v(i,1,7,e) &
                        + Dz(k,8) * v(i,1,8,e) &
                        + Dz(k,9) * v(i,1,9,e) &
                        + Dz(k,10) * v(i,1,10,e) &
                        + Dz(k,11) * v(i,1,11,e)

             wwt(i,1,k) = Dz(k,1) * w(i,1,1,e) &
                        + Dz(k,2) * w(i,1,2,e) &
                        + Dz(k,3) * w(i,1,3,e) &
                        + Dz(k,4) * w(i,1,4,e) &
                        + Dz(k,5) * w(i,1,5,e) &
                        + Dz(k,6) * w(i,1,6,e) &
                        + Dz(k,7) * w(i,1,7,e) &
                        + Dz(k,8) * w(i,1,8,e) &
                        + Dz(k,9) * w(i,1,9,e) &
                        + Dz(k,10) * w(i,1,10,e) &
                        + Dz(k,11) * w(i,1,11,e)
          end do
       end do

       do i = 1, lx*lx*lx

          u1 = wur(i,1,1) * drdx(i,1,1,e) + wus(i,1,1) * dsdx(i,1,1,e) &
                                          + wut(i,1,1) * dtdx(i,1,1,e)
          u2 = wur(i,1,1) * drdy(i,1,1,e) + wus(i,1,1) * dsdy(i,1,1,e) &
                                          + wut(i,1,1) * dtdy(i,1,1,e)
          u3 = wur(i,1,1) * drdz(i,1,1,e) + wus(i,1,1) * dsdz(i,1,1,e) &
                                          + wut(i,1,1) * dtdz(i,1,1,e)

          v1 = wvr(i,1,1) * drdx(i,1,1,e) + wvs(i,1,1) * dsdx(i,1,1,e) &
                                          + wvt(i,1,1) * dtdx(i,1,1,e)
          v2 = wvr(i,1,1) * drdy(i,1,1,e) + wvs(i,1,1) * dsdy(i,1,1,e) &
                                          + wvt(i,1,1) * dtdy(i,1,1,e)
          v3 = wvr(i,1,1) * drdz(i,1,1,e) + wvs(i,1,1) * dsdz(i,1,1,e) &
                                          + wvt(i,1,1) * dtdz(i,1,1,e)

          w1 = wwr(i,1,1) * drdx(i,1,1,e) + wws(i,1,1) * dsdx(i,1,1,e) &
                                          + wwt(i,1,1) * dtdx(i,1,1,e)
          w2 = wwr(i,1,1) * drdy(i,1,1,e) + wws(i,1,1) * dsdy(i,1,1,e) &
                                          + wwt(i,1,1) * dtdy(i,1,1,e)
          w3 = wwr(i,1,1) * drdz(i,1,1,e) + wws(i,1,1) * dsdz(i,1,1,e) &
                                          + wwt(i,1,1) * dtdz(i,1,1,e)

          dj  = h1(i,1,1,e) * weights3(i,1,1) * jacinv(i,1,1,e)
          s11 = dj*(u1 + u1)
          s12 = dj*(u2 + v1)
          s13 = dj*(u3 + w1)
          s21 = dj*(v1 + u2)
          s22 = dj*(v2 + v2)
          s23 = dj*(v3 + w2)
          s31 = dj*(w1 + u3)
          s32 = dj*(w2 + v3)
          s33 = dj*(w3 + w3)

          wur(i,1,1) = drdx(i,1,1,e)*s11 + drdy(i,1,1,e)*s12 + drdz(i,1,1,e)*s13
          wus(i,1,1) = dsdx(i,1,1,e)*s11 + dsdy(i,1,1,e)*s12 + dsdz(i,1,1,e)*s13
          wut(i,1,1) = dtdx(i,1,1,e)*s11 + dtdy(i,1,1,e)*s12 + dtdz(i,1,1,e)*s13

          wvr(i,1,1) = drdx(i,1,1,e)*s21 + drdy(i,1,1,e)*s22 + drdz(i,1,1,e)*s23
          wvs(i,1,1) = dsdx(i,1,1,e)*s21 + dsdy(i,1,1,e)*s22 + dsdz(i,1,1,e)*s23
          wvt(i,1,1) = dtdx(i,1,1,e)*s21 + dtdy(i,1,1,e)*s22 + dtdz(i,1,1,e)*s23

          wwr(i,1,1) = drdx(i,1,1,e)*s31 + drdy(i,1,1,e)*s32 + drdz(i,1,1,e)*s33
          wws(i,1,1) = dsdx(i,1,1,e)*s31 + dsdy(i,1,1,e)*s32 + dsdz(i,1,1,e)*s33
          wwt(i,1,1) = dtdx(i,1,1,e)*s31 + dtdy(i,1,1,e)*s32 + dtdz(i,1,1,e)*s33
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * wur(1,j,1) &
                         + Dxt(i,2) * wur(2,j,1) &
                         + Dxt(i,3) * wur(3,j,1) &
                         + Dxt(i,4) * wur(4,j,1) &
                         + Dxt(i,5) * wur(5,j,1) &
                         + Dxt(i,6) * wur(6,j,1) &
                         + Dxt(i,7) * wur(7,j,1) &
                         + Dxt(i,8) * wur(8,j,1) &
                         + Dxt(i,9) * wur(9,j,1) &
                         + Dxt(i,10) * wur(10,j,1) &
                         + Dxt(i,11) * wur(11,j,1)

             av(i,j,1,e) = Dxt(i,1) * wvr(1,j,1) &
                         + Dxt(i,2) * wvr(2,j,1) &
                         + Dxt(i,3) * wvr(3,j,1) &
                         + Dxt(i,4) * wvr(4,j,1) &
                         + Dxt(i,5) * wvr(5,j,1) &
                         + Dxt(i,6) * wvr(6,j,1) &
                         + Dxt(i,7) * wvr(7,j,1) &
                         + Dxt(i,8) * wvr(8,j,1) &
                         + Dxt(i,9) * wvr(9,j,1) &
                         + Dxt(i,10) * wvr(10,j,1) &
                         + Dxt(i,11) * wvr(11,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wwr(1,j,1) &
                         + Dxt(i,2) * wwr(2,j,1) &
                         + Dxt(i,3) * wwr(3,j,1) &
                         + Dxt(i,4) * wwr(4,j,1) &
                         + Dxt(i,5) * wwr(5,j,1) &
                         + Dxt(i,6) * wwr(6,j,1) &
                         + Dxt(i,7) * wwr(7,j,1) &
                         + Dxt(i,8) * wwr(8,j,1) &
                         + Dxt(i,9) * wwr(9,j,1) &
                         + Dxt(i,10) * wwr(10,j,1) &
                         + Dxt(i,11) * wwr(11,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * wus(i,1,k) &
                            + Dyt(j,2) * wus(i,2,k) &
                            + Dyt(j,3) * wus(i,3,k) &
                            + Dyt(j,4) * wus(i,4,k) &
                            + Dyt(j,5) * wus(i,5,k) &
                            + Dyt(j,6) * wus(i,6,k) &
                            + Dyt(j,7) * wus(i,7,k) &
                            + Dyt(j,8) * wus(i,8,k) &
                            + Dyt(j,9) * wus(i,9,k) &
                            + Dyt(j,10) * wus(i,10,k) &
                            + Dyt(j,11) * wus(i,11,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * wvs(i,1,k) &
                            + Dyt(j,2) * wvs(i,2,k) &
                            + Dyt(j,3) * wvs(i,3,k) &
                            + Dyt(j,4) * wvs(i,4,k) &
                            + Dyt(j,5) * wvs(i,5,k) &
                            + Dyt(j,6) * wvs(i,6,k) &
                            + Dyt(j,7) * wvs(i,7,k) &
                            + Dyt(j,8) * wvs(i,8,k) &
                            + Dyt(j,9) * wvs(i,9,k) &
                            + Dyt(j,10) * wvs(i,10,k) &
                            + Dyt(j,11) * wvs(i,11,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * wws(i,1,k) &
                            + Dyt(j,2) * wws(i,2,k) &
                            + Dyt(j,3) * wws(i,3,k) &
                            + Dyt(j,4) * wws(i,4,k) &
                            + Dyt(j,5) * wws(i,5,k) &
                            + Dyt(j,6) * wws(i,6,k) &
                            + Dyt(j,7) * wws(i,7,k) &
                            + Dyt(j,8) * wws(i,8,k) &
                            + Dyt(j,9) * wws(i,9,k) &
                            + Dyt(j,10) * wws(i,10,k) &
                            + Dyt(j,11) * wws(i,11,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * wut(i,1,1) &
                         + Dzt(k,2) * wut(i,1,2) &
                         + Dzt(k,3) * wut(i,1,3) &
                         + Dzt(k,4) * wut(i,1,4) &
                         + Dzt(k,5) * wut(i,1,5) &
                         + Dzt(k,6) * wut(i,1,6) &
                         + Dzt(k,7) * wut(i,1,7) &
                         + Dzt(k,8) * wut(i,1,8) &
                         + Dzt(k,9) * wut(i,1,9) &
                         + Dzt(k,10) * wut(i,1,10) &
                         + Dzt(k,11) * wut(i,1,11)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * wvt(i,1,1) &
                         + Dzt(k,2) * wvt(i,1,2) &
                         + Dzt(k,3) * wvt(i,1,3) &
                         + Dzt(k,4) * wvt(i,1,4) &
                         + Dzt(k,5) * wvt(i,1,5) &
                         + Dzt(k,6) * wvt(i,1,6) &
                         + Dzt(k,7) * wvt(i,1,7) &
                         + Dzt(k,8) * wvt(i,1,8) &
                         + Dzt(k,9) * wvt(i,1,9) &
                         + Dzt(k,10) * wvt(i,1,10) &
                         + Dzt(k,11) * wvt(i,1,11)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wwt(i,1,1) &
                         + Dzt(k,2) * wwt(i,1,2) &
                         + Dzt(k,3) * wwt(i,1,3) &
                         + Dzt(k,4) * wwt(i,1,4) &
                         + Dzt(k,5) * wwt(i,1,5) &
                         + Dzt(k,6) * wwt(i,1,6) &
                         + Dzt(k,7) * wwt(i,1,7) &
                         + Dzt(k,8) * wwt(i,1,8) &
                         + Dzt(k,9) * wwt(i,1,9) &
                         + Dzt(k,10) * wwt(i,1,10) &
                         + Dzt(k,11) * wwt(i,1,11)
          end do
       end do

    end do

  end subroutine ax_helm_stress_lx11

  subroutine ax_helm_stress_lx10(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, h2, drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
       jacinv, weights3, n)
    integer, parameter :: lx = 10
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: u(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: v(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdz(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: jacinv(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp), intent(in) :: weights3(lx, lx, lx)

    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    real(kind=rp) :: wvr(lx, lx, lx)
    real(kind=rp) :: wvs(lx, lx, lx)
    real(kind=rp) :: wvt(lx, lx, lx)
    real(kind=rp) :: wwr(lx, lx, lx)
    real(kind=rp) :: wws(lx, lx, lx)
    real(kind=rp) :: wwt(lx, lx, lx)

    integer :: e, i, j, k, l
    real(kind=rp) :: dj
    real(kind=rp) :: s11, s12, s13, s21, s22, s23, s31, s32, s33
    real(kind=rp) :: u1, u2, u3, v1, v2, v3, w1, w2, w3

    do e = 1, n

       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) &
                        + Dx(i,4) * u(4,j,1,e) &
                        + Dx(i,5) * u(5,j,1,e) &
                        + Dx(i,6) * u(6,j,1,e) &
                        + Dx(i,7) * u(7,j,1,e) &
                        + Dx(i,8) * u(8,j,1,e) &
                        + Dx(i,9) * u(9,j,1,e) &
                        + Dx(i,10) * u(10,j,1,e)

             wvr(i,j,1) = Dx(i,1) * v(1,j,1,e) &
                        + Dx(i,2) * v(2,j,1,e) &
                        + Dx(i,3) * v(3,j,1,e) &
                        + Dx(i,4) * v(4,j,1,e) &
                        + Dx(i,5) * v(5,j,1,e) &
                        + Dx(i,6) * v(6,j,1,e) &
                        + Dx(i,7) * v(7,j,1,e) &
                        + Dx(i,8) * v(8,j,1,e) &
                        + Dx(i,9) * v(9,j,1,e) &
                        + Dx(i,10) * v(10,j,1,e)

             wwr(i,j,1) = Dx(i,1) * w(1,j,1,e) &
                        + Dx(i,2) * w(2,j,1,e) &
                        + Dx(i,3) * w(3,j,1,e) &
                        + Dx(i,4) * w(4,j,1,e) &
                        + Dx(i,5) * w(5,j,1,e) &
                        + Dx(i,6) * w(6,j,1,e) &
                        + Dx(i,7) * w(7,j,1,e) &
                        + Dx(i,8) * w(8,j,1,e) &
                        + Dx(i,9) * w(9,j,1,e) &
                        + Dx(i,10) * w(10,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) &
                           + Dy(j,4) * u(i,4,k,e) &
                           + Dy(j,5) * u(i,5,k,e) &
                           + Dy(j,6) * u(i,6,k,e) &
                           + Dy(j,7) * u(i,7,k,e) &
                           + Dy(j,8) * u(i,8,k,e) &
                           + Dy(j,9) * u(i,9,k,e) &
                           + Dy(j,10) * u(i,10,k,e)

                wvs(i,j,k) = Dy(j,1) * v(i,1,k,e) &
                           + Dy(j,2) * v(i,2,k,e) &
                           + Dy(j,3) * v(i,3,k,e) &
                           + Dy(j,4) * v(i,4,k,e) &
                           + Dy(j,5) * v(i,5,k,e) &
                           + Dy(j,6) * v(i,6,k,e) &
                           + Dy(j,7) * v(i,7,k,e) &
                           + Dy(j,8) * v(i,8,k,e) &
                           + Dy(j,9) * v(i,9,k,e) &
                           + Dy(j,10) * v(i,10,k,e)

                wws(i,j,k) = Dy(j,1) * w(i,1,k,e) &
                           + Dy(j,2) * w(i,2,k,e) &
                           + Dy(j,3) * w(i,3,k,e) &
                           + Dy(j,4) * w(i,4,k,e) &
                           + Dy(j,5) * w(i,5,k,e) &
                           + Dy(j,6) * w(i,6,k,e) &
                           + Dy(j,7) * w(i,7,k,e) &
                           + Dy(j,8) * w(i,8,k,e) &
                           + Dy(j,9) * w(i,9,k,e) &
                           + Dy(j,10) * w(i,10,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) &
                        + Dz(k,4) * u(i,1,4,e) &
                        + Dz(k,5) * u(i,1,5,e) &
                        + Dz(k,6) * u(i,1,6,e) &
                        + Dz(k,7) * u(i,1,7,e) &
                        + Dz(k,8) * u(i,1,8,e) &
                        + Dz(k,9) * u(i,1,9,e) &
                        + Dz(k,10) * u(i,1,10,e)

             wvt(i,1,k) = Dz(k,1) * v(i,1,1,e) &
                        + Dz(k,2) * v(i,1,2,e) &
                        + Dz(k,3) * v(i,1,3,e) &
                        + Dz(k,4) * v(i,1,4,e) &
                        + Dz(k,5) * v(i,1,5,e) &
                        + Dz(k,6) * v(i,1,6,e) &
                        + Dz(k,7) * v(i,1,7,e) &
                        + Dz(k,8) * v(i,1,8,e) &
                        + Dz(k,9) * v(i,1,9,e) &
                        + Dz(k,10) * v(i,1,10,e)

             wwt(i,1,k) = Dz(k,1) * w(i,1,1,e) &
                        + Dz(k,2) * w(i,1,2,e) &
                        + Dz(k,3) * w(i,1,3,e) &
                        + Dz(k,4) * w(i,1,4,e) &
                        + Dz(k,5) * w(i,1,5,e) &
                        + Dz(k,6) * w(i,1,6,e) &
                        + Dz(k,7) * w(i,1,7,e) &
                        + Dz(k,8) * w(i,1,8,e) &
                        + Dz(k,9) * w(i,1,9,e) &
                        + Dz(k,10) * w(i,1,10,e)
          end do
       end do

       do i = 1, lx*lx*lx

          u1 = wur(i,1,1) * drdx(i,1,1,e) + wus(i,1,1) * dsdx(i,1,1,e) &
                                          + wut(i,1,1) * dtdx(i,1,1,e)
          u2 = wur(i,1,1) * drdy(i,1,1,e) + wus(i,1,1) * dsdy(i,1,1,e) &
                                          + wut(i,1,1) * dtdy(i,1,1,e)
          u3 = wur(i,1,1) * drdz(i,1,1,e) + wus(i,1,1) * dsdz(i,1,1,e) &
                                          + wut(i,1,1) * dtdz(i,1,1,e)

          v1 = wvr(i,1,1) * drdx(i,1,1,e) + wvs(i,1,1) * dsdx(i,1,1,e) &
                                          + wvt(i,1,1) * dtdx(i,1,1,e)
          v2 = wvr(i,1,1) * drdy(i,1,1,e) + wvs(i,1,1) * dsdy(i,1,1,e) &
                                          + wvt(i,1,1) * dtdy(i,1,1,e)
          v3 = wvr(i,1,1) * drdz(i,1,1,e) + wvs(i,1,1) * dsdz(i,1,1,e) &
                                          + wvt(i,1,1) * dtdz(i,1,1,e)

          w1 = wwr(i,1,1) * drdx(i,1,1,e) + wws(i,1,1) * dsdx(i,1,1,e) &
                                          + wwt(i,1,1) * dtdx(i,1,1,e)
          w2 = wwr(i,1,1) * drdy(i,1,1,e) + wws(i,1,1) * dsdy(i,1,1,e) &
                                          + wwt(i,1,1) * dtdy(i,1,1,e)
          w3 = wwr(i,1,1) * drdz(i,1,1,e) + wws(i,1,1) * dsdz(i,1,1,e) &
                                          + wwt(i,1,1) * dtdz(i,1,1,e)

          dj  = h1(i,1,1,e) * weights3(i,1,1) * jacinv(i,1,1,e)
          s11 = dj*(u1 + u1)
          s12 = dj*(u2 + v1)
          s13 = dj*(u3 + w1)
          s21 = dj*(v1 + u2)
          s22 = dj*(v2 + v2)
          s23 = dj*(v3 + w2)
          s31 = dj*(w1 + u3)
          s32 = dj*(w2 + v3)
          s33 = dj*(w3 + w3)

          wur(i,1,1) = drdx(i,1,1,e)*s11 + drdy(i,1,1,e)*s12 + drdz(i,1,1,e)*s13
          wus(i,1,1) = dsdx(i,1,1,e)*s11 + dsdy(i,1,1,e)*s12 + dsdz(i,1,1,e)*s13
          wut(i,1,1) = dtdx(i,1,1,e)*s11 + dtdy(i,1,1,e)*s12 + dtdz(i,1,1,e)*s13

          wvr(i,1,1) = drdx(i,1,1,e)*s21 + drdy(i,1,1,e)*s22 + drdz(i,1,1,e)*s23
          wvs(i,1,1) = dsdx(i,1,1,e)*s21 + dsdy(i,1,1,e)*s22 + dsdz(i,1,1,e)*s23
          wvt(i,1,1) = dtdx(i,1,1,e)*s21 + dtdy(i,1,1,e)*s22 + dtdz(i,1,1,e)*s23

          wwr(i,1,1) = drdx(i,1,1,e)*s31 + drdy(i,1,1,e)*s32 + drdz(i,1,1,e)*s33
          wws(i,1,1) = dsdx(i,1,1,e)*s31 + dsdy(i,1,1,e)*s32 + dsdz(i,1,1,e)*s33
          wwt(i,1,1) = dtdx(i,1,1,e)*s31 + dtdy(i,1,1,e)*s32 + dtdz(i,1,1,e)*s33
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * wur(1,j,1) &
                         + Dxt(i,2) * wur(2,j,1) &
                         + Dxt(i,3) * wur(3,j,1) &
                         + Dxt(i,4) * wur(4,j,1) &
                         + Dxt(i,5) * wur(5,j,1) &
                         + Dxt(i,6) * wur(6,j,1) &
                         + Dxt(i,7) * wur(7,j,1) &
                         + Dxt(i,8) * wur(8,j,1) &
                         + Dxt(i,9) * wur(9,j,1) &
                         + Dxt(i,10) * wur(10,j,1)

             av(i,j,1,e) = Dxt(i,1) * wvr(1,j,1) &
                         + Dxt(i,2) * wvr(2,j,1) &
                         + Dxt(i,3) * wvr(3,j,1) &
                         + Dxt(i,4) * wvr(4,j,1) &
                         + Dxt(i,5) * wvr(5,j,1) &
                         + Dxt(i,6) * wvr(6,j,1) &
                         + Dxt(i,7) * wvr(7,j,1) &
                         + Dxt(i,8) * wvr(8,j,1) &
                         + Dxt(i,9) * wvr(9,j,1) &
                         + Dxt(i,10) * wvr(10,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wwr(1,j,1) &
                         + Dxt(i,2) * wwr(2,j,1) &
                         + Dxt(i,3) * wwr(3,j,1) &
                         + Dxt(i,4) * wwr(4,j,1) &
                         + Dxt(i,5) * wwr(5,j,1) &
                         + Dxt(i,6) * wwr(6,j,1) &
                         + Dxt(i,7) * wwr(7,j,1) &
                         + Dxt(i,8) * wwr(8,j,1) &
                         + Dxt(i,9) * wwr(9,j,1) &
                         + Dxt(i,10) * wwr(10,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * wus(i,1,k) &
                            + Dyt(j,2) * wus(i,2,k) &
                            + Dyt(j,3) * wus(i,3,k) &
                            + Dyt(j,4) * wus(i,4,k) &
                            + Dyt(j,5) * wus(i,5,k) &
                            + Dyt(j,6) * wus(i,6,k) &
                            + Dyt(j,7) * wus(i,7,k) &
                            + Dyt(j,8) * wus(i,8,k) &
                            + Dyt(j,9) * wus(i,9,k) &
                            + Dyt(j,10) * wus(i,10,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * wvs(i,1,k) &
                            + Dyt(j,2) * wvs(i,2,k) &
                            + Dyt(j,3) * wvs(i,3,k) &
                            + Dyt(j,4) * wvs(i,4,k) &
                            + Dyt(j,5) * wvs(i,5,k) &
                            + Dyt(j,6) * wvs(i,6,k) &
                            + Dyt(j,7) * wvs(i,7,k) &
                            + Dyt(j,8) * wvs(i,8,k) &
                            + Dyt(j,9) * wvs(i,9,k) &
                            + Dyt(j,10) * wvs(i,10,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * wws(i,1,k) &
                            + Dyt(j,2) * wws(i,2,k) &
                            + Dyt(j,3) * wws(i,3,k) &
                            + Dyt(j,4) * wws(i,4,k) &
                            + Dyt(j,5) * wws(i,5,k) &
                            + Dyt(j,6) * wws(i,6,k) &
                            + Dyt(j,7) * wws(i,7,k) &
                            + Dyt(j,8) * wws(i,8,k) &
                            + Dyt(j,9) * wws(i,9,k) &
                            + Dyt(j,10) * wws(i,10,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * wut(i,1,1) &
                         + Dzt(k,2) * wut(i,1,2) &
                         + Dzt(k,3) * wut(i,1,3) &
                         + Dzt(k,4) * wut(i,1,4) &
                         + Dzt(k,5) * wut(i,1,5) &
                         + Dzt(k,6) * wut(i,1,6) &
                         + Dzt(k,7) * wut(i,1,7) &
                         + Dzt(k,8) * wut(i,1,8) &
                         + Dzt(k,9) * wut(i,1,9) &
                         + Dzt(k,10) * wut(i,1,10)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * wvt(i,1,1) &
                         + Dzt(k,2) * wvt(i,1,2) &
                         + Dzt(k,3) * wvt(i,1,3) &
                         + Dzt(k,4) * wvt(i,1,4) &
                         + Dzt(k,5) * wvt(i,1,5) &
                         + Dzt(k,6) * wvt(i,1,6) &
                         + Dzt(k,7) * wvt(i,1,7) &
                         + Dzt(k,8) * wvt(i,1,8) &
                         + Dzt(k,9) * wvt(i,1,9) &
                         + Dzt(k,10) * wvt(i,1,10)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wwt(i,1,1) &
                         + Dzt(k,2) * wwt(i,1,2) &
                         + Dzt(k,3) * wwt(i,1,3) &
                         + Dzt(k,4) * wwt(i,1,4) &
                         + Dzt(k,5) * wwt(i,1,5) &
                         + Dzt(k,6) * wwt(i,1,6) &
                         + Dzt(k,7) * wwt(i,1,7) &
                         + Dzt(k,8) * wwt(i,1,8) &
                         + Dzt(k,9) * wwt(i,1,9) &
                         + Dzt(k,10) * wwt(i,1,10)
          end do
       end do

    end do

  end subroutine ax_helm_stress_lx10

  subroutine ax_helm_stress_lx9(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt,  &
       Dzt, h1, h2, drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
       jacinv, weights3, n)
    integer, parameter :: lx = 9
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: u(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: v(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdz(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: jacinv(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp), intent(in) :: weights3(lx, lx, lx)

    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    real(kind=rp) :: wvr(lx, lx, lx)
    real(kind=rp) :: wvs(lx, lx, lx)
    real(kind=rp) :: wvt(lx, lx, lx)
    real(kind=rp) :: wwr(lx, lx, lx)
    real(kind=rp) :: wws(lx, lx, lx)
    real(kind=rp) :: wwt(lx, lx, lx)

    integer :: e, i, j, k, l
    real(kind=rp) :: dj
    real(kind=rp) :: s11, s12, s13, s21, s22, s23, s31, s32, s33
    real(kind=rp) :: u1, u2, u3, v1, v2, v3, w1, w2, w3

    do e = 1, n

       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) &
                        + Dx(i,4) * u(4,j,1,e) &
                        + Dx(i,5) * u(5,j,1,e) &
                        + Dx(i,6) * u(6,j,1,e) &
                        + Dx(i,7) * u(7,j,1,e) &
                        + Dx(i,8) * u(8,j,1,e) &
                        + Dx(i,9) * u(9,j,1,e)

             wvr(i,j,1) = Dx(i,1) * v(1,j,1,e) &
                        + Dx(i,2) * v(2,j,1,e) &
                        + Dx(i,3) * v(3,j,1,e) &
                        + Dx(i,4) * v(4,j,1,e) &
                        + Dx(i,5) * v(5,j,1,e) &
                        + Dx(i,6) * v(6,j,1,e) &
                        + Dx(i,7) * v(7,j,1,e) &
                        + Dx(i,8) * v(8,j,1,e) &
                        + Dx(i,9) * v(9,j,1,e)

             wwr(i,j,1) = Dx(i,1) * w(1,j,1,e) &
                        + Dx(i,2) * w(2,j,1,e) &
                        + Dx(i,3) * w(3,j,1,e) &
                        + Dx(i,4) * w(4,j,1,e) &
                        + Dx(i,5) * w(5,j,1,e) &
                        + Dx(i,6) * w(6,j,1,e) &
                        + Dx(i,7) * w(7,j,1,e) &
                        + Dx(i,8) * w(8,j,1,e) &
                        + Dx(i,9) * w(9,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) &
                           + Dy(j,4) * u(i,4,k,e) &
                           + Dy(j,5) * u(i,5,k,e) &
                           + Dy(j,6) * u(i,6,k,e) &
                           + Dy(j,7) * u(i,7,k,e) &
                           + Dy(j,8) * u(i,8,k,e) &
                           + Dy(j,9) * u(i,9,k,e)

                wvs(i,j,k) = Dy(j,1) * v(i,1,k,e) &
                           + Dy(j,2) * v(i,2,k,e) &
                           + Dy(j,3) * v(i,3,k,e) &
                           + Dy(j,4) * v(i,4,k,e) &
                           + Dy(j,5) * v(i,5,k,e) &
                           + Dy(j,6) * v(i,6,k,e) &
                           + Dy(j,7) * v(i,7,k,e) &
                           + Dy(j,8) * v(i,8,k,e) &
                           + Dy(j,9) * v(i,9,k,e)

                wws(i,j,k) = Dy(j,1) * w(i,1,k,e) &
                           + Dy(j,2) * w(i,2,k,e) &
                           + Dy(j,3) * w(i,3,k,e) &
                           + Dy(j,4) * w(i,4,k,e) &
                           + Dy(j,5) * w(i,5,k,e) &
                           + Dy(j,6) * w(i,6,k,e) &
                           + Dy(j,7) * w(i,7,k,e) &
                           + Dy(j,8) * w(i,8,k,e) &
                           + Dy(j,9) * w(i,9,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) &
                        + Dz(k,4) * u(i,1,4,e) &
                        + Dz(k,5) * u(i,1,5,e) &
                        + Dz(k,6) * u(i,1,6,e) &
                        + Dz(k,7) * u(i,1,7,e) &
                        + Dz(k,8) * u(i,1,8,e) &
                        + Dz(k,9) * u(i,1,9,e)

             wvt(i,1,k) = Dz(k,1) * v(i,1,1,e) &
                        + Dz(k,2) * v(i,1,2,e) &
                        + Dz(k,3) * v(i,1,3,e) &
                        + Dz(k,4) * v(i,1,4,e) &
                        + Dz(k,5) * v(i,1,5,e) &
                        + Dz(k,6) * v(i,1,6,e) &
                        + Dz(k,7) * v(i,1,7,e) &
                        + Dz(k,8) * v(i,1,8,e) &
                        + Dz(k,9) * v(i,1,9,e)

             wwt(i,1,k) = Dz(k,1) * w(i,1,1,e) &
                        + Dz(k,2) * w(i,1,2,e) &
                        + Dz(k,3) * w(i,1,3,e) &
                        + Dz(k,4) * w(i,1,4,e) &
                        + Dz(k,5) * w(i,1,5,e) &
                        + Dz(k,6) * w(i,1,6,e) &
                        + Dz(k,7) * w(i,1,7,e) &
                        + Dz(k,8) * w(i,1,8,e) &
                        + Dz(k,9) * w(i,1,9,e)
          end do
       end do

       do i = 1, lx*lx*lx

          u1 = wur(i,1,1) * drdx(i,1,1,e) + wus(i,1,1) * dsdx(i,1,1,e) &
                                          + wut(i,1,1) * dtdx(i,1,1,e)
          u2 = wur(i,1,1) * drdy(i,1,1,e) + wus(i,1,1) * dsdy(i,1,1,e) &
                                          + wut(i,1,1) * dtdy(i,1,1,e)
          u3 = wur(i,1,1) * drdz(i,1,1,e) + wus(i,1,1) * dsdz(i,1,1,e) &
                                          + wut(i,1,1) * dtdz(i,1,1,e)

          v1 = wvr(i,1,1) * drdx(i,1,1,e) + wvs(i,1,1) * dsdx(i,1,1,e) &
                                          + wvt(i,1,1) * dtdx(i,1,1,e)
          v2 = wvr(i,1,1) * drdy(i,1,1,e) + wvs(i,1,1) * dsdy(i,1,1,e) &
                                          + wvt(i,1,1) * dtdy(i,1,1,e)
          v3 = wvr(i,1,1) * drdz(i,1,1,e) + wvs(i,1,1) * dsdz(i,1,1,e) &
                                          + wvt(i,1,1) * dtdz(i,1,1,e)

          w1 = wwr(i,1,1) * drdx(i,1,1,e) + wws(i,1,1) * dsdx(i,1,1,e) &
                                          + wwt(i,1,1) * dtdx(i,1,1,e)
          w2 = wwr(i,1,1) * drdy(i,1,1,e) + wws(i,1,1) * dsdy(i,1,1,e) &
                                          + wwt(i,1,1) * dtdy(i,1,1,e)
          w3 = wwr(i,1,1) * drdz(i,1,1,e) + wws(i,1,1) * dsdz(i,1,1,e) &
                                          + wwt(i,1,1) * dtdz(i,1,1,e)

          dj  = h1(i,1,1,e) * weights3(i,1,1) * jacinv(i,1,1,e)
          s11 = dj*(u1 + u1)
          s12 = dj*(u2 + v1)
          s13 = dj*(u3 + w1)
          s21 = dj*(v1 + u2)
          s22 = dj*(v2 + v2)
          s23 = dj*(v3 + w2)
          s31 = dj*(w1 + u3)
          s32 = dj*(w2 + v3)
          s33 = dj*(w3 + w3)

          wur(i,1,1) = drdx(i,1,1,e)*s11 + drdy(i,1,1,e)*s12 + drdz(i,1,1,e)*s13
          wus(i,1,1) = dsdx(i,1,1,e)*s11 + dsdy(i,1,1,e)*s12 + dsdz(i,1,1,e)*s13
          wut(i,1,1) = dtdx(i,1,1,e)*s11 + dtdy(i,1,1,e)*s12 + dtdz(i,1,1,e)*s13

          wvr(i,1,1) = drdx(i,1,1,e)*s21 + drdy(i,1,1,e)*s22 + drdz(i,1,1,e)*s23
          wvs(i,1,1) = dsdx(i,1,1,e)*s21 + dsdy(i,1,1,e)*s22 + dsdz(i,1,1,e)*s23
          wvt(i,1,1) = dtdx(i,1,1,e)*s21 + dtdy(i,1,1,e)*s22 + dtdz(i,1,1,e)*s23

          wwr(i,1,1) = drdx(i,1,1,e)*s31 + drdy(i,1,1,e)*s32 + drdz(i,1,1,e)*s33
          wws(i,1,1) = dsdx(i,1,1,e)*s31 + dsdy(i,1,1,e)*s32 + dsdz(i,1,1,e)*s33
          wwt(i,1,1) = dtdx(i,1,1,e)*s31 + dtdy(i,1,1,e)*s32 + dtdz(i,1,1,e)*s33
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * wur(1,j,1) &
                         + Dxt(i,2) * wur(2,j,1) &
                         + Dxt(i,3) * wur(3,j,1) &
                         + Dxt(i,4) * wur(4,j,1) &
                         + Dxt(i,5) * wur(5,j,1) &
                         + Dxt(i,6) * wur(6,j,1) &
                         + Dxt(i,7) * wur(7,j,1) &
                         + Dxt(i,8) * wur(8,j,1) &
                         + Dxt(i,9) * wur(9,j,1)

             av(i,j,1,e) = Dxt(i,1) * wvr(1,j,1) &
                         + Dxt(i,2) * wvr(2,j,1) &
                         + Dxt(i,3) * wvr(3,j,1) &
                         + Dxt(i,4) * wvr(4,j,1) &
                         + Dxt(i,5) * wvr(5,j,1) &
                         + Dxt(i,6) * wvr(6,j,1) &
                         + Dxt(i,7) * wvr(7,j,1) &
                         + Dxt(i,8) * wvr(8,j,1) &
                         + Dxt(i,9) * wvr(9,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wwr(1,j,1) &
                         + Dxt(i,2) * wwr(2,j,1) &
                         + Dxt(i,3) * wwr(3,j,1) &
                         + Dxt(i,4) * wwr(4,j,1) &
                         + Dxt(i,5) * wwr(5,j,1) &
                         + Dxt(i,6) * wwr(6,j,1) &
                         + Dxt(i,7) * wwr(7,j,1) &
                         + Dxt(i,8) * wwr(8,j,1) &
                         + Dxt(i,9) * wwr(9,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * wus(i,1,k) &
                            + Dyt(j,2) * wus(i,2,k) &
                            + Dyt(j,3) * wus(i,3,k) &
                            + Dyt(j,4) * wus(i,4,k) &
                            + Dyt(j,5) * wus(i,5,k) &
                            + Dyt(j,6) * wus(i,6,k) &
                            + Dyt(j,7) * wus(i,7,k) &
                            + Dyt(j,8) * wus(i,8,k) &
                            + Dyt(j,9) * wus(i,9,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * wvs(i,1,k) &
                            + Dyt(j,2) * wvs(i,2,k) &
                            + Dyt(j,3) * wvs(i,3,k) &
                            + Dyt(j,4) * wvs(i,4,k) &
                            + Dyt(j,5) * wvs(i,5,k) &
                            + Dyt(j,6) * wvs(i,6,k) &
                            + Dyt(j,7) * wvs(i,7,k) &
                            + Dyt(j,8) * wvs(i,8,k) &
                            + Dyt(j,9) * wvs(i,9,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * wws(i,1,k) &
                            + Dyt(j,2) * wws(i,2,k) &
                            + Dyt(j,3) * wws(i,3,k) &
                            + Dyt(j,4) * wws(i,4,k) &
                            + Dyt(j,5) * wws(i,5,k) &
                            + Dyt(j,6) * wws(i,6,k) &
                            + Dyt(j,7) * wws(i,7,k) &
                            + Dyt(j,8) * wws(i,8,k) &
                            + Dyt(j,9) * wws(i,9,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * wut(i,1,1) &
                         + Dzt(k,2) * wut(i,1,2) &
                         + Dzt(k,3) * wut(i,1,3) &
                         + Dzt(k,4) * wut(i,1,4) &
                         + Dzt(k,5) * wut(i,1,5) &
                         + Dzt(k,6) * wut(i,1,6) &
                         + Dzt(k,7) * wut(i,1,7) &
                         + Dzt(k,8) * wut(i,1,8) &
                         + Dzt(k,9) * wut(i,1,9)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * wvt(i,1,1) &
                         + Dzt(k,2) * wvt(i,1,2) &
                         + Dzt(k,3) * wvt(i,1,3) &
                         + Dzt(k,4) * wvt(i,1,4) &
                         + Dzt(k,5) * wvt(i,1,5) &
                         + Dzt(k,6) * wvt(i,1,6) &
                         + Dzt(k,7) * wvt(i,1,7) &
                         + Dzt(k,8) * wvt(i,1,8) &
                         + Dzt(k,9) * wvt(i,1,9)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wwt(i,1,1) &
                         + Dzt(k,2) * wwt(i,1,2) &
                         + Dzt(k,3) * wwt(i,1,3) &
                         + Dzt(k,4) * wwt(i,1,4) &
                         + Dzt(k,5) * wwt(i,1,5) &
                         + Dzt(k,6) * wwt(i,1,6) &
                         + Dzt(k,7) * wwt(i,1,7) &
                         + Dzt(k,8) * wwt(i,1,8) &
                         + Dzt(k,9) * wwt(i,1,9)
          end do
       end do

    end do

  end subroutine ax_helm_stress_lx9

  subroutine ax_helm_stress_lx8(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, &
       Dzt, h1, h2, drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
       jacinv, weights3, n)
    integer, parameter :: lx = 8
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: u(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: v(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdz(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: jacinv(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp), intent(in) :: weights3(lx, lx, lx)

    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    real(kind=rp) :: wvr(lx, lx, lx)
    real(kind=rp) :: wvs(lx, lx, lx)
    real(kind=rp) :: wvt(lx, lx, lx)
    real(kind=rp) :: wwr(lx, lx, lx)
    real(kind=rp) :: wws(lx, lx, lx)
    real(kind=rp) :: wwt(lx, lx, lx)

    integer :: e, i, j, k, l
    real(kind=rp) :: dj
    real(kind=rp) :: s11, s12, s13, s21, s22, s23, s31, s32, s33
    real(kind=rp) :: u1, u2, u3, v1, v2, v3, w1, w2, w3

    do e = 1, n

       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) &
                        + Dx(i,4) * u(4,j,1,e) &
                        + Dx(i,5) * u(5,j,1,e) &
                        + Dx(i,6) * u(6,j,1,e) &
                        + Dx(i,7) * u(7,j,1,e) &
                        + Dx(i,8) * u(8,j,1,e)

             wvr(i,j,1) = Dx(i,1) * v(1,j,1,e) &
                        + Dx(i,2) * v(2,j,1,e) &
                        + Dx(i,3) * v(3,j,1,e) &
                        + Dx(i,4) * v(4,j,1,e) &
                        + Dx(i,5) * v(5,j,1,e) &
                        + Dx(i,6) * v(6,j,1,e) &
                        + Dx(i,7) * v(7,j,1,e) &
                        + Dx(i,8) * v(8,j,1,e)

             wwr(i,j,1) = Dx(i,1) * w(1,j,1,e) &
                        + Dx(i,2) * w(2,j,1,e) &
                        + Dx(i,3) * w(3,j,1,e) &
                        + Dx(i,4) * w(4,j,1,e) &
                        + Dx(i,5) * w(5,j,1,e) &
                        + Dx(i,6) * w(6,j,1,e) &
                        + Dx(i,7) * w(7,j,1,e) &
                        + Dx(i,8) * w(8,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) &
                           + Dy(j,4) * u(i,4,k,e) &
                           + Dy(j,5) * u(i,5,k,e) &
                           + Dy(j,6) * u(i,6,k,e) &
                           + Dy(j,7) * u(i,7,k,e) &
                           + Dy(j,8) * u(i,8,k,e)

                wvs(i,j,k) = Dy(j,1) * v(i,1,k,e) &
                           + Dy(j,2) * v(i,2,k,e) &
                           + Dy(j,3) * v(i,3,k,e) &
                           + Dy(j,4) * v(i,4,k,e) &
                           + Dy(j,5) * v(i,5,k,e) &
                           + Dy(j,6) * v(i,6,k,e) &
                           + Dy(j,7) * v(i,7,k,e) &
                           + Dy(j,8) * v(i,8,k,e)

                wws(i,j,k) = Dy(j,1) * w(i,1,k,e) &
                           + Dy(j,2) * w(i,2,k,e) &
                           + Dy(j,3) * w(i,3,k,e) &
                           + Dy(j,4) * w(i,4,k,e) &
                           + Dy(j,5) * w(i,5,k,e) &
                           + Dy(j,6) * w(i,6,k,e) &
                           + Dy(j,7) * w(i,7,k,e) &
                           + Dy(j,8) * w(i,8,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) &
                        + Dz(k,4) * u(i,1,4,e) &
                        + Dz(k,5) * u(i,1,5,e) &
                        + Dz(k,6) * u(i,1,6,e) &
                        + Dz(k,7) * u(i,1,7,e) &
                        + Dz(k,8) * u(i,1,8,e)

             wvt(i,1,k) = Dz(k,1) * v(i,1,1,e) &
                        + Dz(k,2) * v(i,1,2,e) &
                        + Dz(k,3) * v(i,1,3,e) &
                        + Dz(k,4) * v(i,1,4,e) &
                        + Dz(k,5) * v(i,1,5,e) &
                        + Dz(k,6) * v(i,1,6,e) &
                        + Dz(k,7) * v(i,1,7,e) &
                        + Dz(k,8) * v(i,1,8,e)

             wwt(i,1,k) = Dz(k,1) * w(i,1,1,e) &
                        + Dz(k,2) * w(i,1,2,e) &
                        + Dz(k,3) * w(i,1,3,e) &
                        + Dz(k,4) * w(i,1,4,e) &
                        + Dz(k,5) * w(i,1,5,e) &
                        + Dz(k,6) * w(i,1,6,e) &
                        + Dz(k,7) * w(i,1,7,e) &
                        + Dz(k,8) * w(i,1,8,e)
          end do
       end do

       do i = 1, lx*lx*lx

          u1 = wur(i,1,1) * drdx(i,1,1,e) + wus(i,1,1) * dsdx(i,1,1,e) &
                                          + wut(i,1,1) * dtdx(i,1,1,e)
          u2 = wur(i,1,1) * drdy(i,1,1,e) + wus(i,1,1) * dsdy(i,1,1,e) &
                                          + wut(i,1,1) * dtdy(i,1,1,e)
          u3 = wur(i,1,1) * drdz(i,1,1,e) + wus(i,1,1) * dsdz(i,1,1,e) &
                                          + wut(i,1,1) * dtdz(i,1,1,e)

          v1 = wvr(i,1,1) * drdx(i,1,1,e) + wvs(i,1,1) * dsdx(i,1,1,e) &
                                          + wvt(i,1,1) * dtdx(i,1,1,e)
          v2 = wvr(i,1,1) * drdy(i,1,1,e) + wvs(i,1,1) * dsdy(i,1,1,e) &
                                          + wvt(i,1,1) * dtdy(i,1,1,e)
          v3 = wvr(i,1,1) * drdz(i,1,1,e) + wvs(i,1,1) * dsdz(i,1,1,e) &
                                          + wvt(i,1,1) * dtdz(i,1,1,e)

          w1 = wwr(i,1,1) * drdx(i,1,1,e) + wws(i,1,1) * dsdx(i,1,1,e) &
                                          + wwt(i,1,1) * dtdx(i,1,1,e)
          w2 = wwr(i,1,1) * drdy(i,1,1,e) + wws(i,1,1) * dsdy(i,1,1,e) &
                                          + wwt(i,1,1) * dtdy(i,1,1,e)
          w3 = wwr(i,1,1) * drdz(i,1,1,e) + wws(i,1,1) * dsdz(i,1,1,e) &
                                          + wwt(i,1,1) * dtdz(i,1,1,e)

          dj  = h1(i,1,1,e) * weights3(i,1,1) * jacinv(i,1,1,e)
          s11 = dj*(u1 + u1)
          s12 = dj*(u2 + v1)
          s13 = dj*(u3 + w1)
          s21 = dj*(v1 + u2)
          s22 = dj*(v2 + v2)
          s23 = dj*(v3 + w2)
          s31 = dj*(w1 + u3)
          s32 = dj*(w2 + v3)
          s33 = dj*(w3 + w3)

          wur(i,1,1) = drdx(i,1,1,e)*s11 + drdy(i,1,1,e)*s12 + drdz(i,1,1,e)*s13
          wus(i,1,1) = dsdx(i,1,1,e)*s11 + dsdy(i,1,1,e)*s12 + dsdz(i,1,1,e)*s13
          wut(i,1,1) = dtdx(i,1,1,e)*s11 + dtdy(i,1,1,e)*s12 + dtdz(i,1,1,e)*s13

          wvr(i,1,1) = drdx(i,1,1,e)*s21 + drdy(i,1,1,e)*s22 + drdz(i,1,1,e)*s23
          wvs(i,1,1) = dsdx(i,1,1,e)*s21 + dsdy(i,1,1,e)*s22 + dsdz(i,1,1,e)*s23
          wvt(i,1,1) = dtdx(i,1,1,e)*s21 + dtdy(i,1,1,e)*s22 + dtdz(i,1,1,e)*s23

          wwr(i,1,1) = drdx(i,1,1,e)*s31 + drdy(i,1,1,e)*s32 + drdz(i,1,1,e)*s33
          wws(i,1,1) = dsdx(i,1,1,e)*s31 + dsdy(i,1,1,e)*s32 + dsdz(i,1,1,e)*s33
          wwt(i,1,1) = dtdx(i,1,1,e)*s31 + dtdy(i,1,1,e)*s32 + dtdz(i,1,1,e)*s33
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * wur(1,j,1) &
                         + Dxt(i,2) * wur(2,j,1) &
                         + Dxt(i,3) * wur(3,j,1) &
                         + Dxt(i,4) * wur(4,j,1) &
                         + Dxt(i,5) * wur(5,j,1) &
                         + Dxt(i,6) * wur(6,j,1) &
                         + Dxt(i,7) * wur(7,j,1) &
                         + Dxt(i,8) * wur(8,j,1)

             av(i,j,1,e) = Dxt(i,1) * wvr(1,j,1) &
                         + Dxt(i,2) * wvr(2,j,1) &
                         + Dxt(i,3) * wvr(3,j,1) &
                         + Dxt(i,4) * wvr(4,j,1) &
                         + Dxt(i,5) * wvr(5,j,1) &
                         + Dxt(i,6) * wvr(6,j,1) &
                         + Dxt(i,7) * wvr(7,j,1) &
                         + Dxt(i,8) * wvr(8,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wwr(1,j,1) &
                         + Dxt(i,2) * wwr(2,j,1) &
                         + Dxt(i,3) * wwr(3,j,1) &
                         + Dxt(i,4) * wwr(4,j,1) &
                         + Dxt(i,5) * wwr(5,j,1) &
                         + Dxt(i,6) * wwr(6,j,1) &
                         + Dxt(i,7) * wwr(7,j,1) &
                         + Dxt(i,8) * wwr(8,j,1)

          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * wus(i,1,k) &
                            + Dyt(j,2) * wus(i,2,k) &
                            + Dyt(j,3) * wus(i,3,k) &
                            + Dyt(j,4) * wus(i,4,k) &
                            + Dyt(j,5) * wus(i,5,k) &
                            + Dyt(j,6) * wus(i,6,k) &
                            + Dyt(j,7) * wus(i,7,k) &
                            + Dyt(j,8) * wus(i,8,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * wvs(i,1,k) &
                            + Dyt(j,2) * wvs(i,2,k) &
                            + Dyt(j,3) * wvs(i,3,k) &
                            + Dyt(j,4) * wvs(i,4,k) &
                            + Dyt(j,5) * wvs(i,5,k) &
                            + Dyt(j,6) * wvs(i,6,k) &
                            + Dyt(j,7) * wvs(i,7,k) &
                            + Dyt(j,8) * wvs(i,8,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * wws(i,1,k) &
                            + Dyt(j,2) * wws(i,2,k) &
                            + Dyt(j,3) * wws(i,3,k) &
                            + Dyt(j,4) * wws(i,4,k) &
                            + Dyt(j,5) * wws(i,5,k) &
                            + Dyt(j,6) * wws(i,6,k) &
                            + Dyt(j,7) * wws(i,7,k) &
                            + Dyt(j,8) * wws(i,8,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * wut(i,1,1) &
                         + Dzt(k,2) * wut(i,1,2) &
                         + Dzt(k,3) * wut(i,1,3) &
                         + Dzt(k,4) * wut(i,1,4) &
                         + Dzt(k,5) * wut(i,1,5) &
                         + Dzt(k,6) * wut(i,1,6) &
                         + Dzt(k,7) * wut(i,1,7) &
                         + Dzt(k,8) * wut(i,1,8)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * wvt(i,1,1) &
                         + Dzt(k,2) * wvt(i,1,2) &
                         + Dzt(k,3) * wvt(i,1,3) &
                         + Dzt(k,4) * wvt(i,1,4) &
                         + Dzt(k,5) * wvt(i,1,5) &
                         + Dzt(k,6) * wvt(i,1,6) &
                         + Dzt(k,7) * wvt(i,1,7) &
                         + Dzt(k,8) * wvt(i,1,8)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wwt(i,1,1) &
                         + Dzt(k,2) * wwt(i,1,2) &
                         + Dzt(k,3) * wwt(i,1,3) &
                         + Dzt(k,4) * wwt(i,1,4) &
                         + Dzt(k,5) * wwt(i,1,5) &
                         + Dzt(k,6) * wwt(i,1,6) &
                         + Dzt(k,7) * wwt(i,1,7) &
                         + Dzt(k,8) * wwt(i,1,8)
          end do
       end do

    end do

  end subroutine ax_helm_stress_lx8

  subroutine ax_helm_stress_lx7(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, &
       Dzt, h1, h2, drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
       jacinv, weights3, n)
    integer, parameter :: lx = 7
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: u(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: v(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdz(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: jacinv(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp), intent(in) :: weights3(lx, lx, lx)

    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    real(kind=rp) :: wvr(lx, lx, lx)
    real(kind=rp) :: wvs(lx, lx, lx)
    real(kind=rp) :: wvt(lx, lx, lx)
    real(kind=rp) :: wwr(lx, lx, lx)
    real(kind=rp) :: wws(lx, lx, lx)
    real(kind=rp) :: wwt(lx, lx, lx)

    integer :: e, i, j, k, l
    real(kind=rp) :: dj
    real(kind=rp) :: s11, s12, s13, s21, s22, s23, s31, s32, s33
    real(kind=rp) :: u1, u2, u3, v1, v2, v3, w1, w2, w3

    do e = 1, n

       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) &
                        + Dx(i,4) * u(4,j,1,e) &
                        + Dx(i,5) * u(5,j,1,e) &
                        + Dx(i,6) * u(6,j,1,e) &
                        + Dx(i,7) * u(7,j,1,e)

             wvr(i,j,1) = Dx(i,1) * v(1,j,1,e) &
                        + Dx(i,2) * v(2,j,1,e) &
                        + Dx(i,3) * v(3,j,1,e) &
                        + Dx(i,4) * v(4,j,1,e) &
                        + Dx(i,5) * v(5,j,1,e) &
                        + Dx(i,6) * v(6,j,1,e) &
                        + Dx(i,7) * v(7,j,1,e)

             wwr(i,j,1) = Dx(i,1) * w(1,j,1,e) &
                        + Dx(i,2) * w(2,j,1,e) &
                        + Dx(i,3) * w(3,j,1,e) &
                        + Dx(i,4) * w(4,j,1,e) &
                        + Dx(i,5) * w(5,j,1,e) &
                        + Dx(i,6) * w(6,j,1,e) &
                        + Dx(i,7) * w(7,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) &
                           + Dy(j,4) * u(i,4,k,e) &
                           + Dy(j,5) * u(i,5,k,e) &
                           + Dy(j,6) * u(i,6,k,e) &
                           + Dy(j,7) * u(i,7,k,e)

                wvs(i,j,k) = Dy(j,1) * v(i,1,k,e) &
                           + Dy(j,2) * v(i,2,k,e) &
                           + Dy(j,3) * v(i,3,k,e) &
                           + Dy(j,4) * v(i,4,k,e) &
                           + Dy(j,5) * v(i,5,k,e) &
                           + Dy(j,6) * v(i,6,k,e) &
                           + Dy(j,7) * v(i,7,k,e)

                wws(i,j,k) = Dy(j,1) * w(i,1,k,e) &
                           + Dy(j,2) * w(i,2,k,e) &
                           + Dy(j,3) * w(i,3,k,e) &
                           + Dy(j,4) * w(i,4,k,e) &
                           + Dy(j,5) * w(i,5,k,e) &
                           + Dy(j,6) * w(i,6,k,e) &
                           + Dy(j,7) * w(i,7,k,e)

             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) &
                        + Dz(k,4) * u(i,1,4,e) &
                        + Dz(k,5) * u(i,1,5,e) &
                        + Dz(k,6) * u(i,1,6,e) &
                        + Dz(k,7) * u(i,1,7,e)

             wvt(i,1,k) = Dz(k,1) * v(i,1,1,e) &
                        + Dz(k,2) * v(i,1,2,e) &
                        + Dz(k,3) * v(i,1,3,e) &
                        + Dz(k,4) * v(i,1,4,e) &
                        + Dz(k,5) * v(i,1,5,e) &
                        + Dz(k,6) * v(i,1,6,e) &
                        + Dz(k,7) * v(i,1,7,e)

             wwt(i,1,k) = Dz(k,1) * w(i,1,1,e) &
                        + Dz(k,2) * w(i,1,2,e) &
                        + Dz(k,3) * w(i,1,3,e) &
                        + Dz(k,4) * w(i,1,4,e) &
                        + Dz(k,5) * w(i,1,5,e) &
                        + Dz(k,6) * w(i,1,6,e) &
                        + Dz(k,7) * w(i,1,7,e)
          end do
       end do

       do i = 1, lx*lx*lx

          u1 = wur(i,1,1) * drdx(i,1,1,e) + wus(i,1,1) * dsdx(i,1,1,e) &
                                          + wut(i,1,1) * dtdx(i,1,1,e)
          u2 = wur(i,1,1) * drdy(i,1,1,e) + wus(i,1,1) * dsdy(i,1,1,e) &
                                          + wut(i,1,1) * dtdy(i,1,1,e)
          u3 = wur(i,1,1) * drdz(i,1,1,e) + wus(i,1,1) * dsdz(i,1,1,e) &
                                          + wut(i,1,1) * dtdz(i,1,1,e)

          v1 = wvr(i,1,1) * drdx(i,1,1,e) + wvs(i,1,1) * dsdx(i,1,1,e) &
                                          + wvt(i,1,1) * dtdx(i,1,1,e)
          v2 = wvr(i,1,1) * drdy(i,1,1,e) + wvs(i,1,1) * dsdy(i,1,1,e) &
                                          + wvt(i,1,1) * dtdy(i,1,1,e)
          v3 = wvr(i,1,1) * drdz(i,1,1,e) + wvs(i,1,1) * dsdz(i,1,1,e) &
                                          + wvt(i,1,1) * dtdz(i,1,1,e)

          w1 = wwr(i,1,1) * drdx(i,1,1,e) + wws(i,1,1) * dsdx(i,1,1,e) &
                                          + wwt(i,1,1) * dtdx(i,1,1,e)
          w2 = wwr(i,1,1) * drdy(i,1,1,e) + wws(i,1,1) * dsdy(i,1,1,e) &
                                          + wwt(i,1,1) * dtdy(i,1,1,e)
          w3 = wwr(i,1,1) * drdz(i,1,1,e) + wws(i,1,1) * dsdz(i,1,1,e) &
                                          + wwt(i,1,1) * dtdz(i,1,1,e)

          dj  = h1(i,1,1,e) * weights3(i,1,1) * jacinv(i,1,1,e)
          s11 = dj*(u1 + u1)
          s12 = dj*(u2 + v1)
          s13 = dj*(u3 + w1)
          s21 = dj*(v1 + u2)
          s22 = dj*(v2 + v2)
          s23 = dj*(v3 + w2)
          s31 = dj*(w1 + u3)
          s32 = dj*(w2 + v3)
          s33 = dj*(w3 + w3)

          wur(i,1,1) = drdx(i,1,1,e)*s11 + drdy(i,1,1,e)*s12 + drdz(i,1,1,e)*s13
          wus(i,1,1) = dsdx(i,1,1,e)*s11 + dsdy(i,1,1,e)*s12 + dsdz(i,1,1,e)*s13
          wut(i,1,1) = dtdx(i,1,1,e)*s11 + dtdy(i,1,1,e)*s12 + dtdz(i,1,1,e)*s13

          wvr(i,1,1) = drdx(i,1,1,e)*s21 + drdy(i,1,1,e)*s22 + drdz(i,1,1,e)*s23
          wvs(i,1,1) = dsdx(i,1,1,e)*s21 + dsdy(i,1,1,e)*s22 + dsdz(i,1,1,e)*s23
          wvt(i,1,1) = dtdx(i,1,1,e)*s21 + dtdy(i,1,1,e)*s22 + dtdz(i,1,1,e)*s23

          wwr(i,1,1) = drdx(i,1,1,e)*s31 + drdy(i,1,1,e)*s32 + drdz(i,1,1,e)*s33
          wws(i,1,1) = dsdx(i,1,1,e)*s31 + dsdy(i,1,1,e)*s32 + dsdz(i,1,1,e)*s33
          wwt(i,1,1) = dtdx(i,1,1,e)*s31 + dtdy(i,1,1,e)*s32 + dtdz(i,1,1,e)*s33
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * wur(1,j,1) &
                         + Dxt(i,2) * wur(2,j,1) &
                         + Dxt(i,3) * wur(3,j,1) &
                         + Dxt(i,4) * wur(4,j,1) &
                         + Dxt(i,5) * wur(5,j,1) &
                         + Dxt(i,6) * wur(6,j,1) &
                         + Dxt(i,7) * wur(7,j,1)

             av(i,j,1,e) = Dxt(i,1) * wvr(1,j,1) &
                         + Dxt(i,2) * wvr(2,j,1) &
                         + Dxt(i,3) * wvr(3,j,1) &
                         + Dxt(i,4) * wvr(4,j,1) &
                         + Dxt(i,5) * wvr(5,j,1) &
                         + Dxt(i,6) * wvr(6,j,1) &
                         + Dxt(i,7) * wvr(7,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wwr(1,j,1) &
                         + Dxt(i,2) * wwr(2,j,1) &
                         + Dxt(i,3) * wwr(3,j,1) &
                         + Dxt(i,4) * wwr(4,j,1) &
                         + Dxt(i,5) * wwr(5,j,1) &
                         + Dxt(i,6) * wwr(6,j,1) &
                         + Dxt(i,7) * wwr(7,j,1)

          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * wus(i,1,k) &
                            + Dyt(j,2) * wus(i,2,k) &
                            + Dyt(j,3) * wus(i,3,k) &
                            + Dyt(j,4) * wus(i,4,k) &
                            + Dyt(j,5) * wus(i,5,k) &
                            + Dyt(j,6) * wus(i,6,k) &
                            + Dyt(j,7) * wus(i,7,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * wvs(i,1,k) &
                            + Dyt(j,2) * wvs(i,2,k) &
                            + Dyt(j,3) * wvs(i,3,k) &
                            + Dyt(j,4) * wvs(i,4,k) &
                            + Dyt(j,5) * wvs(i,5,k) &
                            + Dyt(j,6) * wvs(i,6,k) &
                            + Dyt(j,7) * wvs(i,7,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * wws(i,1,k) &
                            + Dyt(j,2) * wws(i,2,k) &
                            + Dyt(j,3) * wws(i,3,k) &
                            + Dyt(j,4) * wws(i,4,k) &
                            + Dyt(j,5) * wws(i,5,k) &
                            + Dyt(j,6) * wws(i,6,k) &
                            + Dyt(j,7) * wws(i,7,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * wut(i,1,1) &
                         + Dzt(k,2) * wut(i,1,2) &
                         + Dzt(k,3) * wut(i,1,3) &
                         + Dzt(k,4) * wut(i,1,4) &
                         + Dzt(k,5) * wut(i,1,5) &
                         + Dzt(k,6) * wut(i,1,6) &
                         + Dzt(k,7) * wut(i,1,7)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * wvt(i,1,1) &
                         + Dzt(k,2) * wvt(i,1,2) &
                         + Dzt(k,3) * wvt(i,1,3) &
                         + Dzt(k,4) * wvt(i,1,4) &
                         + Dzt(k,5) * wvt(i,1,5) &
                         + Dzt(k,6) * wvt(i,1,6) &
                         + Dzt(k,7) * wvt(i,1,7)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wwt(i,1,1) &
                         + Dzt(k,2) * wwt(i,1,2) &
                         + Dzt(k,3) * wwt(i,1,3) &
                         + Dzt(k,4) * wwt(i,1,4) &
                         + Dzt(k,5) * wwt(i,1,5) &
                         + Dzt(k,6) * wwt(i,1,6) &
                         + Dzt(k,7) * wwt(i,1,7)
          end do
       end do

    end do

  end subroutine ax_helm_stress_lx7

  subroutine ax_helm_stress_lx6(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, &
       Dzt, h1, h2, drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
       jacinv, weights3, n)
    integer, parameter :: lx = 6
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: u(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: v(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdz(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: jacinv(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp), intent(in) :: weights3(lx, lx, lx)

    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    real(kind=rp) :: wvr(lx, lx, lx)
    real(kind=rp) :: wvs(lx, lx, lx)
    real(kind=rp) :: wvt(lx, lx, lx)
    real(kind=rp) :: wwr(lx, lx, lx)
    real(kind=rp) :: wws(lx, lx, lx)
    real(kind=rp) :: wwt(lx, lx, lx)

    integer :: e, i, j, k, l
    real(kind=rp) :: dj
    real(kind=rp) :: s11, s12, s13, s21, s22, s23, s31, s32, s33
    real(kind=rp) :: u1, u2, u3, v1, v2, v3, w1, w2, w3

    do e = 1, n

       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) &
                        + Dx(i,4) * u(4,j,1,e) &
                        + Dx(i,5) * u(5,j,1,e) &
                        + Dx(i,6) * u(6,j,1,e)

             wvr(i,j,1) = Dx(i,1) * v(1,j,1,e) &
                        + Dx(i,2) * v(2,j,1,e) &
                        + Dx(i,3) * v(3,j,1,e) &
                        + Dx(i,4) * v(4,j,1,e) &
                        + Dx(i,5) * v(5,j,1,e) &
                        + Dx(i,6) * v(6,j,1,e)

             wwr(i,j,1) = Dx(i,1) * w(1,j,1,e) &
                        + Dx(i,2) * w(2,j,1,e) &
                        + Dx(i,3) * w(3,j,1,e) &
                        + Dx(i,4) * w(4,j,1,e) &
                        + Dx(i,5) * w(5,j,1,e) &
                        + Dx(i,6) * w(6,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) &
                           + Dy(j,4) * u(i,4,k,e) &
                           + Dy(j,5) * u(i,5,k,e) &
                           + Dy(j,6) * u(i,6,k,e)

                wvs(i,j,k) = Dy(j,1) * v(i,1,k,e) &
                           + Dy(j,2) * v(i,2,k,e) &
                           + Dy(j,3) * v(i,3,k,e) &
                           + Dy(j,4) * v(i,4,k,e) &
                           + Dy(j,5) * v(i,5,k,e) &
                           + Dy(j,6) * v(i,6,k,e)

                wws(i,j,k) = Dy(j,1) * w(i,1,k,e) &
                           + Dy(j,2) * w(i,2,k,e) &
                           + Dy(j,3) * w(i,3,k,e) &
                           + Dy(j,4) * w(i,4,k,e) &
                           + Dy(j,5) * w(i,5,k,e) &
                           + Dy(j,6) * w(i,6,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) &
                        + Dz(k,4) * u(i,1,4,e) &
                        + Dz(k,5) * u(i,1,5,e) &
                        + Dz(k,6) * u(i,1,6,e)

             wvt(i,1,k) = Dz(k,1) * v(i,1,1,e) &
                        + Dz(k,2) * v(i,1,2,e) &
                        + Dz(k,3) * v(i,1,3,e) &
                        + Dz(k,4) * v(i,1,4,e) &
                        + Dz(k,5) * v(i,1,5,e) &
                        + Dz(k,6) * v(i,1,6,e)

             wwt(i,1,k) = Dz(k,1) * w(i,1,1,e) &
                        + Dz(k,2) * w(i,1,2,e) &
                        + Dz(k,3) * w(i,1,3,e) &
                        + Dz(k,4) * w(i,1,4,e) &
                        + Dz(k,5) * w(i,1,5,e) &
                        + Dz(k,6) * w(i,1,6,e)
          end do
       end do

       do i = 1, lx*lx*lx

          u1 = wur(i,1,1) * drdx(i,1,1,e) + wus(i,1,1) * dsdx(i,1,1,e) &
                                          + wut(i,1,1) * dtdx(i,1,1,e)
          u2 = wur(i,1,1) * drdy(i,1,1,e) + wus(i,1,1) * dsdy(i,1,1,e) &
                                          + wut(i,1,1) * dtdy(i,1,1,e)
          u3 = wur(i,1,1) * drdz(i,1,1,e) + wus(i,1,1) * dsdz(i,1,1,e) &
                                          + wut(i,1,1) * dtdz(i,1,1,e)

          v1 = wvr(i,1,1) * drdx(i,1,1,e) + wvs(i,1,1) * dsdx(i,1,1,e) &
                                          + wvt(i,1,1) * dtdx(i,1,1,e)
          v2 = wvr(i,1,1) * drdy(i,1,1,e) + wvs(i,1,1) * dsdy(i,1,1,e) &
                                          + wvt(i,1,1) * dtdy(i,1,1,e)
          v3 = wvr(i,1,1) * drdz(i,1,1,e) + wvs(i,1,1) * dsdz(i,1,1,e) &
                                          + wvt(i,1,1) * dtdz(i,1,1,e)

          w1 = wwr(i,1,1) * drdx(i,1,1,e) + wws(i,1,1) * dsdx(i,1,1,e) &
                                          + wwt(i,1,1) * dtdx(i,1,1,e)
          w2 = wwr(i,1,1) * drdy(i,1,1,e) + wws(i,1,1) * dsdy(i,1,1,e) &
                                          + wwt(i,1,1) * dtdy(i,1,1,e)
          w3 = wwr(i,1,1) * drdz(i,1,1,e) + wws(i,1,1) * dsdz(i,1,1,e) &
                                          + wwt(i,1,1) * dtdz(i,1,1,e)

          dj  = h1(i,1,1,e) * weights3(i,1,1) * jacinv(i,1,1,e)
          s11 = dj*(u1 + u1)
          s12 = dj*(u2 + v1)
          s13 = dj*(u3 + w1)
          s21 = dj*(v1 + u2)
          s22 = dj*(v2 + v2)
          s23 = dj*(v3 + w2)
          s31 = dj*(w1 + u3)
          s32 = dj*(w2 + v3)
          s33 = dj*(w3 + w3)

          wur(i,1,1) = drdx(i,1,1,e)*s11 + drdy(i,1,1,e)*s12 + drdz(i,1,1,e)*s13
          wus(i,1,1) = dsdx(i,1,1,e)*s11 + dsdy(i,1,1,e)*s12 + dsdz(i,1,1,e)*s13
          wut(i,1,1) = dtdx(i,1,1,e)*s11 + dtdy(i,1,1,e)*s12 + dtdz(i,1,1,e)*s13

          wvr(i,1,1) = drdx(i,1,1,e)*s21 + drdy(i,1,1,e)*s22 + drdz(i,1,1,e)*s23
          wvs(i,1,1) = dsdx(i,1,1,e)*s21 + dsdy(i,1,1,e)*s22 + dsdz(i,1,1,e)*s23
          wvt(i,1,1) = dtdx(i,1,1,e)*s21 + dtdy(i,1,1,e)*s22 + dtdz(i,1,1,e)*s23

          wwr(i,1,1) = drdx(i,1,1,e)*s31 + drdy(i,1,1,e)*s32 + drdz(i,1,1,e)*s33
          wws(i,1,1) = dsdx(i,1,1,e)*s31 + dsdy(i,1,1,e)*s32 + dsdz(i,1,1,e)*s33
          wwt(i,1,1) = dtdx(i,1,1,e)*s31 + dtdy(i,1,1,e)*s32 + dtdz(i,1,1,e)*s33
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * wur(1,j,1) &
                         + Dxt(i,2) * wur(2,j,1) &
                         + Dxt(i,3) * wur(3,j,1) &
                         + Dxt(i,4) * wur(4,j,1) &
                         + Dxt(i,5) * wur(5,j,1) &
                         + Dxt(i,6) * wur(6,j,1)

             av(i,j,1,e) = Dxt(i,1) * wvr(1,j,1) &
                         + Dxt(i,2) * wvr(2,j,1) &
                         + Dxt(i,3) * wvr(3,j,1) &
                         + Dxt(i,4) * wvr(4,j,1) &
                         + Dxt(i,5) * wvr(5,j,1) &
                         + Dxt(i,6) * wvr(6,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wwr(1,j,1) &
                         + Dxt(i,2) * wwr(2,j,1) &
                         + Dxt(i,3) * wwr(3,j,1) &
                         + Dxt(i,4) * wwr(4,j,1) &
                         + Dxt(i,5) * wwr(5,j,1) &
                         + Dxt(i,6) * wwr(6,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * wus(i,1,k) &
                            + Dyt(j,2) * wus(i,2,k) &
                            + Dyt(j,3) * wus(i,3,k) &
                            + Dyt(j,4) * wus(i,4,k) &
                            + Dyt(j,5) * wus(i,5,k) &
                            + Dyt(j,6) * wus(i,6,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * wvs(i,1,k) &
                            + Dyt(j,2) * wvs(i,2,k) &
                            + Dyt(j,3) * wvs(i,3,k) &
                            + Dyt(j,4) * wvs(i,4,k) &
                            + Dyt(j,5) * wvs(i,5,k) &
                            + Dyt(j,6) * wvs(i,6,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * wws(i,1,k) &
                            + Dyt(j,2) * wws(i,2,k) &
                            + Dyt(j,3) * wws(i,3,k) &
                            + Dyt(j,4) * wws(i,4,k) &
                            + Dyt(j,5) * wws(i,5,k) &
                            + Dyt(j,6) * wws(i,6,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * wut(i,1,1) &
                         + Dzt(k,2) * wut(i,1,2) &
                         + Dzt(k,3) * wut(i,1,3) &
                         + Dzt(k,4) * wut(i,1,4) &
                         + Dzt(k,5) * wut(i,1,5) &
                         + Dzt(k,6) * wut(i,1,6)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * wvt(i,1,1) &
                         + Dzt(k,2) * wvt(i,1,2) &
                         + Dzt(k,3) * wvt(i,1,3) &
                         + Dzt(k,4) * wvt(i,1,4) &
                         + Dzt(k,5) * wvt(i,1,5) &
                         + Dzt(k,6) * wvt(i,1,6)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wwt(i,1,1) &
                         + Dzt(k,2) * wwt(i,1,2) &
                         + Dzt(k,3) * wwt(i,1,3) &
                         + Dzt(k,4) * wwt(i,1,4) &
                         + Dzt(k,5) * wwt(i,1,5) &
                         + Dzt(k,6) * wwt(i,1,6)
          end do
       end do

    end do

  end subroutine ax_helm_stress_lx6

  subroutine ax_helm_stress_lx5(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, &
       Dzt, h1, h2, drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
       jacinv, weights3, n)
    integer, parameter :: lx = 5
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: u(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: v(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdz(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: jacinv(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp), intent(in) :: weights3(lx, lx, lx)

    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    real(kind=rp) :: wvr(lx, lx, lx)
    real(kind=rp) :: wvs(lx, lx, lx)
    real(kind=rp) :: wvt(lx, lx, lx)
    real(kind=rp) :: wwr(lx, lx, lx)
    real(kind=rp) :: wws(lx, lx, lx)
    real(kind=rp) :: wwt(lx, lx, lx)

    integer :: e, i, j, k, l
    real(kind=rp) :: dj
    real(kind=rp) :: s11, s12, s13, s21, s22, s23, s31, s32, s33
    real(kind=rp) :: u1, u2, u3, v1, v2, v3, w1, w2, w3

    do e = 1, n

       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) &
                        + Dx(i,4) * u(4,j,1,e) &
                        + Dx(i,5) * u(5,j,1,e)

             wvr(i,j,1) = Dx(i,1) * v(1,j,1,e) &
                        + Dx(i,2) * v(2,j,1,e) &
                        + Dx(i,3) * v(3,j,1,e) &
                        + Dx(i,4) * v(4,j,1,e) &
                        + Dx(i,5) * v(5,j,1,e)

             wwr(i,j,1) = Dx(i,1) * w(1,j,1,e) &
                        + Dx(i,2) * w(2,j,1,e) &
                        + Dx(i,3) * w(3,j,1,e) &
                        + Dx(i,4) * w(4,j,1,e) &
                        + Dx(i,5) * w(5,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) &
                           + Dy(j,4) * u(i,4,k,e) &
                           + Dy(j,5) * u(i,5,k,e)

                wvs(i,j,k) = Dy(j,1) * v(i,1,k,e) &
                           + Dy(j,2) * v(i,2,k,e) &
                           + Dy(j,3) * v(i,3,k,e) &
                           + Dy(j,4) * v(i,4,k,e) &
                           + Dy(j,5) * v(i,5,k,e)

                wws(i,j,k) = Dy(j,1) * w(i,1,k,e) &
                           + Dy(j,2) * w(i,2,k,e) &
                           + Dy(j,3) * w(i,3,k,e) &
                           + Dy(j,4) * w(i,4,k,e) &
                           + Dy(j,5) * w(i,5,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) &
                        + Dz(k,4) * u(i,1,4,e) &
                        + Dz(k,5) * u(i,1,5,e)

             wvt(i,1,k) = Dz(k,1) * v(i,1,1,e) &
                        + Dz(k,2) * v(i,1,2,e) &
                        + Dz(k,3) * v(i,1,3,e) &
                        + Dz(k,4) * v(i,1,4,e) &
                        + Dz(k,5) * v(i,1,5,e)

             wwt(i,1,k) = Dz(k,1) * w(i,1,1,e) &
                        + Dz(k,2) * w(i,1,2,e) &
                        + Dz(k,3) * w(i,1,3,e) &
                        + Dz(k,4) * w(i,1,4,e) &
                        + Dz(k,5) * w(i,1,5,e)
          end do
       end do

       do i = 1, lx*lx*lx

          u1 = wur(i,1,1) * drdx(i,1,1,e) + wus(i,1,1) * dsdx(i,1,1,e) &
                                          + wut(i,1,1) * dtdx(i,1,1,e)
          u2 = wur(i,1,1) * drdy(i,1,1,e) + wus(i,1,1) * dsdy(i,1,1,e) &
                                          + wut(i,1,1) * dtdy(i,1,1,e)
          u3 = wur(i,1,1) * drdz(i,1,1,e) + wus(i,1,1) * dsdz(i,1,1,e) &
                                          + wut(i,1,1) * dtdz(i,1,1,e)

          v1 = wvr(i,1,1) * drdx(i,1,1,e) + wvs(i,1,1) * dsdx(i,1,1,e) &
                                          + wvt(i,1,1) * dtdx(i,1,1,e)
          v2 = wvr(i,1,1) * drdy(i,1,1,e) + wvs(i,1,1) * dsdy(i,1,1,e) &
                                          + wvt(i,1,1) * dtdy(i,1,1,e)
          v3 = wvr(i,1,1) * drdz(i,1,1,e) + wvs(i,1,1) * dsdz(i,1,1,e) &
                                          + wvt(i,1,1) * dtdz(i,1,1,e)

          w1 = wwr(i,1,1) * drdx(i,1,1,e) + wws(i,1,1) * dsdx(i,1,1,e) &
                                          + wwt(i,1,1) * dtdx(i,1,1,e)
          w2 = wwr(i,1,1) * drdy(i,1,1,e) + wws(i,1,1) * dsdy(i,1,1,e) &
                                          + wwt(i,1,1) * dtdy(i,1,1,e)
          w3 = wwr(i,1,1) * drdz(i,1,1,e) + wws(i,1,1) * dsdz(i,1,1,e) &
                                          + wwt(i,1,1) * dtdz(i,1,1,e)

          dj  = h1(i,1,1,e) * weights3(i,1,1) * jacinv(i,1,1,e)
          s11 = dj*(u1 + u1)
          s12 = dj*(u2 + v1)
          s13 = dj*(u3 + w1)
          s21 = dj*(v1 + u2)
          s22 = dj*(v2 + v2)
          s23 = dj*(v3 + w2)
          s31 = dj*(w1 + u3)
          s32 = dj*(w2 + v3)
          s33 = dj*(w3 + w3)

          wur(i,1,1) = drdx(i,1,1,e)*s11 + drdy(i,1,1,e)*s12 + drdz(i,1,1,e)*s13
          wus(i,1,1) = dsdx(i,1,1,e)*s11 + dsdy(i,1,1,e)*s12 + dsdz(i,1,1,e)*s13
          wut(i,1,1) = dtdx(i,1,1,e)*s11 + dtdy(i,1,1,e)*s12 + dtdz(i,1,1,e)*s13

          wvr(i,1,1) = drdx(i,1,1,e)*s21 + drdy(i,1,1,e)*s22 + drdz(i,1,1,e)*s23
          wvs(i,1,1) = dsdx(i,1,1,e)*s21 + dsdy(i,1,1,e)*s22 + dsdz(i,1,1,e)*s23
          wvt(i,1,1) = dtdx(i,1,1,e)*s21 + dtdy(i,1,1,e)*s22 + dtdz(i,1,1,e)*s23

          wwr(i,1,1) = drdx(i,1,1,e)*s31 + drdy(i,1,1,e)*s32 + drdz(i,1,1,e)*s33
          wws(i,1,1) = dsdx(i,1,1,e)*s31 + dsdy(i,1,1,e)*s32 + dsdz(i,1,1,e)*s33
          wwt(i,1,1) = dtdx(i,1,1,e)*s31 + dtdy(i,1,1,e)*s32 + dtdz(i,1,1,e)*s33
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * wur(1,j,1) &
                         + Dxt(i,2) * wur(2,j,1) &
                         + Dxt(i,3) * wur(3,j,1) &
                         + Dxt(i,4) * wur(4,j,1) &
                         + Dxt(i,5) * wur(5,j,1)

             av(i,j,1,e) = Dxt(i,1) * wvr(1,j,1) &
                         + Dxt(i,2) * wvr(2,j,1) &
                         + Dxt(i,3) * wvr(3,j,1) &
                         + Dxt(i,4) * wvr(4,j,1) &
                         + Dxt(i,5) * wvr(5,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wwr(1,j,1) &
                         + Dxt(i,2) * wwr(2,j,1) &
                         + Dxt(i,3) * wwr(3,j,1) &
                         + Dxt(i,4) * wwr(4,j,1) &
                         + Dxt(i,5) * wwr(5,j,1)

          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * wus(i,1,k) &
                            + Dyt(j,2) * wus(i,2,k) &
                            + Dyt(j,3) * wus(i,3,k) &
                            + Dyt(j,4) * wus(i,4,k) &
                            + Dyt(j,5) * wus(i,5,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * wvs(i,1,k) &
                            + Dyt(j,2) * wvs(i,2,k) &
                            + Dyt(j,3) * wvs(i,3,k) &
                            + Dyt(j,4) * wvs(i,4,k) &
                            + Dyt(j,5) * wvs(i,5,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * wws(i,1,k) &
                            + Dyt(j,2) * wws(i,2,k) &
                            + Dyt(j,3) * wws(i,3,k) &
                            + Dyt(j,4) * wws(i,4,k) &
                            + Dyt(j,5) * wws(i,5,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * wut(i,1,1) &
                         + Dzt(k,2) * wut(i,1,2) &
                         + Dzt(k,3) * wut(i,1,3) &
                         + Dzt(k,4) * wut(i,1,4) &
                         + Dzt(k,5) * wut(i,1,5)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * wvt(i,1,1) &
                         + Dzt(k,2) * wvt(i,1,2) &
                         + Dzt(k,3) * wvt(i,1,3) &
                         + Dzt(k,4) * wvt(i,1,4) &
                         + Dzt(k,5) * wvt(i,1,5)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wwt(i,1,1) &
                         + Dzt(k,2) * wwt(i,1,2) &
                         + Dzt(k,3) * wwt(i,1,3) &
                         + Dzt(k,4) * wwt(i,1,4) &
                         + Dzt(k,5) * wwt(i,1,5)
          end do
       end do

    end do

  end subroutine ax_helm_stress_lx5

  subroutine ax_helm_stress_lx4(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, &
       Dzt, h1, h2, drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
       jacinv, weights3, n)
    integer, parameter :: lx = 4
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: u(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: v(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdz(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: jacinv(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp), intent(in) :: weights3(lx, lx, lx)

    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    real(kind=rp) :: wvr(lx, lx, lx)
    real(kind=rp) :: wvs(lx, lx, lx)
    real(kind=rp) :: wvt(lx, lx, lx)
    real(kind=rp) :: wwr(lx, lx, lx)
    real(kind=rp) :: wws(lx, lx, lx)
    real(kind=rp) :: wwt(lx, lx, lx)

    integer :: e, i, j, k, l
    real(kind=rp) :: dj
    real(kind=rp) :: s11, s12, s13, s21, s22, s23, s31, s32, s33
    real(kind=rp) :: u1, u2, u3, v1, v2, v3, w1, w2, w3

    do e = 1, n

       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e) &
                        + Dx(i,4) * u(4,j,1,e)

             wvr(i,j,1) = Dx(i,1) * v(1,j,1,e) &
                        + Dx(i,2) * v(2,j,1,e) &
                        + Dx(i,3) * v(3,j,1,e) &
                        + Dx(i,4) * v(4,j,1,e)

             wwr(i,j,1) = Dx(i,1) * w(1,j,1,e) &
                        + Dx(i,2) * w(2,j,1,e) &
                        + Dx(i,3) * w(3,j,1,e) &
                        + Dx(i,4) * w(4,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e) &
                           + Dy(j,4) * u(i,4,k,e)

                wvs(i,j,k) = Dy(j,1) * v(i,1,k,e) &
                           + Dy(j,2) * v(i,2,k,e) &
                           + Dy(j,3) * v(i,3,k,e) &
                           + Dy(j,4) * v(i,4,k,e)

                wws(i,j,k) = Dy(j,1) * w(i,1,k,e) &
                           + Dy(j,2) * w(i,2,k,e) &
                           + Dy(j,3) * w(i,3,k,e) &
                           + Dy(j,4) * w(i,4,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e) &
                        + Dz(k,4) * u(i,1,4,e)

             wvt(i,1,k) = Dz(k,1) * v(i,1,1,e) &
                        + Dz(k,2) * v(i,1,2,e) &
                        + Dz(k,3) * v(i,1,3,e) &
                        + Dz(k,4) * v(i,1,4,e)

             wwt(i,1,k) = Dz(k,1) * w(i,1,1,e) &
                        + Dz(k,2) * w(i,1,2,e) &
                        + Dz(k,3) * w(i,1,3,e) &
                        + Dz(k,4) * w(i,1,4,e)
          end do
       end do

       do i = 1, lx*lx*lx

          u1 = wur(i,1,1) * drdx(i,1,1,e) + wus(i,1,1) * dsdx(i,1,1,e) &
                                          + wut(i,1,1) * dtdx(i,1,1,e)
          u2 = wur(i,1,1) * drdy(i,1,1,e) + wus(i,1,1) * dsdy(i,1,1,e) &
                                          + wut(i,1,1) * dtdy(i,1,1,e)
          u3 = wur(i,1,1) * drdz(i,1,1,e) + wus(i,1,1) * dsdz(i,1,1,e) &
                                          + wut(i,1,1) * dtdz(i,1,1,e)

          v1 = wvr(i,1,1) * drdx(i,1,1,e) + wvs(i,1,1) * dsdx(i,1,1,e) &
                                          + wvt(i,1,1) * dtdx(i,1,1,e)
          v2 = wvr(i,1,1) * drdy(i,1,1,e) + wvs(i,1,1) * dsdy(i,1,1,e) &
                                          + wvt(i,1,1) * dtdy(i,1,1,e)
          v3 = wvr(i,1,1) * drdz(i,1,1,e) + wvs(i,1,1) * dsdz(i,1,1,e) &
                                          + wvt(i,1,1) * dtdz(i,1,1,e)

          w1 = wwr(i,1,1) * drdx(i,1,1,e) + wws(i,1,1) * dsdx(i,1,1,e) &
                                          + wwt(i,1,1) * dtdx(i,1,1,e)
          w2 = wwr(i,1,1) * drdy(i,1,1,e) + wws(i,1,1) * dsdy(i,1,1,e) &
                                          + wwt(i,1,1) * dtdy(i,1,1,e)
          w3 = wwr(i,1,1) * drdz(i,1,1,e) + wws(i,1,1) * dsdz(i,1,1,e) &
                                          + wwt(i,1,1) * dtdz(i,1,1,e)

          dj  = h1(i,1,1,e) * weights3(i,1,1) * jacinv(i,1,1,e)
          s11 = dj*(u1 + u1)
          s12 = dj*(u2 + v1)
          s13 = dj*(u3 + w1)
          s21 = dj*(v1 + u2)
          s22 = dj*(v2 + v2)
          s23 = dj*(v3 + w2)
          s31 = dj*(w1 + u3)
          s32 = dj*(w2 + v3)
          s33 = dj*(w3 + w3)

          wur(i,1,1) = drdx(i,1,1,e)*s11 + drdy(i,1,1,e)*s12 + drdz(i,1,1,e)*s13
          wus(i,1,1) = dsdx(i,1,1,e)*s11 + dsdy(i,1,1,e)*s12 + dsdz(i,1,1,e)*s13
          wut(i,1,1) = dtdx(i,1,1,e)*s11 + dtdy(i,1,1,e)*s12 + dtdz(i,1,1,e)*s13

          wvr(i,1,1) = drdx(i,1,1,e)*s21 + drdy(i,1,1,e)*s22 + drdz(i,1,1,e)*s23
          wvs(i,1,1) = dsdx(i,1,1,e)*s21 + dsdy(i,1,1,e)*s22 + dsdz(i,1,1,e)*s23
          wvt(i,1,1) = dtdx(i,1,1,e)*s21 + dtdy(i,1,1,e)*s22 + dtdz(i,1,1,e)*s23

          wwr(i,1,1) = drdx(i,1,1,e)*s31 + drdy(i,1,1,e)*s32 + drdz(i,1,1,e)*s33
          wws(i,1,1) = dsdx(i,1,1,e)*s31 + dsdy(i,1,1,e)*s32 + dsdz(i,1,1,e)*s33
          wwt(i,1,1) = dtdx(i,1,1,e)*s31 + dtdy(i,1,1,e)*s32 + dtdz(i,1,1,e)*s33
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * wur(1,j,1) &
                         + Dxt(i,2) * wur(2,j,1) &
                         + Dxt(i,3) * wur(3,j,1) &
                         + Dxt(i,4) * wur(4,j,1)

             av(i,j,1,e) = Dxt(i,1) * wvr(1,j,1) &
                         + Dxt(i,2) * wvr(2,j,1) &
                         + Dxt(i,3) * wvr(3,j,1) &
                         + Dxt(i,4) * wvr(4,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wwr(1,j,1) &
                         + Dxt(i,2) * wwr(2,j,1) &
                         + Dxt(i,3) * wwr(3,j,1) &
                         + Dxt(i,4) * wwr(4,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * wus(i,1,k) &
                            + Dyt(j,2) * wus(i,2,k) &
                            + Dyt(j,3) * wus(i,3,k) &
                            + Dyt(j,4) * wus(i,4,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * wvs(i,1,k) &
                            + Dyt(j,2) * wvs(i,2,k) &
                            + Dyt(j,3) * wvs(i,3,k) &
                            + Dyt(j,4) * wvs(i,4,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * wws(i,1,k) &
                            + Dyt(j,2) * wws(i,2,k) &
                            + Dyt(j,3) * wws(i,3,k) &
                            + Dyt(j,4) * wws(i,4,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * wut(i,1,1) &
                         + Dzt(k,2) * wut(i,1,2) &
                         + Dzt(k,3) * wut(i,1,3) &
                         + Dzt(k,4) * wut(i,1,4)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * wvt(i,1,1) &
                         + Dzt(k,2) * wvt(i,1,2) &
                         + Dzt(k,3) * wvt(i,1,3) &
                         + Dzt(k,4) * wvt(i,1,4)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wwt(i,1,1) &
                         + Dzt(k,2) * wwt(i,1,2) &
                         + Dzt(k,3) * wwt(i,1,3) &
                         + Dzt(k,4) * wwt(i,1,4)
          end do
       end do

    end do

  end subroutine ax_helm_stress_lx4

  subroutine ax_helm_stress_lx3(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, &
       Dzt, h1, h2, drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
       jacinv, weights3, n)
    integer, parameter :: lx = 3
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: u(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: v(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdz(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: jacinv(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp), intent(in) :: weights3(lx, lx, lx)

    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    real(kind=rp) :: wvr(lx, lx, lx)
    real(kind=rp) :: wvs(lx, lx, lx)
    real(kind=rp) :: wvt(lx, lx, lx)
    real(kind=rp) :: wwr(lx, lx, lx)
    real(kind=rp) :: wws(lx, lx, lx)
    real(kind=rp) :: wwt(lx, lx, lx)

    integer :: e, i, j, k, l
    real(kind=rp) :: dj
    real(kind=rp) :: s11, s12, s13, s21, s22, s23, s31, s32, s33
    real(kind=rp) :: u1, u2, u3, v1, v2, v3, w1, w2, w3

    do e = 1, n

       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e) &
                        + Dx(i,3) * u(3,j,1,e)

             wvr(i,j,1) = Dx(i,1) * v(1,j,1,e) &
                        + Dx(i,2) * v(2,j,1,e) &
                        + Dx(i,3) * v(3,j,1,e)

             wwr(i,j,1) = Dx(i,1) * w(1,j,1,e) &
                        + Dx(i,2) * w(2,j,1,e) &
                        + Dx(i,3) * w(3,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e) &
                           + Dy(j,3) * u(i,3,k,e)

                wvs(i,j,k) = Dy(j,1) * v(i,1,k,e) &
                           + Dy(j,2) * v(i,2,k,e) &
                           + Dy(j,3) * v(i,3,k,e)

                wws(i,j,k) = Dy(j,1) * w(i,1,k,e) &
                           + Dy(j,2) * w(i,2,k,e) &
                           + Dy(j,3) * w(i,3,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e) &
                        + Dz(k,3) * u(i,1,3,e)

             wvt(i,1,k) = Dz(k,1) * v(i,1,1,e) &
                        + Dz(k,2) * v(i,1,2,e) &
                        + Dz(k,3) * v(i,1,3,e)

             wwt(i,1,k) = Dz(k,1) * w(i,1,1,e) &
                        + Dz(k,2) * w(i,1,2,e) &
                        + Dz(k,3) * w(i,1,3,e)
          end do
       end do

       do i = 1, lx*lx*lx

          u1 = wur(i,1,1) * drdx(i,1,1,e) + wus(i,1,1) * dsdx(i,1,1,e) &
                                          + wut(i,1,1) * dtdx(i,1,1,e)
          u2 = wur(i,1,1) * drdy(i,1,1,e) + wus(i,1,1) * dsdy(i,1,1,e) &
                                          + wut(i,1,1) * dtdy(i,1,1,e)
          u3 = wur(i,1,1) * drdz(i,1,1,e) + wus(i,1,1) * dsdz(i,1,1,e) &
                                          + wut(i,1,1) * dtdz(i,1,1,e)

          v1 = wvr(i,1,1) * drdx(i,1,1,e) + wvs(i,1,1) * dsdx(i,1,1,e) &
                                          + wvt(i,1,1) * dtdx(i,1,1,e)
          v2 = wvr(i,1,1) * drdy(i,1,1,e) + wvs(i,1,1) * dsdy(i,1,1,e) &
                                          + wvt(i,1,1) * dtdy(i,1,1,e)
          v3 = wvr(i,1,1) * drdz(i,1,1,e) + wvs(i,1,1) * dsdz(i,1,1,e) &
                                          + wvt(i,1,1) * dtdz(i,1,1,e)

          w1 = wwr(i,1,1) * drdx(i,1,1,e) + wws(i,1,1) * dsdx(i,1,1,e) &
                                          + wwt(i,1,1) * dtdx(i,1,1,e)
          w2 = wwr(i,1,1) * drdy(i,1,1,e) + wws(i,1,1) * dsdy(i,1,1,e) &
                                          + wwt(i,1,1) * dtdy(i,1,1,e)
          w3 = wwr(i,1,1) * drdz(i,1,1,e) + wws(i,1,1) * dsdz(i,1,1,e) &
                                          + wwt(i,1,1) * dtdz(i,1,1,e)

          dj  = h1(i,1,1,e) * weights3(i,1,1) * jacinv(i,1,1,e)
          s11 = dj*(u1 + u1)
          s12 = dj*(u2 + v1)
          s13 = dj*(u3 + w1)
          s21 = dj*(v1 + u2)
          s22 = dj*(v2 + v2)
          s23 = dj*(v3 + w2)
          s31 = dj*(w1 + u3)
          s32 = dj*(w2 + v3)
          s33 = dj*(w3 + w3)

          wur(i,1,1) = drdx(i,1,1,e)*s11 + drdy(i,1,1,e)*s12 + drdz(i,1,1,e)*s13
          wus(i,1,1) = dsdx(i,1,1,e)*s11 + dsdy(i,1,1,e)*s12 + dsdz(i,1,1,e)*s13
          wut(i,1,1) = dtdx(i,1,1,e)*s11 + dtdy(i,1,1,e)*s12 + dtdz(i,1,1,e)*s13

          wvr(i,1,1) = drdx(i,1,1,e)*s21 + drdy(i,1,1,e)*s22 + drdz(i,1,1,e)*s23
          wvs(i,1,1) = dsdx(i,1,1,e)*s21 + dsdy(i,1,1,e)*s22 + dsdz(i,1,1,e)*s23
          wvt(i,1,1) = dtdx(i,1,1,e)*s21 + dtdy(i,1,1,e)*s22 + dtdz(i,1,1,e)*s23

          wwr(i,1,1) = drdx(i,1,1,e)*s31 + drdy(i,1,1,e)*s32 + drdz(i,1,1,e)*s33
          wws(i,1,1) = dsdx(i,1,1,e)*s31 + dsdy(i,1,1,e)*s32 + dsdz(i,1,1,e)*s33
          wwt(i,1,1) = dtdx(i,1,1,e)*s31 + dtdy(i,1,1,e)*s32 + dtdz(i,1,1,e)*s33
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * wur(1,j,1) &
                         + Dxt(i,2) * wur(2,j,1) &
                         + Dxt(i,3) * wur(3,j,1)

             av(i,j,1,e) = Dxt(i,1) * wvr(1,j,1) &
                         + Dxt(i,2) * wvr(2,j,1) &
                         + Dxt(i,3) * wvr(3,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wwr(1,j,1) &
                         + Dxt(i,2) * wwr(2,j,1) &
                         + Dxt(i,3) * wwr(3,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * wus(i,1,k) &
                            + Dyt(j,2) * wus(i,2,k) &
                            + Dyt(j,3) * wus(i,3,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * wvs(i,1,k) &
                            + Dyt(j,2) * wvs(i,2,k) &
                            + Dyt(j,3) * wvs(i,3,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * wws(i,1,k) &
                            + Dyt(j,2) * wws(i,2,k) &
                            + Dyt(j,3) * wws(i,3,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * wut(i,1,1) &
                         + Dzt(k,2) * wut(i,1,2) &
                         + Dzt(k,3) * wut(i,1,3)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * wvt(i,1,1) &
                         + Dzt(k,2) * wvt(i,1,2) &
                         + Dzt(k,3) * wvt(i,1,3)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wwt(i,1,1) &
                         + Dzt(k,2) * wwt(i,1,2) &
                         + Dzt(k,3) * wwt(i,1,3)
          end do
       end do

    end do

  end subroutine ax_helm_stress_lx3

  subroutine ax_helm_stress_lx2(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, &
       Dzt, h1, h2, drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
       jacinv, weights3, n)
    integer, parameter :: lx = 2
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: u(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: v(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdz(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: jacinv(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp), intent(in) :: weights3(lx, lx, lx)

    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    real(kind=rp) :: wvr(lx, lx, lx)
    real(kind=rp) :: wvs(lx, lx, lx)
    real(kind=rp) :: wvt(lx, lx, lx)
    real(kind=rp) :: wwr(lx, lx, lx)
    real(kind=rp) :: wws(lx, lx, lx)
    real(kind=rp) :: wwt(lx, lx, lx)

    integer :: e, i, j, k
    real(kind=rp) :: dj
    real(kind=rp) :: s11, s12, s13, s21, s22, s23, s31, s32, s33
    real(kind=rp) :: u1, u2, u3, v1, v2, v3, w1, w2, w3

    do e = 1, n

       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = Dx(i,1) * u(1,j,1,e) &
                        + Dx(i,2) * u(2,j,1,e)

             wvr(i,j,1) = Dx(i,1) * v(1,j,1,e) &
                        + Dx(i,2) * v(2,j,1,e)

             wwr(i,j,1) = Dx(i,1) * w(1,j,1,e) &
                        + Dx(i,2) * w(2,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = Dy(j,1) * u(i,1,k,e) &
                           + Dy(j,2) * u(i,2,k,e)

                wvs(i,j,k) = Dy(j,1) * v(i,1,k,e) &
                           + Dy(j,2) * v(i,2,k,e)

                wws(i,j,k) = Dy(j,1) * w(i,1,k,e) &
                           + Dy(j,2) * w(i,2,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = Dz(k,1) * u(i,1,1,e) &
                        + Dz(k,2) * u(i,1,2,e)

             wvt(i,1,k) = Dz(k,1) * v(i,1,1,e) &
                        + Dz(k,2) * v(i,1,2,e)

             wwt(i,1,k) = Dz(k,1) * w(i,1,1,e) &
                        + Dz(k,2) * w(i,1,2,e)
          end do
       end do

       do i = 1, lx*lx*lx

          u1 = wur(i,1,1) * drdx(i,1,1,e) + wus(i,1,1) * dsdx(i,1,1,e) &
                                          + wut(i,1,1) * dtdx(i,1,1,e)
          u2 = wur(i,1,1) * drdy(i,1,1,e) + wus(i,1,1) * dsdy(i,1,1,e) &
                                          + wut(i,1,1) * dtdy(i,1,1,e)
          u3 = wur(i,1,1) * drdz(i,1,1,e) + wus(i,1,1) * dsdz(i,1,1,e) &
                                          + wut(i,1,1) * dtdz(i,1,1,e)

          v1 = wvr(i,1,1) * drdx(i,1,1,e) + wvs(i,1,1) * dsdx(i,1,1,e) &
                                          + wvt(i,1,1) * dtdx(i,1,1,e)
          v2 = wvr(i,1,1) * drdy(i,1,1,e) + wvs(i,1,1) * dsdy(i,1,1,e) &
                                          + wvt(i,1,1) * dtdy(i,1,1,e)
          v3 = wvr(i,1,1) * drdz(i,1,1,e) + wvs(i,1,1) * dsdz(i,1,1,e) &
                                          + wvt(i,1,1) * dtdz(i,1,1,e)

          w1 = wwr(i,1,1) * drdx(i,1,1,e) + wws(i,1,1) * dsdx(i,1,1,e) &
                                          + wwt(i,1,1) * dtdx(i,1,1,e)
          w2 = wwr(i,1,1) * drdy(i,1,1,e) + wws(i,1,1) * dsdy(i,1,1,e) &
                                          + wwt(i,1,1) * dtdy(i,1,1,e)
          w3 = wwr(i,1,1) * drdz(i,1,1,e) + wws(i,1,1) * dsdz(i,1,1,e) &
                                          + wwt(i,1,1) * dtdz(i,1,1,e)

          dj  = h1(i,1,1,e) * weights3(i,1,1) * jacinv(i,1,1,e)
          s11 = dj*(u1 + u1)
          s12 = dj*(u2 + v1)
          s13 = dj*(u3 + w1)
          s21 = dj*(v1 + u2)
          s22 = dj*(v2 + v2)
          s23 = dj*(v3 + w2)
          s31 = dj*(w1 + u3)
          s32 = dj*(w2 + v3)
          s33 = dj*(w3 + w3)

          wur(i,1,1) = drdx(i,1,1,e)*s11 + drdy(i,1,1,e)*s12 + drdz(i,1,1,e)*s13
          wus(i,1,1) = dsdx(i,1,1,e)*s11 + dsdy(i,1,1,e)*s12 + dsdz(i,1,1,e)*s13
          wut(i,1,1) = dtdx(i,1,1,e)*s11 + dtdy(i,1,1,e)*s12 + dtdz(i,1,1,e)*s13

          wvr(i,1,1) = drdx(i,1,1,e)*s21 + drdy(i,1,1,e)*s22 + drdz(i,1,1,e)*s23
          wvs(i,1,1) = dsdx(i,1,1,e)*s21 + dsdy(i,1,1,e)*s22 + dsdz(i,1,1,e)*s23
          wvt(i,1,1) = dtdx(i,1,1,e)*s21 + dtdy(i,1,1,e)*s22 + dtdz(i,1,1,e)*s23

          wwr(i,1,1) = drdx(i,1,1,e)*s31 + drdy(i,1,1,e)*s32 + drdz(i,1,1,e)*s33
          wws(i,1,1) = dsdx(i,1,1,e)*s31 + dsdy(i,1,1,e)*s32 + dsdz(i,1,1,e)*s33
          wwt(i,1,1) = dtdx(i,1,1,e)*s31 + dtdy(i,1,1,e)*s32 + dtdz(i,1,1,e)*s33
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * wur(1,j,1) &
                         + Dxt(i,2) * wur(2,j,1)

             av(i,j,1,e) = Dxt(i,1) * wvr(1,j,1) &
                         + Dxt(i,2) * wvr(2,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wwr(1,j,1) &
                         + Dxt(i,2) * wwr(2,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * wus(i,1,k) &
                            + Dyt(j,2) * wus(i,2,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * wvs(i,1,k) &
                            + Dyt(j,2) * wvs(i,2,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * wws(i,1,k) &
                            + Dyt(j,2) * wws(i,2,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * wut(i,1,1) &
                         + Dzt(k,2) * wut(i,1,2)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * wvt(i,1,1) &
                         + Dzt(k,2) * wvt(i,1,2)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wwt(i,1,1) &
                         + Dzt(k,2) * wwt(i,1,2)
          end do
       end do

    end do

  end subroutine ax_helm_stress_lx2

end module
