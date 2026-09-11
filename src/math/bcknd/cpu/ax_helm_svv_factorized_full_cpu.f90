! Copyright (c) 2025-2026, The Neko Authors
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
!> CPU implementation of the factorized full-stress SVV Helmholtz operator.
module ax_helm_svv_factorized_full_cpu
  use ax_helm_svv_full, only : ax_helm_svv_full_t
  use num_types, only : rp
  use coefs, only : coef_t
  use space, only : space_t
  use mesh, only : mesh_t
  use math, only : addcol4
  implicit none
  private

  !> CPU matrix-vector product for a Helmholtz problem.
  type, public, extends(ax_helm_svv_full_t) :: ax_helm_svv_factorized_full_cpu_t
   contains
     !> Compute the product.
     procedure, pass(this) :: compute_vector => &
          ax_helm_svv_factorized_full_compute_vector
  end type ax_helm_svv_factorized_full_cpu_t

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
  subroutine ax_helm_svv_factorized_full_compute_vector(this, au, av, aw, &
       u, v, w, coef, msh, Xh)
    class(ax_helm_svv_factorized_full_cpu_t), intent(in) :: this
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: v(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: au(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: av(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: aw(Xh%lx, Xh%ly, Xh%lz, msh%nelv)

    call ax_helm_svv_factorized_full_lx(au, av, aw, u, v, w, &
            Xh%dx, Xh%dy, Xh%dz, Xh%dxt, Xh%dyt, Xh%dzt, &
            coef%h1, coef%drdx, coef%drdy, coef%drdz, coef%dsdx, coef%dsdy, &
            coef%dsdz, coef%dtdx, coef%dtdy, coef%dtdz, &
            coef%jacinv, Xh%w3, this%svv%h1, &
            this%svv%Br, this%svv%Bs, this%svv%Bt, msh%nelv, Xh%lx)

    if (coef%ifh2) then
       call addcol4 (au, coef%h2, coef%B, u, coef%dof%size())
       call addcol4 (av, coef%h2, coef%B, v, coef%dof%size())
       call addcol4 (aw, coef%h2, coef%B, w, coef%dof%size())
    end if


  end subroutine ax_helm_svv_factorized_full_compute_vector

  !> Generic CPU kernel for the factorized full-stress SVV product.
  !! @param au Result for the first vector component.
  !! @param av Result for the second vector component.
  !! @param aw Result for the third vector component.
  !! @param u First input vector component.
  !! @param v Second input vector component.
  !! @param w Third input vector component.
  !! @param Dx Derivative in the first reference direction.
  !! @param Dy Derivative in the second reference direction.
  !! @param Dz Derivative in the third reference direction.
  !! @param Dxt Transpose derivative in the first reference direction.
  !! @param Dyt Transpose derivative in the second reference direction.
  !! @param Dzt Transpose derivative in the third reference direction.
  !! @param h1 Ordinary diffusion coefficient.
  !! @param drdx Geometric derivative factor.
  !! @param drdy Geometric derivative factor.
  !! @param drdz Geometric derivative factor.
  !! @param dsdx Geometric derivative factor.
  !! @param dsdy Geometric derivative factor.
  !! @param dsdz Geometric derivative factor.
  !! @param dtdx Geometric derivative factor.
  !! @param dtdy Geometric derivative factor.
  !! @param dtdz Geometric derivative factor.
  !! @param jacinv Inverse Jacobian determinant.
  !! @param weights3 Tensor-product quadrature weights.
  !! @param svv_h1 Density-weighted SVV viscosity.
  !! @param Br Complementary derivative in the first reference direction.
  !! @param Bs Complementary derivative in the second reference direction.
  !! @param Bt Complementary derivative in the third reference direction.
  !! @param n Number of elements.
  !! @param lx Number of points in one element direction.
  subroutine ax_helm_svv_factorized_full_lx(au, av, aw, u, v, w, &
       Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
       jacinv, weights3, svv_h1, Br, Bs, Bt, n, lx)
    integer, intent(in) :: n, lx
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: drdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dsdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdx(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdy(lx, lx, lx, n)
    real(kind=rp), intent(in) :: dtdz(lx, lx, lx, n)
    real(kind=rp), intent(in) :: jacinv(lx, lx, lx, n)
    real(kind=rp), intent(in) :: weights3(lx, lx, lx)
    real(kind=rp), intent(in) :: Dx(lx,lx)
    real(kind=rp), intent(in) :: Dy(lx,lx)
    real(kind=rp), intent(in) :: Dz(lx,lx)
    real(kind=rp), intent(in) :: Dxt(lx,lx)
    real(kind=rp), intent(in) :: Dyt(lx,lx)
    real(kind=rp), intent(in) :: Dzt(lx,lx)
    real(kind=rp), intent(in) :: svv_h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Br(lx, lx), Bs(lx, lx), Bt(lx, lx)

    real(kind=rp) :: s11_h(lx, lx, lx), s22_h(lx, lx, lx)
    real(kind=rp) :: s33_h(lx, lx, lx), s12_h(lx, lx, lx)
    real(kind=rp) :: s13_h(lx, lx, lx), s23_h(lx, lx, lx)
    
    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    real(kind=rp) :: wvr(lx, lx, lx)
    real(kind=rp) :: wvs(lx, lx, lx)
    real(kind=rp) :: wvt(lx, lx, lx)
    real(kind=rp) :: wwr(lx, lx, lx)
    real(kind=rp) :: wws(lx, lx, lx)
    real(kind=rp) :: wwt(lx, lx, lx)
    real(kind=rp) :: ur_svv(lx, lx, lx)
    real(kind=rp) :: us_svv(lx, lx, lx)
    real(kind=rp) :: ut_svv(lx, lx, lx)
    real(kind=rp) :: vr_svv(lx, lx, lx)
    real(kind=rp) :: vs_svv(lx, lx, lx)
    real(kind=rp) :: vt_svv(lx, lx, lx)
    real(kind=rp) :: wr_svv(lx, lx, lx)
    real(kind=rp) :: ws_svv(lx, lx, lx)
    real(kind=rp) :: wt_svv(lx, lx, lx)

    integer :: e, i, j, k, l

    real(kind=rp) :: t1, t2, t3, t1_svv, t2_svv, t3_svv
    real(kind=rp) :: s11(lx, lx, lx)
    real(kind=rp) :: s22(lx, lx, lx)
    real(kind=rp) :: s33(lx, lx, lx)
    real(kind=rp) :: s12(lx, lx, lx)
    real(kind=rp) :: s13(lx, lx, lx)
    real(kind=rp) :: s23(lx, lx, lx)
    real(kind=rp) :: u1, u2, u3, v1, v2, v3, w1, w2, w3

    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             t1 = 0.0_rp
             t2 = 0.0_rp
             t3 = 0.0_rp
             t1_svv = 0.0_rp
             t2_svv = 0.0_rp
             t3_svv = 0.0_rp
             do k = 1, lx
                t1 = t1 + Dx(i,k) * u(k,j,1,e)
                t2 = t2 + Dx(i,k) * v(k,j,1,e)
                t3 = t3 + Dx(i,k) * w(k,j,1,e)
                t1_svv = t1_svv + Br(i,k) * u(k,j,1,e)
                t2_svv = t2_svv + Br(i,k) * v(k,j,1,e)
                t3_svv = t3_svv + Br(i,k) * w(k,j,1,e)
             end do
             wur(i,j,1) = t1
             wvr(i,j,1) = t2
             wwr(i,j,1) = t3
             ur_svv(i,j,1) = t1_svv
             vr_svv(i,j,1) = t2_svv
             wr_svv(i,j,1) = t3_svv
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                t1 = 0.0_rp
                t2 = 0.0_rp
                t3 = 0.0_rp
                t1_svv = 0.0_rp
                t2_svv = 0.0_rp
                t3_svv = 0.0_rp
                do l = 1, lx
                   t1 = t1 + Dy(j,l) * u(i,l,k,e)
                   t2 = t2 + Dy(j,l) * v(i,l,k,e)
                   t3 = t3 + Dy(j,l) * w(i,l,k,e)
                   t1_svv = t1_svv + Bs(j,l) * u(i,l,k,e)
                   t2_svv = t2_svv + Bs(j,l) * v(i,l,k,e)
                   t3_svv = t3_svv + Bs(j,l) * w(i,l,k,e)
                end do
                wus(i,j,k) = t1
                wvs(i,j,k) = t2
                wws(i,j,k) = t3
                us_svv(i,j,k) = t1_svv
                vs_svv(i,j,k) = t2_svv
                ws_svv(i,j,k) = t3_svv
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             t1 = 0.0_rp
             t2 = 0.0_rp
             t3 = 0.0_rp
             t1_svv = 0.0_rp
             t2_svv = 0.0_rp
             t3_svv = 0.0_rp
             do l = 1, lx
                t1 = t1 + Dz(k,l) * u(i,1,l,e)
                t2 = t2 + Dz(k,l) * v(i,1,l,e)
                t3 = t3 + Dz(k,l) * w(i,1,l,e)
                t1_svv = t1_svv + Bt(k,l) * u(i,1,l,e)
                t2_svv = t2_svv + Bt(k,l) * v(i,1,l,e)
                t3_svv = t3_svv + Bt(k,l) * w(i,1,l,e)
             end do
             wut(i,1,k) = t1
             wvt(i,1,k) = t2
             wwt(i,1,k) = t3
             ut_svv(i,1,k) = t1_svv
             vt_svv(i,1,k) = t2_svv
             wt_svv(i,1,k) = t3_svv
          end do
       end do

       do i = 1, lx*lx*lx
          u1 = (drdx(i,1,1,e) * wur(i,1,1) &
              + dsdx(i,1,1,e) * wus(i,1,1) &
              + dtdx(i,1,1,e) * wut(i,1,1)) * jacinv(i,1,1,e)
          u2 = (drdy(i,1,1,e) * wur(i,1,1) &
              + dsdy(i,1,1,e) * wus(i,1,1) &
              + dtdy(i,1,1,e) * wut(i,1,1)) * jacinv(i,1,1,e)
          u3 = (drdz(i,1,1,e) * wur(i,1,1) &
              + dsdz(i,1,1,e) * wus(i,1,1) &
              + dtdz(i,1,1,e) * wut(i,1,1)) * jacinv(i,1,1,e)
          v1 = (drdx(i,1,1,e) * wvr(i,1,1) &
              + dsdx(i,1,1,e) * wvs(i,1,1) &
              + dtdx(i,1,1,e) * wvt(i,1,1)) * jacinv(i,1,1,e)
          v2 = (drdy(i,1,1,e) * wvr(i,1,1) &
              + dsdy(i,1,1,e) * wvs(i,1,1) &
              + dtdy(i,1,1,e) * wvt(i,1,1)) * jacinv(i,1,1,e)
          v3 = (drdz(i,1,1,e) * wvr(i,1,1) &
              + dsdz(i,1,1,e) * wvs(i,1,1) &
              + dtdz(i,1,1,e) * wvt(i,1,1)) * jacinv(i,1,1,e)
          w1 = (drdx(i,1,1,e) * wwr(i,1,1) &
              + dsdx(i,1,1,e) * wws(i,1,1) &
              + dtdx(i,1,1,e) * wwt(i,1,1)) * jacinv(i,1,1,e)
          w2 = (drdy(i,1,1,e) * wwr(i,1,1) &
              + dsdy(i,1,1,e) * wws(i,1,1) &
              + dtdy(i,1,1,e) * wwt(i,1,1)) * jacinv(i,1,1,e)
          w3 = (drdz(i,1,1,e) * wwr(i,1,1) &
              + dsdz(i,1,1,e) * wws(i,1,1) &
              + dtdz(i,1,1,e) * wwt(i,1,1)) * jacinv(i,1,1,e)
          s11(i,1,1) = u1 + u1
          s22(i,1,1) = v2 + v2
          s33(i,1,1) = w3 + w3
          s12(i,1,1) = u2 + v1
          s13(i,1,1) = u3 + w1
          s23(i,1,1) = v3 + w2
       end do

       do i = 1, lx*lx*lx
          ! Standard unfiltered stress and reference-space flux.
          s11_h(i,1,1) = h1(i,1,1,e) * s11(i,1,1) * weights3(i,1,1)
          s22_h(i,1,1) = h1(i,1,1,e) * s22(i,1,1) * weights3(i,1,1)
          s33_h(i,1,1) = h1(i,1,1,e) * s33(i,1,1) * weights3(i,1,1)
          s12_h(i,1,1) = h1(i,1,1,e) * s12(i,1,1) * weights3(i,1,1)
          s13_h(i,1,1) = h1(i,1,1,e) * s13(i,1,1) * weights3(i,1,1)
          s23_h(i,1,1) = h1(i,1,1,e) * s23(i,1,1) * weights3(i,1,1)
          wur(i,1,1) = drdx(i,1,1,e) * s11_h(i,1,1) &
                     + drdy(i,1,1,e) * s12_h(i,1,1) &
                     + drdz(i,1,1,e) * s13_h(i,1,1)
          wus(i,1,1) = dsdx(i,1,1,e) * s11_h(i,1,1) &
                     + dsdy(i,1,1,e) * s12_h(i,1,1) &
                     + dsdz(i,1,1,e) * s13_h(i,1,1)
          wut(i,1,1) = dtdx(i,1,1,e) * s11_h(i,1,1) &
                     + dtdy(i,1,1,e) * s12_h(i,1,1) &
                     + dtdz(i,1,1,e) * s13_h(i,1,1)
          wvr(i,1,1) = drdx(i,1,1,e) * s12_h(i,1,1) &
                     + drdy(i,1,1,e) * s22_h(i,1,1) &
                     + drdz(i,1,1,e) * s23_h(i,1,1)
          wvs(i,1,1) = dsdx(i,1,1,e) * s12_h(i,1,1) &
                     + dsdy(i,1,1,e) * s22_h(i,1,1) &
                     + dsdz(i,1,1,e) * s23_h(i,1,1)
          wvt(i,1,1) = dtdx(i,1,1,e) * s12_h(i,1,1) &
                     + dtdy(i,1,1,e) * s22_h(i,1,1) &
                     + dtdz(i,1,1,e) * s23_h(i,1,1)
          wwr(i,1,1) = drdx(i,1,1,e) * s13_h(i,1,1) &
                     + drdy(i,1,1,e) * s23_h(i,1,1) &
                     + drdz(i,1,1,e) * s33_h(i,1,1)
          wws(i,1,1) = dsdx(i,1,1,e) * s13_h(i,1,1) &
                     + dsdy(i,1,1,e) * s23_h(i,1,1) &
                     + dsdz(i,1,1,e) * s33_h(i,1,1)
          wwt(i,1,1) = dtdx(i,1,1,e) * s13_h(i,1,1) &
                     + dtdy(i,1,1,e) * s23_h(i,1,1) &
                     + dtdz(i,1,1,e) * s33_h(i,1,1)

          ! Map the filtered reference derivatives to physical space and
          ! form the six factorized SVV stress components.
          u1 = (drdx(i,1,1,e) * ur_svv(i,1,1) &
              + dsdx(i,1,1,e) * us_svv(i,1,1) &
              + dtdx(i,1,1,e) * ut_svv(i,1,1)) * jacinv(i,1,1,e)
          u2 = (drdy(i,1,1,e) * ur_svv(i,1,1) &
              + dsdy(i,1,1,e) * us_svv(i,1,1) &
              + dtdy(i,1,1,e) * ut_svv(i,1,1)) * jacinv(i,1,1,e)
          u3 = (drdz(i,1,1,e) * ur_svv(i,1,1) &
              + dsdz(i,1,1,e) * us_svv(i,1,1) &
              + dtdz(i,1,1,e) * ut_svv(i,1,1)) * jacinv(i,1,1,e)
          v1 = (drdx(i,1,1,e) * vr_svv(i,1,1) &
              + dsdx(i,1,1,e) * vs_svv(i,1,1) &
              + dtdx(i,1,1,e) * vt_svv(i,1,1)) * jacinv(i,1,1,e)
          v2 = (drdy(i,1,1,e) * vr_svv(i,1,1) &
              + dsdy(i,1,1,e) * vs_svv(i,1,1) &
              + dtdy(i,1,1,e) * vt_svv(i,1,1)) * jacinv(i,1,1,e)
          v3 = (drdz(i,1,1,e) * vr_svv(i,1,1) &
              + dsdz(i,1,1,e) * vs_svv(i,1,1) &
              + dtdz(i,1,1,e) * vt_svv(i,1,1)) * jacinv(i,1,1,e)
          w1 = (drdx(i,1,1,e) * wr_svv(i,1,1) &
              + dsdx(i,1,1,e) * ws_svv(i,1,1) &
              + dtdx(i,1,1,e) * wt_svv(i,1,1)) * jacinv(i,1,1,e)
          w2 = (drdy(i,1,1,e) * wr_svv(i,1,1) &
              + dsdy(i,1,1,e) * ws_svv(i,1,1) &
              + dtdy(i,1,1,e) * wt_svv(i,1,1)) * jacinv(i,1,1,e)
          w3 = (drdz(i,1,1,e) * wr_svv(i,1,1) &
              + dsdz(i,1,1,e) * ws_svv(i,1,1) &
              + dtdz(i,1,1,e) * wt_svv(i,1,1)) * jacinv(i,1,1,e)

          s11_h(i,1,1) = svv_h1(i,1,1,e) * (u1 + u1) * weights3(i,1,1)
          s22_h(i,1,1) = svv_h1(i,1,1,e) * (v2 + v2) * weights3(i,1,1)
          s33_h(i,1,1) = svv_h1(i,1,1,e) * (w3 + w3) * weights3(i,1,1)
          s12_h(i,1,1) = svv_h1(i,1,1,e) * (u2 + v1) * weights3(i,1,1)
          s13_h(i,1,1) = svv_h1(i,1,1,e) * (u3 + w1) * weights3(i,1,1)
          s23_h(i,1,1) = svv_h1(i,1,1,e) * (v3 + w2) * weights3(i,1,1)

          ur_svv(i,1,1) = drdx(i,1,1,e) * s11_h(i,1,1) &
                        + drdy(i,1,1,e) * s12_h(i,1,1) &
                        + drdz(i,1,1,e) * s13_h(i,1,1)
          us_svv(i,1,1) = dsdx(i,1,1,e) * s11_h(i,1,1) &
                        + dsdy(i,1,1,e) * s12_h(i,1,1) &
                        + dsdz(i,1,1,e) * s13_h(i,1,1)
          ut_svv(i,1,1) = dtdx(i,1,1,e) * s11_h(i,1,1) &
                        + dtdy(i,1,1,e) * s12_h(i,1,1) &
                        + dtdz(i,1,1,e) * s13_h(i,1,1)
          vr_svv(i,1,1) = drdx(i,1,1,e) * s12_h(i,1,1) &
                        + drdy(i,1,1,e) * s22_h(i,1,1) &
                        + drdz(i,1,1,e) * s23_h(i,1,1)
          vs_svv(i,1,1) = dsdx(i,1,1,e) * s12_h(i,1,1) &
                        + dsdy(i,1,1,e) * s22_h(i,1,1) &
                        + dsdz(i,1,1,e) * s23_h(i,1,1)
          vt_svv(i,1,1) = dtdx(i,1,1,e) * s12_h(i,1,1) &
                        + dtdy(i,1,1,e) * s22_h(i,1,1) &
                        + dtdz(i,1,1,e) * s23_h(i,1,1)
          wr_svv(i,1,1) = drdx(i,1,1,e) * s13_h(i,1,1) &
                        + drdy(i,1,1,e) * s23_h(i,1,1) &
                        + drdz(i,1,1,e) * s33_h(i,1,1)
          ws_svv(i,1,1) = dsdx(i,1,1,e) * s13_h(i,1,1) &
                        + dsdy(i,1,1,e) * s23_h(i,1,1) &
                        + dsdz(i,1,1,e) * s33_h(i,1,1)
          wt_svv(i,1,1) = dtdx(i,1,1,e) * s13_h(i,1,1) &
                        + dtdy(i,1,1,e) * s23_h(i,1,1) &
                        + dtdz(i,1,1,e) * s33_h(i,1,1)
       end do

       do j = 1, lx*lx
          do i = 1, lx
             t1 = 0.0_rp
             t2 = 0.0_rp
             t3 = 0.0_rp
             do k = 1, lx
                t1 = t1 + Dxt(i,k) * wur(k,j,1) &
                     + Br(k,i) * ur_svv(k,j,1)
                t2 = t2 + Dxt(i,k) * wvr(k,j,1) &
                     + Br(k,i) * vr_svv(k,j,1)
                t3 = t3 + Dxt(i,k) * wwr(k,j,1) &
                     + Br(k,i) * wr_svv(k,j,1)
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
                   t1 = t1 + Dyt(j,l) * wus(i,l,k) &
                        + Bs(l,j) * us_svv(i,l,k)
                   t2 = t2 + Dyt(j,l) * wvs(i,l,k) &
                        + Bs(l,j) * vs_svv(i,l,k)
                   t3 = t3 + Dyt(j,l) * wws(i,l,k) &
                        + Bs(l,j) * ws_svv(i,l,k)
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
                t1 = t1 + Dzt(k,l) * wut(i,1,l) &
                     + Bt(l,k) * ut_svv(i,1,l)
                t2 = t2 + Dzt(k,l) * wvt(i,1,l) &
                     + Bt(l,k) * vt_svv(i,1,l)
                t3 = t3 + Dzt(k,l) * wwt(i,1,l) &
                     + Bt(l,k) * wt_svv(i,1,l)
             end do
             au(i,1,k,e) = au(i,1,k,e) + t1
             av(i,1,k,e) = av(i,1,k,e) + t2
             aw(i,1,k,e) = aw(i,1,k,e) + t3
          end do
       end do

    end do
  end subroutine ax_helm_svv_factorized_full_lx

end module ax_helm_svv_factorized_full_cpu
