! Copyright (c) 2026, The Neko Authors
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
!     copyright notice, this list of conditions and the following disclaimer
!     in the documentation and/or other materials provided with the
!     distribution.
!
!   * Neither the name of the authors nor the names of its contributors may
!     be used to endorse or promote products derived from this software
!     without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> CPU implementation of the Kirby-Sherwin SVV Helmholtz operator.
module ax_helm_svv_KS_cpu
  use ax_helm_svv_KS, only : ax_helm_svv_KS_t
  use num_types, only : rp
  use coefs, only : coef_t
  use space, only : space_t
  use mesh, only : mesh_t
  use math, only : addcol4
  use tensor, only : tnsr3d_el
  implicit none
  private

  !> CPU Kirby-Sherwin SVV Helmholtz operator.
  type, public, extends(ax_helm_svv_KS_t) :: ax_helm_svv_KS_cpu_t
   contains
     procedure, pass(this) :: compute => ax_helm_svv_KS_compute
  end type ax_helm_svv_KS_cpu_t

contains

  !> Compute the Kirby-Sherwin SVV Helmholtz product.
  !! @param this CPU SVV operator.
  !! @param w Result.
  !! @param u Input field.
  !! @param coef SEM coefficients.
  !! @param msh Mesh.
  !! @param Xh Function space.
  subroutine ax_helm_svv_KS_compute(this, w, u, coef, msh, Xh)
    class(ax_helm_svv_KS_cpu_t), intent(in) :: this
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)

    call ax_helm_svv_KS_lx(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, Xh%dyt, &
         Xh%dzt, coef%h1, coef%drdx, coef%drdy, coef%drdz, coef%dsdx, &
         coef%dsdy, coef%dsdz, coef%dtdx, coef%dtdy, coef%dtdz, &
         coef%jacinv, Xh%w3, this%svv%h1, this%svv%filter%fh, &
         this%svv%filter%fht, this%svv%direction, this%svv%ident, &
         msh%nelv, Xh%lx)

    if (coef%ifh2) then
       call addcol4(w, coef%h2, coef%B, u, coef%dof%size())
    end if
  end subroutine ax_helm_svv_KS_compute

  !> Generic CPU kernel for the Kirby-Sherwin SVV Helmholtz product.
  !! @param w Result.
  !! @param u Input field.
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
  !! @param svv_F One-dimensional low-pass filter.
  !! @param svv_Ft Transpose one-dimensional low-pass filter.
  !! @param svv_direction Filtered reference directions.
  !! @param ident Identity matrix.
  !! @param n Number of elements.
  !! @param lx Number of points in one element direction.
  subroutine ax_helm_svv_KS_lx(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
       jacinv, weights3, svv_h1, svv_F, svv_Ft, svv_direction, ident, n, lx)
    integer, intent(in) :: n, lx
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
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
    real(kind=rp), intent(in) :: Dx(lx,lx), Dy(lx,lx), Dz(lx,lx)
    real(kind=rp), intent(in) :: Dxt(lx,lx), Dyt(lx,lx), Dzt(lx,lx)
    real(kind=rp), intent(in) :: svv_h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: svv_F(lx, lx), svv_Ft(lx, lx)
    real(kind=rp), intent(in) :: ident(lx, lx)
    character(len=*), intent(in) :: svv_direction
    real(kind=rp) :: u1(lx, lx, lx), u2(lx, lx, lx), u3(lx, lx, lx)
    real(kind=rp) :: u1_svv(lx, lx, lx), u2_svv(lx, lx, lx)
    real(kind=rp) :: u3_svv(lx, lx, lx)
    real(kind=rp) :: wur(lx, lx, lx), wus(lx, lx, lx), wut(lx, lx, lx)
    real(kind=rp) :: filter_r(lx, lx), filter_s(lx, lx), filter_t(lx, lx)
    real(kind=rp) :: ur_h, us_h, ut_h, tmp
    integer :: e, i, j, k, l

    if (index(svv_direction, "r") > 0) then
       filter_r = svv_F
    else
       filter_r = ident
    end if
    if (index(svv_direction, "s") > 0) then
       filter_s = svv_Ft
    else
       filter_s = ident
    end if
    if (index(svv_direction, "t") > 0) then
       filter_t = svv_Ft
    else
       filter_t = ident
    end if

    !$omp parallel do private(e, i, j, k, l, tmp, ur_h, us_h, ut_h) &
    !$omp& private(u1, u2, u3, u1_svv, u2_svv, u3_svv, wur, wus, wut)
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + Dx(i,k) * u(k,j,1,e)
             end do
             wur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + Dy(j,l) * u(i,l,k,e)
                end do
                wus(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx * lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + Dz(k,l) * u(i,1,l,e)
             end do
             wut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          u1(i,1,1) = (drdx(i,1,1,e) * wur(i,1,1) + &
               dsdx(i,1,1,e) * wus(i,1,1) + &
               dtdx(i,1,1,e) * wut(i,1,1)) * jacinv(i,1,1,e)
          u2(i,1,1) = (drdy(i,1,1,e) * wur(i,1,1) + &
               dsdy(i,1,1,e) * wus(i,1,1) + &
               dtdy(i,1,1,e) * wut(i,1,1)) * jacinv(i,1,1,e)
          u3(i,1,1) = (drdz(i,1,1,e) * wur(i,1,1) + &
               dsdz(i,1,1,e) * wus(i,1,1) + &
               dtdz(i,1,1,e) * wut(i,1,1)) * jacinv(i,1,1,e)
       end do

       call tnsr3d_el(u1_svv, lx, u1, lx, filter_r, filter_s, filter_t)
       call tnsr3d_el(u2_svv, lx, u2, lx, filter_r, filter_s, filter_t)
       call tnsr3d_el(u3_svv, lx, u3, lx, filter_r, filter_s, filter_t)

       do i = 1, lx * lx * lx
          u1_svv(i,1,1) = u1(i,1,1) - u1_svv(i,1,1)
          u2_svv(i,1,1) = u2(i,1,1) - u2_svv(i,1,1)
          u3_svv(i,1,1) = u3(i,1,1) - u3_svv(i,1,1)

          ur_h = (svv_h1(i,1,1,e) * u1_svv(i,1,1) + &
               h1(i,1,1,e) * u1(i,1,1)) * weights3(i,1,1)
          us_h = (svv_h1(i,1,1,e) * u2_svv(i,1,1) + &
               h1(i,1,1,e) * u2(i,1,1)) * weights3(i,1,1)
          ut_h = (svv_h1(i,1,1,e) * u3_svv(i,1,1) + &
               h1(i,1,1,e) * u3(i,1,1)) * weights3(i,1,1)

          wur(i,1,1) = drdx(i,1,1,e) * ur_h + &
               drdy(i,1,1,e) * us_h + drdz(i,1,1,e) * ut_h
          wus(i,1,1) = dsdx(i,1,1,e) * ur_h + &
               dsdy(i,1,1,e) * us_h + dsdz(i,1,1,e) * ut_h
          wut(i,1,1) = dtdx(i,1,1,e) * ur_h + &
               dtdy(i,1,1,e) * us_h + dtdz(i,1,1,e) * ut_h
       end do

       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + Dxt(i,k) * wur(k,j,1)
             end do
             w(i,j,1,e) = tmp
          end do
       end do
       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + Dyt(j,l) * wus(i,l,k)
                end do
                w(i,j,k,e) = w(i,j,k,e) + tmp
             end do
          end do
       end do
       do k = 1, lx
          do i = 1, lx * lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + Dzt(k,l) * wut(i,1,l)
             end do
             w(i,1,k,e) = w(i,1,k,e) + tmp
          end do
       end do
    end do
    !$omp end parallel do
  end subroutine ax_helm_svv_KS_lx

end module ax_helm_svv_KS_cpu
