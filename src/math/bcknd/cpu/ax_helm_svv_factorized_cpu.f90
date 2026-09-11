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
!> CPU implementation of the factorized SVV Helmholtz operator.
module ax_helm_svv_factorized_cpu
  use ax_helm_svv, only : ax_helm_svv_t
  use num_types, only : rp
  use coefs, only : coef_t
  use space, only : space_t
  use mesh, only : mesh_t
  use math, only : addcol4
  implicit none
  private

  !> CPU matrix-vector product for a Helmholtz problem.
  type, public, extends(ax_helm_svv_t) :: ax_helm_svv_factorized_cpu_t
   contains
     !> Compute the product.
     procedure, pass(this) :: compute => ax_helm_svv_factorized_compute
  end type ax_helm_svv_factorized_cpu_t

contains

  !> Compute the product.
  !! @param w Vector of size @a (lx,ly,lz,nelv).
  !! @param u Vector of size @a (lx,ly,lz,nelv).
  !! @param coef Coefficients.
  !! @param msh Mesh.
  !! @param Xh Function space \f$ X_h \f$.
  !! @note Since this is a performance-crtical routine, it is implemented in
  !! several kernels corresponding to different polynmial orders.
  subroutine ax_helm_svv_factorized_compute(this, w, u, coef, msh, Xh)
    class(ax_helm_svv_factorized_cpu_t), intent(in) :: this
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)

    call ax_helm_svv_factorized_lx(w, u, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, &
            this%svv%Br, this%svv%Bs, this%svv%Bt, coef%h1, this%svv%h1, &
            coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
            msh%nelv, Xh%lx)

    if (coef%ifh2) call addcol4 (w,coef%h2,coef%B,u,coef%dof%size())


  end subroutine ax_helm_svv_factorized_compute

  !> Generic CPU kernel for the Helmholz matrix-vector product.
  !! @param w Result.
  !! @param u Velocity field.
  !! @param Dx Derivative operator in first dimension.
  !! @param Dy Derivative operator in second dimension.
  !! @param Dz Derivative operator in third dimension.
  !! @param Dxt Derivative operator transpose in first dimension.
  !! @param Dyt Derivative operator transpose in second dimension.
  !! @param Dzt Derivative operator transpose in third dimension.
  !! @param Br Complementary derivative operator in first dimension.
  !! @param Bs Complementary derivative operator in second dimension.
  !! @param Bt Complementary derivative operator in third dimension.
  !! @param h1 Ordinary viscosity coefficient.
  !! @param svv_h1 SVV viscosity coefficient.
  !! @param G11 Geometric factor \f$G_{11}\f$.
  !! @param G22 Geometric factor \f$G_{22}\f$.
  !! @param G33 Geometric factor \f$G_{33}\f$.
  !! @param G12 Geometric factor \f$G_{12}\f$.
  !! @param G13 Geometric factor \f$G_{13}\f$.
  !! @param G23 Geometric factor \f$G_{23}\f$.
  !! @param n Number of elements.
  !! @param lx Polynomial order.
  subroutine ax_helm_svv_factorized_lx(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       Br, Bs, Bt, h1, svv_h1, G11, G22, G33, G12, G13, G23, n, lx)
    integer, intent(in) :: n, lx
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: svv_h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx), Dy(lx, lx), Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx), Dyt(lx, lx), Dzt(lx, lx)
    real(kind=rp), intent(in) :: Br(lx, lx), Bs(lx, lx), Bt(lx, lx)
    real(kind=rp) :: ur_h
    real(kind=rp) :: us_h
    real(kind=rp) :: ut_h
    real(kind=rp) :: u1_svv(lx, lx, lx)
    real(kind=rp) :: u2_svv(lx, lx, lx)
    real(kind=rp) :: u3_svv(lx, lx, lx)
    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    real(kind=rp) :: tmp, tmp_svv
    integer :: e, i, j, k, l

    do e = 1, n
       ! Ordinary and complementary reference-space derivatives, D u and B u.
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             tmp_svv = 0.0_rp
             do k = 1, lx
                tmp = tmp + Dx(i,k) * u(k,j,1,e)
                tmp_svv = tmp_svv + Br(i,k) * u(k,j,1,e)
             end do
             wur(i,j,1) = tmp
             u1_svv(i,j,1) = tmp_svv
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                tmp_svv = 0.0_rp
                do l = 1, lx
                   tmp = tmp + Dy(j,l) * u(i,l,k,e)
                   tmp_svv = tmp_svv + Bs(j,l) * u(i,l,k,e)
                end do
                wus(i,j,k) = tmp
                u2_svv(i,j,k) = tmp_svv
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             tmp_svv = 0.0_rp
             do l = 1, lx
                tmp = tmp + Dz(k,l) * u(i,1,l,e)
                tmp_svv = tmp_svv + Bt(k,l) * u(i,1,l,e)
             end do
             wut(i,1,k) = tmp
             u3_svv(i,1,k) = tmp_svv
          end do
       end do

       ! Construct the ordinary and SVV reference-space fluxes.
       do i = 1, lx*lx*lx
          ur_h = h1(i,1,1,e) * (G11(i,1,1,e) * wur(i,1,1) &
               + G12(i,1,1,e) * wus(i,1,1) &
               + G13(i,1,1,e) * wut(i,1,1))
          us_h = h1(i,1,1,e) * (G12(i,1,1,e) * wur(i,1,1) &
               + G22(i,1,1,e) * wus(i,1,1) &
               + G23(i,1,1,e) * wut(i,1,1))
          ut_h = h1(i,1,1,e) * (G13(i,1,1,e) * wur(i,1,1) &
               + G23(i,1,1,e) * wus(i,1,1) &
               + G33(i,1,1,e) * wut(i,1,1))
          wur(i,1,1) = ur_h
          wus(i,1,1) = us_h
          wut(i,1,1) = ut_h

          ur_h = svv_h1(i,1,1,e) * (G11(i,1,1,e) * u1_svv(i,1,1) &
               + G12(i,1,1,e) * u2_svv(i,1,1) &
               + G13(i,1,1,e) * u3_svv(i,1,1))
          us_h = svv_h1(i,1,1,e) * (G12(i,1,1,e) * u1_svv(i,1,1) &
               + G22(i,1,1,e) * u2_svv(i,1,1) &
               + G23(i,1,1,e) * u3_svv(i,1,1))
          ut_h = svv_h1(i,1,1,e) * (G13(i,1,1,e) * u1_svv(i,1,1) &
               + G23(i,1,1,e) * u2_svv(i,1,1) &
               + G33(i,1,1,e) * u3_svv(i,1,1))
          u1_svv(i,1,1) = ur_h
          u2_svv(i,1,1) = us_h
          u3_svv(i,1,1) = ut_h
       end do

       do j = 1, lx*lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + Dxt(i,k) * wur(k,j,1) &
                     + Br(k,i) * u1_svv(k,j,1)
             end do
             w(i,j,1,e) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + Dyt(j,l) * wus(i,l,k) &
                        + Bs(l,j) * u2_svv(i,l,k)
                end do
                w(i,j,k,e) = w(i,j,k,e) + tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + Dzt(k,l) * wut(i,1,l) &
                     + Bt(l,k) * u3_svv(i,1,l)
             end do
             w(i,1,k,e) = w(i,1,k,e) + tmp
          end do
       end do

    end do
  end subroutine ax_helm_svv_factorized_lx

end module ax_helm_svv_factorized_cpu
