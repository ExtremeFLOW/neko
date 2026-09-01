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
submodule(ax_helm_cpu) ax_helm_vector_cpu
  implicit none

contains

  !> Compute \f$ Ax \f$ product, taking three
  !! components of a vector field in an uncoupled manner.
  !! @param au Result for the first component of the vector.
  !! @param av Result for the first component of the vector.
  !! @param aw Result for the first component of the vector.
  !! @param u The first component of the vector.
  !! @param v The second component of the vector.
  !! @param w The third component of the vector.
  !! @param coef Coefficients.
  !! @param msh Mesh.
  !! @param Xh Function space \f$ X_h \f$.
  !! @note Since this is a performance-crtical routine, it is implemented in
  !! several kernels corresponding to different polynmial orders.
  module subroutine ax_helm_compute_vector(this, au, av, aw, &
       u, v, w, coef, msh, Xh)
    class(ax_helm_cpu_t), intent(in) :: this
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: au(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: av(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: aw(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: v(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    integer :: i

    !$omp parallel
    select case(Xh%lx)
    case (14)
       call ax_helm_vector_lx14(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%G11, coef%G22, coef%G33, &
            coef%G12, coef%G13, coef%G23, msh%nelv)
    case (13)
       call ax_helm_vector_lx13(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%G11, coef%G22, coef%G33, &
            coef%G12, coef%G13, coef%G23, msh%nelv)
    case (12)
       call ax_helm_vector_lx12(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%G11, coef%G22, coef%G33, &
            coef%G12, coef%G13, coef%G23, msh%nelv)
    case (11)
       call ax_helm_vector_lx11(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%G11, coef%G22, coef%G33, &
            coef%G12, coef%G13, coef%G23, msh%nelv)
    case (10)
       call ax_helm_vector_lx10(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%G11, coef%G22, coef%G33, &
            coef%G12, coef%G13, coef%G23, msh%nelv)
    case (9)
       call ax_helm_vector_lx9(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%G11, coef%G22, coef%G33, &
            coef%G12, coef%G13, coef%G23, msh%nelv)
    case (8)
       call ax_helm_vector_lx8(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%G11, coef%G22, coef%G33, &
            coef%G12, coef%G13, coef%G23, msh%nelv)
    case (7)
       call ax_helm_vector_lx7(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%G11, coef%G22, coef%G33, &
            coef%G12, coef%G13, coef%G23, msh%nelv)
    case (6)
       call ax_helm_vector_lx6(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%G11, coef%G22, coef%G33, &
            coef%G12, coef%G13, coef%G23, msh%nelv)
    case (5)
       call ax_helm_vector_lx5(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%G11, coef%G22, coef%G33, &
            coef%G12, coef%G13, coef%G23, msh%nelv)
    case (4)
       call ax_helm_vector_lx4(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%G11, coef%G22, coef%G33, &
            coef%G12, coef%G13, coef%G23, msh%nelv)
    case (3)
       call ax_helm_vector_lx3(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%G11, coef%G22, coef%G33, &
            coef%G12, coef%G13, coef%G23, msh%nelv)
    case (2)
       call ax_helm_vector_lx2(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%G11, coef%G22, coef%G33, &
            coef%G12, coef%G13, coef%G23, msh%nelv)
    case default
       call ax_helm_vector_lx(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%G11, coef%G22, coef%G33, &
            coef%G12, coef%G13, coef%G23, msh%nelv, Xh%lx)
    end select

    if (coef%ifh2) then
       !$omp do private(i)
       do i = 1, coef%dof%size()
          au(i,1,1,1) = au(i,1,1,1) + &
               coef%h2(i,1,1,1) * coef%B(i,1,1,1) * u(i,1,1,1)
          av(i,1,1,1) = av(i,1,1,1) + &
               coef%h2(i,1,1,1) * coef%B(i,1,1,1) * v(i,1,1,1)
          aw(i,1,1,1) = aw(i,1,1,1) + &
               coef%h2(i,1,1,1) * coef%B(i,1,1,1) * w(i,1,1,1)
       end do
       !$omp end do
    end if
    !$omp end parallel

  end subroutine ax_helm_compute_vector

  !> Generic CPU kernel for the Helmholz matrix-vector product.
  !! @param w Result.
  !! @param u Velocity field.
  !! @param Dx Derivative operator in first dimension.
  !! @param Dy Derivative operator in second dimension.
  !! @param Dz Derivative operator in third dimension.
  !! @param Dxt Derivative operator transpose in first dimension.
  !! @param Dyt Derivative operator transpose in second dimension.
  !! @param Dzt Derivative operator transpose in third dimension.
  !! @param G11 Geometric factor.
  !! @param n Number of elements.
  !! @param lx Polynomial order.
  subroutine ax_helm_vector_lx(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n, lx)
    integer, intent(in) :: n, lx
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
    real(kind=rp) :: wur(lx, lx, lx)
    real(kind=rp) :: wvr(lx, lx, lx)
    real(kind=rp) :: wwr(lx, lx, lx)
    real(kind=rp) :: wus(lx, lx, lx)
    real(kind=rp) :: wvs(lx, lx, lx)
    real(kind=rp) :: wws(lx, lx, lx)
    real(kind=rp) :: wut(lx, lx, lx)
    real(kind=rp) :: wvt(lx, lx, lx)
    real(kind=rp) :: wwt(lx, lx, lx)
    real(kind=rp) :: t1, t2, t3
    integer :: e, i, j, k, l

    !$omp do
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
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )

          vr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wvr(i,1,1) &
                      + G12(i,1,1,e) * wvs(i,1,1) &
                      + G13(i,1,1,e) * wvt(i,1,1) )
          vs(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wvr(i,1,1) &
                      + G22(i,1,1,e) * wvs(i,1,1) &
                      + G23(i,1,1,e) * wvt(i,1,1) )
          vt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wvr(i,1,1) &
                      + G23(i,1,1,e) * wvs(i,1,1) &
                      + G33(i,1,1,e) * wvt(i,1,1) )

          wr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wwr(i,1,1) &
                      + G12(i,1,1,e) * wws(i,1,1) &
                      + G13(i,1,1,e) * wwt(i,1,1) )
          ws(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wwr(i,1,1) &
                      + G22(i,1,1,e) * wws(i,1,1) &
                      + G23(i,1,1,e) * wwt(i,1,1) )
          wt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wwr(i,1,1) &
                      + G23(i,1,1,e) * wws(i,1,1) &
                      + G33(i,1,1,e) * wwt(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             t1 = 0.0_rp
             t2 = 0.0_rp
             t3 = 0.0_rp
             do k = 1, lx
                t1 = t1 + Dxt(i,k) * ur(k,j,1)
                t2 = t2 + Dxt(i,k) * vr(k,j,1)
                t3 = t3 + Dxt(i,k) * wr(k,j,1)
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
                   t1 = t1 + Dyt(j,l) * us(i,l,k)
                   t2 = t2 + Dyt(j,l) * vs(i,l,k)
                   t3 = t3 + Dyt(j,l) * ws(i,l,k)
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
                t1 = t1 + Dzt(k,l) * ut(i,1,l)
                t2 = t2 + Dzt(k,l) * vt(i,1,l)
                t3 = t3 + Dzt(k,l) * wt(i,1,l)
             end do
             au(i,1,k,e) = au(i,1,k,e) + t1
             av(i,1,k,e) = av(i,1,k,e) + t2
             aw(i,1,k,e) = aw(i,1,k,e) + t3
          end do
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx

  subroutine ax_helm_vector_lx14(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 14
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
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

    !$omp do
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
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )

          vr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wvr(i,1,1) &
                      + G12(i,1,1,e) * wvs(i,1,1) &
                      + G13(i,1,1,e) * wvt(i,1,1) )
          vs(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wvr(i,1,1) &
                      + G22(i,1,1,e) * wvs(i,1,1) &
                      + G23(i,1,1,e) * wvt(i,1,1) )
          vt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wvr(i,1,1) &
                      + G23(i,1,1,e) * wvs(i,1,1) &
                      + G33(i,1,1,e) * wvt(i,1,1) )

          wr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wwr(i,1,1) &
                      + G12(i,1,1,e) * wws(i,1,1) &
                      + G13(i,1,1,e) * wwt(i,1,1) )
          ws(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wwr(i,1,1) &
                      + G22(i,1,1,e) * wws(i,1,1) &
                      + G23(i,1,1,e) * wwt(i,1,1) )
          wt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wwr(i,1,1) &
                      + G23(i,1,1,e) * wws(i,1,1) &
                      + G33(i,1,1,e) * wwt(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                         + Dxt(i,2) * ur(2,j,1) &
                         + Dxt(i,3) * ur(3,j,1) &
                         + Dxt(i,4) * ur(4,j,1) &
                         + Dxt(i,5) * ur(5,j,1) &
                         + Dxt(i,6) * ur(6,j,1) &
                         + Dxt(i,7) * ur(7,j,1) &
                         + Dxt(i,8) * ur(8,j,1) &
                         + Dxt(i,9) * ur(9,j,1) &
                         + Dxt(i,10) * ur(10,j,1) &
                         + Dxt(i,11) * ur(11,j,1) &
                         + Dxt(i,12) * ur(12,j,1) &
                         + Dxt(i,13) * ur(13,j,1) &
                         + Dxt(i,14) * ur(14,j,1)

             av(i,j,1,e) = Dxt(i,1) * vr(1,j,1) &
                         + Dxt(i,2) * vr(2,j,1) &
                         + Dxt(i,3) * vr(3,j,1) &
                         + Dxt(i,4) * vr(4,j,1) &
                         + Dxt(i,5) * vr(5,j,1) &
                         + Dxt(i,6) * vr(6,j,1) &
                         + Dxt(i,7) * vr(7,j,1) &
                         + Dxt(i,8) * vr(8,j,1) &
                         + Dxt(i,9) * vr(9,j,1) &
                         + Dxt(i,10) * vr(10,j,1) &
                         + Dxt(i,11) * vr(11,j,1) &
                         + Dxt(i,12) * vr(12,j,1) &
                         + Dxt(i,13) * vr(13,j,1) &
                         + Dxt(i,14) * vr(14,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wr(1,j,1) &
                         + Dxt(i,2) * wr(2,j,1) &
                         + Dxt(i,3) * wr(3,j,1) &
                         + Dxt(i,4) * wr(4,j,1) &
                         + Dxt(i,5) * wr(5,j,1) &
                         + Dxt(i,6) * wr(6,j,1) &
                         + Dxt(i,7) * wr(7,j,1) &
                         + Dxt(i,8) * wr(8,j,1) &
                         + Dxt(i,9) * wr(9,j,1) &
                         + Dxt(i,10) * wr(10,j,1) &
                         + Dxt(i,11) * wr(11,j,1) &
                         + Dxt(i,12) * wr(12,j,1) &
                         + Dxt(i,13) * wr(13,j,1) &
                         + Dxt(i,14) * wr(14,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * us(i,1,k) &
                            + Dyt(j,2) * us(i,2,k) &
                            + Dyt(j,3) * us(i,3,k) &
                            + Dyt(j,4) * us(i,4,k) &
                            + Dyt(j,5) * us(i,5,k) &
                            + Dyt(j,6) * us(i,6,k) &
                            + Dyt(j,7) * us(i,7,k) &
                            + Dyt(j,8) * us(i,8,k) &
                            + Dyt(j,9) * us(i,9,k) &
                            + Dyt(j,10) * us(i,10,k) &
                            + Dyt(j,11) * us(i,11,k) &
                            + Dyt(j,12) * us(i,12,k) &
                            + Dyt(j,13) * us(i,13,k) &
                            + Dyt(j,14) * us(i,14,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * vs(i,1,k) &
                            + Dyt(j,2) * vs(i,2,k) &
                            + Dyt(j,3) * vs(i,3,k) &
                            + Dyt(j,4) * vs(i,4,k) &
                            + Dyt(j,5) * vs(i,5,k) &
                            + Dyt(j,6) * vs(i,6,k) &
                            + Dyt(j,7) * vs(i,7,k) &
                            + Dyt(j,8) * vs(i,8,k) &
                            + Dyt(j,9) * vs(i,9,k) &
                            + Dyt(j,10) * vs(i,10,k) &
                            + Dyt(j,11) * vs(i,11,k) &
                            + Dyt(j,12) * vs(i,12,k) &
                            + Dyt(j,13) * vs(i,13,k) &
                            + Dyt(j,14) * vs(i,14,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * ws(i,1,k) &
                            + Dyt(j,2) * ws(i,2,k) &
                            + Dyt(j,3) * ws(i,3,k) &
                            + Dyt(j,4) * ws(i,4,k) &
                            + Dyt(j,5) * ws(i,5,k) &
                            + Dyt(j,6) * ws(i,6,k) &
                            + Dyt(j,7) * ws(i,7,k) &
                            + Dyt(j,8) * ws(i,8,k) &
                            + Dyt(j,9) * ws(i,9,k) &
                            + Dyt(j,10) * ws(i,10,k) &
                            + Dyt(j,11) * ws(i,11,k) &
                            + Dyt(j,12) * ws(i,12,k) &
                            + Dyt(j,13) * ws(i,13,k) &
                            + Dyt(j,14) * ws(i,14,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * ut(i,1,1) &
                         + Dzt(k,2) * ut(i,1,2) &
                         + Dzt(k,3) * ut(i,1,3) &
                         + Dzt(k,4) * ut(i,1,4) &
                         + Dzt(k,5) * ut(i,1,5) &
                         + Dzt(k,6) * ut(i,1,6) &
                         + Dzt(k,7) * ut(i,1,7) &
                         + Dzt(k,8) * ut(i,1,8) &
                         + Dzt(k,9) * ut(i,1,9) &
                         + Dzt(k,10) * ut(i,1,10) &
                         + Dzt(k,11) * ut(i,1,11) &
                         + Dzt(k,12) * ut(i,1,12) &
                         + Dzt(k,13) * ut(i,1,13) &
                         + Dzt(k,14) * ut(i,1,14)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * vt(i,1,1) &
                         + Dzt(k,2) * vt(i,1,2) &
                         + Dzt(k,3) * vt(i,1,3) &
                         + Dzt(k,4) * vt(i,1,4) &
                         + Dzt(k,5) * vt(i,1,5) &
                         + Dzt(k,6) * vt(i,1,6) &
                         + Dzt(k,7) * vt(i,1,7) &
                         + Dzt(k,8) * vt(i,1,8) &
                         + Dzt(k,9) * vt(i,1,9) &
                         + Dzt(k,10) * vt(i,1,10) &
                         + Dzt(k,11) * vt(i,1,11) &
                         + Dzt(k,12) * vt(i,1,12) &
                         + Dzt(k,13) * vt(i,1,13) &
                         + Dzt(k,14) * vt(i,1,14)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wt(i,1,1) &
                         + Dzt(k,2) * wt(i,1,2) &
                         + Dzt(k,3) * wt(i,1,3) &
                         + Dzt(k,4) * wt(i,1,4) &
                         + Dzt(k,5) * wt(i,1,5) &
                         + Dzt(k,6) * wt(i,1,6) &
                         + Dzt(k,7) * wt(i,1,7) &
                         + Dzt(k,8) * wt(i,1,8) &
                         + Dzt(k,9) * wt(i,1,9) &
                         + Dzt(k,10) * wt(i,1,10) &
                         + Dzt(k,11) * wt(i,1,11) &
                         + Dzt(k,12) * wt(i,1,12) &
                         + Dzt(k,13) * wt(i,1,13) &
                         + Dzt(k,14) * wt(i,1,14)
          end do
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx14

  subroutine ax_helm_vector_lx13(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 13
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
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

    !$omp do
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
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )

          vr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wvr(i,1,1) &
                      + G12(i,1,1,e) * wvs(i,1,1) &
                      + G13(i,1,1,e) * wvt(i,1,1) )
          vs(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wvr(i,1,1) &
                      + G22(i,1,1,e) * wvs(i,1,1) &
                      + G23(i,1,1,e) * wvt(i,1,1) )
          vt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wvr(i,1,1) &
                      + G23(i,1,1,e) * wvs(i,1,1) &
                      + G33(i,1,1,e) * wvt(i,1,1) )

          wr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wwr(i,1,1) &
                      + G12(i,1,1,e) * wws(i,1,1) &
                      + G13(i,1,1,e) * wwt(i,1,1) )
          ws(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wwr(i,1,1) &
                      + G22(i,1,1,e) * wws(i,1,1) &
                      + G23(i,1,1,e) * wwt(i,1,1) )
          wt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wwr(i,1,1) &
                      + G23(i,1,1,e) * wws(i,1,1) &
                      + G33(i,1,1,e) * wwt(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                         + Dxt(i,2) * ur(2,j,1) &
                         + Dxt(i,3) * ur(3,j,1) &
                         + Dxt(i,4) * ur(4,j,1) &
                         + Dxt(i,5) * ur(5,j,1) &
                         + Dxt(i,6) * ur(6,j,1) &
                         + Dxt(i,7) * ur(7,j,1) &
                         + Dxt(i,8) * ur(8,j,1) &
                         + Dxt(i,9) * ur(9,j,1) &
                         + Dxt(i,10) * ur(10,j,1) &
                         + Dxt(i,11) * ur(11,j,1) &
                         + Dxt(i,12) * ur(12,j,1) &
                         + Dxt(i,13) * ur(13,j,1)

             av(i,j,1,e) = Dxt(i,1) * vr(1,j,1) &
                         + Dxt(i,2) * vr(2,j,1) &
                         + Dxt(i,3) * vr(3,j,1) &
                         + Dxt(i,4) * vr(4,j,1) &
                         + Dxt(i,5) * vr(5,j,1) &
                         + Dxt(i,6) * vr(6,j,1) &
                         + Dxt(i,7) * vr(7,j,1) &
                         + Dxt(i,8) * vr(8,j,1) &
                         + Dxt(i,9) * vr(9,j,1) &
                         + Dxt(i,10) * vr(10,j,1) &
                         + Dxt(i,11) * vr(11,j,1) &
                         + Dxt(i,12) * vr(12,j,1) &
                         + Dxt(i,13) * vr(13,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wr(1,j,1) &
                         + Dxt(i,2) * wr(2,j,1) &
                         + Dxt(i,3) * wr(3,j,1) &
                         + Dxt(i,4) * wr(4,j,1) &
                         + Dxt(i,5) * wr(5,j,1) &
                         + Dxt(i,6) * wr(6,j,1) &
                         + Dxt(i,7) * wr(7,j,1) &
                         + Dxt(i,8) * wr(8,j,1) &
                         + Dxt(i,9) * wr(9,j,1) &
                         + Dxt(i,10) * wr(10,j,1) &
                         + Dxt(i,11) * wr(11,j,1) &
                         + Dxt(i,12) * wr(12,j,1) &
                         + Dxt(i,13) * wr(13,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * us(i,1,k) &
                            + Dyt(j,2) * us(i,2,k) &
                            + Dyt(j,3) * us(i,3,k) &
                            + Dyt(j,4) * us(i,4,k) &
                            + Dyt(j,5) * us(i,5,k) &
                            + Dyt(j,6) * us(i,6,k) &
                            + Dyt(j,7) * us(i,7,k) &
                            + Dyt(j,8) * us(i,8,k) &
                            + Dyt(j,9) * us(i,9,k) &
                            + Dyt(j,10) * us(i,10,k) &
                            + Dyt(j,11) * us(i,11,k) &
                            + Dyt(j,12) * us(i,12,k) &
                            + Dyt(j,13) * us(i,13,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * vs(i,1,k) &
                            + Dyt(j,2) * vs(i,2,k) &
                            + Dyt(j,3) * vs(i,3,k) &
                            + Dyt(j,4) * vs(i,4,k) &
                            + Dyt(j,5) * vs(i,5,k) &
                            + Dyt(j,6) * vs(i,6,k) &
                            + Dyt(j,7) * vs(i,7,k) &
                            + Dyt(j,8) * vs(i,8,k) &
                            + Dyt(j,9) * vs(i,9,k) &
                            + Dyt(j,10) * vs(i,10,k) &
                            + Dyt(j,11) * vs(i,11,k) &
                            + Dyt(j,12) * vs(i,12,k) &
                            + Dyt(j,13) * vs(i,13,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * ws(i,1,k) &
                            + Dyt(j,2) * ws(i,2,k) &
                            + Dyt(j,3) * ws(i,3,k) &
                            + Dyt(j,4) * ws(i,4,k) &
                            + Dyt(j,5) * ws(i,5,k) &
                            + Dyt(j,6) * ws(i,6,k) &
                            + Dyt(j,7) * ws(i,7,k) &
                            + Dyt(j,8) * ws(i,8,k) &
                            + Dyt(j,9) * ws(i,9,k) &
                            + Dyt(j,10) * ws(i,10,k) &
                            + Dyt(j,11) * ws(i,11,k) &
                            + Dyt(j,12) * ws(i,12,k) &
                            + Dyt(j,13) * ws(i,13,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * ut(i,1,1) &
                         + Dzt(k,2) * ut(i,1,2) &
                         + Dzt(k,3) * ut(i,1,3) &
                         + Dzt(k,4) * ut(i,1,4) &
                         + Dzt(k,5) * ut(i,1,5) &
                         + Dzt(k,6) * ut(i,1,6) &
                         + Dzt(k,7) * ut(i,1,7) &
                         + Dzt(k,8) * ut(i,1,8) &
                         + Dzt(k,9) * ut(i,1,9) &
                         + Dzt(k,10) * ut(i,1,10) &
                         + Dzt(k,11) * ut(i,1,11) &
                         + Dzt(k,12) * ut(i,1,12) &
                         + Dzt(k,13) * ut(i,1,13)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * vt(i,1,1) &
                         + Dzt(k,2) * vt(i,1,2) &
                         + Dzt(k,3) * vt(i,1,3) &
                         + Dzt(k,4) * vt(i,1,4) &
                         + Dzt(k,5) * vt(i,1,5) &
                         + Dzt(k,6) * vt(i,1,6) &
                         + Dzt(k,7) * vt(i,1,7) &
                         + Dzt(k,8) * vt(i,1,8) &
                         + Dzt(k,9) * vt(i,1,9) &
                         + Dzt(k,10) * vt(i,1,10) &
                         + Dzt(k,11) * vt(i,1,11) &
                         + Dzt(k,12) * vt(i,1,12) &
                         + Dzt(k,13) * vt(i,1,13)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wt(i,1,1) &
                         + Dzt(k,2) * wt(i,1,2) &
                         + Dzt(k,3) * wt(i,1,3) &
                         + Dzt(k,4) * wt(i,1,4) &
                         + Dzt(k,5) * wt(i,1,5) &
                         + Dzt(k,6) * wt(i,1,6) &
                         + Dzt(k,7) * wt(i,1,7) &
                         + Dzt(k,8) * wt(i,1,8) &
                         + Dzt(k,9) * wt(i,1,9) &
                         + Dzt(k,10) * wt(i,1,10) &
                         + Dzt(k,11) * wt(i,1,11) &
                         + Dzt(k,12) * wt(i,1,12) &
                         + Dzt(k,13) * wt(i,1,13)
          end do
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx13

  subroutine ax_helm_vector_lx12(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 12
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
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

    !$omp do
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
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )

          vr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wvr(i,1,1) &
                      + G12(i,1,1,e) * wvs(i,1,1) &
                      + G13(i,1,1,e) * wvt(i,1,1) )
          vs(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wvr(i,1,1) &
                      + G22(i,1,1,e) * wvs(i,1,1) &
                      + G23(i,1,1,e) * wvt(i,1,1) )
          vt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wvr(i,1,1) &
                      + G23(i,1,1,e) * wvs(i,1,1) &
                      + G33(i,1,1,e) * wvt(i,1,1) )

          wr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wwr(i,1,1) &
                      + G12(i,1,1,e) * wws(i,1,1) &
                      + G13(i,1,1,e) * wwt(i,1,1) )
          ws(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wwr(i,1,1) &
                      + G22(i,1,1,e) * wws(i,1,1) &
                      + G23(i,1,1,e) * wwt(i,1,1) )
          wt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wwr(i,1,1) &
                      + G23(i,1,1,e) * wws(i,1,1) &
                      + G33(i,1,1,e) * wwt(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                         + Dxt(i,2) * ur(2,j,1) &
                         + Dxt(i,3) * ur(3,j,1) &
                         + Dxt(i,4) * ur(4,j,1) &
                         + Dxt(i,5) * ur(5,j,1) &
                         + Dxt(i,6) * ur(6,j,1) &
                         + Dxt(i,7) * ur(7,j,1) &
                         + Dxt(i,8) * ur(8,j,1) &
                         + Dxt(i,9) * ur(9,j,1) &
                         + Dxt(i,10) * ur(10,j,1) &
                         + Dxt(i,11) * ur(11,j,1) &
                         + Dxt(i,12) * ur(12,j,1)

             av(i,j,1,e) = Dxt(i,1) * vr(1,j,1) &
                         + Dxt(i,2) * vr(2,j,1) &
                         + Dxt(i,3) * vr(3,j,1) &
                         + Dxt(i,4) * vr(4,j,1) &
                         + Dxt(i,5) * vr(5,j,1) &
                         + Dxt(i,6) * vr(6,j,1) &
                         + Dxt(i,7) * vr(7,j,1) &
                         + Dxt(i,8) * vr(8,j,1) &
                         + Dxt(i,9) * vr(9,j,1) &
                         + Dxt(i,10) * vr(10,j,1) &
                         + Dxt(i,11) * vr(11,j,1) &
                         + Dxt(i,12) * vr(12,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wr(1,j,1) &
                         + Dxt(i,2) * wr(2,j,1) &
                         + Dxt(i,3) * wr(3,j,1) &
                         + Dxt(i,4) * wr(4,j,1) &
                         + Dxt(i,5) * wr(5,j,1) &
                         + Dxt(i,6) * wr(6,j,1) &
                         + Dxt(i,7) * wr(7,j,1) &
                         + Dxt(i,8) * wr(8,j,1) &
                         + Dxt(i,9) * wr(9,j,1) &
                         + Dxt(i,10) * wr(10,j,1) &
                         + Dxt(i,11) * wr(11,j,1) &
                         + Dxt(i,12) * wr(12,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * us(i,1,k) &
                            + Dyt(j,2) * us(i,2,k) &
                            + Dyt(j,3) * us(i,3,k) &
                            + Dyt(j,4) * us(i,4,k) &
                            + Dyt(j,5) * us(i,5,k) &
                            + Dyt(j,6) * us(i,6,k) &
                            + Dyt(j,7) * us(i,7,k) &
                            + Dyt(j,8) * us(i,8,k) &
                            + Dyt(j,9) * us(i,9,k) &
                            + Dyt(j,10) * us(i,10,k) &
                            + Dyt(j,11) * us(i,11,k) &
                            + Dyt(j,12) * us(i,12,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * vs(i,1,k) &
                            + Dyt(j,2) * vs(i,2,k) &
                            + Dyt(j,3) * vs(i,3,k) &
                            + Dyt(j,4) * vs(i,4,k) &
                            + Dyt(j,5) * vs(i,5,k) &
                            + Dyt(j,6) * vs(i,6,k) &
                            + Dyt(j,7) * vs(i,7,k) &
                            + Dyt(j,8) * vs(i,8,k) &
                            + Dyt(j,9) * vs(i,9,k) &
                            + Dyt(j,10) * vs(i,10,k) &
                            + Dyt(j,11) * vs(i,11,k) &
                            + Dyt(j,12) * vs(i,12,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * ws(i,1,k) &
                            + Dyt(j,2) * ws(i,2,k) &
                            + Dyt(j,3) * ws(i,3,k) &
                            + Dyt(j,4) * ws(i,4,k) &
                            + Dyt(j,5) * ws(i,5,k) &
                            + Dyt(j,6) * ws(i,6,k) &
                            + Dyt(j,7) * ws(i,7,k) &
                            + Dyt(j,8) * ws(i,8,k) &
                            + Dyt(j,9) * ws(i,9,k) &
                            + Dyt(j,10) * ws(i,10,k) &
                            + Dyt(j,11) * ws(i,11,k) &
                            + Dyt(j,12) * ws(i,12,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * ut(i,1,1) &
                         + Dzt(k,2) * ut(i,1,2) &
                         + Dzt(k,3) * ut(i,1,3) &
                         + Dzt(k,4) * ut(i,1,4) &
                         + Dzt(k,5) * ut(i,1,5) &
                         + Dzt(k,6) * ut(i,1,6) &
                         + Dzt(k,7) * ut(i,1,7) &
                         + Dzt(k,8) * ut(i,1,8) &
                         + Dzt(k,9) * ut(i,1,9) &
                         + Dzt(k,10) * ut(i,1,10) &
                         + Dzt(k,11) * ut(i,1,11) &
                         + Dzt(k,12) * ut(i,1,12)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * vt(i,1,1) &
                         + Dzt(k,2) * vt(i,1,2) &
                         + Dzt(k,3) * vt(i,1,3) &
                         + Dzt(k,4) * vt(i,1,4) &
                         + Dzt(k,5) * vt(i,1,5) &
                         + Dzt(k,6) * vt(i,1,6) &
                         + Dzt(k,7) * vt(i,1,7) &
                         + Dzt(k,8) * vt(i,1,8) &
                         + Dzt(k,9) * vt(i,1,9) &
                         + Dzt(k,10) * vt(i,1,10) &
                         + Dzt(k,11) * vt(i,1,11) &
                         + Dzt(k,12) * vt(i,1,12)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wt(i,1,1) &
                         + Dzt(k,2) * wt(i,1,2) &
                         + Dzt(k,3) * wt(i,1,3) &
                         + Dzt(k,4) * wt(i,1,4) &
                         + Dzt(k,5) * wt(i,1,5) &
                         + Dzt(k,6) * wt(i,1,6) &
                         + Dzt(k,7) * wt(i,1,7) &
                         + Dzt(k,8) * wt(i,1,8) &
                         + Dzt(k,9) * wt(i,1,9) &
                         + Dzt(k,10) * wt(i,1,10) &
                         + Dzt(k,11) * wt(i,1,11) &
                         + Dzt(k,12) * wt(i,1,12)
          end do
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx12

  subroutine ax_helm_vector_lx11(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 11
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
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

    !$omp do
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
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )

          vr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wvr(i,1,1) &
                      + G12(i,1,1,e) * wvs(i,1,1) &
                      + G13(i,1,1,e) * wvt(i,1,1) )
          vs(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wvr(i,1,1) &
                      + G22(i,1,1,e) * wvs(i,1,1) &
                      + G23(i,1,1,e) * wvt(i,1,1) )
          vt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wvr(i,1,1) &
                      + G23(i,1,1,e) * wvs(i,1,1) &
                      + G33(i,1,1,e) * wvt(i,1,1) )

          wr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wwr(i,1,1) &
                      + G12(i,1,1,e) * wws(i,1,1) &
                      + G13(i,1,1,e) * wwt(i,1,1) )
          ws(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wwr(i,1,1) &
                      + G22(i,1,1,e) * wws(i,1,1) &
                      + G23(i,1,1,e) * wwt(i,1,1) )
          wt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wwr(i,1,1) &
                      + G23(i,1,1,e) * wws(i,1,1) &
                      + G33(i,1,1,e) * wwt(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                         + Dxt(i,2) * ur(2,j,1) &
                         + Dxt(i,3) * ur(3,j,1) &
                         + Dxt(i,4) * ur(4,j,1) &
                         + Dxt(i,5) * ur(5,j,1) &
                         + Dxt(i,6) * ur(6,j,1) &
                         + Dxt(i,7) * ur(7,j,1) &
                         + Dxt(i,8) * ur(8,j,1) &
                         + Dxt(i,9) * ur(9,j,1) &
                         + Dxt(i,10) * ur(10,j,1) &
                         + Dxt(i,11) * ur(11,j,1)

             av(i,j,1,e) = Dxt(i,1) * vr(1,j,1) &
                         + Dxt(i,2) * vr(2,j,1) &
                         + Dxt(i,3) * vr(3,j,1) &
                         + Dxt(i,4) * vr(4,j,1) &
                         + Dxt(i,5) * vr(5,j,1) &
                         + Dxt(i,6) * vr(6,j,1) &
                         + Dxt(i,7) * vr(7,j,1) &
                         + Dxt(i,8) * vr(8,j,1) &
                         + Dxt(i,9) * vr(9,j,1) &
                         + Dxt(i,10) * vr(10,j,1) &
                         + Dxt(i,11) * vr(11,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wr(1,j,1) &
                         + Dxt(i,2) * wr(2,j,1) &
                         + Dxt(i,3) * wr(3,j,1) &
                         + Dxt(i,4) * wr(4,j,1) &
                         + Dxt(i,5) * wr(5,j,1) &
                         + Dxt(i,6) * wr(6,j,1) &
                         + Dxt(i,7) * wr(7,j,1) &
                         + Dxt(i,8) * wr(8,j,1) &
                         + Dxt(i,9) * wr(9,j,1) &
                         + Dxt(i,10) * wr(10,j,1) &
                         + Dxt(i,11) * wr(11,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * us(i,1,k) &
                            + Dyt(j,2) * us(i,2,k) &
                            + Dyt(j,3) * us(i,3,k) &
                            + Dyt(j,4) * us(i,4,k) &
                            + Dyt(j,5) * us(i,5,k) &
                            + Dyt(j,6) * us(i,6,k) &
                            + Dyt(j,7) * us(i,7,k) &
                            + Dyt(j,8) * us(i,8,k) &
                            + Dyt(j,9) * us(i,9,k) &
                            + Dyt(j,10) * us(i,10,k) &
                            + Dyt(j,11) * us(i,11,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * vs(i,1,k) &
                            + Dyt(j,2) * vs(i,2,k) &
                            + Dyt(j,3) * vs(i,3,k) &
                            + Dyt(j,4) * vs(i,4,k) &
                            + Dyt(j,5) * vs(i,5,k) &
                            + Dyt(j,6) * vs(i,6,k) &
                            + Dyt(j,7) * vs(i,7,k) &
                            + Dyt(j,8) * vs(i,8,k) &
                            + Dyt(j,9) * vs(i,9,k) &
                            + Dyt(j,10) * vs(i,10,k) &
                            + Dyt(j,11) * vs(i,11,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * ws(i,1,k) &
                            + Dyt(j,2) * ws(i,2,k) &
                            + Dyt(j,3) * ws(i,3,k) &
                            + Dyt(j,4) * ws(i,4,k) &
                            + Dyt(j,5) * ws(i,5,k) &
                            + Dyt(j,6) * ws(i,6,k) &
                            + Dyt(j,7) * ws(i,7,k) &
                            + Dyt(j,8) * ws(i,8,k) &
                            + Dyt(j,9) * ws(i,9,k) &
                            + Dyt(j,10) * ws(i,10,k) &
                            + Dyt(j,11) * ws(i,11,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * ut(i,1,1) &
                         + Dzt(k,2) * ut(i,1,2) &
                         + Dzt(k,3) * ut(i,1,3) &
                         + Dzt(k,4) * ut(i,1,4) &
                         + Dzt(k,5) * ut(i,1,5) &
                         + Dzt(k,6) * ut(i,1,6) &
                         + Dzt(k,7) * ut(i,1,7) &
                         + Dzt(k,8) * ut(i,1,8) &
                         + Dzt(k,9) * ut(i,1,9) &
                         + Dzt(k,10) * ut(i,1,10) &
                         + Dzt(k,11) * ut(i,1,11)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * vt(i,1,1) &
                         + Dzt(k,2) * vt(i,1,2) &
                         + Dzt(k,3) * vt(i,1,3) &
                         + Dzt(k,4) * vt(i,1,4) &
                         + Dzt(k,5) * vt(i,1,5) &
                         + Dzt(k,6) * vt(i,1,6) &
                         + Dzt(k,7) * vt(i,1,7) &
                         + Dzt(k,8) * vt(i,1,8) &
                         + Dzt(k,9) * vt(i,1,9) &
                         + Dzt(k,10) * vt(i,1,10) &
                         + Dzt(k,11) * vt(i,1,11)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wt(i,1,1) &
                         + Dzt(k,2) * wt(i,1,2) &
                         + Dzt(k,3) * wt(i,1,3) &
                         + Dzt(k,4) * wt(i,1,4) &
                         + Dzt(k,5) * wt(i,1,5) &
                         + Dzt(k,6) * wt(i,1,6) &
                         + Dzt(k,7) * wt(i,1,7) &
                         + Dzt(k,8) * wt(i,1,8) &
                         + Dzt(k,9) * wt(i,1,9) &
                         + Dzt(k,10) * wt(i,1,10) &
                         + Dzt(k,11) * wt(i,1,11)
          end do
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx11

  subroutine ax_helm_vector_lx10(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 10
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
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

    !$omp do
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
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )

          vr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wvr(i,1,1) &
                      + G12(i,1,1,e) * wvs(i,1,1) &
                      + G13(i,1,1,e) * wvt(i,1,1) )
          vs(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wvr(i,1,1) &
                      + G22(i,1,1,e) * wvs(i,1,1) &
                      + G23(i,1,1,e) * wvt(i,1,1) )
          vt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wvr(i,1,1) &
                      + G23(i,1,1,e) * wvs(i,1,1) &
                      + G33(i,1,1,e) * wvt(i,1,1) )

          wr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wwr(i,1,1) &
                      + G12(i,1,1,e) * wws(i,1,1) &
                      + G13(i,1,1,e) * wwt(i,1,1) )
          ws(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wwr(i,1,1) &
                      + G22(i,1,1,e) * wws(i,1,1) &
                      + G23(i,1,1,e) * wwt(i,1,1) )
          wt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wwr(i,1,1) &
                      + G23(i,1,1,e) * wws(i,1,1) &
                      + G33(i,1,1,e) * wwt(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                         + Dxt(i,2) * ur(2,j,1) &
                         + Dxt(i,3) * ur(3,j,1) &
                         + Dxt(i,4) * ur(4,j,1) &
                         + Dxt(i,5) * ur(5,j,1) &
                         + Dxt(i,6) * ur(6,j,1) &
                         + Dxt(i,7) * ur(7,j,1) &
                         + Dxt(i,8) * ur(8,j,1) &
                         + Dxt(i,9) * ur(9,j,1) &
                         + Dxt(i,10) * ur(10,j,1)

             av(i,j,1,e) = Dxt(i,1) * vr(1,j,1) &
                         + Dxt(i,2) * vr(2,j,1) &
                         + Dxt(i,3) * vr(3,j,1) &
                         + Dxt(i,4) * vr(4,j,1) &
                         + Dxt(i,5) * vr(5,j,1) &
                         + Dxt(i,6) * vr(6,j,1) &
                         + Dxt(i,7) * vr(7,j,1) &
                         + Dxt(i,8) * vr(8,j,1) &
                         + Dxt(i,9) * vr(9,j,1) &
                         + Dxt(i,10) * vr(10,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wr(1,j,1) &
                         + Dxt(i,2) * wr(2,j,1) &
                         + Dxt(i,3) * wr(3,j,1) &
                         + Dxt(i,4) * wr(4,j,1) &
                         + Dxt(i,5) * wr(5,j,1) &
                         + Dxt(i,6) * wr(6,j,1) &
                         + Dxt(i,7) * wr(7,j,1) &
                         + Dxt(i,8) * wr(8,j,1) &
                         + Dxt(i,9) * wr(9,j,1) &
                         + Dxt(i,10) * wr(10,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * us(i,1,k) &
                            + Dyt(j,2) * us(i,2,k) &
                            + Dyt(j,3) * us(i,3,k) &
                            + Dyt(j,4) * us(i,4,k) &
                            + Dyt(j,5) * us(i,5,k) &
                            + Dyt(j,6) * us(i,6,k) &
                            + Dyt(j,7) * us(i,7,k) &
                            + Dyt(j,8) * us(i,8,k) &
                            + Dyt(j,9) * us(i,9,k) &
                            + Dyt(j,10) * us(i,10,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * vs(i,1,k) &
                            + Dyt(j,2) * vs(i,2,k) &
                            + Dyt(j,3) * vs(i,3,k) &
                            + Dyt(j,4) * vs(i,4,k) &
                            + Dyt(j,5) * vs(i,5,k) &
                            + Dyt(j,6) * vs(i,6,k) &
                            + Dyt(j,7) * vs(i,7,k) &
                            + Dyt(j,8) * vs(i,8,k) &
                            + Dyt(j,9) * vs(i,9,k) &
                            + Dyt(j,10) * vs(i,10,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * ws(i,1,k) &
                            + Dyt(j,2) * ws(i,2,k) &
                            + Dyt(j,3) * ws(i,3,k) &
                            + Dyt(j,4) * ws(i,4,k) &
                            + Dyt(j,5) * ws(i,5,k) &
                            + Dyt(j,6) * ws(i,6,k) &
                            + Dyt(j,7) * ws(i,7,k) &
                            + Dyt(j,8) * ws(i,8,k) &
                            + Dyt(j,9) * ws(i,9,k) &
                            + Dyt(j,10) * ws(i,10,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * ut(i,1,1) &
                         + Dzt(k,2) * ut(i,1,2) &
                         + Dzt(k,3) * ut(i,1,3) &
                         + Dzt(k,4) * ut(i,1,4) &
                         + Dzt(k,5) * ut(i,1,5) &
                         + Dzt(k,6) * ut(i,1,6) &
                         + Dzt(k,7) * ut(i,1,7) &
                         + Dzt(k,8) * ut(i,1,8) &
                         + Dzt(k,9) * ut(i,1,9) &
                         + Dzt(k,10) * ut(i,1,10)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * vt(i,1,1) &
                         + Dzt(k,2) * vt(i,1,2) &
                         + Dzt(k,3) * vt(i,1,3) &
                         + Dzt(k,4) * vt(i,1,4) &
                         + Dzt(k,5) * vt(i,1,5) &
                         + Dzt(k,6) * vt(i,1,6) &
                         + Dzt(k,7) * vt(i,1,7) &
                         + Dzt(k,8) * vt(i,1,8) &
                         + Dzt(k,9) * vt(i,1,9) &
                         + Dzt(k,10) * vt(i,1,10)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wt(i,1,1) &
                         + Dzt(k,2) * wt(i,1,2) &
                         + Dzt(k,3) * wt(i,1,3) &
                         + Dzt(k,4) * wt(i,1,4) &
                         + Dzt(k,5) * wt(i,1,5) &
                         + Dzt(k,6) * wt(i,1,6) &
                         + Dzt(k,7) * wt(i,1,7) &
                         + Dzt(k,8) * wt(i,1,8) &
                         + Dzt(k,9) * wt(i,1,9) &
                         + Dzt(k,10) * wt(i,1,10)
          end do
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx10

  subroutine ax_helm_vector_lx9(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 9
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
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

    !$omp do
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
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )

          vr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wvr(i,1,1) &
                      + G12(i,1,1,e) * wvs(i,1,1) &
                      + G13(i,1,1,e) * wvt(i,1,1) )
          vs(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wvr(i,1,1) &
                      + G22(i,1,1,e) * wvs(i,1,1) &
                      + G23(i,1,1,e) * wvt(i,1,1) )
          vt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wvr(i,1,1) &
                      + G23(i,1,1,e) * wvs(i,1,1) &
                      + G33(i,1,1,e) * wvt(i,1,1) )

          wr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wwr(i,1,1) &
                      + G12(i,1,1,e) * wws(i,1,1) &
                      + G13(i,1,1,e) * wwt(i,1,1) )
          ws(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wwr(i,1,1) &
                      + G22(i,1,1,e) * wws(i,1,1) &
                      + G23(i,1,1,e) * wwt(i,1,1) )
          wt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wwr(i,1,1) &
                      + G23(i,1,1,e) * wws(i,1,1) &
                      + G33(i,1,1,e) * wwt(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                         + Dxt(i,2) * ur(2,j,1) &
                         + Dxt(i,3) * ur(3,j,1) &
                         + Dxt(i,4) * ur(4,j,1) &
                         + Dxt(i,5) * ur(5,j,1) &
                         + Dxt(i,6) * ur(6,j,1) &
                         + Dxt(i,7) * ur(7,j,1) &
                         + Dxt(i,8) * ur(8,j,1) &
                         + Dxt(i,9) * ur(9,j,1)

             av(i,j,1,e) = Dxt(i,1) * vr(1,j,1) &
                         + Dxt(i,2) * vr(2,j,1) &
                         + Dxt(i,3) * vr(3,j,1) &
                         + Dxt(i,4) * vr(4,j,1) &
                         + Dxt(i,5) * vr(5,j,1) &
                         + Dxt(i,6) * vr(6,j,1) &
                         + Dxt(i,7) * vr(7,j,1) &
                         + Dxt(i,8) * vr(8,j,1) &
                         + Dxt(i,9) * vr(9,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wr(1,j,1) &
                         + Dxt(i,2) * wr(2,j,1) &
                         + Dxt(i,3) * wr(3,j,1) &
                         + Dxt(i,4) * wr(4,j,1) &
                         + Dxt(i,5) * wr(5,j,1) &
                         + Dxt(i,6) * wr(6,j,1) &
                         + Dxt(i,7) * wr(7,j,1) &
                         + Dxt(i,8) * wr(8,j,1) &
                         + Dxt(i,9) * wr(9,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * us(i,1,k) &
                            + Dyt(j,2) * us(i,2,k) &
                            + Dyt(j,3) * us(i,3,k) &
                            + Dyt(j,4) * us(i,4,k) &
                            + Dyt(j,5) * us(i,5,k) &
                            + Dyt(j,6) * us(i,6,k) &
                            + Dyt(j,7) * us(i,7,k) &
                            + Dyt(j,8) * us(i,8,k) &
                            + Dyt(j,9) * us(i,9,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * vs(i,1,k) &
                            + Dyt(j,2) * vs(i,2,k) &
                            + Dyt(j,3) * vs(i,3,k) &
                            + Dyt(j,4) * vs(i,4,k) &
                            + Dyt(j,5) * vs(i,5,k) &
                            + Dyt(j,6) * vs(i,6,k) &
                            + Dyt(j,7) * vs(i,7,k) &
                            + Dyt(j,8) * vs(i,8,k) &
                            + Dyt(j,9) * vs(i,9,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * ws(i,1,k) &
                            + Dyt(j,2) * ws(i,2,k) &
                            + Dyt(j,3) * ws(i,3,k) &
                            + Dyt(j,4) * ws(i,4,k) &
                            + Dyt(j,5) * ws(i,5,k) &
                            + Dyt(j,6) * ws(i,6,k) &
                            + Dyt(j,7) * ws(i,7,k) &
                            + Dyt(j,8) * ws(i,8,k) &
                            + Dyt(j,9) * ws(i,9,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * ut(i,1,1) &
                         + Dzt(k,2) * ut(i,1,2) &
                         + Dzt(k,3) * ut(i,1,3) &
                         + Dzt(k,4) * ut(i,1,4) &
                         + Dzt(k,5) * ut(i,1,5) &
                         + Dzt(k,6) * ut(i,1,6) &
                         + Dzt(k,7) * ut(i,1,7) &
                         + Dzt(k,8) * ut(i,1,8) &
                         + Dzt(k,9) * ut(i,1,9)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * vt(i,1,1) &
                         + Dzt(k,2) * vt(i,1,2) &
                         + Dzt(k,3) * vt(i,1,3) &
                         + Dzt(k,4) * vt(i,1,4) &
                         + Dzt(k,5) * vt(i,1,5) &
                         + Dzt(k,6) * vt(i,1,6) &
                         + Dzt(k,7) * vt(i,1,7) &
                         + Dzt(k,8) * vt(i,1,8) &
                         + Dzt(k,9) * vt(i,1,9)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wt(i,1,1) &
                         + Dzt(k,2) * wt(i,1,2) &
                         + Dzt(k,3) * wt(i,1,3) &
                         + Dzt(k,4) * wt(i,1,4) &
                         + Dzt(k,5) * wt(i,1,5) &
                         + Dzt(k,6) * wt(i,1,6) &
                         + Dzt(k,7) * wt(i,1,7) &
                         + Dzt(k,8) * wt(i,1,8) &
                         + Dzt(k,9) * wt(i,1,9)
          end do
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx9

  subroutine ax_helm_vector_lx8(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 8
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
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

    !$omp do
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
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )

          vr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wvr(i,1,1) &
                      + G12(i,1,1,e) * wvs(i,1,1) &
                      + G13(i,1,1,e) * wvt(i,1,1) )
          vs(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wvr(i,1,1) &
                      + G22(i,1,1,e) * wvs(i,1,1) &
                      + G23(i,1,1,e) * wvt(i,1,1) )
          vt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wvr(i,1,1) &
                      + G23(i,1,1,e) * wvs(i,1,1) &
                      + G33(i,1,1,e) * wvt(i,1,1) )

          wr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wwr(i,1,1) &
                      + G12(i,1,1,e) * wws(i,1,1) &
                      + G13(i,1,1,e) * wwt(i,1,1) )
          ws(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wwr(i,1,1) &
                      + G22(i,1,1,e) * wws(i,1,1) &
                      + G23(i,1,1,e) * wwt(i,1,1) )
          wt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wwr(i,1,1) &
                      + G23(i,1,1,e) * wws(i,1,1) &
                      + G33(i,1,1,e) * wwt(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                         + Dxt(i,2) * ur(2,j,1) &
                         + Dxt(i,3) * ur(3,j,1) &
                         + Dxt(i,4) * ur(4,j,1) &
                         + Dxt(i,5) * ur(5,j,1) &
                         + Dxt(i,6) * ur(6,j,1) &
                         + Dxt(i,7) * ur(7,j,1) &
                         + Dxt(i,8) * ur(8,j,1)

             av(i,j,1,e) = Dxt(i,1) * vr(1,j,1) &
                         + Dxt(i,2) * vr(2,j,1) &
                         + Dxt(i,3) * vr(3,j,1) &
                         + Dxt(i,4) * vr(4,j,1) &
                         + Dxt(i,5) * vr(5,j,1) &
                         + Dxt(i,6) * vr(6,j,1) &
                         + Dxt(i,7) * vr(7,j,1) &
                         + Dxt(i,8) * vr(8,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wr(1,j,1) &
                         + Dxt(i,2) * wr(2,j,1) &
                         + Dxt(i,3) * wr(3,j,1) &
                         + Dxt(i,4) * wr(4,j,1) &
                         + Dxt(i,5) * wr(5,j,1) &
                         + Dxt(i,6) * wr(6,j,1) &
                         + Dxt(i,7) * wr(7,j,1) &
                         + Dxt(i,8) * wr(8,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * us(i,1,k) &
                            + Dyt(j,2) * us(i,2,k) &
                            + Dyt(j,3) * us(i,3,k) &
                            + Dyt(j,4) * us(i,4,k) &
                            + Dyt(j,5) * us(i,5,k) &
                            + Dyt(j,6) * us(i,6,k) &
                            + Dyt(j,7) * us(i,7,k) &
                            + Dyt(j,8) * us(i,8,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * vs(i,1,k) &
                            + Dyt(j,2) * vs(i,2,k) &
                            + Dyt(j,3) * vs(i,3,k) &
                            + Dyt(j,4) * vs(i,4,k) &
                            + Dyt(j,5) * vs(i,5,k) &
                            + Dyt(j,6) * vs(i,6,k) &
                            + Dyt(j,7) * vs(i,7,k) &
                            + Dyt(j,8) * vs(i,8,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * ws(i,1,k) &
                            + Dyt(j,2) * ws(i,2,k) &
                            + Dyt(j,3) * ws(i,3,k) &
                            + Dyt(j,4) * ws(i,4,k) &
                            + Dyt(j,5) * ws(i,5,k) &
                            + Dyt(j,6) * ws(i,6,k) &
                            + Dyt(j,7) * ws(i,7,k) &
                            + Dyt(j,8) * ws(i,8,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * ut(i,1,1) &
                         + Dzt(k,2) * ut(i,1,2) &
                         + Dzt(k,3) * ut(i,1,3) &
                         + Dzt(k,4) * ut(i,1,4) &
                         + Dzt(k,5) * ut(i,1,5) &
                         + Dzt(k,6) * ut(i,1,6) &
                         + Dzt(k,7) * ut(i,1,7) &
                         + Dzt(k,8) * ut(i,1,8)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * vt(i,1,1) &
                         + Dzt(k,2) * vt(i,1,2) &
                         + Dzt(k,3) * vt(i,1,3) &
                         + Dzt(k,4) * vt(i,1,4) &
                         + Dzt(k,5) * vt(i,1,5) &
                         + Dzt(k,6) * vt(i,1,6) &
                         + Dzt(k,7) * vt(i,1,7) &
                         + Dzt(k,8) * vt(i,1,8)

             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * wt(i,1,1) &
                         + Dzt(k,2) * wt(i,1,2) &
                         + Dzt(k,3) * wt(i,1,3) &
                         + Dzt(k,4) * wt(i,1,4) &
                         + Dzt(k,5) * wt(i,1,5) &
                         + Dzt(k,6) * wt(i,1,6) &
                         + Dzt(k,7) * wt(i,1,7) &
                         + Dzt(k,8) * wt(i,1,8)
          end do
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx8

  subroutine ax_helm_vector_lx7(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 7
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
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

    !$omp do
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
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )

          vr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wvr(i,1,1) &
                      + G12(i,1,1,e) * wvs(i,1,1) &
                      + G13(i,1,1,e) * wvt(i,1,1) )
          vs(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wvr(i,1,1) &
                      + G22(i,1,1,e) * wvs(i,1,1) &
                      + G23(i,1,1,e) * wvt(i,1,1) )
          vt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wvr(i,1,1) &
                      + G23(i,1,1,e) * wvs(i,1,1) &
                      + G33(i,1,1,e) * wvt(i,1,1) )

          wr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wwr(i,1,1) &
                      + G12(i,1,1,e) * wws(i,1,1) &
                      + G13(i,1,1,e) * wwt(i,1,1) )
          ws(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wwr(i,1,1) &
                      + G22(i,1,1,e) * wws(i,1,1) &
                      + G23(i,1,1,e) * wwt(i,1,1) )
          wt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wwr(i,1,1) &
                      + G23(i,1,1,e) * wws(i,1,1) &
                      + G33(i,1,1,e) * wwt(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                         + Dxt(i,2) * ur(2,j,1) &
                         + Dxt(i,3) * ur(3,j,1) &
                         + Dxt(i,4) * ur(4,j,1) &
                         + Dxt(i,5) * ur(5,j,1) &
                         + Dxt(i,6) * ur(6,j,1) &
                         + Dxt(i,7) * ur(7,j,1)

             av(i,j,1,e) = Dxt(i,1) * vr(1,j,1) &
                         + Dxt(i,2) * vr(2,j,1) &
                         + Dxt(i,3) * vr(3,j,1) &
                         + Dxt(i,4) * vr(4,j,1) &
                         + Dxt(i,5) * vr(5,j,1) &
                         + Dxt(i,6) * vr(6,j,1) &
                         + Dxt(i,7) * vr(7,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wr(1,j,1) &
                         + Dxt(i,2) * wr(2,j,1) &
                         + Dxt(i,3) * wr(3,j,1) &
                         + Dxt(i,4) * wr(4,j,1) &
                         + Dxt(i,5) * wr(5,j,1) &
                         + Dxt(i,6) * wr(6,j,1) &
                         + Dxt(i,7) * wr(7,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * us(i,1,k) &
                            + Dyt(j,2) * us(i,2,k) &
                            + Dyt(j,3) * us(i,3,k) &
                            + Dyt(j,4) * us(i,4,k) &
                            + Dyt(j,5) * us(i,5,k) &
                            + Dyt(j,6) * us(i,6,k) &
                            + Dyt(j,7) * us(i,7,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * vs(i,1,k) &
                            + Dyt(j,2) * vs(i,2,k) &
                            + Dyt(j,3) * vs(i,3,k) &
                            + Dyt(j,4) * vs(i,4,k) &
                            + Dyt(j,5) * vs(i,5,k) &
                            + Dyt(j,6) * vs(i,6,k) &
                            + Dyt(j,7) * vs(i,7,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * ws(i,1,k) &
                            + Dyt(j,2) * ws(i,2,k) &
                            + Dyt(j,3) * ws(i,3,k) &
                            + Dyt(j,4) * ws(i,4,k) &
                            + Dyt(j,5) * ws(i,5,k) &
                            + Dyt(j,6) * ws(i,6,k) &
                            + Dyt(j,7) * ws(i,7,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * ut(i,1,1) &
                         + Dzt(k,2) * ut(i,1,2) &
                         + Dzt(k,3) * ut(i,1,3) &
                         + Dzt(k,4) * ut(i,1,4) &
                         + Dzt(k,5) * ut(i,1,5) &
                         + Dzt(k,6) * ut(i,1,6) &
                         + Dzt(k,7) * ut(i,1,7)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * vt(i,1,1) &
                         + Dzt(k,2) * vt(i,1,2) &
                         + Dzt(k,3) * vt(i,1,3) &
                         + Dzt(k,4) * vt(i,1,4) &
                         + Dzt(k,5) * vt(i,1,5) &
                         + Dzt(k,6) * vt(i,1,6) &
                         + Dzt(k,7) * vt(i,1,7)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wt(i,1,1) &
                         + Dzt(k,2) * wt(i,1,2) &
                         + Dzt(k,3) * wt(i,1,3) &
                         + Dzt(k,4) * wt(i,1,4) &
                         + Dzt(k,5) * wt(i,1,5) &
                         + Dzt(k,6) * wt(i,1,6) &
                         + Dzt(k,7) * wt(i,1,7)
          end do
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx7

  subroutine ax_helm_vector_lx6(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 6
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
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

    !$omp do
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
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )

          vr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wvr(i,1,1) &
                      + G12(i,1,1,e) * wvs(i,1,1) &
                      + G13(i,1,1,e) * wvt(i,1,1) )
          vs(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wvr(i,1,1) &
                      + G22(i,1,1,e) * wvs(i,1,1) &
                      + G23(i,1,1,e) * wvt(i,1,1) )
          vt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wvr(i,1,1) &
                      + G23(i,1,1,e) * wvs(i,1,1) &
                      + G33(i,1,1,e) * wvt(i,1,1) )

          wr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wwr(i,1,1) &
                      + G12(i,1,1,e) * wws(i,1,1) &
                      + G13(i,1,1,e) * wwt(i,1,1) )
          ws(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wwr(i,1,1) &
                      + G22(i,1,1,e) * wws(i,1,1) &
                      + G23(i,1,1,e) * wwt(i,1,1) )
          wt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wwr(i,1,1) &
                      + G23(i,1,1,e) * wws(i,1,1) &
                      + G33(i,1,1,e) * wwt(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                         + Dxt(i,2) * ur(2,j,1) &
                         + Dxt(i,3) * ur(3,j,1) &
                         + Dxt(i,4) * ur(4,j,1) &
                         + Dxt(i,5) * ur(5,j,1) &
                         + Dxt(i,6) * ur(6,j,1)

             av(i,j,1,e) = Dxt(i,1) * vr(1,j,1) &
                         + Dxt(i,2) * vr(2,j,1) &
                         + Dxt(i,3) * vr(3,j,1) &
                         + Dxt(i,4) * vr(4,j,1) &
                         + Dxt(i,5) * vr(5,j,1) &
                         + Dxt(i,6) * vr(6,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wr(1,j,1) &
                         + Dxt(i,2) * wr(2,j,1) &
                         + Dxt(i,3) * wr(3,j,1) &
                         + Dxt(i,4) * wr(4,j,1) &
                         + Dxt(i,5) * wr(5,j,1) &
                         + Dxt(i,6) * wr(6,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * us(i,1,k) &
                            + Dyt(j,2) * us(i,2,k) &
                            + Dyt(j,3) * us(i,3,k) &
                            + Dyt(j,4) * us(i,4,k) &
                            + Dyt(j,5) * us(i,5,k) &
                            + Dyt(j,6) * us(i,6,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * vs(i,1,k) &
                            + Dyt(j,2) * vs(i,2,k) &
                            + Dyt(j,3) * vs(i,3,k) &
                            + Dyt(j,4) * vs(i,4,k) &
                            + Dyt(j,5) * vs(i,5,k) &
                            + Dyt(j,6) * vs(i,6,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * ws(i,1,k) &
                            + Dyt(j,2) * ws(i,2,k) &
                            + Dyt(j,3) * ws(i,3,k) &
                            + Dyt(j,4) * ws(i,4,k) &
                            + Dyt(j,5) * ws(i,5,k) &
                            + Dyt(j,6) * ws(i,6,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * ut(i,1,1) &
                         + Dzt(k,2) * ut(i,1,2) &
                         + Dzt(k,3) * ut(i,1,3) &
                         + Dzt(k,4) * ut(i,1,4) &
                         + Dzt(k,5) * ut(i,1,5) &
                         + Dzt(k,6) * ut(i,1,6)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * vt(i,1,1) &
                         + Dzt(k,2) * vt(i,1,2) &
                         + Dzt(k,3) * vt(i,1,3) &
                         + Dzt(k,4) * vt(i,1,4) &
                         + Dzt(k,5) * vt(i,1,5) &
                         + Dzt(k,6) * vt(i,1,6)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wt(i,1,1) &
                         + Dzt(k,2) * wt(i,1,2) &
                         + Dzt(k,3) * wt(i,1,3) &
                         + Dzt(k,4) * wt(i,1,4) &
                         + Dzt(k,5) * wt(i,1,5) &
                         + Dzt(k,6) * wt(i,1,6)
          end do
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx6

  subroutine ax_helm_vector_lx5(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 5
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
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

    !$omp do
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
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )

          vr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wvr(i,1,1) &
                      + G12(i,1,1,e) * wvs(i,1,1) &
                      + G13(i,1,1,e) * wvt(i,1,1) )
          vs(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wvr(i,1,1) &
                      + G22(i,1,1,e) * wvs(i,1,1) &
                      + G23(i,1,1,e) * wvt(i,1,1) )
          vt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wvr(i,1,1) &
                      + G23(i,1,1,e) * wvs(i,1,1) &
                      + G33(i,1,1,e) * wvt(i,1,1) )

          wr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wwr(i,1,1) &
                      + G12(i,1,1,e) * wws(i,1,1) &
                      + G13(i,1,1,e) * wwt(i,1,1) )
          ws(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wwr(i,1,1) &
                      + G22(i,1,1,e) * wws(i,1,1) &
                      + G23(i,1,1,e) * wwt(i,1,1) )
          wt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wwr(i,1,1) &
                      + G23(i,1,1,e) * wws(i,1,1) &
                      + G33(i,1,1,e) * wwt(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                         + Dxt(i,2) * ur(2,j,1) &
                         + Dxt(i,3) * ur(3,j,1) &
                         + Dxt(i,4) * ur(4,j,1) &
                         + Dxt(i,5) * ur(5,j,1)

             av(i,j,1,e) = Dxt(i,1) * vr(1,j,1) &
                         + Dxt(i,2) * vr(2,j,1) &
                         + Dxt(i,3) * vr(3,j,1) &
                         + Dxt(i,4) * vr(4,j,1) &
                         + Dxt(i,5) * vr(5,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wr(1,j,1) &
                         + Dxt(i,2) * wr(2,j,1) &
                         + Dxt(i,3) * wr(3,j,1) &
                         + Dxt(i,4) * wr(4,j,1) &
                         + Dxt(i,5) * wr(5,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * us(i,1,k) &
                            + Dyt(j,2) * us(i,2,k) &
                            + Dyt(j,3) * us(i,3,k) &
                            + Dyt(j,4) * us(i,4,k) &
                            + Dyt(j,5) * us(i,5,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * vs(i,1,k) &
                            + Dyt(j,2) * vs(i,2,k) &
                            + Dyt(j,3) * vs(i,3,k) &
                            + Dyt(j,4) * vs(i,4,k) &
                            + Dyt(j,5) * vs(i,5,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * ws(i,1,k) &
                            + Dyt(j,2) * ws(i,2,k) &
                            + Dyt(j,3) * ws(i,3,k) &
                            + Dyt(j,4) * ws(i,4,k) &
                            + Dyt(j,5) * ws(i,5,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * ut(i,1,1) &
                         + Dzt(k,2) * ut(i,1,2) &
                         + Dzt(k,3) * ut(i,1,3) &
                         + Dzt(k,4) * ut(i,1,4) &
                         + Dzt(k,5) * ut(i,1,5)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * vt(i,1,1) &
                         + Dzt(k,2) * vt(i,1,2) &
                         + Dzt(k,3) * vt(i,1,3) &
                         + Dzt(k,4) * vt(i,1,4) &
                         + Dzt(k,5) * vt(i,1,5)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wt(i,1,1) &
                         + Dzt(k,2) * wt(i,1,2) &
                         + Dzt(k,3) * wt(i,1,3) &
                         + Dzt(k,4) * wt(i,1,4) &
                         + Dzt(k,5) * wt(i,1,5)
          end do
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx5

  subroutine ax_helm_vector_lx4(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 4
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
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

    !$omp do
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
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )

          vr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wvr(i,1,1) &
                      + G12(i,1,1,e) * wvs(i,1,1) &
                      + G13(i,1,1,e) * wvt(i,1,1) )
          vs(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wvr(i,1,1) &
                      + G22(i,1,1,e) * wvs(i,1,1) &
                      + G23(i,1,1,e) * wvt(i,1,1) )
          vt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wvr(i,1,1) &
                      + G23(i,1,1,e) * wvs(i,1,1) &
                      + G33(i,1,1,e) * wvt(i,1,1) )

          wr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wwr(i,1,1) &
                      + G12(i,1,1,e) * wws(i,1,1) &
                      + G13(i,1,1,e) * wwt(i,1,1) )
          ws(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wwr(i,1,1) &
                      + G22(i,1,1,e) * wws(i,1,1) &
                      + G23(i,1,1,e) * wwt(i,1,1) )
          wt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wwr(i,1,1) &
                      + G23(i,1,1,e) * wws(i,1,1) &
                      + G33(i,1,1,e) * wwt(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                         + Dxt(i,2) * ur(2,j,1) &
                         + Dxt(i,3) * ur(3,j,1) &
                         + Dxt(i,4) * ur(4,j,1)

             av(i,j,1,e) = Dxt(i,1) * vr(1,j,1) &
                         + Dxt(i,2) * vr(2,j,1) &
                         + Dxt(i,3) * vr(3,j,1) &
                         + Dxt(i,4) * vr(4,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wr(1,j,1) &
                         + Dxt(i,2) * wr(2,j,1) &
                         + Dxt(i,3) * wr(3,j,1) &
                         + Dxt(i,4) * wr(4,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * us(i,1,k) &
                            + Dyt(j,2) * us(i,2,k) &
                            + Dyt(j,3) * us(i,3,k) &
                            + Dyt(j,4) * us(i,4,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * vs(i,1,k) &
                            + Dyt(j,2) * vs(i,2,k) &
                            + Dyt(j,3) * vs(i,3,k) &
                            + Dyt(j,4) * vs(i,4,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * ws(i,1,k) &
                            + Dyt(j,2) * ws(i,2,k) &
                            + Dyt(j,3) * ws(i,3,k) &
                            + Dyt(j,4) * ws(i,4,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * ut(i,1,1) &
                         + Dzt(k,2) * ut(i,1,2) &
                         + Dzt(k,3) * ut(i,1,3) &
                         + Dzt(k,4) * ut(i,1,4)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * vt(i,1,1) &
                         + Dzt(k,2) * vt(i,1,2) &
                         + Dzt(k,3) * vt(i,1,3) &
                         + Dzt(k,4) * vt(i,1,4)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wt(i,1,1) &
                         + Dzt(k,2) * wt(i,1,2) &
                         + Dzt(k,3) * wt(i,1,3) &
                         + Dzt(k,4) * wt(i,1,4)
          end do
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx4

  subroutine ax_helm_vector_lx3(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 3
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
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

    !$omp do
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
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )

          vr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wvr(i,1,1) &
                      + G12(i,1,1,e) * wvs(i,1,1) &
                      + G13(i,1,1,e) * wvt(i,1,1) )
          vs(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wvr(i,1,1) &
                      + G22(i,1,1,e) * wvs(i,1,1) &
                      + G23(i,1,1,e) * wvt(i,1,1) )
          vt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wwr(i,1,1) &
                      + G23(i,1,1,e) * wws(i,1,1) &
                      + G33(i,1,1,e) * wwt(i,1,1) )

          wr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wwr(i,1,1) &
                      + G12(i,1,1,e) * wws(i,1,1) &
                      + G13(i,1,1,e) * wwt(i,1,1) )
          ws(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wwr(i,1,1) &
                      + G22(i,1,1,e) * wws(i,1,1) &
                      + G23(i,1,1,e) * wwt(i,1,1) )
          wt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wwr(i,1,1) &
                      + G23(i,1,1,e) * wws(i,1,1) &
                      + G33(i,1,1,e) * wwt(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                         + Dxt(i,2) * ur(2,j,1) &
                         + Dxt(i,3) * ur(3,j,1)

             av(i,j,1,e) = Dxt(i,1) * vr(1,j,1) &
                         + Dxt(i,2) * vr(2,j,1) &
                         + Dxt(i,3) * vr(3,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wr(1,j,1) &
                         + Dxt(i,2) * wr(2,j,1) &
                         + Dxt(i,3) * wr(3,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * us(i,1,k) &
                            + Dyt(j,2) * us(i,2,k) &
                            + Dyt(j,3) * us(i,3,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * vs(i,1,k) &
                            + Dyt(j,2) * vs(i,2,k) &
                            + Dyt(j,3) * vs(i,3,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * ws(i,1,k) &
                            + Dyt(j,2) * ws(i,2,k) &
                            + Dyt(j,3) * ws(i,3,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * ut(i,1,1) &
                         + Dzt(k,2) * ut(i,1,2) &
                         + Dzt(k,3) * ut(i,1,3)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * vt(i,1,1) &
                         + Dzt(k,2) * vt(i,1,2) &
                         + Dzt(k,3) * vt(i,1,3)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wt(i,1,1) &
                         + Dzt(k,2) * wt(i,1,2) &
                         + Dzt(k,3) * wt(i,1,3)
          end do
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx3

  subroutine ax_helm_vector_lx2(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 2
    integer, intent(in) :: n
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=rp) :: ur(lx, lx, lx)
    real(kind=rp) :: us(lx, lx, lx)
    real(kind=rp) :: ut(lx, lx, lx)
    real(kind=rp) :: vr(lx, lx, lx)
    real(kind=rp) :: vs(lx, lx, lx)
    real(kind=rp) :: vt(lx, lx, lx)
    real(kind=rp) :: wr(lx, lx, lx)
    real(kind=rp) :: ws(lx, lx, lx)
    real(kind=rp) :: wt(lx, lx, lx)
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

    !$omp do
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
          ur(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wur(i,1,1) &
                      + G12(i,1,1,e) * wus(i,1,1) &
                      + G13(i,1,1,e) * wut(i,1,1) )
          us(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wur(i,1,1) &
                      + G22(i,1,1,e) * wus(i,1,1) &
                      + G23(i,1,1,e) * wut(i,1,1) )
          ut(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wur(i,1,1) &
                      + G23(i,1,1,e) * wus(i,1,1) &
                      + G33(i,1,1,e) * wut(i,1,1) )

          vr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wvr(i,1,1) &
                      + G12(i,1,1,e) * wvs(i,1,1) &
                      + G13(i,1,1,e) * wvt(i,1,1) )
          vs(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wvr(i,1,1) &
                      + G22(i,1,1,e) * wvs(i,1,1) &
                      + G23(i,1,1,e) * wvt(i,1,1) )
          vt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wvr(i,1,1) &
                      + G23(i,1,1,e) * wvs(i,1,1) &
                      + G33(i,1,1,e) * wvt(i,1,1) )

          wr(i,1,1) = h1(i,1,1,e) &
                    * ( G11(i,1,1,e) * wwr(i,1,1) &
                      + G12(i,1,1,e) * wws(i,1,1) &
                      + G13(i,1,1,e) * wwt(i,1,1) )
          ws(i,1,1) = h1(i,1,1,e) &
                    * ( G12(i,1,1,e) * wwr(i,1,1) &
                      + G22(i,1,1,e) * wws(i,1,1) &
                      + G23(i,1,1,e) * wwt(i,1,1) )
          wt(i,1,1) = h1(i,1,1,e) &
                    * ( G13(i,1,1,e) * wwr(i,1,1) &
                      + G23(i,1,1,e) * wws(i,1,1) &
                      + G33(i,1,1,e) * wwt(i,1,1) )
       end do

       do j = 1, lx*lx
          do i = 1, lx
             au(i,j,1,e) = Dxt(i,1) * ur(1,j,1) &
                         + Dxt(i,2) * ur(2,j,1)

             av(i,j,1,e) = Dxt(i,1) * vr(1,j,1) &
                         + Dxt(i,2) * vr(2,j,1)

             aw(i,j,1,e) = Dxt(i,1) * wr(1,j,1) &
                         + Dxt(i,2) * wr(2,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                au(i,j,k,e) = au(i,j,k,e) &
                            + Dyt(j,1) * us(i,1,k) &
                            + Dyt(j,2) * us(i,2,k)

                av(i,j,k,e) = av(i,j,k,e) &
                            + Dyt(j,1) * vs(i,1,k) &
                            + Dyt(j,2) * vs(i,2,k)

                aw(i,j,k,e) = aw(i,j,k,e) &
                            + Dyt(j,1) * ws(i,1,k) &
                            + Dyt(j,2) * ws(i,2,k)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             au(i,1,k,e) = au(i,1,k,e) &
                         + Dzt(k,1) * ut(i,1,1) &
                         + Dzt(k,2) * ut(i,1,2)

             av(i,1,k,e) = av(i,1,k,e) &
                         + Dzt(k,1) * vt(i,1,1) &
                         + Dzt(k,2) * vt(i,1,2)

             aw(i,1,k,e) = aw(i,1,k,e) &
                         + Dzt(k,1) * wt(i,1,1) &
                         + Dzt(k,2) * wt(i,1,2)
          end do
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx2

end submodule ax_helm_vector_cpu
