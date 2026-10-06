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
  !! @note The mass term \f$ h_2 B u \f$ is applied inside the element
  !! kernels (when @a coef%ifh2 is set), while the element is still in cache.
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

    !$omp parallel
    select case(Xh%lx)
    case (14)
       call ax_helm_vector_lx14(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, coef%B, coef%ifh2, &
            coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
            msh%nelv)
    case (13)
       call ax_helm_vector_lx13(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, coef%B, coef%ifh2, &
            coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
            msh%nelv)
    case (12)
       call ax_helm_vector_lx12(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, coef%B, coef%ifh2, &
            coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
            msh%nelv)
    case (11)
       call ax_helm_vector_lx11(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, coef%B, coef%ifh2, &
            coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
            msh%nelv)
    case (10)
       call ax_helm_vector_lx10(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, coef%B, coef%ifh2, &
            coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
            msh%nelv)
    case (9)
       call ax_helm_vector_lx9(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, coef%B, coef%ifh2, &
            coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
            msh%nelv)
    case (8)
       call ax_helm_vector_lx8(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, coef%B, coef%ifh2, &
            coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
            msh%nelv)
    case (7)
       call ax_helm_vector_lx7(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, coef%B, coef%ifh2, &
            coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
            msh%nelv)
    case (6)
       call ax_helm_vector_lx6(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, coef%B, coef%ifh2, &
            coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
            msh%nelv)
    case (5)
       call ax_helm_vector_lx5(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, coef%B, coef%ifh2, &
            coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
            msh%nelv)
    case (4)
       call ax_helm_vector_lx4(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, coef%B, coef%ifh2, &
            coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
            msh%nelv)
    case (3)
       call ax_helm_vector_lx3(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, coef%B, coef%ifh2, &
            coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
            msh%nelv)
    case (2)
       call ax_helm_vector_lx2(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, coef%B, coef%ifh2, &
            coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
            msh%nelv)
    case default
       call ax_helm_vector_lx(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, &
            Xh%dxt, Xh%dyt, Xh%dzt, coef%h1, coef%h2, coef%B, coef%ifh2, &
            coef%G11, coef%G22, coef%G33, coef%G12, coef%G13, coef%G23, &
            msh%nelv, Xh%lx)
    end select
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
  !! @param h1 Stiffness scaling.
  !! @param h2 Mass scaling.
  !! @param B Mass matrix.
  !! @param ifh2 True if the mass term \f$ h_2 B u \f$ should be added.
  !! @param G11 Geometric factor.
  !! @param n Number of elements.
  !! @param lx Polynomial order.
  subroutine ax_helm_vector_lx(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, h2, B, ifh2, G11, G22, G33, G12, G13, G23, n, lx)
    integer, intent(in) :: n, lx
    logical, intent(in) :: ifh2
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: B(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=xp) :: ur(lx, lx, lx)
    real(kind=xp) :: us(lx, lx, lx)
    real(kind=xp) :: ut(lx, lx, lx)
    real(kind=xp) :: vr(lx, lx, lx)
    real(kind=xp) :: vs(lx, lx, lx)
    real(kind=xp) :: vt(lx, lx, lx)
    real(kind=xp) :: wr(lx, lx, lx)
    real(kind=xp) :: ws(lx, lx, lx)
    real(kind=xp) :: wt(lx, lx, lx)
    real(kind=xp) :: wur(lx, lx, lx)
    real(kind=xp) :: wvr(lx, lx, lx)
    real(kind=xp) :: wwr(lx, lx, lx)
    real(kind=xp) :: wus(lx, lx, lx)
    real(kind=xp) :: wvs(lx, lx, lx)
    real(kind=xp) :: wws(lx, lx, lx)
    real(kind=xp) :: wut(lx, lx, lx)
    real(kind=xp) :: wvt(lx, lx, lx)
    real(kind=xp) :: wwt(lx, lx, lx)
    real(kind=xp) :: aud(lx, lx, lx), avd(lx, lx, lx), awd(lx, lx, lx)
    real(kind=xp) :: t1, t2, t3
    integer :: e, i, j, k, l

    !$omp do
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             t1 = 0.0_xp
             t2 = 0.0_xp
             t3 = 0.0_xp
             do k = 1, lx
                t1 = t1 + real(Dx(i,k), xp) * u(k,j,1,e)
                t2 = t2 + real(Dx(i,k), xp) * v(k,j,1,e)
                t3 = t3 + real(Dx(i,k), xp) * w(k,j,1,e)
             end do
             wur(i,j,1) = t1
             wvr(i,j,1) = t2
             wwr(i,j,1) = t3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                t1 = 0.0_xp
                t2 = 0.0_xp
                t3 = 0.0_xp
                do l = 1, lx
                   t1 = t1 + real(Dy(j,l), xp) * u(i,l,k,e)
                   t2 = t2 + real(Dy(j,l), xp) * v(i,l,k,e)
                   t3 = t3 + real(Dy(j,l), xp) * w(i,l,k,e)
                end do
                wus(i,j,k) = t1
                wvs(i,j,k) = t2
                wws(i,j,k) = t3
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             t1 = 0.0_xp
             t2 = 0.0_xp
             t3 = 0.0_xp
             do l = 1, lx
                t1 = t1 + real(Dz(k,l), xp) * u(i,1,l,e)
                t2 = t2 + real(Dz(k,l), xp) * v(i,1,l,e)
                t3 = t3 + real(Dz(k,l), xp) * w(i,1,l,e)
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
             t1 = 0.0_xp
             t2 = 0.0_xp
             t3 = 0.0_xp
             do k = 1, lx
                t1 = t1 + Dxt(i,k) * ur(k,j,1)
                t2 = t2 + Dxt(i,k) * vr(k,j,1)
                t3 = t3 + Dxt(i,k) * wr(k,j,1)
             end do
             aud(i,j,1) = t1
             avd(i,j,1) = t2
             awd(i,j,1) = t3
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                t1 = 0.0_xp
                t2 = 0.0_xp
                t3 = 0.0_xp
                do l = 1, lx
                   t1 = t1 + Dyt(j,l) * us(i,l,k)
                   t2 = t2 + Dyt(j,l) * vs(i,l,k)
                   t3 = t3 + Dyt(j,l) * ws(i,l,k)
                end do
                aud(i,j,k) = aud(i,j,k) + t1
                avd(i,j,k) = avd(i,j,k) + t2
                awd(i,j,k) = awd(i,j,k) + t3
             end do
          end do
       end do

       if (ifh2) then
          do k = 1, lx
             do i = 1, lx*lx
                t1 = 0.0_xp
                t2 = 0.0_xp
                t3 = 0.0_xp
                do l = 1, lx
                   t1 = t1 + Dzt(k,l) * ut(i,1,l)
                   t2 = t2 + Dzt(k,l) * vt(i,1,l)
                   t3 = t3 + Dzt(k,l) * wt(i,1,l)
                end do
                aud(i,1,k) = aud(i,1,k) + t1 &
                     + h2(i,1,k,e) * B(i,1,k,e) * u(i,1,k,e)
                avd(i,1,k) = avd(i,1,k) + t2 &
                     + h2(i,1,k,e) * B(i,1,k,e) * v(i,1,k,e)
                awd(i,1,k) = awd(i,1,k) + t3 &
                     + h2(i,1,k,e) * B(i,1,k,e) * w(i,1,k,e)
             end do
          end do
       else
          do k = 1, lx
             do i = 1, lx*lx
                t1 = 0.0_xp
                t2 = 0.0_xp
                t3 = 0.0_xp
                do l = 1, lx
                   t1 = t1 + Dzt(k,l) * ut(i,1,l)
                   t2 = t2 + Dzt(k,l) * vt(i,1,l)
                   t3 = t3 + Dzt(k,l) * wt(i,1,l)
                end do
                aud(i,1,k) = aud(i,1,k) + t1
                avd(i,1,k) = avd(i,1,k) + t2
                awd(i,1,k) = awd(i,1,k) + t3
             end do
          end do
       end if


       ! Single truncation of the dp-accumulated operator
       do i = 1, lx*lx*lx
          au(i,1,1,e) = aud(i,1,1)
          av(i,1,1,e) = avd(i,1,1)
          aw(i,1,1,e) = awd(i,1,1)
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx

  subroutine ax_helm_vector_lx14(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, h2, B, ifh2, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 14
    integer, intent(in) :: n
    logical, intent(in) :: ifh2
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: B(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=xp) :: ur(lx, lx, lx)
    real(kind=xp) :: us(lx, lx, lx)
    real(kind=xp) :: ut(lx, lx, lx)
    real(kind=xp) :: vr(lx, lx, lx)
    real(kind=xp) :: vs(lx, lx, lx)
    real(kind=xp) :: vt(lx, lx, lx)
    real(kind=xp) :: wr(lx, lx, lx)
    real(kind=xp) :: ws(lx, lx, lx)
    real(kind=xp) :: wt(lx, lx, lx)
    real(kind=xp) :: wur(lx, lx, lx)
    real(kind=xp) :: wus(lx, lx, lx)
    real(kind=xp) :: wut(lx, lx, lx)
    real(kind=xp) :: wvr(lx, lx, lx)
    real(kind=xp) :: wvs(lx, lx, lx)
    real(kind=xp) :: wvt(lx, lx, lx)
    real(kind=xp) :: wwr(lx, lx, lx)
    real(kind=xp) :: wws(lx, lx, lx)
    real(kind=xp) :: wwt(lx, lx, lx)
    real(kind=xp) :: aud(lx, lx, lx), avd(lx, lx, lx), awd(lx, lx, lx)
    integer :: e, i, j, k

    !$omp do
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = real(Dx(i,1), xp) * u(1,j,1,e) &
                  + real(Dx(i,2), xp) * u(2,j,1,e) &
                  + real(Dx(i,3), xp) * u(3,j,1,e) &
                  + real(Dx(i,4), xp) * u(4,j,1,e) &
                  + real(Dx(i,5), xp) * u(5,j,1,e) &
                  + real(Dx(i,6), xp) * u(6,j,1,e) &
                  + real(Dx(i,7), xp) * u(7,j,1,e) &
                  + real(Dx(i,8), xp) * u(8,j,1,e) &
                  + real(Dx(i,9), xp) * u(9,j,1,e) &
                  + real(Dx(i,10), xp) * u(10,j,1,e) &
                  + real(Dx(i,11), xp) * u(11,j,1,e) &
                  + real(Dx(i,12), xp) * u(12,j,1,e) &
                  + real(Dx(i,13), xp) * u(13,j,1,e) &
                  + real(Dx(i,14), xp) * u(14,j,1,e)

             wvr(i,j,1) = real(Dx(i,1), xp) * v(1,j,1,e) &
                  + real(Dx(i,2), xp) * v(2,j,1,e) &
                  + real(Dx(i,3), xp) * v(3,j,1,e) &
                  + real(Dx(i,4), xp) * v(4,j,1,e) &
                  + real(Dx(i,5), xp) * v(5,j,1,e) &
                  + real(Dx(i,6), xp) * v(6,j,1,e) &
                  + real(Dx(i,7), xp) * v(7,j,1,e) &
                  + real(Dx(i,8), xp) * v(8,j,1,e) &
                  + real(Dx(i,9), xp) * v(9,j,1,e) &
                  + real(Dx(i,10), xp) * v(10,j,1,e) &
                  + real(Dx(i,11), xp) * v(11,j,1,e) &
                  + real(Dx(i,12), xp) * v(12,j,1,e) &
                  + real(Dx(i,13), xp) * v(13,j,1,e) &
                  + real(Dx(i,14), xp) * v(14,j,1,e)

             wwr(i,j,1) = real(Dx(i,1), xp) * w(1,j,1,e) &
                  + real(Dx(i,2), xp) * w(2,j,1,e) &
                  + real(Dx(i,3), xp) * w(3,j,1,e) &
                  + real(Dx(i,4), xp) * w(4,j,1,e) &
                  + real(Dx(i,5), xp) * w(5,j,1,e) &
                  + real(Dx(i,6), xp) * w(6,j,1,e) &
                  + real(Dx(i,7), xp) * w(7,j,1,e) &
                  + real(Dx(i,8), xp) * w(8,j,1,e) &
                  + real(Dx(i,9), xp) * w(9,j,1,e) &
                  + real(Dx(i,10), xp) * w(10,j,1,e) &
                  + real(Dx(i,11), xp) * w(11,j,1,e) &
                  + real(Dx(i,12), xp) * w(12,j,1,e) &
                  + real(Dx(i,13), xp) * w(13,j,1,e) &
                  + real(Dx(i,14), xp) * w(14,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = real(Dy(j,1), xp) * u(i,1,k,e) &
                     + real(Dy(j,2), xp) * u(i,2,k,e) &
                     + real(Dy(j,3), xp) * u(i,3,k,e) &
                     + real(Dy(j,4), xp) * u(i,4,k,e) &
                     + real(Dy(j,5), xp) * u(i,5,k,e) &
                     + real(Dy(j,6), xp) * u(i,6,k,e) &
                     + real(Dy(j,7), xp) * u(i,7,k,e) &
                     + real(Dy(j,8), xp) * u(i,8,k,e) &
                     + real(Dy(j,9), xp) * u(i,9,k,e) &
                     + real(Dy(j,10), xp) * u(i,10,k,e) &
                     + real(Dy(j,11), xp) * u(i,11,k,e) &
                     + real(Dy(j,12), xp) * u(i,12,k,e) &
                     + real(Dy(j,13), xp) * u(i,13,k,e) &
                     + real(Dy(j,14), xp) * u(i,14,k,e)

                wvs(i,j,k) = real(Dy(j,1), xp) * v(i,1,k,e) &
                     + real(Dy(j,2), xp) * v(i,2,k,e) &
                     + real(Dy(j,3), xp) * v(i,3,k,e) &
                     + real(Dy(j,4), xp) * v(i,4,k,e) &
                     + real(Dy(j,5), xp) * v(i,5,k,e) &
                     + real(Dy(j,6), xp) * v(i,6,k,e) &
                     + real(Dy(j,7), xp) * v(i,7,k,e) &
                     + real(Dy(j,8), xp) * v(i,8,k,e) &
                     + real(Dy(j,9), xp) * v(i,9,k,e) &
                     + real(Dy(j,10), xp) * v(i,10,k,e) &
                     + real(Dy(j,11), xp) * v(i,11,k,e) &
                     + real(Dy(j,12), xp) * v(i,12,k,e) &
                     + real(Dy(j,13), xp) * v(i,13,k,e) &
                     + real(Dy(j,14), xp) * v(i,14,k,e)

                wws(i,j,k) = real(Dy(j,1), xp) * w(i,1,k,e) &
                     + real(Dy(j,2), xp) * w(i,2,k,e) &
                     + real(Dy(j,3), xp) * w(i,3,k,e) &
                     + real(Dy(j,4), xp) * w(i,4,k,e) &
                     + real(Dy(j,5), xp) * w(i,5,k,e) &
                     + real(Dy(j,6), xp) * w(i,6,k,e) &
                     + real(Dy(j,7), xp) * w(i,7,k,e) &
                     + real(Dy(j,8), xp) * w(i,8,k,e) &
                     + real(Dy(j,9), xp) * w(i,9,k,e) &
                     + real(Dy(j,10), xp) * w(i,10,k,e) &
                     + real(Dy(j,11), xp) * w(i,11,k,e) &
                     + real(Dy(j,12), xp) * w(i,12,k,e) &
                     + real(Dy(j,13), xp) * w(i,13,k,e) &
                     + real(Dy(j,14), xp) * w(i,14,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = real(Dz(k,1), xp) * u(i,1,1,e) &
                  + real(Dz(k,2), xp) * u(i,1,2,e) &
                  + real(Dz(k,3), xp) * u(i,1,3,e) &
                  + real(Dz(k,4), xp) * u(i,1,4,e) &
                  + real(Dz(k,5), xp) * u(i,1,5,e) &
                  + real(Dz(k,6), xp) * u(i,1,6,e) &
                  + real(Dz(k,7), xp) * u(i,1,7,e) &
                  + real(Dz(k,8), xp) * u(i,1,8,e) &
                  + real(Dz(k,9), xp) * u(i,1,9,e) &
                  + real(Dz(k,10), xp) * u(i,1,10,e) &
                  + real(Dz(k,11), xp) * u(i,1,11,e) &
                  + real(Dz(k,12), xp) * u(i,1,12,e) &
                  + real(Dz(k,13), xp) * u(i,1,13,e) &
                  + real(Dz(k,14), xp) * u(i,1,14,e)

             wvt(i,1,k) = real(Dz(k,1), xp) * v(i,1,1,e) &
                  + real(Dz(k,2), xp) * v(i,1,2,e) &
                  + real(Dz(k,3), xp) * v(i,1,3,e) &
                  + real(Dz(k,4), xp) * v(i,1,4,e) &
                  + real(Dz(k,5), xp) * v(i,1,5,e) &
                  + real(Dz(k,6), xp) * v(i,1,6,e) &
                  + real(Dz(k,7), xp) * v(i,1,7,e) &
                  + real(Dz(k,8), xp) * v(i,1,8,e) &
                  + real(Dz(k,9), xp) * v(i,1,9,e) &
                  + real(Dz(k,10), xp) * v(i,1,10,e) &
                  + real(Dz(k,11), xp) * v(i,1,11,e) &
                  + real(Dz(k,12), xp) * v(i,1,12,e) &
                  + real(Dz(k,13), xp) * v(i,1,13,e) &
                  + real(Dz(k,14), xp) * v(i,1,14,e)

             wwt(i,1,k) = real(Dz(k,1), xp) * w(i,1,1,e) &
                  + real(Dz(k,2), xp) * w(i,1,2,e) &
                  + real(Dz(k,3), xp) * w(i,1,3,e) &
                  + real(Dz(k,4), xp) * w(i,1,4,e) &
                  + real(Dz(k,5), xp) * w(i,1,5,e) &
                  + real(Dz(k,6), xp) * w(i,1,6,e) &
                  + real(Dz(k,7), xp) * w(i,1,7,e) &
                  + real(Dz(k,8), xp) * w(i,1,8,e) &
                  + real(Dz(k,9), xp) * w(i,1,9,e) &
                  + real(Dz(k,10), xp) * w(i,1,10,e) &
                  + real(Dz(k,11), xp) * w(i,1,11,e) &
                  + real(Dz(k,12), xp) * w(i,1,12,e) &
                  + real(Dz(k,13), xp) * w(i,1,13,e) &
                  + real(Dz(k,14), xp) * w(i,1,14,e)
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
             aud(i,j,1) = Dxt(i,1) * ur(1,j,1) &
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

             avd(i,j,1) = Dxt(i,1) * vr(1,j,1) &
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

             awd(i,j,1) = Dxt(i,1) * wr(1,j,1) &
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
                aud(i,j,k) = aud(i,j,k) &
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

                avd(i,j,k) = avd(i,j,k) &
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

                awd(i,j,k) = awd(i,j,k) &
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

       if (ifh2) then
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
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
                     + Dzt(k,14) * ut(i,1,14) &
                     + h2(i,1,k,e) * B(i,1,k,e) * u(i,1,k,e)

                avd(i,1,k) = avd(i,1,k) &
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
                     + Dzt(k,14) * vt(i,1,14) &
                     + h2(i,1,k,e) * B(i,1,k,e) * v(i,1,k,e)

                awd(i,1,k) = awd(i,1,k) &
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
                     + Dzt(k,14) * wt(i,1,14) &
                     + h2(i,1,k,e) * B(i,1,k,e) * w(i,1,k,e)
             end do
          end do
       else
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
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

                avd(i,1,k) = avd(i,1,k) &
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

                awd(i,1,k) = awd(i,1,k) &
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
       end if


       ! Single truncation of the dp-accumulated operator
       do i = 1, lx*lx*lx
          au(i,1,1,e) = aud(i,1,1)
          av(i,1,1,e) = avd(i,1,1)
          aw(i,1,1,e) = awd(i,1,1)
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx14

  subroutine ax_helm_vector_lx13(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, h2, B, ifh2, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 13
    integer, intent(in) :: n
    logical, intent(in) :: ifh2
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: B(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=xp) :: ur(lx, lx, lx)
    real(kind=xp) :: us(lx, lx, lx)
    real(kind=xp) :: ut(lx, lx, lx)
    real(kind=xp) :: vr(lx, lx, lx)
    real(kind=xp) :: vs(lx, lx, lx)
    real(kind=xp) :: vt(lx, lx, lx)
    real(kind=xp) :: wr(lx, lx, lx)
    real(kind=xp) :: ws(lx, lx, lx)
    real(kind=xp) :: wt(lx, lx, lx)
    real(kind=xp) :: wur(lx, lx, lx)
    real(kind=xp) :: wus(lx, lx, lx)
    real(kind=xp) :: wut(lx, lx, lx)
    real(kind=xp) :: wvr(lx, lx, lx)
    real(kind=xp) :: wvs(lx, lx, lx)
    real(kind=xp) :: wvt(lx, lx, lx)
    real(kind=xp) :: wwr(lx, lx, lx)
    real(kind=xp) :: wws(lx, lx, lx)
    real(kind=xp) :: wwt(lx, lx, lx)
    real(kind=xp) :: aud(lx, lx, lx), avd(lx, lx, lx), awd(lx, lx, lx)
    integer :: e, i, j, k

    !$omp do
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = real(Dx(i,1), xp) * u(1,j,1,e) &
                  + real(Dx(i,2), xp) * u(2,j,1,e) &
                  + real(Dx(i,3), xp) * u(3,j,1,e) &
                  + real(Dx(i,4), xp) * u(4,j,1,e) &
                  + real(Dx(i,5), xp) * u(5,j,1,e) &
                  + real(Dx(i,6), xp) * u(6,j,1,e) &
                  + real(Dx(i,7), xp) * u(7,j,1,e) &
                  + real(Dx(i,8), xp) * u(8,j,1,e) &
                  + real(Dx(i,9), xp) * u(9,j,1,e) &
                  + real(Dx(i,10), xp) * u(10,j,1,e) &
                  + real(Dx(i,11), xp) * u(11,j,1,e) &
                  + real(Dx(i,12), xp) * u(12,j,1,e) &
                  + real(Dx(i,13), xp) * u(13,j,1,e)

             wvr(i,j,1) = real(Dx(i,1), xp) * v(1,j,1,e) &
                  + real(Dx(i,2), xp) * v(2,j,1,e) &
                  + real(Dx(i,3), xp) * v(3,j,1,e) &
                  + real(Dx(i,4), xp) * v(4,j,1,e) &
                  + real(Dx(i,5), xp) * v(5,j,1,e) &
                  + real(Dx(i,6), xp) * v(6,j,1,e) &
                  + real(Dx(i,7), xp) * v(7,j,1,e) &
                  + real(Dx(i,8), xp) * v(8,j,1,e) &
                  + real(Dx(i,9), xp) * v(9,j,1,e) &
                  + real(Dx(i,10), xp) * v(10,j,1,e) &
                  + real(Dx(i,11), xp) * v(11,j,1,e) &
                  + real(Dx(i,12), xp) * v(12,j,1,e) &
                  + real(Dx(i,13), xp) * v(13,j,1,e)

             wwr(i,j,1) = real(Dx(i,1), xp) * w(1,j,1,e) &
                  + real(Dx(i,2), xp) * w(2,j,1,e) &
                  + real(Dx(i,3), xp) * w(3,j,1,e) &
                  + real(Dx(i,4), xp) * w(4,j,1,e) &
                  + real(Dx(i,5), xp) * w(5,j,1,e) &
                  + real(Dx(i,6), xp) * w(6,j,1,e) &
                  + real(Dx(i,7), xp) * w(7,j,1,e) &
                  + real(Dx(i,8), xp) * w(8,j,1,e) &
                  + real(Dx(i,9), xp) * w(9,j,1,e) &
                  + real(Dx(i,10), xp) * w(10,j,1,e) &
                  + real(Dx(i,11), xp) * w(11,j,1,e) &
                  + real(Dx(i,12), xp) * w(12,j,1,e) &
                  + real(Dx(i,13), xp) * w(13,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = real(Dy(j,1), xp) * u(i,1,k,e) &
                     + real(Dy(j,2), xp) * u(i,2,k,e) &
                     + real(Dy(j,3), xp) * u(i,3,k,e) &
                     + real(Dy(j,4), xp) * u(i,4,k,e) &
                     + real(Dy(j,5), xp) * u(i,5,k,e) &
                     + real(Dy(j,6), xp) * u(i,6,k,e) &
                     + real(Dy(j,7), xp) * u(i,7,k,e) &
                     + real(Dy(j,8), xp) * u(i,8,k,e) &
                     + real(Dy(j,9), xp) * u(i,9,k,e) &
                     + real(Dy(j,10), xp) * u(i,10,k,e) &
                     + real(Dy(j,11), xp) * u(i,11,k,e) &
                     + real(Dy(j,12), xp) * u(i,12,k,e) &
                     + real(Dy(j,13), xp) * u(i,13,k,e)

                wvs(i,j,k) = real(Dy(j,1), xp) * v(i,1,k,e) &
                     + real(Dy(j,2), xp) * v(i,2,k,e) &
                     + real(Dy(j,3), xp) * v(i,3,k,e) &
                     + real(Dy(j,4), xp) * v(i,4,k,e) &
                     + real(Dy(j,5), xp) * v(i,5,k,e) &
                     + real(Dy(j,6), xp) * v(i,6,k,e) &
                     + real(Dy(j,7), xp) * v(i,7,k,e) &
                     + real(Dy(j,8), xp) * v(i,8,k,e) &
                     + real(Dy(j,9), xp) * v(i,9,k,e) &
                     + real(Dy(j,10), xp) * v(i,10,k,e) &
                     + real(Dy(j,11), xp) * v(i,11,k,e) &
                     + real(Dy(j,12), xp) * v(i,12,k,e) &
                     + real(Dy(j,13), xp) * v(i,13,k,e)

                wws(i,j,k) = real(Dy(j,1), xp) * w(i,1,k,e) &
                     + real(Dy(j,2), xp) * w(i,2,k,e) &
                     + real(Dy(j,3), xp) * w(i,3,k,e) &
                     + real(Dy(j,4), xp) * w(i,4,k,e) &
                     + real(Dy(j,5), xp) * w(i,5,k,e) &
                     + real(Dy(j,6), xp) * w(i,6,k,e) &
                     + real(Dy(j,7), xp) * w(i,7,k,e) &
                     + real(Dy(j,8), xp) * w(i,8,k,e) &
                     + real(Dy(j,9), xp) * w(i,9,k,e) &
                     + real(Dy(j,10), xp) * w(i,10,k,e) &
                     + real(Dy(j,11), xp) * w(i,11,k,e) &
                     + real(Dy(j,12), xp) * w(i,12,k,e) &
                     + real(Dy(j,13), xp) * w(i,13,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = real(Dz(k,1), xp) * u(i,1,1,e) &
                  + real(Dz(k,2), xp) * u(i,1,2,e) &
                  + real(Dz(k,3), xp) * u(i,1,3,e) &
                  + real(Dz(k,4), xp) * u(i,1,4,e) &
                  + real(Dz(k,5), xp) * u(i,1,5,e) &
                  + real(Dz(k,6), xp) * u(i,1,6,e) &
                  + real(Dz(k,7), xp) * u(i,1,7,e) &
                  + real(Dz(k,8), xp) * u(i,1,8,e) &
                  + real(Dz(k,9), xp) * u(i,1,9,e) &
                  + real(Dz(k,10), xp) * u(i,1,10,e) &
                  + real(Dz(k,11), xp) * u(i,1,11,e) &
                  + real(Dz(k,12), xp) * u(i,1,12,e) &
                  + real(Dz(k,13), xp) * u(i,1,13,e)

             wvt(i,1,k) = real(Dz(k,1), xp) * v(i,1,1,e) &
                  + real(Dz(k,2), xp) * v(i,1,2,e) &
                  + real(Dz(k,3), xp) * v(i,1,3,e) &
                  + real(Dz(k,4), xp) * v(i,1,4,e) &
                  + real(Dz(k,5), xp) * v(i,1,5,e) &
                  + real(Dz(k,6), xp) * v(i,1,6,e) &
                  + real(Dz(k,7), xp) * v(i,1,7,e) &
                  + real(Dz(k,8), xp) * v(i,1,8,e) &
                  + real(Dz(k,9), xp) * v(i,1,9,e) &
                  + real(Dz(k,10), xp) * v(i,1,10,e) &
                  + real(Dz(k,11), xp) * v(i,1,11,e) &
                  + real(Dz(k,12), xp) * v(i,1,12,e) &
                  + real(Dz(k,13), xp) * v(i,1,13,e)

             wwt(i,1,k) = real(Dz(k,1), xp) * w(i,1,1,e) &
                  + real(Dz(k,2), xp) * w(i,1,2,e) &
                  + real(Dz(k,3), xp) * w(i,1,3,e) &
                  + real(Dz(k,4), xp) * w(i,1,4,e) &
                  + real(Dz(k,5), xp) * w(i,1,5,e) &
                  + real(Dz(k,6), xp) * w(i,1,6,e) &
                  + real(Dz(k,7), xp) * w(i,1,7,e) &
                  + real(Dz(k,8), xp) * w(i,1,8,e) &
                  + real(Dz(k,9), xp) * w(i,1,9,e) &
                  + real(Dz(k,10), xp) * w(i,1,10,e) &
                  + real(Dz(k,11), xp) * w(i,1,11,e) &
                  + real(Dz(k,12), xp) * w(i,1,12,e) &
                  + real(Dz(k,13), xp) * w(i,1,13,e)
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
             aud(i,j,1) = Dxt(i,1) * ur(1,j,1) &
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

             avd(i,j,1) = Dxt(i,1) * vr(1,j,1) &
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

             awd(i,j,1) = Dxt(i,1) * wr(1,j,1) &
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
                aud(i,j,k) = aud(i,j,k) &
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

                avd(i,j,k) = avd(i,j,k) &
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

                awd(i,j,k) = awd(i,j,k) &
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

       if (ifh2) then
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
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
                     + h2(i,1,k,e) * B(i,1,k,e) * u(i,1,k,e)

                avd(i,1,k) = avd(i,1,k) &
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
                     + h2(i,1,k,e) * B(i,1,k,e) * v(i,1,k,e)

                awd(i,1,k) = awd(i,1,k) &
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
                     + h2(i,1,k,e) * B(i,1,k,e) * w(i,1,k,e)
             end do
          end do
       else
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
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

                avd(i,1,k) = avd(i,1,k) &
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

                awd(i,1,k) = awd(i,1,k) &
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
       end if


       ! Single truncation of the dp-accumulated operator
       do i = 1, lx*lx*lx
          au(i,1,1,e) = aud(i,1,1)
          av(i,1,1,e) = avd(i,1,1)
          aw(i,1,1,e) = awd(i,1,1)
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx13

  subroutine ax_helm_vector_lx12(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, h2, B, ifh2, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 12
    integer, intent(in) :: n
    logical, intent(in) :: ifh2
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: B(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=xp) :: ur(lx, lx, lx)
    real(kind=xp) :: us(lx, lx, lx)
    real(kind=xp) :: ut(lx, lx, lx)
    real(kind=xp) :: vr(lx, lx, lx)
    real(kind=xp) :: vs(lx, lx, lx)
    real(kind=xp) :: vt(lx, lx, lx)
    real(kind=xp) :: wr(lx, lx, lx)
    real(kind=xp) :: ws(lx, lx, lx)
    real(kind=xp) :: wt(lx, lx, lx)
    real(kind=xp) :: wur(lx, lx, lx)
    real(kind=xp) :: wus(lx, lx, lx)
    real(kind=xp) :: wut(lx, lx, lx)
    real(kind=xp) :: wvr(lx, lx, lx)
    real(kind=xp) :: wvs(lx, lx, lx)
    real(kind=xp) :: wvt(lx, lx, lx)
    real(kind=xp) :: wwr(lx, lx, lx)
    real(kind=xp) :: wws(lx, lx, lx)
    real(kind=xp) :: wwt(lx, lx, lx)
    real(kind=xp) :: aud(lx, lx, lx), avd(lx, lx, lx), awd(lx, lx, lx)
    integer :: e, i, j, k

    !$omp do
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = real(Dx(i,1), xp) * u(1,j,1,e) &
                  + real(Dx(i,2), xp) * u(2,j,1,e) &
                  + real(Dx(i,3), xp) * u(3,j,1,e) &
                  + real(Dx(i,4), xp) * u(4,j,1,e) &
                  + real(Dx(i,5), xp) * u(5,j,1,e) &
                  + real(Dx(i,6), xp) * u(6,j,1,e) &
                  + real(Dx(i,7), xp) * u(7,j,1,e) &
                  + real(Dx(i,8), xp) * u(8,j,1,e) &
                  + real(Dx(i,9), xp) * u(9,j,1,e) &
                  + real(Dx(i,10), xp) * u(10,j,1,e) &
                  + real(Dx(i,11), xp) * u(11,j,1,e) &
                  + real(Dx(i,12), xp) * u(12,j,1,e)

             wvr(i,j,1) = real(Dx(i,1), xp) * v(1,j,1,e) &
                  + real(Dx(i,2), xp) * v(2,j,1,e) &
                  + real(Dx(i,3), xp) * v(3,j,1,e) &
                  + real(Dx(i,4), xp) * v(4,j,1,e) &
                  + real(Dx(i,5), xp) * v(5,j,1,e) &
                  + real(Dx(i,6), xp) * v(6,j,1,e) &
                  + real(Dx(i,7), xp) * v(7,j,1,e) &
                  + real(Dx(i,8), xp) * v(8,j,1,e) &
                  + real(Dx(i,9), xp) * v(9,j,1,e) &
                  + real(Dx(i,10), xp) * v(10,j,1,e) &
                  + real(Dx(i,11), xp) * v(11,j,1,e) &
                  + real(Dx(i,12), xp) * v(12,j,1,e)

             wwr(i,j,1) = real(Dx(i,1), xp) * w(1,j,1,e) &
                  + real(Dx(i,2), xp) * w(2,j,1,e) &
                  + real(Dx(i,3), xp) * w(3,j,1,e) &
                  + real(Dx(i,4), xp) * w(4,j,1,e) &
                  + real(Dx(i,5), xp) * w(5,j,1,e) &
                  + real(Dx(i,6), xp) * w(6,j,1,e) &
                  + real(Dx(i,7), xp) * w(7,j,1,e) &
                  + real(Dx(i,8), xp) * w(8,j,1,e) &
                  + real(Dx(i,9), xp) * w(9,j,1,e) &
                  + real(Dx(i,10), xp) * w(10,j,1,e) &
                  + real(Dx(i,11), xp) * w(11,j,1,e) &
                  + real(Dx(i,12), xp) * w(12,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = real(Dy(j,1), xp) * u(i,1,k,e) &
                     + real(Dy(j,2), xp) * u(i,2,k,e) &
                     + real(Dy(j,3), xp) * u(i,3,k,e) &
                     + real(Dy(j,4), xp) * u(i,4,k,e) &
                     + real(Dy(j,5), xp) * u(i,5,k,e) &
                     + real(Dy(j,6), xp) * u(i,6,k,e) &
                     + real(Dy(j,7), xp) * u(i,7,k,e) &
                     + real(Dy(j,8), xp) * u(i,8,k,e) &
                     + real(Dy(j,9), xp) * u(i,9,k,e) &
                     + real(Dy(j,10), xp) * u(i,10,k,e) &
                     + real(Dy(j,11), xp) * u(i,11,k,e) &
                     + real(Dy(j,12), xp) * u(i,12,k,e)

                wvs(i,j,k) = real(Dy(j,1), xp) * v(i,1,k,e) &
                     + real(Dy(j,2), xp) * v(i,2,k,e) &
                     + real(Dy(j,3), xp) * v(i,3,k,e) &
                     + real(Dy(j,4), xp) * v(i,4,k,e) &
                     + real(Dy(j,5), xp) * v(i,5,k,e) &
                     + real(Dy(j,6), xp) * v(i,6,k,e) &
                     + real(Dy(j,7), xp) * v(i,7,k,e) &
                     + real(Dy(j,8), xp) * v(i,8,k,e) &
                     + real(Dy(j,9), xp) * v(i,9,k,e) &
                     + real(Dy(j,10), xp) * v(i,10,k,e) &
                     + real(Dy(j,11), xp) * v(i,11,k,e) &
                     + real(Dy(j,12), xp) * v(i,12,k,e)

                wws(i,j,k) = real(Dy(j,1), xp) * w(i,1,k,e) &
                     + real(Dy(j,2), xp) * w(i,2,k,e) &
                     + real(Dy(j,3), xp) * w(i,3,k,e) &
                     + real(Dy(j,4), xp) * w(i,4,k,e) &
                     + real(Dy(j,5), xp) * w(i,5,k,e) &
                     + real(Dy(j,6), xp) * w(i,6,k,e) &
                     + real(Dy(j,7), xp) * w(i,7,k,e) &
                     + real(Dy(j,8), xp) * w(i,8,k,e) &
                     + real(Dy(j,9), xp) * w(i,9,k,e) &
                     + real(Dy(j,10), xp) * w(i,10,k,e) &
                     + real(Dy(j,11), xp) * w(i,11,k,e) &
                     + real(Dy(j,12), xp) * w(i,12,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = real(Dz(k,1), xp) * u(i,1,1,e) &
                  + real(Dz(k,2), xp) * u(i,1,2,e) &
                  + real(Dz(k,3), xp) * u(i,1,3,e) &
                  + real(Dz(k,4), xp) * u(i,1,4,e) &
                  + real(Dz(k,5), xp) * u(i,1,5,e) &
                  + real(Dz(k,6), xp) * u(i,1,6,e) &
                  + real(Dz(k,7), xp) * u(i,1,7,e) &
                  + real(Dz(k,8), xp) * u(i,1,8,e) &
                  + real(Dz(k,9), xp) * u(i,1,9,e) &
                  + real(Dz(k,10), xp) * u(i,1,10,e) &
                  + real(Dz(k,11), xp) * u(i,1,11,e) &
                  + real(Dz(k,12), xp) * u(i,1,12,e)

             wvt(i,1,k) = real(Dz(k,1), xp) * v(i,1,1,e) &
                  + real(Dz(k,2), xp) * v(i,1,2,e) &
                  + real(Dz(k,3), xp) * v(i,1,3,e) &
                  + real(Dz(k,4), xp) * v(i,1,4,e) &
                  + real(Dz(k,5), xp) * v(i,1,5,e) &
                  + real(Dz(k,6), xp) * v(i,1,6,e) &
                  + real(Dz(k,7), xp) * v(i,1,7,e) &
                  + real(Dz(k,8), xp) * v(i,1,8,e) &
                  + real(Dz(k,9), xp) * v(i,1,9,e) &
                  + real(Dz(k,10), xp) * v(i,1,10,e) &
                  + real(Dz(k,11), xp) * v(i,1,11,e) &
                  + real(Dz(k,12), xp) * v(i,1,12,e)

             wwt(i,1,k) = real(Dz(k,1), xp) * w(i,1,1,e) &
                  + real(Dz(k,2), xp) * w(i,1,2,e) &
                  + real(Dz(k,3), xp) * w(i,1,3,e) &
                  + real(Dz(k,4), xp) * w(i,1,4,e) &
                  + real(Dz(k,5), xp) * w(i,1,5,e) &
                  + real(Dz(k,6), xp) * w(i,1,6,e) &
                  + real(Dz(k,7), xp) * w(i,1,7,e) &
                  + real(Dz(k,8), xp) * w(i,1,8,e) &
                  + real(Dz(k,9), xp) * w(i,1,9,e) &
                  + real(Dz(k,10), xp) * w(i,1,10,e) &
                  + real(Dz(k,11), xp) * w(i,1,11,e) &
                  + real(Dz(k,12), xp) * w(i,1,12,e)
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
             aud(i,j,1) = Dxt(i,1) * ur(1,j,1) &
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

             avd(i,j,1) = Dxt(i,1) * vr(1,j,1) &
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

             awd(i,j,1) = Dxt(i,1) * wr(1,j,1) &
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
                aud(i,j,k) = aud(i,j,k) &
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

                avd(i,j,k) = avd(i,j,k) &
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

                awd(i,j,k) = awd(i,j,k) &
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

       if (ifh2) then
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
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
                     + h2(i,1,k,e) * B(i,1,k,e) * u(i,1,k,e)

                avd(i,1,k) = avd(i,1,k) &
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
                     + h2(i,1,k,e) * B(i,1,k,e) * v(i,1,k,e)

                awd(i,1,k) = awd(i,1,k) &
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
                     + h2(i,1,k,e) * B(i,1,k,e) * w(i,1,k,e)
             end do
          end do
       else
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
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

                avd(i,1,k) = avd(i,1,k) &
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

                awd(i,1,k) = awd(i,1,k) &
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
       end if


       ! Single truncation of the dp-accumulated operator
       do i = 1, lx*lx*lx
          au(i,1,1,e) = aud(i,1,1)
          av(i,1,1,e) = avd(i,1,1)
          aw(i,1,1,e) = awd(i,1,1)
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx12

  subroutine ax_helm_vector_lx11(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, h2, B, ifh2, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 11
    integer, intent(in) :: n
    logical, intent(in) :: ifh2
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: B(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=xp) :: ur(lx, lx, lx)
    real(kind=xp) :: us(lx, lx, lx)
    real(kind=xp) :: ut(lx, lx, lx)
    real(kind=xp) :: vr(lx, lx, lx)
    real(kind=xp) :: vs(lx, lx, lx)
    real(kind=xp) :: vt(lx, lx, lx)
    real(kind=xp) :: wr(lx, lx, lx)
    real(kind=xp) :: ws(lx, lx, lx)
    real(kind=xp) :: wt(lx, lx, lx)
    real(kind=xp) :: wur(lx, lx, lx)
    real(kind=xp) :: wus(lx, lx, lx)
    real(kind=xp) :: wut(lx, lx, lx)
    real(kind=xp) :: wvr(lx, lx, lx)
    real(kind=xp) :: wvs(lx, lx, lx)
    real(kind=xp) :: wvt(lx, lx, lx)
    real(kind=xp) :: wwr(lx, lx, lx)
    real(kind=xp) :: wws(lx, lx, lx)
    real(kind=xp) :: wwt(lx, lx, lx)
    real(kind=xp) :: aud(lx, lx, lx), avd(lx, lx, lx), awd(lx, lx, lx)
    integer :: e, i, j, k

    !$omp do
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = real(Dx(i,1), xp) * u(1,j,1,e) &
                  + real(Dx(i,2), xp) * u(2,j,1,e) &
                  + real(Dx(i,3), xp) * u(3,j,1,e) &
                  + real(Dx(i,4), xp) * u(4,j,1,e) &
                  + real(Dx(i,5), xp) * u(5,j,1,e) &
                  + real(Dx(i,6), xp) * u(6,j,1,e) &
                  + real(Dx(i,7), xp) * u(7,j,1,e) &
                  + real(Dx(i,8), xp) * u(8,j,1,e) &
                  + real(Dx(i,9), xp) * u(9,j,1,e) &
                  + real(Dx(i,10), xp) * u(10,j,1,e) &
                  + real(Dx(i,11), xp) * u(11,j,1,e)

             wvr(i,j,1) = real(Dx(i,1), xp) * v(1,j,1,e) &
                  + real(Dx(i,2), xp) * v(2,j,1,e) &
                  + real(Dx(i,3), xp) * v(3,j,1,e) &
                  + real(Dx(i,4), xp) * v(4,j,1,e) &
                  + real(Dx(i,5), xp) * v(5,j,1,e) &
                  + real(Dx(i,6), xp) * v(6,j,1,e) &
                  + real(Dx(i,7), xp) * v(7,j,1,e) &
                  + real(Dx(i,8), xp) * v(8,j,1,e) &
                  + real(Dx(i,9), xp) * v(9,j,1,e) &
                  + real(Dx(i,10), xp) * v(10,j,1,e) &
                  + real(Dx(i,11), xp) * v(11,j,1,e)

             wwr(i,j,1) = real(Dx(i,1), xp) * w(1,j,1,e) &
                  + real(Dx(i,2), xp) * w(2,j,1,e) &
                  + real(Dx(i,3), xp) * w(3,j,1,e) &
                  + real(Dx(i,4), xp) * w(4,j,1,e) &
                  + real(Dx(i,5), xp) * w(5,j,1,e) &
                  + real(Dx(i,6), xp) * w(6,j,1,e) &
                  + real(Dx(i,7), xp) * w(7,j,1,e) &
                  + real(Dx(i,8), xp) * w(8,j,1,e) &
                  + real(Dx(i,9), xp) * w(9,j,1,e) &
                  + real(Dx(i,10), xp) * w(10,j,1,e) &
                  + real(Dx(i,11), xp) * w(11,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = real(Dy(j,1), xp) * u(i,1,k,e) &
                     + real(Dy(j,2), xp) * u(i,2,k,e) &
                     + real(Dy(j,3), xp) * u(i,3,k,e) &
                     + real(Dy(j,4), xp) * u(i,4,k,e) &
                     + real(Dy(j,5), xp) * u(i,5,k,e) &
                     + real(Dy(j,6), xp) * u(i,6,k,e) &
                     + real(Dy(j,7), xp) * u(i,7,k,e) &
                     + real(Dy(j,8), xp) * u(i,8,k,e) &
                     + real(Dy(j,9), xp) * u(i,9,k,e) &
                     + real(Dy(j,10), xp) * u(i,10,k,e) &
                     + real(Dy(j,11), xp) * u(i,11,k,e)

                wvs(i,j,k) = real(Dy(j,1), xp) * v(i,1,k,e) &
                     + real(Dy(j,2), xp) * v(i,2,k,e) &
                     + real(Dy(j,3), xp) * v(i,3,k,e) &
                     + real(Dy(j,4), xp) * v(i,4,k,e) &
                     + real(Dy(j,5), xp) * v(i,5,k,e) &
                     + real(Dy(j,6), xp) * v(i,6,k,e) &
                     + real(Dy(j,7), xp) * v(i,7,k,e) &
                     + real(Dy(j,8), xp) * v(i,8,k,e) &
                     + real(Dy(j,9), xp) * v(i,9,k,e) &
                     + real(Dy(j,10), xp) * v(i,10,k,e) &
                     + real(Dy(j,11), xp) * v(i,11,k,e)

                wws(i,j,k) = real(Dy(j,1), xp) * w(i,1,k,e) &
                     + real(Dy(j,2), xp) * w(i,2,k,e) &
                     + real(Dy(j,3), xp) * w(i,3,k,e) &
                     + real(Dy(j,4), xp) * w(i,4,k,e) &
                     + real(Dy(j,5), xp) * w(i,5,k,e) &
                     + real(Dy(j,6), xp) * w(i,6,k,e) &
                     + real(Dy(j,7), xp) * w(i,7,k,e) &
                     + real(Dy(j,8), xp) * w(i,8,k,e) &
                     + real(Dy(j,9), xp) * w(i,9,k,e) &
                     + real(Dy(j,10), xp) * w(i,10,k,e) &
                     + real(Dy(j,11), xp) * w(i,11,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = real(Dz(k,1), xp) * u(i,1,1,e) &
                  + real(Dz(k,2), xp) * u(i,1,2,e) &
                  + real(Dz(k,3), xp) * u(i,1,3,e) &
                  + real(Dz(k,4), xp) * u(i,1,4,e) &
                  + real(Dz(k,5), xp) * u(i,1,5,e) &
                  + real(Dz(k,6), xp) * u(i,1,6,e) &
                  + real(Dz(k,7), xp) * u(i,1,7,e) &
                  + real(Dz(k,8), xp) * u(i,1,8,e) &
                  + real(Dz(k,9), xp) * u(i,1,9,e) &
                  + real(Dz(k,10), xp) * u(i,1,10,e) &
                  + real(Dz(k,11), xp) * u(i,1,11,e)

             wvt(i,1,k) = real(Dz(k,1), xp) * v(i,1,1,e) &
                  + real(Dz(k,2), xp) * v(i,1,2,e) &
                  + real(Dz(k,3), xp) * v(i,1,3,e) &
                  + real(Dz(k,4), xp) * v(i,1,4,e) &
                  + real(Dz(k,5), xp) * v(i,1,5,e) &
                  + real(Dz(k,6), xp) * v(i,1,6,e) &
                  + real(Dz(k,7), xp) * v(i,1,7,e) &
                  + real(Dz(k,8), xp) * v(i,1,8,e) &
                  + real(Dz(k,9), xp) * v(i,1,9,e) &
                  + real(Dz(k,10), xp) * v(i,1,10,e) &
                  + real(Dz(k,11), xp) * v(i,1,11,e)

             wwt(i,1,k) = real(Dz(k,1), xp) * w(i,1,1,e) &
                  + real(Dz(k,2), xp) * w(i,1,2,e) &
                  + real(Dz(k,3), xp) * w(i,1,3,e) &
                  + real(Dz(k,4), xp) * w(i,1,4,e) &
                  + real(Dz(k,5), xp) * w(i,1,5,e) &
                  + real(Dz(k,6), xp) * w(i,1,6,e) &
                  + real(Dz(k,7), xp) * w(i,1,7,e) &
                  + real(Dz(k,8), xp) * w(i,1,8,e) &
                  + real(Dz(k,9), xp) * w(i,1,9,e) &
                  + real(Dz(k,10), xp) * w(i,1,10,e) &
                  + real(Dz(k,11), xp) * w(i,1,11,e)
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
             aud(i,j,1) = Dxt(i,1) * ur(1,j,1) &
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

             avd(i,j,1) = Dxt(i,1) * vr(1,j,1) &
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

             awd(i,j,1) = Dxt(i,1) * wr(1,j,1) &
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
                aud(i,j,k) = aud(i,j,k) &
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

                avd(i,j,k) = avd(i,j,k) &
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

                awd(i,j,k) = awd(i,j,k) &
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

       if (ifh2) then
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
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
                     + h2(i,1,k,e) * B(i,1,k,e) * u(i,1,k,e)

                avd(i,1,k) = avd(i,1,k) &
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
                     + h2(i,1,k,e) * B(i,1,k,e) * v(i,1,k,e)

                awd(i,1,k) = awd(i,1,k) &
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
                     + h2(i,1,k,e) * B(i,1,k,e) * w(i,1,k,e)
             end do
          end do
       else
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
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

                avd(i,1,k) = avd(i,1,k) &
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

                awd(i,1,k) = awd(i,1,k) &
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
       end if


       ! Single truncation of the dp-accumulated operator
       do i = 1, lx*lx*lx
          au(i,1,1,e) = aud(i,1,1)
          av(i,1,1,e) = avd(i,1,1)
          aw(i,1,1,e) = awd(i,1,1)
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx11

  subroutine ax_helm_vector_lx10(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, h2, B, ifh2, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 10
    integer, intent(in) :: n
    logical, intent(in) :: ifh2
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: B(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=xp) :: ur(lx, lx, lx)
    real(kind=xp) :: us(lx, lx, lx)
    real(kind=xp) :: ut(lx, lx, lx)
    real(kind=xp) :: vr(lx, lx, lx)
    real(kind=xp) :: vs(lx, lx, lx)
    real(kind=xp) :: vt(lx, lx, lx)
    real(kind=xp) :: wr(lx, lx, lx)
    real(kind=xp) :: ws(lx, lx, lx)
    real(kind=xp) :: wt(lx, lx, lx)
    real(kind=xp) :: wur(lx, lx, lx)
    real(kind=xp) :: wus(lx, lx, lx)
    real(kind=xp) :: wut(lx, lx, lx)
    real(kind=xp) :: wvr(lx, lx, lx)
    real(kind=xp) :: wvs(lx, lx, lx)
    real(kind=xp) :: wvt(lx, lx, lx)
    real(kind=xp) :: wwr(lx, lx, lx)
    real(kind=xp) :: wws(lx, lx, lx)
    real(kind=xp) :: wwt(lx, lx, lx)
    real(kind=xp) :: aud(lx, lx, lx), avd(lx, lx, lx), awd(lx, lx, lx)
    integer :: e, i, j, k

    !$omp do
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = real(Dx(i,1), xp) * u(1,j,1,e) &
                  + real(Dx(i,2), xp) * u(2,j,1,e) &
                  + real(Dx(i,3), xp) * u(3,j,1,e) &
                  + real(Dx(i,4), xp) * u(4,j,1,e) &
                  + real(Dx(i,5), xp) * u(5,j,1,e) &
                  + real(Dx(i,6), xp) * u(6,j,1,e) &
                  + real(Dx(i,7), xp) * u(7,j,1,e) &
                  + real(Dx(i,8), xp) * u(8,j,1,e) &
                  + real(Dx(i,9), xp) * u(9,j,1,e) &
                  + real(Dx(i,10), xp) * u(10,j,1,e)

             wvr(i,j,1) = real(Dx(i,1), xp) * v(1,j,1,e) &
                  + real(Dx(i,2), xp) * v(2,j,1,e) &
                  + real(Dx(i,3), xp) * v(3,j,1,e) &
                  + real(Dx(i,4), xp) * v(4,j,1,e) &
                  + real(Dx(i,5), xp) * v(5,j,1,e) &
                  + real(Dx(i,6), xp) * v(6,j,1,e) &
                  + real(Dx(i,7), xp) * v(7,j,1,e) &
                  + real(Dx(i,8), xp) * v(8,j,1,e) &
                  + real(Dx(i,9), xp) * v(9,j,1,e) &
                  + real(Dx(i,10), xp) * v(10,j,1,e)

             wwr(i,j,1) = real(Dx(i,1), xp) * w(1,j,1,e) &
                  + real(Dx(i,2), xp) * w(2,j,1,e) &
                  + real(Dx(i,3), xp) * w(3,j,1,e) &
                  + real(Dx(i,4), xp) * w(4,j,1,e) &
                  + real(Dx(i,5), xp) * w(5,j,1,e) &
                  + real(Dx(i,6), xp) * w(6,j,1,e) &
                  + real(Dx(i,7), xp) * w(7,j,1,e) &
                  + real(Dx(i,8), xp) * w(8,j,1,e) &
                  + real(Dx(i,9), xp) * w(9,j,1,e) &
                  + real(Dx(i,10), xp) * w(10,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = real(Dy(j,1), xp) * u(i,1,k,e) &
                     + real(Dy(j,2), xp) * u(i,2,k,e) &
                     + real(Dy(j,3), xp) * u(i,3,k,e) &
                     + real(Dy(j,4), xp) * u(i,4,k,e) &
                     + real(Dy(j,5), xp) * u(i,5,k,e) &
                     + real(Dy(j,6), xp) * u(i,6,k,e) &
                     + real(Dy(j,7), xp) * u(i,7,k,e) &
                     + real(Dy(j,8), xp) * u(i,8,k,e) &
                     + real(Dy(j,9), xp) * u(i,9,k,e) &
                     + real(Dy(j,10), xp) * u(i,10,k,e)

                wvs(i,j,k) = real(Dy(j,1), xp) * v(i,1,k,e) &
                     + real(Dy(j,2), xp) * v(i,2,k,e) &
                     + real(Dy(j,3), xp) * v(i,3,k,e) &
                     + real(Dy(j,4), xp) * v(i,4,k,e) &
                     + real(Dy(j,5), xp) * v(i,5,k,e) &
                     + real(Dy(j,6), xp) * v(i,6,k,e) &
                     + real(Dy(j,7), xp) * v(i,7,k,e) &
                     + real(Dy(j,8), xp) * v(i,8,k,e) &
                     + real(Dy(j,9), xp) * v(i,9,k,e) &
                     + real(Dy(j,10), xp) * v(i,10,k,e)

                wws(i,j,k) = real(Dy(j,1), xp) * w(i,1,k,e) &
                     + real(Dy(j,2), xp) * w(i,2,k,e) &
                     + real(Dy(j,3), xp) * w(i,3,k,e) &
                     + real(Dy(j,4), xp) * w(i,4,k,e) &
                     + real(Dy(j,5), xp) * w(i,5,k,e) &
                     + real(Dy(j,6), xp) * w(i,6,k,e) &
                     + real(Dy(j,7), xp) * w(i,7,k,e) &
                     + real(Dy(j,8), xp) * w(i,8,k,e) &
                     + real(Dy(j,9), xp) * w(i,9,k,e) &
                     + real(Dy(j,10), xp) * w(i,10,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = real(Dz(k,1), xp) * u(i,1,1,e) &
                  + real(Dz(k,2), xp) * u(i,1,2,e) &
                  + real(Dz(k,3), xp) * u(i,1,3,e) &
                  + real(Dz(k,4), xp) * u(i,1,4,e) &
                  + real(Dz(k,5), xp) * u(i,1,5,e) &
                  + real(Dz(k,6), xp) * u(i,1,6,e) &
                  + real(Dz(k,7), xp) * u(i,1,7,e) &
                  + real(Dz(k,8), xp) * u(i,1,8,e) &
                  + real(Dz(k,9), xp) * u(i,1,9,e) &
                  + real(Dz(k,10), xp) * u(i,1,10,e)

             wvt(i,1,k) = real(Dz(k,1), xp) * v(i,1,1,e) &
                  + real(Dz(k,2), xp) * v(i,1,2,e) &
                  + real(Dz(k,3), xp) * v(i,1,3,e) &
                  + real(Dz(k,4), xp) * v(i,1,4,e) &
                  + real(Dz(k,5), xp) * v(i,1,5,e) &
                  + real(Dz(k,6), xp) * v(i,1,6,e) &
                  + real(Dz(k,7), xp) * v(i,1,7,e) &
                  + real(Dz(k,8), xp) * v(i,1,8,e) &
                  + real(Dz(k,9), xp) * v(i,1,9,e) &
                  + real(Dz(k,10), xp) * v(i,1,10,e)

             wwt(i,1,k) = real(Dz(k,1), xp) * w(i,1,1,e) &
                  + real(Dz(k,2), xp) * w(i,1,2,e) &
                  + real(Dz(k,3), xp) * w(i,1,3,e) &
                  + real(Dz(k,4), xp) * w(i,1,4,e) &
                  + real(Dz(k,5), xp) * w(i,1,5,e) &
                  + real(Dz(k,6), xp) * w(i,1,6,e) &
                  + real(Dz(k,7), xp) * w(i,1,7,e) &
                  + real(Dz(k,8), xp) * w(i,1,8,e) &
                  + real(Dz(k,9), xp) * w(i,1,9,e) &
                  + real(Dz(k,10), xp) * w(i,1,10,e)
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
             aud(i,j,1) = Dxt(i,1) * ur(1,j,1) &
                  + Dxt(i,2) * ur(2,j,1) &
                  + Dxt(i,3) * ur(3,j,1) &
                  + Dxt(i,4) * ur(4,j,1) &
                  + Dxt(i,5) * ur(5,j,1) &
                  + Dxt(i,6) * ur(6,j,1) &
                  + Dxt(i,7) * ur(7,j,1) &
                  + Dxt(i,8) * ur(8,j,1) &
                  + Dxt(i,9) * ur(9,j,1) &
                  + Dxt(i,10) * ur(10,j,1)

             avd(i,j,1) = Dxt(i,1) * vr(1,j,1) &
                  + Dxt(i,2) * vr(2,j,1) &
                  + Dxt(i,3) * vr(3,j,1) &
                  + Dxt(i,4) * vr(4,j,1) &
                  + Dxt(i,5) * vr(5,j,1) &
                  + Dxt(i,6) * vr(6,j,1) &
                  + Dxt(i,7) * vr(7,j,1) &
                  + Dxt(i,8) * vr(8,j,1) &
                  + Dxt(i,9) * vr(9,j,1) &
                  + Dxt(i,10) * vr(10,j,1)

             awd(i,j,1) = Dxt(i,1) * wr(1,j,1) &
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
                aud(i,j,k) = aud(i,j,k) &
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

                avd(i,j,k) = avd(i,j,k) &
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

                awd(i,j,k) = awd(i,j,k) &
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

       if (ifh2) then
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
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
                     + h2(i,1,k,e) * B(i,1,k,e) * u(i,1,k,e)

                avd(i,1,k) = avd(i,1,k) &
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
                     + h2(i,1,k,e) * B(i,1,k,e) * v(i,1,k,e)

                awd(i,1,k) = awd(i,1,k) &
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
                     + h2(i,1,k,e) * B(i,1,k,e) * w(i,1,k,e)
             end do
          end do
       else
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
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

                avd(i,1,k) = avd(i,1,k) &
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

                awd(i,1,k) = awd(i,1,k) &
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
       end if


       ! Single truncation of the dp-accumulated operator
       do i = 1, lx*lx*lx
          au(i,1,1,e) = aud(i,1,1)
          av(i,1,1,e) = avd(i,1,1)
          aw(i,1,1,e) = awd(i,1,1)
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx10

  subroutine ax_helm_vector_lx9(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, h2, B, ifh2, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 9
    integer, intent(in) :: n
    logical, intent(in) :: ifh2
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: B(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=xp) :: ur(lx, lx, lx)
    real(kind=xp) :: us(lx, lx, lx)
    real(kind=xp) :: ut(lx, lx, lx)
    real(kind=xp) :: vr(lx, lx, lx)
    real(kind=xp) :: vs(lx, lx, lx)
    real(kind=xp) :: vt(lx, lx, lx)
    real(kind=xp) :: wr(lx, lx, lx)
    real(kind=xp) :: ws(lx, lx, lx)
    real(kind=xp) :: wt(lx, lx, lx)
    real(kind=xp) :: wur(lx, lx, lx)
    real(kind=xp) :: wus(lx, lx, lx)
    real(kind=xp) :: wut(lx, lx, lx)
    real(kind=xp) :: wvr(lx, lx, lx)
    real(kind=xp) :: wvs(lx, lx, lx)
    real(kind=xp) :: wvt(lx, lx, lx)
    real(kind=xp) :: wwr(lx, lx, lx)
    real(kind=xp) :: wws(lx, lx, lx)
    real(kind=xp) :: wwt(lx, lx, lx)
    real(kind=xp) :: aud(lx, lx, lx), avd(lx, lx, lx), awd(lx, lx, lx)
    integer :: e, i, j, k

    !$omp do
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = real(Dx(i,1), xp) * u(1,j,1,e) &
                  + real(Dx(i,2), xp) * u(2,j,1,e) &
                  + real(Dx(i,3), xp) * u(3,j,1,e) &
                  + real(Dx(i,4), xp) * u(4,j,1,e) &
                  + real(Dx(i,5), xp) * u(5,j,1,e) &
                  + real(Dx(i,6), xp) * u(6,j,1,e) &
                  + real(Dx(i,7), xp) * u(7,j,1,e) &
                  + real(Dx(i,8), xp) * u(8,j,1,e) &
                  + real(Dx(i,9), xp) * u(9,j,1,e)

             wvr(i,j,1) = real(Dx(i,1), xp) * v(1,j,1,e) &
                  + real(Dx(i,2), xp) * v(2,j,1,e) &
                  + real(Dx(i,3), xp) * v(3,j,1,e) &
                  + real(Dx(i,4), xp) * v(4,j,1,e) &
                  + real(Dx(i,5), xp) * v(5,j,1,e) &
                  + real(Dx(i,6), xp) * v(6,j,1,e) &
                  + real(Dx(i,7), xp) * v(7,j,1,e) &
                  + real(Dx(i,8), xp) * v(8,j,1,e) &
                  + real(Dx(i,9), xp) * v(9,j,1,e)

             wwr(i,j,1) = real(Dx(i,1), xp) * w(1,j,1,e) &
                  + real(Dx(i,2), xp) * w(2,j,1,e) &
                  + real(Dx(i,3), xp) * w(3,j,1,e) &
                  + real(Dx(i,4), xp) * w(4,j,1,e) &
                  + real(Dx(i,5), xp) * w(5,j,1,e) &
                  + real(Dx(i,6), xp) * w(6,j,1,e) &
                  + real(Dx(i,7), xp) * w(7,j,1,e) &
                  + real(Dx(i,8), xp) * w(8,j,1,e) &
                  + real(Dx(i,9), xp) * w(9,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = real(Dy(j,1), xp) * u(i,1,k,e) &
                     + real(Dy(j,2), xp) * u(i,2,k,e) &
                     + real(Dy(j,3), xp) * u(i,3,k,e) &
                     + real(Dy(j,4), xp) * u(i,4,k,e) &
                     + real(Dy(j,5), xp) * u(i,5,k,e) &
                     + real(Dy(j,6), xp) * u(i,6,k,e) &
                     + real(Dy(j,7), xp) * u(i,7,k,e) &
                     + real(Dy(j,8), xp) * u(i,8,k,e) &
                     + real(Dy(j,9), xp) * u(i,9,k,e)

                wvs(i,j,k) = real(Dy(j,1), xp) * v(i,1,k,e) &
                     + real(Dy(j,2), xp) * v(i,2,k,e) &
                     + real(Dy(j,3), xp) * v(i,3,k,e) &
                     + real(Dy(j,4), xp) * v(i,4,k,e) &
                     + real(Dy(j,5), xp) * v(i,5,k,e) &
                     + real(Dy(j,6), xp) * v(i,6,k,e) &
                     + real(Dy(j,7), xp) * v(i,7,k,e) &
                     + real(Dy(j,8), xp) * v(i,8,k,e) &
                     + real(Dy(j,9), xp) * v(i,9,k,e)

                wws(i,j,k) = real(Dy(j,1), xp) * w(i,1,k,e) &
                     + real(Dy(j,2), xp) * w(i,2,k,e) &
                     + real(Dy(j,3), xp) * w(i,3,k,e) &
                     + real(Dy(j,4), xp) * w(i,4,k,e) &
                     + real(Dy(j,5), xp) * w(i,5,k,e) &
                     + real(Dy(j,6), xp) * w(i,6,k,e) &
                     + real(Dy(j,7), xp) * w(i,7,k,e) &
                     + real(Dy(j,8), xp) * w(i,8,k,e) &
                     + real(Dy(j,9), xp) * w(i,9,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = real(Dz(k,1), xp) * u(i,1,1,e) &
                  + real(Dz(k,2), xp) * u(i,1,2,e) &
                  + real(Dz(k,3), xp) * u(i,1,3,e) &
                  + real(Dz(k,4), xp) * u(i,1,4,e) &
                  + real(Dz(k,5), xp) * u(i,1,5,e) &
                  + real(Dz(k,6), xp) * u(i,1,6,e) &
                  + real(Dz(k,7), xp) * u(i,1,7,e) &
                  + real(Dz(k,8), xp) * u(i,1,8,e) &
                  + real(Dz(k,9), xp) * u(i,1,9,e)

             wvt(i,1,k) = real(Dz(k,1), xp) * v(i,1,1,e) &
                  + real(Dz(k,2), xp) * v(i,1,2,e) &
                  + real(Dz(k,3), xp) * v(i,1,3,e) &
                  + real(Dz(k,4), xp) * v(i,1,4,e) &
                  + real(Dz(k,5), xp) * v(i,1,5,e) &
                  + real(Dz(k,6), xp) * v(i,1,6,e) &
                  + real(Dz(k,7), xp) * v(i,1,7,e) &
                  + real(Dz(k,8), xp) * v(i,1,8,e) &
                  + real(Dz(k,9), xp) * v(i,1,9,e)

             wwt(i,1,k) = real(Dz(k,1), xp) * w(i,1,1,e) &
                  + real(Dz(k,2), xp) * w(i,1,2,e) &
                  + real(Dz(k,3), xp) * w(i,1,3,e) &
                  + real(Dz(k,4), xp) * w(i,1,4,e) &
                  + real(Dz(k,5), xp) * w(i,1,5,e) &
                  + real(Dz(k,6), xp) * w(i,1,6,e) &
                  + real(Dz(k,7), xp) * w(i,1,7,e) &
                  + real(Dz(k,8), xp) * w(i,1,8,e) &
                  + real(Dz(k,9), xp) * w(i,1,9,e)
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
             aud(i,j,1) = Dxt(i,1) * ur(1,j,1) &
                  + Dxt(i,2) * ur(2,j,1) &
                  + Dxt(i,3) * ur(3,j,1) &
                  + Dxt(i,4) * ur(4,j,1) &
                  + Dxt(i,5) * ur(5,j,1) &
                  + Dxt(i,6) * ur(6,j,1) &
                  + Dxt(i,7) * ur(7,j,1) &
                  + Dxt(i,8) * ur(8,j,1) &
                  + Dxt(i,9) * ur(9,j,1)

             avd(i,j,1) = Dxt(i,1) * vr(1,j,1) &
                  + Dxt(i,2) * vr(2,j,1) &
                  + Dxt(i,3) * vr(3,j,1) &
                  + Dxt(i,4) * vr(4,j,1) &
                  + Dxt(i,5) * vr(5,j,1) &
                  + Dxt(i,6) * vr(6,j,1) &
                  + Dxt(i,7) * vr(7,j,1) &
                  + Dxt(i,8) * vr(8,j,1) &
                  + Dxt(i,9) * vr(9,j,1)

             awd(i,j,1) = Dxt(i,1) * wr(1,j,1) &
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
                aud(i,j,k) = aud(i,j,k) &
                     + Dyt(j,1) * us(i,1,k) &
                     + Dyt(j,2) * us(i,2,k) &
                     + Dyt(j,3) * us(i,3,k) &
                     + Dyt(j,4) * us(i,4,k) &
                     + Dyt(j,5) * us(i,5,k) &
                     + Dyt(j,6) * us(i,6,k) &
                     + Dyt(j,7) * us(i,7,k) &
                     + Dyt(j,8) * us(i,8,k) &
                     + Dyt(j,9) * us(i,9,k)

                avd(i,j,k) = avd(i,j,k) &
                     + Dyt(j,1) * vs(i,1,k) &
                     + Dyt(j,2) * vs(i,2,k) &
                     + Dyt(j,3) * vs(i,3,k) &
                     + Dyt(j,4) * vs(i,4,k) &
                     + Dyt(j,5) * vs(i,5,k) &
                     + Dyt(j,6) * vs(i,6,k) &
                     + Dyt(j,7) * vs(i,7,k) &
                     + Dyt(j,8) * vs(i,8,k) &
                     + Dyt(j,9) * vs(i,9,k)

                awd(i,j,k) = awd(i,j,k) &
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

       if (ifh2) then
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
                     + Dzt(k,1) * ut(i,1,1) &
                     + Dzt(k,2) * ut(i,1,2) &
                     + Dzt(k,3) * ut(i,1,3) &
                     + Dzt(k,4) * ut(i,1,4) &
                     + Dzt(k,5) * ut(i,1,5) &
                     + Dzt(k,6) * ut(i,1,6) &
                     + Dzt(k,7) * ut(i,1,7) &
                     + Dzt(k,8) * ut(i,1,8) &
                     + Dzt(k,9) * ut(i,1,9) &
                     + h2(i,1,k,e) * B(i,1,k,e) * u(i,1,k,e)

                avd(i,1,k) = avd(i,1,k) &
                     + Dzt(k,1) * vt(i,1,1) &
                     + Dzt(k,2) * vt(i,1,2) &
                     + Dzt(k,3) * vt(i,1,3) &
                     + Dzt(k,4) * vt(i,1,4) &
                     + Dzt(k,5) * vt(i,1,5) &
                     + Dzt(k,6) * vt(i,1,6) &
                     + Dzt(k,7) * vt(i,1,7) &
                     + Dzt(k,8) * vt(i,1,8) &
                     + Dzt(k,9) * vt(i,1,9) &
                     + h2(i,1,k,e) * B(i,1,k,e) * v(i,1,k,e)

                awd(i,1,k) = awd(i,1,k) &
                     + Dzt(k,1) * wt(i,1,1) &
                     + Dzt(k,2) * wt(i,1,2) &
                     + Dzt(k,3) * wt(i,1,3) &
                     + Dzt(k,4) * wt(i,1,4) &
                     + Dzt(k,5) * wt(i,1,5) &
                     + Dzt(k,6) * wt(i,1,6) &
                     + Dzt(k,7) * wt(i,1,7) &
                     + Dzt(k,8) * wt(i,1,8) &
                     + Dzt(k,9) * wt(i,1,9) &
                     + h2(i,1,k,e) * B(i,1,k,e) * w(i,1,k,e)
             end do
          end do
       else
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
                     + Dzt(k,1) * ut(i,1,1) &
                     + Dzt(k,2) * ut(i,1,2) &
                     + Dzt(k,3) * ut(i,1,3) &
                     + Dzt(k,4) * ut(i,1,4) &
                     + Dzt(k,5) * ut(i,1,5) &
                     + Dzt(k,6) * ut(i,1,6) &
                     + Dzt(k,7) * ut(i,1,7) &
                     + Dzt(k,8) * ut(i,1,8) &
                     + Dzt(k,9) * ut(i,1,9)

                avd(i,1,k) = avd(i,1,k) &
                     + Dzt(k,1) * vt(i,1,1) &
                     + Dzt(k,2) * vt(i,1,2) &
                     + Dzt(k,3) * vt(i,1,3) &
                     + Dzt(k,4) * vt(i,1,4) &
                     + Dzt(k,5) * vt(i,1,5) &
                     + Dzt(k,6) * vt(i,1,6) &
                     + Dzt(k,7) * vt(i,1,7) &
                     + Dzt(k,8) * vt(i,1,8) &
                     + Dzt(k,9) * vt(i,1,9)

                awd(i,1,k) = awd(i,1,k) &
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
       end if


       ! Single truncation of the dp-accumulated operator
       do i = 1, lx*lx*lx
          au(i,1,1,e) = aud(i,1,1)
          av(i,1,1,e) = avd(i,1,1)
          aw(i,1,1,e) = awd(i,1,1)
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx9

  subroutine ax_helm_vector_lx8(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, h2, B, ifh2, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 8
    integer, intent(in) :: n
    logical, intent(in) :: ifh2
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: B(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=xp) :: ur(lx, lx, lx)
    real(kind=xp) :: us(lx, lx, lx)
    real(kind=xp) :: ut(lx, lx, lx)
    real(kind=xp) :: vr(lx, lx, lx)
    real(kind=xp) :: vs(lx, lx, lx)
    real(kind=xp) :: vt(lx, lx, lx)
    real(kind=xp) :: wr(lx, lx, lx)
    real(kind=xp) :: ws(lx, lx, lx)
    real(kind=xp) :: wt(lx, lx, lx)
    real(kind=xp) :: wur(lx, lx, lx)
    real(kind=xp) :: wus(lx, lx, lx)
    real(kind=xp) :: wut(lx, lx, lx)
    real(kind=xp) :: wvr(lx, lx, lx)
    real(kind=xp) :: wvs(lx, lx, lx)
    real(kind=xp) :: wvt(lx, lx, lx)
    real(kind=xp) :: wwr(lx, lx, lx)
    real(kind=xp) :: wws(lx, lx, lx)
    real(kind=xp) :: wwt(lx, lx, lx)
    real(kind=xp) :: aud(lx, lx, lx), avd(lx, lx, lx), awd(lx, lx, lx)
    integer :: e, i, j, k

    !$omp do
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = real(Dx(i,1), xp) * u(1,j,1,e) &
                  + real(Dx(i,2), xp) * u(2,j,1,e) &
                  + real(Dx(i,3), xp) * u(3,j,1,e) &
                  + real(Dx(i,4), xp) * u(4,j,1,e) &
                  + real(Dx(i,5), xp) * u(5,j,1,e) &
                  + real(Dx(i,6), xp) * u(6,j,1,e) &
                  + real(Dx(i,7), xp) * u(7,j,1,e) &
                  + real(Dx(i,8), xp) * u(8,j,1,e)

             wvr(i,j,1) = real(Dx(i,1), xp) * v(1,j,1,e) &
                  + real(Dx(i,2), xp) * v(2,j,1,e) &
                  + real(Dx(i,3), xp) * v(3,j,1,e) &
                  + real(Dx(i,4), xp) * v(4,j,1,e) &
                  + real(Dx(i,5), xp) * v(5,j,1,e) &
                  + real(Dx(i,6), xp) * v(6,j,1,e) &
                  + real(Dx(i,7), xp) * v(7,j,1,e) &
                  + real(Dx(i,8), xp) * v(8,j,1,e)

             wwr(i,j,1) = real(Dx(i,1), xp) * w(1,j,1,e) &
                  + real(Dx(i,2), xp) * w(2,j,1,e) &
                  + real(Dx(i,3), xp) * w(3,j,1,e) &
                  + real(Dx(i,4), xp) * w(4,j,1,e) &
                  + real(Dx(i,5), xp) * w(5,j,1,e) &
                  + real(Dx(i,6), xp) * w(6,j,1,e) &
                  + real(Dx(i,7), xp) * w(7,j,1,e) &
                  + real(Dx(i,8), xp) * w(8,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = real(Dy(j,1), xp) * u(i,1,k,e) &
                     + real(Dy(j,2), xp) * u(i,2,k,e) &
                     + real(Dy(j,3), xp) * u(i,3,k,e) &
                     + real(Dy(j,4), xp) * u(i,4,k,e) &
                     + real(Dy(j,5), xp) * u(i,5,k,e) &
                     + real(Dy(j,6), xp) * u(i,6,k,e) &
                     + real(Dy(j,7), xp) * u(i,7,k,e) &
                     + real(Dy(j,8), xp) * u(i,8,k,e)

                wvs(i,j,k) = real(Dy(j,1), xp) * v(i,1,k,e) &
                     + real(Dy(j,2), xp) * v(i,2,k,e) &
                     + real(Dy(j,3), xp) * v(i,3,k,e) &
                     + real(Dy(j,4), xp) * v(i,4,k,e) &
                     + real(Dy(j,5), xp) * v(i,5,k,e) &
                     + real(Dy(j,6), xp) * v(i,6,k,e) &
                     + real(Dy(j,7), xp) * v(i,7,k,e) &
                     + real(Dy(j,8), xp) * v(i,8,k,e)

                wws(i,j,k) = real(Dy(j,1), xp) * w(i,1,k,e) &
                     + real(Dy(j,2), xp) * w(i,2,k,e) &
                     + real(Dy(j,3), xp) * w(i,3,k,e) &
                     + real(Dy(j,4), xp) * w(i,4,k,e) &
                     + real(Dy(j,5), xp) * w(i,5,k,e) &
                     + real(Dy(j,6), xp) * w(i,6,k,e) &
                     + real(Dy(j,7), xp) * w(i,7,k,e) &
                     + real(Dy(j,8), xp) * w(i,8,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = real(Dz(k,1), xp) * u(i,1,1,e) &
                  + real(Dz(k,2), xp) * u(i,1,2,e) &
                  + real(Dz(k,3), xp) * u(i,1,3,e) &
                  + real(Dz(k,4), xp) * u(i,1,4,e) &
                  + real(Dz(k,5), xp) * u(i,1,5,e) &
                  + real(Dz(k,6), xp) * u(i,1,6,e) &
                  + real(Dz(k,7), xp) * u(i,1,7,e) &
                  + real(Dz(k,8), xp) * u(i,1,8,e)

             wvt(i,1,k) = real(Dz(k,1), xp) * v(i,1,1,e) &
                  + real(Dz(k,2), xp) * v(i,1,2,e) &
                  + real(Dz(k,3), xp) * v(i,1,3,e) &
                  + real(Dz(k,4), xp) * v(i,1,4,e) &
                  + real(Dz(k,5), xp) * v(i,1,5,e) &
                  + real(Dz(k,6), xp) * v(i,1,6,e) &
                  + real(Dz(k,7), xp) * v(i,1,7,e) &
                  + real(Dz(k,8), xp) * v(i,1,8,e)

             wwt(i,1,k) = real(Dz(k,1), xp) * w(i,1,1,e) &
                  + real(Dz(k,2), xp) * w(i,1,2,e) &
                  + real(Dz(k,3), xp) * w(i,1,3,e) &
                  + real(Dz(k,4), xp) * w(i,1,4,e) &
                  + real(Dz(k,5), xp) * w(i,1,5,e) &
                  + real(Dz(k,6), xp) * w(i,1,6,e) &
                  + real(Dz(k,7), xp) * w(i,1,7,e) &
                  + real(Dz(k,8), xp) * w(i,1,8,e)
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
             aud(i,j,1) = Dxt(i,1) * ur(1,j,1) &
                  + Dxt(i,2) * ur(2,j,1) &
                  + Dxt(i,3) * ur(3,j,1) &
                  + Dxt(i,4) * ur(4,j,1) &
                  + Dxt(i,5) * ur(5,j,1) &
                  + Dxt(i,6) * ur(6,j,1) &
                  + Dxt(i,7) * ur(7,j,1) &
                  + Dxt(i,8) * ur(8,j,1)

             avd(i,j,1) = Dxt(i,1) * vr(1,j,1) &
                  + Dxt(i,2) * vr(2,j,1) &
                  + Dxt(i,3) * vr(3,j,1) &
                  + Dxt(i,4) * vr(4,j,1) &
                  + Dxt(i,5) * vr(5,j,1) &
                  + Dxt(i,6) * vr(6,j,1) &
                  + Dxt(i,7) * vr(7,j,1) &
                  + Dxt(i,8) * vr(8,j,1)

             awd(i,j,1) = Dxt(i,1) * wr(1,j,1) &
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
                aud(i,j,k) = aud(i,j,k) &
                     + Dyt(j,1) * us(i,1,k) &
                     + Dyt(j,2) * us(i,2,k) &
                     + Dyt(j,3) * us(i,3,k) &
                     + Dyt(j,4) * us(i,4,k) &
                     + Dyt(j,5) * us(i,5,k) &
                     + Dyt(j,6) * us(i,6,k) &
                     + Dyt(j,7) * us(i,7,k) &
                     + Dyt(j,8) * us(i,8,k)

                avd(i,j,k) = avd(i,j,k) &
                     + Dyt(j,1) * vs(i,1,k) &
                     + Dyt(j,2) * vs(i,2,k) &
                     + Dyt(j,3) * vs(i,3,k) &
                     + Dyt(j,4) * vs(i,4,k) &
                     + Dyt(j,5) * vs(i,5,k) &
                     + Dyt(j,6) * vs(i,6,k) &
                     + Dyt(j,7) * vs(i,7,k) &
                     + Dyt(j,8) * vs(i,8,k)

                awd(i,j,k) = awd(i,j,k) &
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

       if (ifh2) then
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
                     + Dzt(k,1) * ut(i,1,1) &
                     + Dzt(k,2) * ut(i,1,2) &
                     + Dzt(k,3) * ut(i,1,3) &
                     + Dzt(k,4) * ut(i,1,4) &
                     + Dzt(k,5) * ut(i,1,5) &
                     + Dzt(k,6) * ut(i,1,6) &
                     + Dzt(k,7) * ut(i,1,7) &
                     + Dzt(k,8) * ut(i,1,8) &
                     + h2(i,1,k,e) * B(i,1,k,e) * u(i,1,k,e)

                avd(i,1,k) = avd(i,1,k) &
                     + Dzt(k,1) * vt(i,1,1) &
                     + Dzt(k,2) * vt(i,1,2) &
                     + Dzt(k,3) * vt(i,1,3) &
                     + Dzt(k,4) * vt(i,1,4) &
                     + Dzt(k,5) * vt(i,1,5) &
                     + Dzt(k,6) * vt(i,1,6) &
                     + Dzt(k,7) * vt(i,1,7) &
                     + Dzt(k,8) * vt(i,1,8) &
                     + h2(i,1,k,e) * B(i,1,k,e) * v(i,1,k,e)

                awd(i,1,k) = awd(i,1,k) &
                     + Dzt(k,1) * wt(i,1,1) &
                     + Dzt(k,2) * wt(i,1,2) &
                     + Dzt(k,3) * wt(i,1,3) &
                     + Dzt(k,4) * wt(i,1,4) &
                     + Dzt(k,5) * wt(i,1,5) &
                     + Dzt(k,6) * wt(i,1,6) &
                     + Dzt(k,7) * wt(i,1,7) &
                     + Dzt(k,8) * wt(i,1,8) &
                     + h2(i,1,k,e) * B(i,1,k,e) * w(i,1,k,e)
             end do
          end do
       else
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
                     + Dzt(k,1) * ut(i,1,1) &
                     + Dzt(k,2) * ut(i,1,2) &
                     + Dzt(k,3) * ut(i,1,3) &
                     + Dzt(k,4) * ut(i,1,4) &
                     + Dzt(k,5) * ut(i,1,5) &
                     + Dzt(k,6) * ut(i,1,6) &
                     + Dzt(k,7) * ut(i,1,7) &
                     + Dzt(k,8) * ut(i,1,8)

                avd(i,1,k) = avd(i,1,k) &
                     + Dzt(k,1) * vt(i,1,1) &
                     + Dzt(k,2) * vt(i,1,2) &
                     + Dzt(k,3) * vt(i,1,3) &
                     + Dzt(k,4) * vt(i,1,4) &
                     + Dzt(k,5) * vt(i,1,5) &
                     + Dzt(k,6) * vt(i,1,6) &
                     + Dzt(k,7) * vt(i,1,7) &
                     + Dzt(k,8) * vt(i,1,8)

                awd(i,1,k) = awd(i,1,k) &
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
       end if


       ! Single truncation of the dp-accumulated operator
       do i = 1, lx*lx*lx
          au(i,1,1,e) = aud(i,1,1)
          av(i,1,1,e) = avd(i,1,1)
          aw(i,1,1,e) = awd(i,1,1)
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx8

  subroutine ax_helm_vector_lx7(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, h2, B, ifh2, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 7
    integer, intent(in) :: n
    logical, intent(in) :: ifh2
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: B(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=xp) :: ur(lx, lx, lx)
    real(kind=xp) :: us(lx, lx, lx)
    real(kind=xp) :: ut(lx, lx, lx)
    real(kind=xp) :: vr(lx, lx, lx)
    real(kind=xp) :: vs(lx, lx, lx)
    real(kind=xp) :: vt(lx, lx, lx)
    real(kind=xp) :: wr(lx, lx, lx)
    real(kind=xp) :: ws(lx, lx, lx)
    real(kind=xp) :: wt(lx, lx, lx)
    real(kind=xp) :: wur(lx, lx, lx)
    real(kind=xp) :: wus(lx, lx, lx)
    real(kind=xp) :: wut(lx, lx, lx)
    real(kind=xp) :: wvr(lx, lx, lx)
    real(kind=xp) :: wvs(lx, lx, lx)
    real(kind=xp) :: wvt(lx, lx, lx)
    real(kind=xp) :: wwr(lx, lx, lx)
    real(kind=xp) :: wws(lx, lx, lx)
    real(kind=xp) :: wwt(lx, lx, lx)
    real(kind=xp) :: aud(lx, lx, lx), avd(lx, lx, lx), awd(lx, lx, lx)
    integer :: e, i, j, k

    !$omp do
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = real(Dx(i,1), xp) * u(1,j,1,e) &
                  + real(Dx(i,2), xp) * u(2,j,1,e) &
                  + real(Dx(i,3), xp) * u(3,j,1,e) &
                  + real(Dx(i,4), xp) * u(4,j,1,e) &
                  + real(Dx(i,5), xp) * u(5,j,1,e) &
                  + real(Dx(i,6), xp) * u(6,j,1,e) &
                  + real(Dx(i,7), xp) * u(7,j,1,e)

             wvr(i,j,1) = real(Dx(i,1), xp) * v(1,j,1,e) &
                  + real(Dx(i,2), xp) * v(2,j,1,e) &
                  + real(Dx(i,3), xp) * v(3,j,1,e) &
                  + real(Dx(i,4), xp) * v(4,j,1,e) &
                  + real(Dx(i,5), xp) * v(5,j,1,e) &
                  + real(Dx(i,6), xp) * v(6,j,1,e) &
                  + real(Dx(i,7), xp) * v(7,j,1,e)

             wwr(i,j,1) = real(Dx(i,1), xp) * w(1,j,1,e) &
                  + real(Dx(i,2), xp) * w(2,j,1,e) &
                  + real(Dx(i,3), xp) * w(3,j,1,e) &
                  + real(Dx(i,4), xp) * w(4,j,1,e) &
                  + real(Dx(i,5), xp) * w(5,j,1,e) &
                  + real(Dx(i,6), xp) * w(6,j,1,e) &
                  + real(Dx(i,7), xp) * w(7,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = real(Dy(j,1), xp) * u(i,1,k,e) &
                     + real(Dy(j,2), xp) * u(i,2,k,e) &
                     + real(Dy(j,3), xp) * u(i,3,k,e) &
                     + real(Dy(j,4), xp) * u(i,4,k,e) &
                     + real(Dy(j,5), xp) * u(i,5,k,e) &
                     + real(Dy(j,6), xp) * u(i,6,k,e) &
                     + real(Dy(j,7), xp) * u(i,7,k,e)

                wvs(i,j,k) = real(Dy(j,1), xp) * v(i,1,k,e) &
                     + real(Dy(j,2), xp) * v(i,2,k,e) &
                     + real(Dy(j,3), xp) * v(i,3,k,e) &
                     + real(Dy(j,4), xp) * v(i,4,k,e) &
                     + real(Dy(j,5), xp) * v(i,5,k,e) &
                     + real(Dy(j,6), xp) * v(i,6,k,e) &
                     + real(Dy(j,7), xp) * v(i,7,k,e)

                wws(i,j,k) = real(Dy(j,1), xp) * w(i,1,k,e) &
                     + real(Dy(j,2), xp) * w(i,2,k,e) &
                     + real(Dy(j,3), xp) * w(i,3,k,e) &
                     + real(Dy(j,4), xp) * w(i,4,k,e) &
                     + real(Dy(j,5), xp) * w(i,5,k,e) &
                     + real(Dy(j,6), xp) * w(i,6,k,e) &
                     + real(Dy(j,7), xp) * w(i,7,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = real(Dz(k,1), xp) * u(i,1,1,e) &
                  + real(Dz(k,2), xp) * u(i,1,2,e) &
                  + real(Dz(k,3), xp) * u(i,1,3,e) &
                  + real(Dz(k,4), xp) * u(i,1,4,e) &
                  + real(Dz(k,5), xp) * u(i,1,5,e) &
                  + real(Dz(k,6), xp) * u(i,1,6,e) &
                  + real(Dz(k,7), xp) * u(i,1,7,e)

             wvt(i,1,k) = real(Dz(k,1), xp) * v(i,1,1,e) &
                  + real(Dz(k,2), xp) * v(i,1,2,e) &
                  + real(Dz(k,3), xp) * v(i,1,3,e) &
                  + real(Dz(k,4), xp) * v(i,1,4,e) &
                  + real(Dz(k,5), xp) * v(i,1,5,e) &
                  + real(Dz(k,6), xp) * v(i,1,6,e) &
                  + real(Dz(k,7), xp) * v(i,1,7,e)

             wwt(i,1,k) = real(Dz(k,1), xp) * w(i,1,1,e) &
                  + real(Dz(k,2), xp) * w(i,1,2,e) &
                  + real(Dz(k,3), xp) * w(i,1,3,e) &
                  + real(Dz(k,4), xp) * w(i,1,4,e) &
                  + real(Dz(k,5), xp) * w(i,1,5,e) &
                  + real(Dz(k,6), xp) * w(i,1,6,e) &
                  + real(Dz(k,7), xp) * w(i,1,7,e)
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
             aud(i,j,1) = Dxt(i,1) * ur(1,j,1) &
                  + Dxt(i,2) * ur(2,j,1) &
                  + Dxt(i,3) * ur(3,j,1) &
                  + Dxt(i,4) * ur(4,j,1) &
                  + Dxt(i,5) * ur(5,j,1) &
                  + Dxt(i,6) * ur(6,j,1) &
                  + Dxt(i,7) * ur(7,j,1)

             avd(i,j,1) = Dxt(i,1) * vr(1,j,1) &
                  + Dxt(i,2) * vr(2,j,1) &
                  + Dxt(i,3) * vr(3,j,1) &
                  + Dxt(i,4) * vr(4,j,1) &
                  + Dxt(i,5) * vr(5,j,1) &
                  + Dxt(i,6) * vr(6,j,1) &
                  + Dxt(i,7) * vr(7,j,1)

             awd(i,j,1) = Dxt(i,1) * wr(1,j,1) &
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
                aud(i,j,k) = aud(i,j,k) &
                     + Dyt(j,1) * us(i,1,k) &
                     + Dyt(j,2) * us(i,2,k) &
                     + Dyt(j,3) * us(i,3,k) &
                     + Dyt(j,4) * us(i,4,k) &
                     + Dyt(j,5) * us(i,5,k) &
                     + Dyt(j,6) * us(i,6,k) &
                     + Dyt(j,7) * us(i,7,k)

                avd(i,j,k) = avd(i,j,k) &
                     + Dyt(j,1) * vs(i,1,k) &
                     + Dyt(j,2) * vs(i,2,k) &
                     + Dyt(j,3) * vs(i,3,k) &
                     + Dyt(j,4) * vs(i,4,k) &
                     + Dyt(j,5) * vs(i,5,k) &
                     + Dyt(j,6) * vs(i,6,k) &
                     + Dyt(j,7) * vs(i,7,k)

                awd(i,j,k) = awd(i,j,k) &
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

       if (ifh2) then
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
                     + Dzt(k,1) * ut(i,1,1) &
                     + Dzt(k,2) * ut(i,1,2) &
                     + Dzt(k,3) * ut(i,1,3) &
                     + Dzt(k,4) * ut(i,1,4) &
                     + Dzt(k,5) * ut(i,1,5) &
                     + Dzt(k,6) * ut(i,1,6) &
                     + Dzt(k,7) * ut(i,1,7) &
                     + h2(i,1,k,e) * B(i,1,k,e) * u(i,1,k,e)

                avd(i,1,k) = avd(i,1,k) &
                     + Dzt(k,1) * vt(i,1,1) &
                     + Dzt(k,2) * vt(i,1,2) &
                     + Dzt(k,3) * vt(i,1,3) &
                     + Dzt(k,4) * vt(i,1,4) &
                     + Dzt(k,5) * vt(i,1,5) &
                     + Dzt(k,6) * vt(i,1,6) &
                     + Dzt(k,7) * vt(i,1,7) &
                     + h2(i,1,k,e) * B(i,1,k,e) * v(i,1,k,e)

                awd(i,1,k) = awd(i,1,k) &
                     + Dzt(k,1) * wt(i,1,1) &
                     + Dzt(k,2) * wt(i,1,2) &
                     + Dzt(k,3) * wt(i,1,3) &
                     + Dzt(k,4) * wt(i,1,4) &
                     + Dzt(k,5) * wt(i,1,5) &
                     + Dzt(k,6) * wt(i,1,6) &
                     + Dzt(k,7) * wt(i,1,7) &
                     + h2(i,1,k,e) * B(i,1,k,e) * w(i,1,k,e)
             end do
          end do
       else
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
                     + Dzt(k,1) * ut(i,1,1) &
                     + Dzt(k,2) * ut(i,1,2) &
                     + Dzt(k,3) * ut(i,1,3) &
                     + Dzt(k,4) * ut(i,1,4) &
                     + Dzt(k,5) * ut(i,1,5) &
                     + Dzt(k,6) * ut(i,1,6) &
                     + Dzt(k,7) * ut(i,1,7)

                avd(i,1,k) = avd(i,1,k) &
                     + Dzt(k,1) * vt(i,1,1) &
                     + Dzt(k,2) * vt(i,1,2) &
                     + Dzt(k,3) * vt(i,1,3) &
                     + Dzt(k,4) * vt(i,1,4) &
                     + Dzt(k,5) * vt(i,1,5) &
                     + Dzt(k,6) * vt(i,1,6) &
                     + Dzt(k,7) * vt(i,1,7)

                awd(i,1,k) = awd(i,1,k) &
                     + Dzt(k,1) * wt(i,1,1) &
                     + Dzt(k,2) * wt(i,1,2) &
                     + Dzt(k,3) * wt(i,1,3) &
                     + Dzt(k,4) * wt(i,1,4) &
                     + Dzt(k,5) * wt(i,1,5) &
                     + Dzt(k,6) * wt(i,1,6) &
                     + Dzt(k,7) * wt(i,1,7)
             end do
          end do
       end if


       ! Single truncation of the dp-accumulated operator
       do i = 1, lx*lx*lx
          au(i,1,1,e) = aud(i,1,1)
          av(i,1,1,e) = avd(i,1,1)
          aw(i,1,1,e) = awd(i,1,1)
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx7

  subroutine ax_helm_vector_lx6(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, h2, B, ifh2, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 6
    integer, intent(in) :: n
    logical, intent(in) :: ifh2
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: B(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=xp) :: ur(lx, lx, lx)
    real(kind=xp) :: us(lx, lx, lx)
    real(kind=xp) :: ut(lx, lx, lx)
    real(kind=xp) :: vr(lx, lx, lx)
    real(kind=xp) :: vs(lx, lx, lx)
    real(kind=xp) :: vt(lx, lx, lx)
    real(kind=xp) :: wr(lx, lx, lx)
    real(kind=xp) :: ws(lx, lx, lx)
    real(kind=xp) :: wt(lx, lx, lx)
    real(kind=xp) :: wur(lx, lx, lx)
    real(kind=xp) :: wus(lx, lx, lx)
    real(kind=xp) :: wut(lx, lx, lx)
    real(kind=xp) :: wvr(lx, lx, lx)
    real(kind=xp) :: wvs(lx, lx, lx)
    real(kind=xp) :: wvt(lx, lx, lx)
    real(kind=xp) :: wwr(lx, lx, lx)
    real(kind=xp) :: wws(lx, lx, lx)
    real(kind=xp) :: wwt(lx, lx, lx)
    real(kind=xp) :: aud(lx, lx, lx), avd(lx, lx, lx), awd(lx, lx, lx)
    integer :: e, i, j, k

    !$omp do
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = real(Dx(i,1), xp) * u(1,j,1,e) &
                  + real(Dx(i,2), xp) * u(2,j,1,e) &
                  + real(Dx(i,3), xp) * u(3,j,1,e) &
                  + real(Dx(i,4), xp) * u(4,j,1,e) &
                  + real(Dx(i,5), xp) * u(5,j,1,e) &
                  + real(Dx(i,6), xp) * u(6,j,1,e)

             wvr(i,j,1) = real(Dx(i,1), xp) * v(1,j,1,e) &
                  + real(Dx(i,2), xp) * v(2,j,1,e) &
                  + real(Dx(i,3), xp) * v(3,j,1,e) &
                  + real(Dx(i,4), xp) * v(4,j,1,e) &
                  + real(Dx(i,5), xp) * v(5,j,1,e) &
                  + real(Dx(i,6), xp) * v(6,j,1,e)

             wwr(i,j,1) = real(Dx(i,1), xp) * w(1,j,1,e) &
                  + real(Dx(i,2), xp) * w(2,j,1,e) &
                  + real(Dx(i,3), xp) * w(3,j,1,e) &
                  + real(Dx(i,4), xp) * w(4,j,1,e) &
                  + real(Dx(i,5), xp) * w(5,j,1,e) &
                  + real(Dx(i,6), xp) * w(6,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = real(Dy(j,1), xp) * u(i,1,k,e) &
                     + real(Dy(j,2), xp) * u(i,2,k,e) &
                     + real(Dy(j,3), xp) * u(i,3,k,e) &
                     + real(Dy(j,4), xp) * u(i,4,k,e) &
                     + real(Dy(j,5), xp) * u(i,5,k,e) &
                     + real(Dy(j,6), xp) * u(i,6,k,e)

                wvs(i,j,k) = real(Dy(j,1), xp) * v(i,1,k,e) &
                     + real(Dy(j,2), xp) * v(i,2,k,e) &
                     + real(Dy(j,3), xp) * v(i,3,k,e) &
                     + real(Dy(j,4), xp) * v(i,4,k,e) &
                     + real(Dy(j,5), xp) * v(i,5,k,e) &
                     + real(Dy(j,6), xp) * v(i,6,k,e)

                wws(i,j,k) = real(Dy(j,1), xp) * w(i,1,k,e) &
                     + real(Dy(j,2), xp) * w(i,2,k,e) &
                     + real(Dy(j,3), xp) * w(i,3,k,e) &
                     + real(Dy(j,4), xp) * w(i,4,k,e) &
                     + real(Dy(j,5), xp) * w(i,5,k,e) &
                     + real(Dy(j,6), xp) * w(i,6,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = real(Dz(k,1), xp) * u(i,1,1,e) &
                  + real(Dz(k,2), xp) * u(i,1,2,e) &
                  + real(Dz(k,3), xp) * u(i,1,3,e) &
                  + real(Dz(k,4), xp) * u(i,1,4,e) &
                  + real(Dz(k,5), xp) * u(i,1,5,e) &
                  + real(Dz(k,6), xp) * u(i,1,6,e)

             wvt(i,1,k) = real(Dz(k,1), xp) * v(i,1,1,e) &
                  + real(Dz(k,2), xp) * v(i,1,2,e) &
                  + real(Dz(k,3), xp) * v(i,1,3,e) &
                  + real(Dz(k,4), xp) * v(i,1,4,e) &
                  + real(Dz(k,5), xp) * v(i,1,5,e) &
                  + real(Dz(k,6), xp) * v(i,1,6,e)

             wwt(i,1,k) = real(Dz(k,1), xp) * w(i,1,1,e) &
                  + real(Dz(k,2), xp) * w(i,1,2,e) &
                  + real(Dz(k,3), xp) * w(i,1,3,e) &
                  + real(Dz(k,4), xp) * w(i,1,4,e) &
                  + real(Dz(k,5), xp) * w(i,1,5,e) &
                  + real(Dz(k,6), xp) * w(i,1,6,e)
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
             aud(i,j,1) = Dxt(i,1) * ur(1,j,1) &
                  + Dxt(i,2) * ur(2,j,1) &
                  + Dxt(i,3) * ur(3,j,1) &
                  + Dxt(i,4) * ur(4,j,1) &
                  + Dxt(i,5) * ur(5,j,1) &
                  + Dxt(i,6) * ur(6,j,1)

             avd(i,j,1) = Dxt(i,1) * vr(1,j,1) &
                  + Dxt(i,2) * vr(2,j,1) &
                  + Dxt(i,3) * vr(3,j,1) &
                  + Dxt(i,4) * vr(4,j,1) &
                  + Dxt(i,5) * vr(5,j,1) &
                  + Dxt(i,6) * vr(6,j,1)

             awd(i,j,1) = Dxt(i,1) * wr(1,j,1) &
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
                aud(i,j,k) = aud(i,j,k) &
                     + Dyt(j,1) * us(i,1,k) &
                     + Dyt(j,2) * us(i,2,k) &
                     + Dyt(j,3) * us(i,3,k) &
                     + Dyt(j,4) * us(i,4,k) &
                     + Dyt(j,5) * us(i,5,k) &
                     + Dyt(j,6) * us(i,6,k)

                avd(i,j,k) = avd(i,j,k) &
                     + Dyt(j,1) * vs(i,1,k) &
                     + Dyt(j,2) * vs(i,2,k) &
                     + Dyt(j,3) * vs(i,3,k) &
                     + Dyt(j,4) * vs(i,4,k) &
                     + Dyt(j,5) * vs(i,5,k) &
                     + Dyt(j,6) * vs(i,6,k)

                awd(i,j,k) = awd(i,j,k) &
                     + Dyt(j,1) * ws(i,1,k) &
                     + Dyt(j,2) * ws(i,2,k) &
                     + Dyt(j,3) * ws(i,3,k) &
                     + Dyt(j,4) * ws(i,4,k) &
                     + Dyt(j,5) * ws(i,5,k) &
                     + Dyt(j,6) * ws(i,6,k)
             end do
          end do
       end do

       if (ifh2) then
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
                     + Dzt(k,1) * ut(i,1,1) &
                     + Dzt(k,2) * ut(i,1,2) &
                     + Dzt(k,3) * ut(i,1,3) &
                     + Dzt(k,4) * ut(i,1,4) &
                     + Dzt(k,5) * ut(i,1,5) &
                     + Dzt(k,6) * ut(i,1,6) &
                     + h2(i,1,k,e) * B(i,1,k,e) * u(i,1,k,e)

                avd(i,1,k) = avd(i,1,k) &
                     + Dzt(k,1) * vt(i,1,1) &
                     + Dzt(k,2) * vt(i,1,2) &
                     + Dzt(k,3) * vt(i,1,3) &
                     + Dzt(k,4) * vt(i,1,4) &
                     + Dzt(k,5) * vt(i,1,5) &
                     + Dzt(k,6) * vt(i,1,6) &
                     + h2(i,1,k,e) * B(i,1,k,e) * v(i,1,k,e)

                awd(i,1,k) = awd(i,1,k) &
                     + Dzt(k,1) * wt(i,1,1) &
                     + Dzt(k,2) * wt(i,1,2) &
                     + Dzt(k,3) * wt(i,1,3) &
                     + Dzt(k,4) * wt(i,1,4) &
                     + Dzt(k,5) * wt(i,1,5) &
                     + Dzt(k,6) * wt(i,1,6) &
                     + h2(i,1,k,e) * B(i,1,k,e) * w(i,1,k,e)
             end do
          end do
       else
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
                     + Dzt(k,1) * ut(i,1,1) &
                     + Dzt(k,2) * ut(i,1,2) &
                     + Dzt(k,3) * ut(i,1,3) &
                     + Dzt(k,4) * ut(i,1,4) &
                     + Dzt(k,5) * ut(i,1,5) &
                     + Dzt(k,6) * ut(i,1,6)

                avd(i,1,k) = avd(i,1,k) &
                     + Dzt(k,1) * vt(i,1,1) &
                     + Dzt(k,2) * vt(i,1,2) &
                     + Dzt(k,3) * vt(i,1,3) &
                     + Dzt(k,4) * vt(i,1,4) &
                     + Dzt(k,5) * vt(i,1,5) &
                     + Dzt(k,6) * vt(i,1,6)

                awd(i,1,k) = awd(i,1,k) &
                     + Dzt(k,1) * wt(i,1,1) &
                     + Dzt(k,2) * wt(i,1,2) &
                     + Dzt(k,3) * wt(i,1,3) &
                     + Dzt(k,4) * wt(i,1,4) &
                     + Dzt(k,5) * wt(i,1,5) &
                     + Dzt(k,6) * wt(i,1,6)
             end do
          end do
       end if


       ! Single truncation of the dp-accumulated operator
       do i = 1, lx*lx*lx
          au(i,1,1,e) = aud(i,1,1)
          av(i,1,1,e) = avd(i,1,1)
          aw(i,1,1,e) = awd(i,1,1)
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx6

  subroutine ax_helm_vector_lx5(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, h2, B, ifh2, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 5
    integer, intent(in) :: n
    logical, intent(in) :: ifh2
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: B(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=xp) :: ur(lx, lx, lx)
    real(kind=xp) :: us(lx, lx, lx)
    real(kind=xp) :: ut(lx, lx, lx)
    real(kind=xp) :: vr(lx, lx, lx)
    real(kind=xp) :: vs(lx, lx, lx)
    real(kind=xp) :: vt(lx, lx, lx)
    real(kind=xp) :: wr(lx, lx, lx)
    real(kind=xp) :: ws(lx, lx, lx)
    real(kind=xp) :: wt(lx, lx, lx)
    real(kind=xp) :: wur(lx, lx, lx)
    real(kind=xp) :: wus(lx, lx, lx)
    real(kind=xp) :: wut(lx, lx, lx)
    real(kind=xp) :: wvr(lx, lx, lx)
    real(kind=xp) :: wvs(lx, lx, lx)
    real(kind=xp) :: wvt(lx, lx, lx)
    real(kind=xp) :: wwr(lx, lx, lx)
    real(kind=xp) :: wws(lx, lx, lx)
    real(kind=xp) :: wwt(lx, lx, lx)
    real(kind=xp) :: aud(lx, lx, lx), avd(lx, lx, lx), awd(lx, lx, lx)
    integer :: e, i, j, k

    !$omp do
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = real(Dx(i,1), xp) * u(1,j,1,e) &
                  + real(Dx(i,2), xp) * u(2,j,1,e) &
                  + real(Dx(i,3), xp) * u(3,j,1,e) &
                  + real(Dx(i,4), xp) * u(4,j,1,e) &
                  + real(Dx(i,5), xp) * u(5,j,1,e)

             wvr(i,j,1) = real(Dx(i,1), xp) * v(1,j,1,e) &
                  + real(Dx(i,2), xp) * v(2,j,1,e) &
                  + real(Dx(i,3), xp) * v(3,j,1,e) &
                  + real(Dx(i,4), xp) * v(4,j,1,e) &
                  + real(Dx(i,5), xp) * v(5,j,1,e)

             wwr(i,j,1) = real(Dx(i,1), xp) * w(1,j,1,e) &
                  + real(Dx(i,2), xp) * w(2,j,1,e) &
                  + real(Dx(i,3), xp) * w(3,j,1,e) &
                  + real(Dx(i,4), xp) * w(4,j,1,e) &
                  + real(Dx(i,5), xp) * w(5,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = real(Dy(j,1), xp) * u(i,1,k,e) &
                     + real(Dy(j,2), xp) * u(i,2,k,e) &
                     + real(Dy(j,3), xp) * u(i,3,k,e) &
                     + real(Dy(j,4), xp) * u(i,4,k,e) &
                     + real(Dy(j,5), xp) * u(i,5,k,e)

                wvs(i,j,k) = real(Dy(j,1), xp) * v(i,1,k,e) &
                     + real(Dy(j,2), xp) * v(i,2,k,e) &
                     + real(Dy(j,3), xp) * v(i,3,k,e) &
                     + real(Dy(j,4), xp) * v(i,4,k,e) &
                     + real(Dy(j,5), xp) * v(i,5,k,e)

                wws(i,j,k) = real(Dy(j,1), xp) * w(i,1,k,e) &
                     + real(Dy(j,2), xp) * w(i,2,k,e) &
                     + real(Dy(j,3), xp) * w(i,3,k,e) &
                     + real(Dy(j,4), xp) * w(i,4,k,e) &
                     + real(Dy(j,5), xp) * w(i,5,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = real(Dz(k,1), xp) * u(i,1,1,e) &
                  + real(Dz(k,2), xp) * u(i,1,2,e) &
                  + real(Dz(k,3), xp) * u(i,1,3,e) &
                  + real(Dz(k,4), xp) * u(i,1,4,e) &
                  + real(Dz(k,5), xp) * u(i,1,5,e)

             wvt(i,1,k) = real(Dz(k,1), xp) * v(i,1,1,e) &
                  + real(Dz(k,2), xp) * v(i,1,2,e) &
                  + real(Dz(k,3), xp) * v(i,1,3,e) &
                  + real(Dz(k,4), xp) * v(i,1,4,e) &
                  + real(Dz(k,5), xp) * v(i,1,5,e)

             wwt(i,1,k) = real(Dz(k,1), xp) * w(i,1,1,e) &
                  + real(Dz(k,2), xp) * w(i,1,2,e) &
                  + real(Dz(k,3), xp) * w(i,1,3,e) &
                  + real(Dz(k,4), xp) * w(i,1,4,e) &
                  + real(Dz(k,5), xp) * w(i,1,5,e)
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
             aud(i,j,1) = Dxt(i,1) * ur(1,j,1) &
                  + Dxt(i,2) * ur(2,j,1) &
                  + Dxt(i,3) * ur(3,j,1) &
                  + Dxt(i,4) * ur(4,j,1) &
                  + Dxt(i,5) * ur(5,j,1)

             avd(i,j,1) = Dxt(i,1) * vr(1,j,1) &
                  + Dxt(i,2) * vr(2,j,1) &
                  + Dxt(i,3) * vr(3,j,1) &
                  + Dxt(i,4) * vr(4,j,1) &
                  + Dxt(i,5) * vr(5,j,1)

             awd(i,j,1) = Dxt(i,1) * wr(1,j,1) &
                  + Dxt(i,2) * wr(2,j,1) &
                  + Dxt(i,3) * wr(3,j,1) &
                  + Dxt(i,4) * wr(4,j,1) &
                  + Dxt(i,5) * wr(5,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                aud(i,j,k) = aud(i,j,k) &
                     + Dyt(j,1) * us(i,1,k) &
                     + Dyt(j,2) * us(i,2,k) &
                     + Dyt(j,3) * us(i,3,k) &
                     + Dyt(j,4) * us(i,4,k) &
                     + Dyt(j,5) * us(i,5,k)

                avd(i,j,k) = avd(i,j,k) &
                     + Dyt(j,1) * vs(i,1,k) &
                     + Dyt(j,2) * vs(i,2,k) &
                     + Dyt(j,3) * vs(i,3,k) &
                     + Dyt(j,4) * vs(i,4,k) &
                     + Dyt(j,5) * vs(i,5,k)

                awd(i,j,k) = awd(i,j,k) &
                     + Dyt(j,1) * ws(i,1,k) &
                     + Dyt(j,2) * ws(i,2,k) &
                     + Dyt(j,3) * ws(i,3,k) &
                     + Dyt(j,4) * ws(i,4,k) &
                     + Dyt(j,5) * ws(i,5,k)
             end do
          end do
       end do

       if (ifh2) then
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
                     + Dzt(k,1) * ut(i,1,1) &
                     + Dzt(k,2) * ut(i,1,2) &
                     + Dzt(k,3) * ut(i,1,3) &
                     + Dzt(k,4) * ut(i,1,4) &
                     + Dzt(k,5) * ut(i,1,5) &
                     + h2(i,1,k,e) * B(i,1,k,e) * u(i,1,k,e)

                avd(i,1,k) = avd(i,1,k) &
                     + Dzt(k,1) * vt(i,1,1) &
                     + Dzt(k,2) * vt(i,1,2) &
                     + Dzt(k,3) * vt(i,1,3) &
                     + Dzt(k,4) * vt(i,1,4) &
                     + Dzt(k,5) * vt(i,1,5) &
                     + h2(i,1,k,e) * B(i,1,k,e) * v(i,1,k,e)

                awd(i,1,k) = awd(i,1,k) &
                     + Dzt(k,1) * wt(i,1,1) &
                     + Dzt(k,2) * wt(i,1,2) &
                     + Dzt(k,3) * wt(i,1,3) &
                     + Dzt(k,4) * wt(i,1,4) &
                     + Dzt(k,5) * wt(i,1,5) &
                     + h2(i,1,k,e) * B(i,1,k,e) * w(i,1,k,e)
             end do
          end do
       else
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
                     + Dzt(k,1) * ut(i,1,1) &
                     + Dzt(k,2) * ut(i,1,2) &
                     + Dzt(k,3) * ut(i,1,3) &
                     + Dzt(k,4) * ut(i,1,4) &
                     + Dzt(k,5) * ut(i,1,5)

                avd(i,1,k) = avd(i,1,k) &
                     + Dzt(k,1) * vt(i,1,1) &
                     + Dzt(k,2) * vt(i,1,2) &
                     + Dzt(k,3) * vt(i,1,3) &
                     + Dzt(k,4) * vt(i,1,4) &
                     + Dzt(k,5) * vt(i,1,5)

                awd(i,1,k) = awd(i,1,k) &
                     + Dzt(k,1) * wt(i,1,1) &
                     + Dzt(k,2) * wt(i,1,2) &
                     + Dzt(k,3) * wt(i,1,3) &
                     + Dzt(k,4) * wt(i,1,4) &
                     + Dzt(k,5) * wt(i,1,5)
             end do
          end do
       end if


       ! Single truncation of the dp-accumulated operator
       do i = 1, lx*lx*lx
          au(i,1,1,e) = aud(i,1,1)
          av(i,1,1,e) = avd(i,1,1)
          aw(i,1,1,e) = awd(i,1,1)
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx5

  subroutine ax_helm_vector_lx4(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, h2, B, ifh2, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 4
    integer, intent(in) :: n
    logical, intent(in) :: ifh2
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: B(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=xp) :: ur(lx, lx, lx)
    real(kind=xp) :: us(lx, lx, lx)
    real(kind=xp) :: ut(lx, lx, lx)
    real(kind=xp) :: vr(lx, lx, lx)
    real(kind=xp) :: vs(lx, lx, lx)
    real(kind=xp) :: vt(lx, lx, lx)
    real(kind=xp) :: wr(lx, lx, lx)
    real(kind=xp) :: ws(lx, lx, lx)
    real(kind=xp) :: wt(lx, lx, lx)
    real(kind=xp) :: wur(lx, lx, lx)
    real(kind=xp) :: wus(lx, lx, lx)
    real(kind=xp) :: wut(lx, lx, lx)
    real(kind=xp) :: wvr(lx, lx, lx)
    real(kind=xp) :: wvs(lx, lx, lx)
    real(kind=xp) :: wvt(lx, lx, lx)
    real(kind=xp) :: wwr(lx, lx, lx)
    real(kind=xp) :: wws(lx, lx, lx)
    real(kind=xp) :: wwt(lx, lx, lx)
    real(kind=xp) :: aud(lx, lx, lx), avd(lx, lx, lx), awd(lx, lx, lx)
    integer :: e, i, j, k

    !$omp do
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = real(Dx(i,1), xp) * u(1,j,1,e) &
                  + real(Dx(i,2), xp) * u(2,j,1,e) &
                  + real(Dx(i,3), xp) * u(3,j,1,e) &
                  + real(Dx(i,4), xp) * u(4,j,1,e)

             wvr(i,j,1) = real(Dx(i,1), xp) * v(1,j,1,e) &
                  + real(Dx(i,2), xp) * v(2,j,1,e) &
                  + real(Dx(i,3), xp) * v(3,j,1,e) &
                  + real(Dx(i,4), xp) * v(4,j,1,e)

             wwr(i,j,1) = real(Dx(i,1), xp) * w(1,j,1,e) &
                  + real(Dx(i,2), xp) * w(2,j,1,e) &
                  + real(Dx(i,3), xp) * w(3,j,1,e) &
                  + real(Dx(i,4), xp) * w(4,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = real(Dy(j,1), xp) * u(i,1,k,e) &
                     + real(Dy(j,2), xp) * u(i,2,k,e) &
                     + real(Dy(j,3), xp) * u(i,3,k,e) &
                     + real(Dy(j,4), xp) * u(i,4,k,e)

                wvs(i,j,k) = real(Dy(j,1), xp) * v(i,1,k,e) &
                     + real(Dy(j,2), xp) * v(i,2,k,e) &
                     + real(Dy(j,3), xp) * v(i,3,k,e) &
                     + real(Dy(j,4), xp) * v(i,4,k,e)

                wws(i,j,k) = real(Dy(j,1), xp) * w(i,1,k,e) &
                     + real(Dy(j,2), xp) * w(i,2,k,e) &
                     + real(Dy(j,3), xp) * w(i,3,k,e) &
                     + real(Dy(j,4), xp) * w(i,4,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = real(Dz(k,1), xp) * u(i,1,1,e) &
                  + real(Dz(k,2), xp) * u(i,1,2,e) &
                  + real(Dz(k,3), xp) * u(i,1,3,e) &
                  + real(Dz(k,4), xp) * u(i,1,4,e)

             wvt(i,1,k) = real(Dz(k,1), xp) * v(i,1,1,e) &
                  + real(Dz(k,2), xp) * v(i,1,2,e) &
                  + real(Dz(k,3), xp) * v(i,1,3,e) &
                  + real(Dz(k,4), xp) * v(i,1,4,e)

             wwt(i,1,k) = real(Dz(k,1), xp) * w(i,1,1,e) &
                  + real(Dz(k,2), xp) * w(i,1,2,e) &
                  + real(Dz(k,3), xp) * w(i,1,3,e) &
                  + real(Dz(k,4), xp) * w(i,1,4,e)
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
             aud(i,j,1) = Dxt(i,1) * ur(1,j,1) &
                  + Dxt(i,2) * ur(2,j,1) &
                  + Dxt(i,3) * ur(3,j,1) &
                  + Dxt(i,4) * ur(4,j,1)

             avd(i,j,1) = Dxt(i,1) * vr(1,j,1) &
                  + Dxt(i,2) * vr(2,j,1) &
                  + Dxt(i,3) * vr(3,j,1) &
                  + Dxt(i,4) * vr(4,j,1)

             awd(i,j,1) = Dxt(i,1) * wr(1,j,1) &
                  + Dxt(i,2) * wr(2,j,1) &
                  + Dxt(i,3) * wr(3,j,1) &
                  + Dxt(i,4) * wr(4,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                aud(i,j,k) = aud(i,j,k) &
                     + Dyt(j,1) * us(i,1,k) &
                     + Dyt(j,2) * us(i,2,k) &
                     + Dyt(j,3) * us(i,3,k) &
                     + Dyt(j,4) * us(i,4,k)

                avd(i,j,k) = avd(i,j,k) &
                     + Dyt(j,1) * vs(i,1,k) &
                     + Dyt(j,2) * vs(i,2,k) &
                     + Dyt(j,3) * vs(i,3,k) &
                     + Dyt(j,4) * vs(i,4,k)

                awd(i,j,k) = awd(i,j,k) &
                     + Dyt(j,1) * ws(i,1,k) &
                     + Dyt(j,2) * ws(i,2,k) &
                     + Dyt(j,3) * ws(i,3,k) &
                     + Dyt(j,4) * ws(i,4,k)
             end do
          end do
       end do

       if (ifh2) then
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
                     + Dzt(k,1) * ut(i,1,1) &
                     + Dzt(k,2) * ut(i,1,2) &
                     + Dzt(k,3) * ut(i,1,3) &
                     + Dzt(k,4) * ut(i,1,4) &
                     + h2(i,1,k,e) * B(i,1,k,e) * u(i,1,k,e)

                avd(i,1,k) = avd(i,1,k) &
                     + Dzt(k,1) * vt(i,1,1) &
                     + Dzt(k,2) * vt(i,1,2) &
                     + Dzt(k,3) * vt(i,1,3) &
                     + Dzt(k,4) * vt(i,1,4) &
                     + h2(i,1,k,e) * B(i,1,k,e) * v(i,1,k,e)

                awd(i,1,k) = awd(i,1,k) &
                     + Dzt(k,1) * wt(i,1,1) &
                     + Dzt(k,2) * wt(i,1,2) &
                     + Dzt(k,3) * wt(i,1,3) &
                     + Dzt(k,4) * wt(i,1,4) &
                     + h2(i,1,k,e) * B(i,1,k,e) * w(i,1,k,e)
             end do
          end do
       else
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
                     + Dzt(k,1) * ut(i,1,1) &
                     + Dzt(k,2) * ut(i,1,2) &
                     + Dzt(k,3) * ut(i,1,3) &
                     + Dzt(k,4) * ut(i,1,4)

                avd(i,1,k) = avd(i,1,k) &
                     + Dzt(k,1) * vt(i,1,1) &
                     + Dzt(k,2) * vt(i,1,2) &
                     + Dzt(k,3) * vt(i,1,3) &
                     + Dzt(k,4) * vt(i,1,4)

                awd(i,1,k) = awd(i,1,k) &
                     + Dzt(k,1) * wt(i,1,1) &
                     + Dzt(k,2) * wt(i,1,2) &
                     + Dzt(k,3) * wt(i,1,3) &
                     + Dzt(k,4) * wt(i,1,4)
             end do
          end do
       end if


       ! Single truncation of the dp-accumulated operator
       do i = 1, lx*lx*lx
          au(i,1,1,e) = aud(i,1,1)
          av(i,1,1,e) = avd(i,1,1)
          aw(i,1,1,e) = awd(i,1,1)
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx4

  subroutine ax_helm_vector_lx3(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, h2, B, ifh2, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 3
    integer, intent(in) :: n
    logical, intent(in) :: ifh2
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: B(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=xp) :: ur(lx, lx, lx)
    real(kind=xp) :: us(lx, lx, lx)
    real(kind=xp) :: ut(lx, lx, lx)
    real(kind=xp) :: vr(lx, lx, lx)
    real(kind=xp) :: vs(lx, lx, lx)
    real(kind=xp) :: vt(lx, lx, lx)
    real(kind=xp) :: wr(lx, lx, lx)
    real(kind=xp) :: ws(lx, lx, lx)
    real(kind=xp) :: wt(lx, lx, lx)
    real(kind=xp) :: wur(lx, lx, lx)
    real(kind=xp) :: wus(lx, lx, lx)
    real(kind=xp) :: wut(lx, lx, lx)
    real(kind=xp) :: wvr(lx, lx, lx)
    real(kind=xp) :: wvs(lx, lx, lx)
    real(kind=xp) :: wvt(lx, lx, lx)
    real(kind=xp) :: wwr(lx, lx, lx)
    real(kind=xp) :: wws(lx, lx, lx)
    real(kind=xp) :: wwt(lx, lx, lx)
    real(kind=xp) :: aud(lx, lx, lx), avd(lx, lx, lx), awd(lx, lx, lx)
    integer :: e, i, j, k

    !$omp do
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = real(Dx(i,1), xp) * u(1,j,1,e) &
                  + real(Dx(i,2), xp) * u(2,j,1,e) &
                  + real(Dx(i,3), xp) * u(3,j,1,e)

             wvr(i,j,1) = real(Dx(i,1), xp) * v(1,j,1,e) &
                  + real(Dx(i,2), xp) * v(2,j,1,e) &
                  + real(Dx(i,3), xp) * v(3,j,1,e)

             wwr(i,j,1) = real(Dx(i,1), xp) * w(1,j,1,e) &
                  + real(Dx(i,2), xp) * w(2,j,1,e) &
                  + real(Dx(i,3), xp) * w(3,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = real(Dy(j,1), xp) * u(i,1,k,e) &
                     + real(Dy(j,2), xp) * u(i,2,k,e) &
                     + real(Dy(j,3), xp) * u(i,3,k,e)

                wvs(i,j,k) = real(Dy(j,1), xp) * v(i,1,k,e) &
                     + real(Dy(j,2), xp) * v(i,2,k,e) &
                     + real(Dy(j,3), xp) * v(i,3,k,e)

                wws(i,j,k) = real(Dy(j,1), xp) * w(i,1,k,e) &
                     + real(Dy(j,2), xp) * w(i,2,k,e) &
                     + real(Dy(j,3), xp) * w(i,3,k,e)
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = real(Dz(k,1), xp) * u(i,1,1,e) &
                  + real(Dz(k,2), xp) * u(i,1,2,e) &
                  + real(Dz(k,3), xp) * u(i,1,3,e)

             wvt(i,1,k) = real(Dz(k,1), xp) * v(i,1,1,e) &
                  + real(Dz(k,2), xp) * v(i,1,2,e) &
                  + real(Dz(k,3), xp) * v(i,1,3,e)

             wwt(i,1,k) = real(Dz(k,1), xp) * w(i,1,1,e) &
                  + real(Dz(k,2), xp) * w(i,1,2,e) &
                  + real(Dz(k,3), xp) * w(i,1,3,e)
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
             aud(i,j,1) = Dxt(i,1) * ur(1,j,1) &
                  + Dxt(i,2) * ur(2,j,1) &
                  + Dxt(i,3) * ur(3,j,1)

             avd(i,j,1) = Dxt(i,1) * vr(1,j,1) &
                  + Dxt(i,2) * vr(2,j,1) &
                  + Dxt(i,3) * vr(3,j,1)

             awd(i,j,1) = Dxt(i,1) * wr(1,j,1) &
                  + Dxt(i,2) * wr(2,j,1) &
                  + Dxt(i,3) * wr(3,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                aud(i,j,k) = aud(i,j,k) &
                     + Dyt(j,1) * us(i,1,k) &
                     + Dyt(j,2) * us(i,2,k) &
                     + Dyt(j,3) * us(i,3,k)

                avd(i,j,k) = avd(i,j,k) &
                     + Dyt(j,1) * vs(i,1,k) &
                     + Dyt(j,2) * vs(i,2,k) &
                     + Dyt(j,3) * vs(i,3,k)

                awd(i,j,k) = awd(i,j,k) &
                     + Dyt(j,1) * ws(i,1,k) &
                     + Dyt(j,2) * ws(i,2,k) &
                     + Dyt(j,3) * ws(i,3,k)
             end do
          end do
       end do

       if (ifh2) then
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
                     + Dzt(k,1) * ut(i,1,1) &
                     + Dzt(k,2) * ut(i,1,2) &
                     + Dzt(k,3) * ut(i,1,3) &
                     + h2(i,1,k,e) * B(i,1,k,e) * u(i,1,k,e)

                avd(i,1,k) = avd(i,1,k) &
                     + Dzt(k,1) * vt(i,1,1) &
                     + Dzt(k,2) * vt(i,1,2) &
                     + Dzt(k,3) * vt(i,1,3) &
                     + h2(i,1,k,e) * B(i,1,k,e) * v(i,1,k,e)

                awd(i,1,k) = awd(i,1,k) &
                     + Dzt(k,1) * wt(i,1,1) &
                     + Dzt(k,2) * wt(i,1,2) &
                     + Dzt(k,3) * wt(i,1,3) &
                     + h2(i,1,k,e) * B(i,1,k,e) * w(i,1,k,e)
             end do
          end do
       else
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
                     + Dzt(k,1) * ut(i,1,1) &
                     + Dzt(k,2) * ut(i,1,2) &
                     + Dzt(k,3) * ut(i,1,3)

                avd(i,1,k) = avd(i,1,k) &
                     + Dzt(k,1) * vt(i,1,1) &
                     + Dzt(k,2) * vt(i,1,2) &
                     + Dzt(k,3) * vt(i,1,3)

                awd(i,1,k) = awd(i,1,k) &
                     + Dzt(k,1) * wt(i,1,1) &
                     + Dzt(k,2) * wt(i,1,2) &
                     + Dzt(k,3) * wt(i,1,3)
             end do
          end do
       end if


       ! Single truncation of the dp-accumulated operator
       do i = 1, lx*lx*lx
          au(i,1,1,e) = aud(i,1,1)
          av(i,1,1,e) = avd(i,1,1)
          aw(i,1,1,e) = awd(i,1,1)
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx3

  subroutine ax_helm_vector_lx2(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, h2, B, ifh2, G11, G22, G33, G12, G13, G23, n)
    integer, parameter :: lx = 2
    integer, intent(in) :: n
    logical, intent(in) :: ifh2
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: v(lx, lx, lx, n)
    real(kind=rp), intent(in) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: B(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=xp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    real(kind=xp) :: ur(lx, lx, lx)
    real(kind=xp) :: us(lx, lx, lx)
    real(kind=xp) :: ut(lx, lx, lx)
    real(kind=xp) :: vr(lx, lx, lx)
    real(kind=xp) :: vs(lx, lx, lx)
    real(kind=xp) :: vt(lx, lx, lx)
    real(kind=xp) :: wr(lx, lx, lx)
    real(kind=xp) :: ws(lx, lx, lx)
    real(kind=xp) :: wt(lx, lx, lx)
    real(kind=xp) :: wur(lx, lx, lx)
    real(kind=xp) :: wus(lx, lx, lx)
    real(kind=xp) :: wut(lx, lx, lx)
    real(kind=xp) :: wvr(lx, lx, lx)
    real(kind=xp) :: wvs(lx, lx, lx)
    real(kind=xp) :: wvt(lx, lx, lx)
    real(kind=xp) :: wwr(lx, lx, lx)
    real(kind=xp) :: wws(lx, lx, lx)
    real(kind=xp) :: wwt(lx, lx, lx)
    real(kind=xp) :: aud(lx, lx, lx), avd(lx, lx, lx), awd(lx, lx, lx)
    integer :: e, i, j, k

    !$omp do
    do e = 1, n
       do j = 1, lx * lx
          do i = 1, lx
             wur(i,j,1) = real(Dx(i,1), xp) * u(1,j,1,e) &
                  + real(Dx(i,2), xp) * u(2,j,1,e)

             wvr(i,j,1) = real(Dx(i,1), xp) * v(1,j,1,e) &
                  + real(Dx(i,2), xp) * v(2,j,1,e)

             wwr(i,j,1) = real(Dx(i,1), xp) * w(1,j,1,e) &
                  + real(Dx(i,2), xp) * w(2,j,1,e)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                wus(i,j,k) = real(Dy(j,1), xp) * u(i,1,k,e) &
                     + real(Dy(j,2), xp) * u(i,2,k,e)

                wvs(i,j,k) = real(Dy(j,1), xp) * v(i,1,k,e) &
                     + real(Dy(j,2), xp) * v(i,2,k,e)

                wws(i,j,k) = real(Dy(j,1), xp) * w(i,1,k,e) &
                     + real(Dy(j,2), xp) * w(i,2,k,e)

             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             wut(i,1,k) = real(Dz(k,1), xp) * u(i,1,1,e) &
                  + real(Dz(k,2), xp) * u(i,1,2,e)

             wvt(i,1,k) = real(Dz(k,1), xp) * v(i,1,1,e) &
                  + real(Dz(k,2), xp) * v(i,1,2,e)

             wwt(i,1,k) = real(Dz(k,1), xp) * w(i,1,1,e) &
                  + real(Dz(k,2), xp) * w(i,1,2,e)
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
             aud(i,j,1) = Dxt(i,1) * ur(1,j,1) &
                  + Dxt(i,2) * ur(2,j,1)

             avd(i,j,1) = Dxt(i,1) * vr(1,j,1) &
                  + Dxt(i,2) * vr(2,j,1)

             awd(i,j,1) = Dxt(i,1) * wr(1,j,1) &
                  + Dxt(i,2) * wr(2,j,1)
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                aud(i,j,k) = aud(i,j,k) &
                     + Dyt(j,1) * us(i,1,k) &
                     + Dyt(j,2) * us(i,2,k)

                avd(i,j,k) = avd(i,j,k) &
                     + Dyt(j,1) * vs(i,1,k) &
                     + Dyt(j,2) * vs(i,2,k)

                awd(i,j,k) = awd(i,j,k) &
                     + Dyt(j,1) * ws(i,1,k) &
                     + Dyt(j,2) * ws(i,2,k)
             end do
          end do
       end do

       if (ifh2) then
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
                     + Dzt(k,1) * ut(i,1,1) &
                     + Dzt(k,2) * ut(i,1,2) &
                     + h2(i,1,k,e) * B(i,1,k,e) * u(i,1,k,e)

                avd(i,1,k) = avd(i,1,k) &
                     + Dzt(k,1) * vt(i,1,1) &
                     + Dzt(k,2) * vt(i,1,2) &
                     + h2(i,1,k,e) * B(i,1,k,e) * v(i,1,k,e)

                awd(i,1,k) = awd(i,1,k) &
                     + Dzt(k,1) * wt(i,1,1) &
                     + Dzt(k,2) * wt(i,1,2) &
                     + h2(i,1,k,e) * B(i,1,k,e) * w(i,1,k,e)
             end do
          end do
       else
          do k = 1, lx
             do i = 1, lx*lx
                aud(i,1,k) = aud(i,1,k) &
                     + Dzt(k,1) * ut(i,1,1) &
                     + Dzt(k,2) * ut(i,1,2)

                avd(i,1,k) = avd(i,1,k) &
                     + Dzt(k,1) * vt(i,1,1) &
                     + Dzt(k,2) * vt(i,1,2)

                awd(i,1,k) = awd(i,1,k) &
                     + Dzt(k,1) * wt(i,1,1) &
                     + Dzt(k,2) * wt(i,1,2)
             end do
          end do
       end if


       ! Single truncation of the dp-accumulated operator
       do i = 1, lx*lx*lx
          au(i,1,1,e) = aud(i,1,1)
          av(i,1,1,e) = avd(i,1,1)
          aw(i,1,1,e) = awd(i,1,1)
       end do

    end do
    !$omp end do
  end subroutine ax_helm_vector_lx2

end submodule ax_helm_vector_cpu
