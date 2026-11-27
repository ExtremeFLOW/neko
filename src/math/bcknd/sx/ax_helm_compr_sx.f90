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
!> SX Helmholtz matrix-vector product on compressed geometric factors
!! @note The kernels are those of ax_helm_sx with the geometric factors read
!! at G(:,:,:,ind(e)), so a change to one belongs in the other.
module ax_helm_compr_sx
  use ax_helm_sx, only : ax_helm_sx_t
  use num_types, only : rp
  use coefs, only : coef_t
  use space, only : space_t
  use mesh, only : mesh_t
  use math, only : addcol4
  use utils, only : neko_error
  implicit none
  private

  !> SX matrix-vector product for a Helmholtz problem, reading the geometric
  !! factors \f$ G_{ij} \f$ from their compressed form in coef_t.
  !! @details Element e reads the copy coef%compression_inds(e) of
  !! coef%G11_compr etc. Without compressed factors, see
  !! coef_t::enable_geo_compression, the product is the one of ax_helm_sx_t.
  type, public, extends(ax_helm_sx_t) :: ax_helm_compr_sx_t
   contains
     procedure, pass(this) :: compute => ax_helm_compr_sx_compute
  end type ax_helm_compr_sx_t

contains

  subroutine ax_helm_compr_sx_compute(this, w, u, coef, msh, Xh)
    class(ax_helm_compr_sx_t), intent(in) :: this
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    integer :: m

    if (.not. allocated(coef%compression_inds)) then
       call this%ax_helm_sx_t%compute(w, u, coef, msh, Xh)
       return
    end if

    if (size(coef%compression_inds) .ne. msh%nelv) then
       call neko_error('Compressed geometric factors do not match the mesh')
    end if
    m = size(coef%G11_compr, 4)

    select case(Xh%lx)
    case(14)
       call sx_ax_helm_compr_lx14(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, &
            Xh%dyt, Xh%dzt, coef%h1, coef%G11_compr, &
            coef%G22_compr, coef%G33_compr, coef%G12_compr, &
            coef%G13_compr, coef%G23_compr, &
            coef%compression_inds, msh%nelv, m)
    case(13)
       call sx_ax_helm_compr_lx13(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, &
            Xh%dyt, Xh%dzt, coef%h1, coef%G11_compr, &
            coef%G22_compr, coef%G33_compr, coef%G12_compr, &
            coef%G13_compr, coef%G23_compr, &
            coef%compression_inds, msh%nelv, m)
    case(12)
       call sx_ax_helm_compr_lx12(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, &
            Xh%dyt, Xh%dzt, coef%h1, coef%G11_compr, &
            coef%G22_compr, coef%G33_compr, coef%G12_compr, &
            coef%G13_compr, coef%G23_compr, &
            coef%compression_inds, msh%nelv, m)
    case(11)
       call sx_ax_helm_compr_lx11(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, &
            Xh%dyt, Xh%dzt, coef%h1, coef%G11_compr, &
            coef%G22_compr, coef%G33_compr, coef%G12_compr, &
            coef%G13_compr, coef%G23_compr, &
            coef%compression_inds, msh%nelv, m)
    case(10)
       call sx_ax_helm_compr_lx10(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, &
            Xh%dyt, Xh%dzt, coef%h1, coef%G11_compr, &
            coef%G22_compr, coef%G33_compr, coef%G12_compr, &
            coef%G13_compr, coef%G23_compr, &
            coef%compression_inds, msh%nelv, m)
    case(9)
       call sx_ax_helm_compr_lx9(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, &
            Xh%dyt, Xh%dzt, coef%h1, coef%G11_compr, &
            coef%G22_compr, coef%G33_compr, coef%G12_compr, &
            coef%G13_compr, coef%G23_compr, &
            coef%compression_inds, msh%nelv, m)
    case(8)
       call sx_ax_helm_compr_lx8(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, &
            Xh%dyt, Xh%dzt, coef%h1, coef%G11_compr, &
            coef%G22_compr, coef%G33_compr, coef%G12_compr, &
            coef%G13_compr, coef%G23_compr, &
            coef%compression_inds, msh%nelv, m)
    case(7)
       call sx_ax_helm_compr_lx7(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, &
            Xh%dyt, Xh%dzt, coef%h1, coef%G11_compr, &
            coef%G22_compr, coef%G33_compr, coef%G12_compr, &
            coef%G13_compr, coef%G23_compr, &
            coef%compression_inds, msh%nelv, m)
    case(6)
       call sx_ax_helm_compr_lx6(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, &
            Xh%dyt, Xh%dzt, coef%h1, coef%G11_compr, &
            coef%G22_compr, coef%G33_compr, coef%G12_compr, &
            coef%G13_compr, coef%G23_compr, &
            coef%compression_inds, msh%nelv, m)
    case(5)
       call sx_ax_helm_compr_lx5(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, &
            Xh%dyt, Xh%dzt, coef%h1, coef%G11_compr, &
            coef%G22_compr, coef%G33_compr, coef%G12_compr, &
            coef%G13_compr, coef%G23_compr, &
            coef%compression_inds, msh%nelv, m)
    case(4)
       call sx_ax_helm_compr_lx4(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, &
            Xh%dyt, Xh%dzt, coef%h1, coef%G11_compr, &
            coef%G22_compr, coef%G33_compr, coef%G12_compr, &
            coef%G13_compr, coef%G23_compr, &
            coef%compression_inds, msh%nelv, m)
    case(3)
       call sx_ax_helm_compr_lx3(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, &
            Xh%dyt, Xh%dzt, coef%h1, coef%G11_compr, &
            coef%G22_compr, coef%G33_compr, coef%G12_compr, &
            coef%G13_compr, coef%G23_compr, &
            coef%compression_inds, msh%nelv, m)
    case(2)
       call sx_ax_helm_compr_lx2(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, &
            Xh%dyt, Xh%dzt, coef%h1, coef%G11_compr, &
            coef%G22_compr, coef%G33_compr, coef%G12_compr, &
            coef%G13_compr, coef%G23_compr, &
            coef%compression_inds, msh%nelv, m)
    case default
       call sx_ax_helm_compr_lx(w, u, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, &
            Xh%dyt, Xh%dzt, coef%h1, coef%G11_compr, &
            coef%G22_compr, coef%G33_compr, coef%G12_compr, &
            coef%G13_compr, coef%G23_compr, &
            coef%compression_inds, msh%nelv, m, Xh%lx)
    end select

    if (coef%ifh2) call addcol4(w, coef%h2, coef%B, u, coef%dof%size())

  end subroutine ax_helm_compr_sx_compute

  subroutine sx_ax_helm_compr_lx(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, ind, n, m, lx)
    integer, intent(in) :: n, m, lx
    integer, intent(in) :: ind(n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, m)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    integer :: e, c, i, j, k, jj, kk
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: uur(lx, lx, lx, n)
    real(kind=rp) :: uus(lx, lx, lx, n)
    real(kind=rp) :: uut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + Dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          c = ind(e)
          uur(i,1,1,e) = h1(i,1,1,e) * &
               ( G11(i,1,1,c) * ur(i,1,1,e) &
               + G12(i,1,1,c) * us(i,1,1,e) &
               + G13(i,1,1,c) * ut(i,1,1,e))

          uus(i,1,1,e) = h1(i,1,1,e) * &
               ( G22(i,1,1,c) * us(i,1,1,e) &
               + G12(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * ut(i,1,1,e) )

          uut(i,1,1,e) = h1(i,1,1,e) * &
               ( G33(i,1,1,c) * ut(i,1,1,e) &
               + G13(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * us(i,1,1,e))
       end do
    end do

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dxt(i,kk) * uur(kk,jj,1,1)
          end do
          w(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dyt(j, kk)*uus(i,kk,k,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dzt(k, kk)*uut(i,j,kk,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + wt
             end do
          end do
       end do
    end do


  end subroutine sx_ax_helm_compr_lx

  subroutine sx_ax_helm_compr_lx14(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, ind, n, m)
    integer, parameter :: lx = 14
    integer, intent(in) :: n, m
    integer, intent(in) :: ind(n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, m)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    integer :: e, c, i, j, k, jj, kk
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: uur(lx, lx, lx, n)
    real(kind=rp) :: uus(lx, lx, lx, n)
    real(kind=rp) :: uut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + Dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          c = ind(e)
          uur(i,1,1,e) = h1(i,1,1,e) * &
               ( G11(i,1,1,c) * ur(i,1,1,e) &
               + G12(i,1,1,c) * us(i,1,1,e) &
               + G13(i,1,1,c) * ut(i,1,1,e))

          uus(i,1,1,e) = h1(i,1,1,e) * &
               ( G22(i,1,1,c) * us(i,1,1,e) &
               + G12(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * ut(i,1,1,e) )

          uut(i,1,1,e) = h1(i,1,1,e) * &
               ( G33(i,1,1,c) * ut(i,1,1,e) &
               + G13(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * us(i,1,1,e))
       end do
    end do

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dxt(i,kk) * uur(kk,jj,1,1)
          end do
          w(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dyt(j, kk)*uus(i,kk,k,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dzt(k, kk)*uut(i,j,kk,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + wt
             end do
          end do
       end do
    end do


  end subroutine sx_ax_helm_compr_lx14

  subroutine sx_ax_helm_compr_lx13(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, ind, n, m)
    integer, parameter :: lx = 13
    integer, intent(in) :: n, m
    integer, intent(in) :: ind(n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, m)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    integer :: e, c, i, j, k, jj, kk
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: uur(lx, lx, lx, n)
    real(kind=rp) :: uus(lx, lx, lx, n)
    real(kind=rp) :: uut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + Dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          c = ind(e)
          uur(i,1,1,e) = h1(i,1,1,e) * &
               ( G11(i,1,1,c) * ur(i,1,1,e) &
               + G12(i,1,1,c) * us(i,1,1,e) &
               + G13(i,1,1,c) * ut(i,1,1,e))

          uus(i,1,1,e) = h1(i,1,1,e) * &
               ( G22(i,1,1,c) * us(i,1,1,e) &
               + G12(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * ut(i,1,1,e) )

          uut(i,1,1,e) = h1(i,1,1,e) * &
               ( G33(i,1,1,c) * ut(i,1,1,e) &
               + G13(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * us(i,1,1,e))
       end do
    end do

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dxt(i,kk) * uur(kk,jj,1,1)
          end do
          w(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dyt(j, kk)*uus(i,kk,k,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dzt(k, kk)*uut(i,j,kk,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + wt
             end do
          end do
       end do
    end do


  end subroutine sx_ax_helm_compr_lx13

  subroutine sx_ax_helm_compr_lx12(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, ind, n, m)
    integer, parameter :: lx = 12
    integer, intent(in) :: n, m
    integer, intent(in) :: ind(n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, m)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    integer :: e, c, i, j, k, jj, kk
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: uur(lx, lx, lx, n)
    real(kind=rp) :: uus(lx, lx, lx, n)
    real(kind=rp) :: uut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + Dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          c = ind(e)
          uur(i,1,1,e) = h1(i,1,1,e) * &
               ( G11(i,1,1,c) * ur(i,1,1,e) &
               + G12(i,1,1,c) * us(i,1,1,e) &
               + G13(i,1,1,c) * ut(i,1,1,e))

          uus(i,1,1,e) = h1(i,1,1,e) * &
               ( G22(i,1,1,c) * us(i,1,1,e) &
               + G12(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * ut(i,1,1,e) )

          uut(i,1,1,e) = h1(i,1,1,e) * &
               ( G33(i,1,1,c) * ut(i,1,1,e) &
               + G13(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * us(i,1,1,e))
       end do
    end do

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dxt(i,kk) * uur(kk,jj,1,1)
          end do
          w(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dyt(j, kk)*uus(i,kk,k,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dzt(k, kk)*uut(i,j,kk,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + wt
             end do
          end do
       end do
    end do


  end subroutine sx_ax_helm_compr_lx12

  subroutine sx_ax_helm_compr_lx11(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, ind, n, m)
    integer, parameter :: lx = 11
    integer, intent(in) :: n, m
    integer, intent(in) :: ind(n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, m)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    integer :: e, c, i, j, k, jj, kk
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: uur(lx, lx, lx, n)
    real(kind=rp) :: uus(lx, lx, lx, n)
    real(kind=rp) :: uut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + Dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          c = ind(e)
          uur(i,1,1,e) = h1(i,1,1,e) * &
               ( G11(i,1,1,c) * ur(i,1,1,e) &
               + G12(i,1,1,c) * us(i,1,1,e) &
               + G13(i,1,1,c) * ut(i,1,1,e))

          uus(i,1,1,e) = h1(i,1,1,e) * &
               ( G22(i,1,1,c) * us(i,1,1,e) &
               + G12(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * ut(i,1,1,e) )

          uut(i,1,1,e) = h1(i,1,1,e) * &
               ( G33(i,1,1,c) * ut(i,1,1,e) &
               + G13(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * us(i,1,1,e))
       end do
    end do

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dxt(i,kk) * uur(kk,jj,1,1)
          end do
          w(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dyt(j, kk)*uus(i,kk,k,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dzt(k, kk)*uut(i,j,kk,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + wt
             end do
          end do
       end do
    end do


  end subroutine sx_ax_helm_compr_lx11

  subroutine sx_ax_helm_compr_lx10(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, ind, n, m)
    integer, parameter :: lx = 10
    integer, intent(in) :: n, m
    integer, intent(in) :: ind(n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, m)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    integer :: e, c, i, j, k, jj, kk
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: uur(lx, lx, lx, n)
    real(kind=rp) :: uus(lx, lx, lx, n)
    real(kind=rp) :: uut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + Dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          c = ind(e)
          uur(i,1,1,e) = h1(i,1,1,e) * &
               ( G11(i,1,1,c) * ur(i,1,1,e) &
               + G12(i,1,1,c) * us(i,1,1,e) &
               + G13(i,1,1,c) * ut(i,1,1,e))

          uus(i,1,1,e) = h1(i,1,1,e) * &
               ( G22(i,1,1,c) * us(i,1,1,e) &
               + G12(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * ut(i,1,1,e) )

          uut(i,1,1,e) = h1(i,1,1,e) * &
               ( G33(i,1,1,c) * ut(i,1,1,e) &
               + G13(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * us(i,1,1,e))
       end do
    end do

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dxt(i,kk) * uur(kk,jj,1,1)
          end do
          w(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dyt(j, kk)*uus(i,kk,k,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dzt(k, kk)*uut(i,j,kk,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + wt
             end do
          end do
       end do
    end do


  end subroutine sx_ax_helm_compr_lx10

  subroutine sx_ax_helm_compr_lx9(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, ind, n, m)
    integer, parameter :: lx = 9
    integer, intent(in) :: n, m
    integer, intent(in) :: ind(n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, m)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    integer :: e, c, i, j, k, jj, kk
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: uur(lx, lx, lx, n)
    real(kind=rp) :: uus(lx, lx, lx, n)
    real(kind=rp) :: uut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + Dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          c = ind(e)
          uur(i,1,1,e) = h1(i,1,1,e) * &
               ( G11(i,1,1,c) * ur(i,1,1,e) &
               + G12(i,1,1,c) * us(i,1,1,e) &
               + G13(i,1,1,c) * ut(i,1,1,e))

          uus(i,1,1,e) = h1(i,1,1,e) * &
               ( G22(i,1,1,c) * us(i,1,1,e) &
               + G12(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * ut(i,1,1,e) )

          uut(i,1,1,e) = h1(i,1,1,e) * &
               ( G33(i,1,1,c) * ut(i,1,1,e) &
               + G13(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * us(i,1,1,e))
       end do
    end do

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dxt(i,kk) * uur(kk,jj,1,1)
          end do
          w(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dyt(j, kk)*uus(i,kk,k,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dzt(k, kk)*uut(i,j,kk,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + wt
             end do
          end do
       end do
    end do


  end subroutine sx_ax_helm_compr_lx9

  subroutine sx_ax_helm_compr_lx8(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, ind, n, m)
    integer, parameter :: lx = 8
    integer, intent(in) :: n, m
    integer, intent(in) :: ind(n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, m)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    integer :: e, c, i, j, k, jj, kk
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: uur(lx, lx, lx, n)
    real(kind=rp) :: uus(lx, lx, lx, n)
    real(kind=rp) :: uut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + Dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          c = ind(e)
          uur(i,1,1,e) = h1(i,1,1,e) * &
               ( G11(i,1,1,c) * ur(i,1,1,e) &
               + G12(i,1,1,c) * us(i,1,1,e) &
               + G13(i,1,1,c) * ut(i,1,1,e))

          uus(i,1,1,e) = h1(i,1,1,e) * &
               ( G22(i,1,1,c) * us(i,1,1,e) &
               + G12(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * ut(i,1,1,e) )

          uut(i,1,1,e) = h1(i,1,1,e) * &
               ( G33(i,1,1,c) * ut(i,1,1,e) &
               + G13(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * us(i,1,1,e))
       end do
    end do

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dxt(i,kk) * uur(kk,jj,1,1)
          end do
          w(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dyt(j, kk)*uus(i,kk,k,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dzt(k, kk)*uut(i,j,kk,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + wt
             end do
          end do
       end do
    end do


  end subroutine sx_ax_helm_compr_lx8

  subroutine sx_ax_helm_compr_lx7(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, ind, n, m)
    integer, parameter :: lx = 7
    integer, intent(in) :: n, m
    integer, intent(in) :: ind(n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, m)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    integer :: e, c, i, j, k, jj, kk
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: uur(lx, lx, lx, n)
    real(kind=rp) :: uus(lx, lx, lx, n)
    real(kind=rp) :: uut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + Dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          c = ind(e)
          uur(i,1,1,e) = h1(i,1,1,e) * &
               ( G11(i,1,1,c) * ur(i,1,1,e) &
               + G12(i,1,1,c) * us(i,1,1,e) &
               + G13(i,1,1,c) * ut(i,1,1,e))

          uus(i,1,1,e) = h1(i,1,1,e) * &
               ( G22(i,1,1,c) * us(i,1,1,e) &
               + G12(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * ut(i,1,1,e) )

          uut(i,1,1,e) = h1(i,1,1,e) * &
               ( G33(i,1,1,c) * ut(i,1,1,e) &
               + G13(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * us(i,1,1,e))
       end do
    end do

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dxt(i,kk) * uur(kk,jj,1,1)
          end do
          w(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dyt(j, kk)*uus(i,kk,k,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dzt(k, kk)*uut(i,j,kk,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + wt
             end do
          end do
       end do
    end do


  end subroutine sx_ax_helm_compr_lx7

  subroutine sx_ax_helm_compr_lx6(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, ind, n, m)
    integer, parameter :: lx = 6
    integer, intent(in) :: n, m
    integer, intent(in) :: ind(n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, m)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    integer :: e, c, i, j, k, jj, kk
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: uur(lx, lx, lx, n)
    real(kind=rp) :: uus(lx, lx, lx, n)
    real(kind=rp) :: uut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + Dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          c = ind(e)
          uur(i,1,1,e) = h1(i,1,1,e) * &
               ( G11(i,1,1,c) * ur(i,1,1,e) &
               + G12(i,1,1,c) * us(i,1,1,e) &
               + G13(i,1,1,c) * ut(i,1,1,e))

          uus(i,1,1,e) = h1(i,1,1,e) * &
               ( G22(i,1,1,c) * us(i,1,1,e) &
               + G12(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * ut(i,1,1,e) )

          uut(i,1,1,e) = h1(i,1,1,e) * &
               ( G33(i,1,1,c) * ut(i,1,1,e) &
               + G13(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * us(i,1,1,e))
       end do
    end do

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dxt(i,kk) * uur(kk,jj,1,1)
          end do
          w(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dyt(j, kk)*uus(i,kk,k,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dzt(k, kk)*uut(i,j,kk,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + wt
             end do
          end do
       end do
    end do


  end subroutine sx_ax_helm_compr_lx6

  subroutine sx_ax_helm_compr_lx5(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, ind, n, m)
    integer, parameter :: lx = 5
    integer, intent(in) :: n, m
    integer, intent(in) :: ind(n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, m)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    integer :: e, c, i, j, k, jj, kk
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: uur(lx, lx, lx, n)
    real(kind=rp) :: uus(lx, lx, lx, n)
    real(kind=rp) :: uut(lx, lx, lx, n)
    real(kind=rp) :: wr, ws, wt

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dx(i,kk)*u(kk,jj,1,1)
          end do
          ur(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dy(j,kk)*u(i,kk,k,e)
                end do
                us(i,j,k,e) = ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + Dz(k,kk)*u(i,j,kk,e)
                end do
                ut(i,j,k,e) = wt
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          c = ind(e)
          uur(i,1,1,e) = h1(i,1,1,e) * &
               ( G11(i,1,1,c) * ur(i,1,1,e) &
               + G12(i,1,1,c) * us(i,1,1,e) &
               + G13(i,1,1,c) * ut(i,1,1,e))

          uus(i,1,1,e) = h1(i,1,1,e) * &
               ( G22(i,1,1,c) * us(i,1,1,e) &
               + G12(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * ut(i,1,1,e) )

          uut(i,1,1,e) = h1(i,1,1,e) * &
               ( G33(i,1,1,c) * ut(i,1,1,e) &
               + G13(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * us(i,1,1,e))
       end do
    end do

    do i = 1, lx
       do jj = 1, lx * lx * n
          wr = 0d0
          do kk = 1, lx
             wr = wr + Dxt(i,kk) * uur(kk,jj,1,1)
          end do
          w(i,jj,1,1) = wr
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                ws = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   ws = ws + Dyt(j, kk)*uus(i,kk,k,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + ws
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                wt = 0d0
                !NEC$ unroll_completely
                do kk = 1, lx
                   wt = wt + dzt(k, kk)*uut(i,j,kk,e)
                end do
                w(i,j,k,e) = w(i,j,k,e) + wt
             end do
          end do
       end do
    end do


  end subroutine sx_ax_helm_compr_lx5

  subroutine sx_ax_helm_compr_lx4(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, ind, n, m)
    integer, parameter :: lx = 4
    integer, intent(in) :: n, m
    integer, intent(in) :: ind(n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, m)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    integer :: e, c, i, j, k, jj
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: uur(lx, lx, lx, n)
    real(kind=rp) :: uus(lx, lx, lx, n)
    real(kind=rp) :: uut(lx, lx, lx, n)

    do i = 1, lx
       do jj = 1, lx * lx * n
          ur(i,jj,1,1) = Dx(i,1)*u(1,jj,1,1) &
               + Dx(i,2)*u(2,jj,1,1) &
               + Dx(i,3)*u(3,jj,1,1) &
               + Dx(i,4)*u(4,jj,1,1)
       end do
    end do


    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                us(i,j,k,e) = Dy(j,1) * u(i,1,k,e) &
                     + Dy(j,2) * u(i,2,k,e) &
                     + Dy(j,3) * u(i,3,k,e) &
                     + Dy(j,4) * u(i,4,k,e)
             end do
          end do
       end do
    end do



    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                ut(i,j,k,e) = Dz(k,1) * u(i,j,1,e) &
                     + Dz(k,2) * u(i,j,2,e) &
                     + Dz(k,3) * u(i,j,3,e) &
                     + Dz(k,4) * u(i,j,4,e)
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          c = ind(e)
          uur(i,1,1,e) = h1(i,1,1,e) * &
               ( G11(i,1,1,c) * ur(i,1,1,e) &
               + G12(i,1,1,c) * us(i,1,1,e) &
               + G13(i,1,1,c) * ut(i,1,1,e))

          uus(i,1,1,e) = h1(i,1,1,e) * &
               ( G22(i,1,1,c) * us(i,1,1,e) &
               + G12(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * ut(i,1,1,e) )

          uut(i,1,1,e) = h1(i,1,1,e) * &
               ( G33(i,1,1,c) * ut(i,1,1,e) &
               + G13(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * us(i,1,1,e))
       end do
    end do

    do i = 1, lx
       do jj = 1, lx * lx * n
          w(i,jj,1,1) = Dxt(i,1) * uur(1,jj,1,1) &
               + Dxt(i,2) * uur(2,jj,1,1) &
               + Dxt(i,3) * uur(3,jj,1,1) &
               + Dxt(i,4) * uur(4,jj,1,1)
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                w(i,j,k,e) = w(i,j,k,e) + Dyt(j,1) * uus(i,1,k,e) &
                     + Dyt(j,2) * uus(i,2,k,e) &
                     + Dyt(j,3) * uus(i,3,k,e) &
                     + Dyt(j,4) * uus(i,4,k,e)
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                w(i,j,k,e) = w(i,j,k,e) + Dzt(k,1) * uut(i,j,1,e) &
                     + Dzt(k,2) * uut(i,j,2,e) &
                     + Dzt(k,3) * uut(i,j,3,e) &
                     + Dzt(k,4) * uut(i,j,4,e)
             end do
          end do
       end do
    end do

  end subroutine sx_ax_helm_compr_lx4

  subroutine sx_ax_helm_compr_lx3(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, ind, n, m)
    integer, parameter :: lx = 3
    integer, intent(in) :: n, m
    integer, intent(in) :: ind(n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, m)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    integer :: e, c, i, j, k, jj
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: uur(lx, lx, lx, n)
    real(kind=rp) :: uus(lx, lx, lx, n)
    real(kind=rp) :: uut(lx, lx, lx, n)

    do i = 1, lx
       do jj = 1, lx * lx * n
          ur(i,jj,1,1) = Dx(i,1)*u(1,jj,1,1) &
               + Dx(i,2)*u(2,jj,1,1) &
               + Dx(i,3)*u(3,jj,1,1)
       end do
    end do


    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                us(i,j,k,e) = Dy(j,1) * u(i,1,k,e) &
                     + Dy(j,2) * u(i,2,k,e) &
                     + Dy(j,3) * u(i,3,k,e)
             end do
          end do
       end do
    end do



    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                ut(i,j,k,e) = Dz(k,1) * u(i,j,1,e) &
                     + Dz(k,2) * u(i,j,2,e) &
                     + Dz(k,3) * u(i,j,3,e)
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          c = ind(e)
          uur(i,1,1,e) = h1(i,1,1,e) * &
               ( G11(i,1,1,c) * ur(i,1,1,e) &
               + G12(i,1,1,c) * us(i,1,1,e) &
               + G13(i,1,1,c) * ut(i,1,1,e))

          uus(i,1,1,e) = h1(i,1,1,e) * &
               ( G22(i,1,1,c) * us(i,1,1,e) &
               + G12(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * ut(i,1,1,e) )

          uut(i,1,1,e) = h1(i,1,1,e) * &
               ( G33(i,1,1,c) * ut(i,1,1,e) &
               + G13(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * us(i,1,1,e))
       end do
    end do

    do i = 1, lx
       do jj = 1, lx * lx * n
          w(i,jj,1,1) = Dxt(i,1) * uur(1,jj,1,1) &
               + Dxt(i,2) * uur(2,jj,1,1) &
               + Dxt(i,3) * uur(3,jj,1,1)
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                w(i,j,k,e) = w(i,j,k,e) + Dyt(j,1) * uus(i,1,k,e) &
                     + Dyt(j,2) * uus(i,2,k,e) &
                     + Dyt(j,3) * uus(i,3,k,e)
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                w(i,j,k,e) = w(i,j,k,e) + Dzt(k,1) * uut(i,j,1,e) &
                     + Dzt(k,2) * uut(i,j,2,e) &
                     + Dzt(k,3) * uut(i,j,3,e)
             end do
          end do
       end do
    end do

  end subroutine sx_ax_helm_compr_lx3

  subroutine sx_ax_helm_compr_lx2(w, u, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
       h1, G11, G22, G33, G12, G13, G23, ind, n, m)
    integer, parameter :: lx = 2
    integer, intent(in) :: n, m
    integer, intent(in) :: ind(n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(in) :: u(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, m)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, m)
    real(kind=rp), intent(in) :: Dx(lx, lx)
    real(kind=rp), intent(in) :: Dy(lx, lx)
    real(kind=rp), intent(in) :: Dz(lx, lx)
    real(kind=rp), intent(in) :: Dxt(lx, lx)
    real(kind=rp), intent(in) :: Dyt(lx, lx)
    real(kind=rp), intent(in) :: Dzt(lx, lx)
    integer :: e, c, i, j, k, jj
    real(kind=rp) :: ur(lx, lx, lx, n)
    real(kind=rp) :: us(lx, lx, lx, n)
    real(kind=rp) :: ut(lx, lx, lx, n)
    real(kind=rp) :: uur(lx, lx, lx, n)
    real(kind=rp) :: uus(lx, lx, lx, n)
    real(kind=rp) :: uut(lx, lx, lx, n)

    do i = 1, lx
       do jj = 1, lx * lx * n
          ur(i,jj,1,1) = Dx(i,1) * u(1,jj,1,1) &
               + Dx(i,2) * u(2,jj,1,1)
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                us(i,j,k,e) = Dy(j,1) * u(i,1,k,e) &
                     + Dy(j,2) * u(i,2,k,e)
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                ut(i,j,k,e) = Dz(k,1) * u(i,j,1,e) &
                     + Dz(k,2) * u(i,j,2,e)
             end do
          end do
       end do
    end do

    do i = 1, lx * lx * lx
       do e = 1, n
          c = ind(e)
          uur(i,1,1,e) = h1(i,1,1,e) * &
               ( G11(i,1,1,c) * ur(i,1,1,e) &
               + G12(i,1,1,c) * us(i,1,1,e) &
               + G13(i,1,1,c) * ut(i,1,1,e))

          uus(i,1,1,e) = h1(i,1,1,e) * &
               ( G22(i,1,1,c) * us(i,1,1,e) &
               + G12(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * ut(i,1,1,e) )

          uut(i,1,1,e) = h1(i,1,1,e) * &
               ( G33(i,1,1,c) * ut(i,1,1,e) &
               + G13(i,1,1,c) * ur(i,1,1,e) &
               + G23(i,1,1,c) * us(i,1,1,e))
       end do
    end do

    do i = 1, lx
       do jj = 1, lx * lx * n
          w(i,jj,1,1) = Dxt(i,1) * uur(1,jj,1,1) &
               + Dxt(i,2) * uur(2,jj,1,1)
       end do
    end do

    do k = 1, lx
       do i = 1, lx
          do j = 1, lx
             do e = 1, n
                w(i,j,k,e) = w(i,j,k,e) + Dyt(j,1) * uus(i,1,k,e) &
                     + Dyt(j,2) * uus(i,2,k,e)
             end do
          end do
       end do
    end do

    do j = 1, lx
       do i = 1, lx
          do k = 1, lx
             do e = 1, n
                w(i,j,k,e) = w(i,j,k,e) + Dzt(k,1) * uut(i,j,1,e) &
                     + Dzt(k,2) * uut(i,j,2,e)
             end do
          end do
       end do
    end do

  end subroutine sx_ax_helm_compr_lx2

end module ax_helm_compr_sx
