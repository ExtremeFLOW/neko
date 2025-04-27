! Copyright (c) 2021-2025, The Neko Authors
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
module ax_poisson
  use ax_product
  use utils, only : neko_error
  use num_types, only : rp
  use coefs, only : coef_t
  use space, only : space_t
  use mesh, only : mesh_t
  use math, only : addcol4
  implicit none
  private

  type, public, extends(ax_t) :: ax_poisson_t
   contains
     procedure, nopass :: compute => ax_poisson_compute
     procedure, pass(this) :: compute_vector => ax_poisson_compute_vector
  end type ax_poisson_t

contains

  subroutine ax_poisson_compute(w, u, coef, msh, Xh)
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp) :: ur(Xh%lx, Xh%lx, Xh%lx)
    real(kind=rp) :: us(Xh%lx, Xh%lx, Xh%lx)
    real(kind=rp) :: ut(Xh%lx, Xh%lx, Xh%lx)
    real(kind=rp) :: wur(Xh%lx, Xh%lx, Xh%lx)
    real(kind=rp) :: wus(Xh%lx, Xh%lx, Xh%lx)
    real(kind=rp) :: wut(Xh%lx, Xh%lx, Xh%lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l

    associate( D => Xh%dx, Dt => Xh%dxt, &
         G11 => coef%G11, G22 => coef%G22, G33 => coef%G33, &
         G12 => coef%G12, G13 => coef%G13, G23 => coef%G23, &
         n => msh%nelv, lx => Xh%lx)

      do e = 1, n
         do j = 1, lx * lx
            do i = 1, lx
               tmp = 0.0_rp
               do k = 1, lx
                  tmp = tmp + D(i,k) * u(k,j,1,e)
               end do
               wur(i,j,1) = tmp
            end do
         end do

         do k = 1, lx
            do j = 1, lx
               do i = 1, lx
                  tmp = 0.0_rp
                  do l = 1, lx
                     tmp = tmp + D(j,l) * u(i,l,k,e)
                  end do
                  wus(i,j,k) = tmp
               end do
            end do
         end do

         do k = 1, lx
            do i = 1, lx*lx
               tmp = 0.0_rp
               do l = 1, lx
                  tmp = tmp + D(k,l) * u(i,1,l,e)
               end do
               wut(i,1,k) = tmp
            end do
         end do

         do i = 1, lx*lx*lx
            ur(i,1,1) = ( G11(i,1,1,e) * wur(i,1,1) &
                        + G12(i,1,1,e) * wus(i,1,1) &
                        + G13(i,1,1,e) * wut(i,1,1) )
            us(i,1,1) = ( G12(i,1,1,e) * wur(i,1,1) &
                        + G22(i,1,1,e) * wus(i,1,1) &
                        + G23(i,1,1,e) * wut(i,1,1) )
            ut(i,1,1) = ( G13(i,1,1,e) * wur(i,1,1) &
                        + G23(i,1,1,e) * wus(i,1,1) &
                        + G33(i,1,1,e) * wut(i,1,1) )
         end do

         do j = 1, lx*lx
            do i = 1, lx
               tmp = 0.0_rp
               do k = 1, lx
                  tmp = tmp + Dt(i,k) * ur(k,j,1)
               end do
               w(i,j,1,e) = tmp
            end do
         end do

         do k = 1, lx
            do j = 1, lx
               do i = 1, lx
                  tmp = 0.0_rp
                  do l = 1, lx
                     tmp = tmp + Dt(j,l) * us(i,l,k)
                  end do
                  w(i,j,k,e) = w(i,j,k,e) + tmp
               end do
            end do
         end do

         do k = 1, lx
            do i = 1, lx*lx
               tmp = 0.0_rp
               do l = 1, lx
                  tmp = tmp + Dt(k,l) * ut(i,1,l)
               end do
               w(i,1,k,e) = w(i,1,k,e) + tmp
            end do
         end do

      end do
    end associate
  end subroutine ax_poisson_compute

  subroutine ax_poisson_compute_vector(this, au, av, aw, u, v, w, coef, msh, Xh)
    class(ax_poisson_t), intent(in) :: this
    type(space_t), intent(in) :: Xh
    type(mesh_t), intent(in) :: msh
    type(coef_t), intent(in) :: coef
    real(kind=rp), intent(inout) :: au(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: av(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: aw(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: v(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)

    call neko_error('Not in use')

  end subroutine ax_poisson_compute_vector

end module ax_poisson
