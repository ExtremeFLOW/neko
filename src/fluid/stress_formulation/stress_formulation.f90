! Copyright (c) 2023, The Neko Authors
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
module stress_formulation
  use num_types, only : rp
  use coefs, only : coef_t
  use mxm_wrapper, only : mxm
  use space, only : space_t
  use mesh, only : mesh_t
  use math, only : addcol4, add2
  implicit none
  private

  public :: ax_helm_stress_compute

contains

  subroutine ax_helm_stress_compute(au, av, aw, u, v, w, coef, msh, Xh)
    type(mesh_t), intent(inout) :: msh
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(inout) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: v(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: au(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: av(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: aw(Xh%lx, Xh%ly, Xh%lz, msh%nelv)

    call ax_helm_stress_lx(au, av, aw, u, v, w, Xh%dx, Xh%dy, Xh%dz, Xh%dxt, &
      Xh%dyt, Xh%dzt, coef%h1, coef%h2, coef%G11, coef%G22, coef%G33,&
      coef%G12, coef%G13, coef%G23, coef%jacinv, Xh%w3, msh%nelv, Xh%lx)

    if (coef%ifh2) then
       call addcol4 (au, coef%h2, coef%B, u, coef%dof%size())
       call addcol4 (av, coef%h2, coef%B, v, coef%dof%size())
       call addcol4 (aw, coef%h2, coef%B, w, coef%dof%size())
    end if

  end subroutine ax_helm_stress_compute

  subroutine ax_helm_stress_lx(au, av, aw, u, v, w, Dx, Dy, Dz, Dxt, Dyt, Dzt, &
    h1, h2, G11, G22, G33, G12, G13, G23, jackinv, weights3, n, lx)

    integer, intent(in) :: n, lx
    real(kind=rp), intent(inout) :: u(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: v(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: w(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: au(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: av(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: aw(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h1(lx, lx, lx, n)
    real(kind=rp), intent(in) :: h2(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G11(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G22(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G33(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G12(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G13(lx, lx, lx, n)
    real(kind=rp), intent(in) :: G23(lx, lx, lx, n)
    real(kind=rp), intent(inout) :: jackinv(lx, lx, lx, n)
    real(kind=rp), intent(in) :: Dx(lx,lx)
    real(kind=rp), intent(in) :: Dy(lx,lx)
    real(kind=rp), intent(in) :: Dz(lx,lx)
    real(kind=rp), intent(in) :: Dxt(lx,lx)
    real(kind=rp), intent(in) :: Dyt(lx,lx)
    real(kind=rp), intent(in) :: Dzt(lx,lx)
    real(kind=rp), intent(in) :: weights3(lx, lx, lx)

    real(kind=rp) :: ur(lx*lx*lx, 3, 3)
    integer :: e, p, i, ngll
    real(kind=rp) :: dj
    real(kind=rp) :: s11, s12, s13, s21, s22, s23, s31, s32, s33
    real(kind=rp) :: u1, u2, u3, v1, v2, v3, w1, w2, w3

    do e=1, n

       p = lx - 1
       ngll = lx*lx*lx

       ! NB: Passing the same derivative matrix everywhere, as in Nek5000.
       call local_grad(ur(1,1,1), ur(1,2,1), ur(1,3,1), u(:,:,:,e), p, 1, Dx, Dxt)
       call local_grad(ur(1,1,2), ur(1,2,2), ur(1,3,2), v(:,:,:,e), p, 1, Dx, Dxt)
       call local_grad(ur(1,1,3), ur(1,2,3), ur(1,3,3), w(:,:,:,e), p, 1, Dx, Dxt)

       do i=1, ngll

          u1 = ur(i,1,1)*G11(i,1,1,e) + ur(i,2,1)*G12(i,1,1,e) &
                                      + ur(i,3,1)*G13(i,1,1,e)
          u2 = ur(i,1,1)*G12(i,1,1,e) + ur(i,2,1)*G22(i,1,1,e) &
                                      + ur(i,3,1)*G23(i,1,1,e)
          u3 = ur(i,1,1)*G13(i,1,1,e) + ur(i,2,1)*G23(i,1,1,e) &
                                      + ur(i,3,1)*G33(i,1,1,e)

          v1 = ur(i,1,2)*G11(i,1,1,e) + ur(i,2,2)*G12(i,1,1,e) &
                                      + ur(i,3,2)*G13(i,1,1,e)
          v2 = ur(i,1,2)*G12(i,1,1,e) + ur(i,2,2)*G22(i,1,1,e) &
                                      + ur(i,3,2)*G23(i,1,1,e)
          v3 = ur(i,1,2)*G13(i,1,1,e) + ur(i,2,2)*G23(i,1,1,e) &
                                      + ur(i,3,2)*G33(i,1,1,e)

          w1 = ur(i,1,3)*G11(i,1,1,e) + ur(i,2,3)*G12(i,1,1,e) &
                                      + ur(i,3,3)*G13(i,1,1,e)
          w2 = ur(i,1,3)*G12(i,1,1,e) + ur(i,2,3)*G22(i,1,1,e) &
                                      + ur(i,3,3)*G23(i,1,1,e)
          w3 = ur(i,1,3)*G13(i,1,1,e) + ur(i,2,3)*G23(i,1,1,e) &
                                      + ur(i,3,3)*G33(i,1,1,e)

          dj  = h1(i,1,1,e) * weights3(i,1,1) * jackinv(i,1,1,e)
          s11 = dj*(u1 + u1)
          s12 = dj*(u2 + v1)
          s13 = dj*(u3 + w1)
          s21 = dj*(v1 + u2)
          s22 = dj*(v2 + v2)
          s23 = dj*(v3 + w2)
          s31 = dj*(w1 + u3)
          s32 = dj*(w2 + v3)
          s33 = dj*(w3 + w3)

          ur(i,1,1) = G11(i,1,1,e)*s11 + G12(i,1,1,e)*s12 + G13(i,1,1,e)*s13
          ur(i,2,1) = G12(i,1,1,e)*s11 + G22(i,1,1,e)*s12 + G23(i,1,1,e)*s13
          ur(i,3,1) = G13(i,1,1,e)*s11 + G23(i,1,1,e)*s12 + G33(i,1,1,e)*s13

          ur(i,1,2) = G11(i,1,1,e)*s21 + G12(i,1,1,e)*s22 + G13(i,1,1,e)*s23
          ur(i,2,2) = G12(i,1,1,e)*s21 + G22(i,1,1,e)*s22 + G23(i,1,1,e)*s23
          ur(i,3,2) = G13(i,1,1,e)*s21 + G23(i,1,1,e)*s22 + G33(i,1,1,e)*s23

          ur(i,1,3) = G11(i,1,1,e)*s31 + G12(i,1,1,e)*s32 + G13(i,1,1,e)*s33
          ur(i,2,3) = G12(i,1,1,e)*s31 + G22(i,1,1,e)*s32 + G23(i,1,1,e)*s33
          ur(i,3,3) = G13(i,1,1,e)*s31 + G23(i,1,1,e)*s32 + G33(i,1,1,e)*s33
       enddo

       call local_grad3_t(au(:,:,:,e), ur(1,1,1), ur(1,2,1), ur(1,3,1), p, 1, Dx, Dxt, av(:,:,:,e))
       call local_grad3_t(av(:,:,:,e), ur(1,1,2), ur(1,2,2), ur(1,3,2), p, 1, Dx, Dxt, ur)
       call local_grad3_t(aw(:,:,:,e), ur(1,1,3), ur(1,2,3), ur(1,3,3), p, 1, Dx, Dxt, ur)
    end do

  end subroutine ax_helm_stress_lx

  !> Compute gradient with respect to parametric coordinates.
  subroutine local_grad(ur, us, ut, u, N, e, D, Dt)
    integer, intent(in) :: N
    real(kind=rp), intent(inout) :: ur(0:N,0:N,0:N)
    real(kind=rp), intent(inout) :: us(0:N,0:N,0:N)
    real(kind=rp), intent(inout) :: ut(0:N,0:N,0:N)
    real(kind=rp), intent(inout) :: u(0:N,0:N,0:N,1)
    real(kind=rp), intent(in) :: D(0:N,0:N), Dt(0:N,0:N)
    integer, intent(in) ::  e
    integer :: m1, m2, k


    m1 = N+1
    m2 = m1*m1

    call mxm(D, m1, u(0,0,0,e), m1, ur, m2)
    do k=0,N
       call mxm(u(0,0,k,e), m1, Dt, m1, us(0,0,k), m1)
    end do
    call mxm(u(0,0,0,e), m2, Dt, m1, ut, m1)
  end subroutine local_grad

  subroutine local_grad3_t(u, ur, us, ut, N, e, D, Dt, w)
    integer, intent(in) :: N
    real(kind=rp), intent(inout) :: u(0:N, 0:N, 0:N, 1)
    real(kind=rp), intent(inout) :: ur(0:N, 0:N, 0:N)
    real(kind=rp), intent(inout) :: us(0:N, 0:N, 0:N)
    real(kind=rp), intent(inout) :: ut(0:N, 0:N, 0:N)
    real(kind=rp), intent(in) :: D(0:N, 0:N), Dt(0:N, 0:N)
    real(kind=rp), intent(inout) :: w(0:N, 0:N, 0:N)
    integer, intent(in) :: e
    integer :: m1, m2, m3, k

    m1 = N+1
    m2 = m1*m1
    m3 = m1*m1*m1

    call mxm(Dt, m1, ur, m1, u(0,0,0,e), m2)

    do k=0, N
       call mxm(us(0,0,k), m1, D, m1, w(0,0,k), m1)
    end do
    call add2(u(0,0,0,e), w, m3)

    call mxm(ut, m2, D, m1, w, m1)
    call add2(u(0,0,0,e), w, m3)
  end

end module