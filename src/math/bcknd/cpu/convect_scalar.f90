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
!>
submodule (opr_cpu) cpu_convect_scalar
  use num_types, only : rp
  use math, only : col2
  implicit none

contains

  module subroutine opr_cpu_convect_scalar(du, u, c, Xh_GLL, Xh_GL, coef_GLL, &
                                    coef_GL, GLL_to_GL)
    type(space_t), intent(in) :: Xh_GL
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(coef_t), intent(in) :: coef_GL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: &
                   du(Xh_GLL%lx, Xh_GLL%ly, Xh_GLL%lz, coef_GL%msh%nelv)
    real(kind=rp), intent(inout) :: &
                   u(Xh_GL%lx, Xh_GL%lx, Xh_GL%lx, coef_GL%msh%nelv)
    real(kind=rp), intent(inout) :: c(Xh_GL%lxyz, coef_GL%msh%nelv, 3)
    associate(dx => Xh_GL%dx, dy => Xh_GL%dy, dz => Xh_GL%dz, &
         lx => Xh_GL%lx, nelv => coef_GL%msh%nelv)

      select case (lx)
      case (18)
         call cpu_convect_scalar_lx18(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case (17)
         call cpu_convect_scalar_lx17(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case (16)
         call cpu_convect_scalar_lx16(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case (15)
         call cpu_convect_scalar_lx15(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case (14)
         call cpu_convect_scalar_lx14(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case (13)
         call cpu_convect_scalar_lx13(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case (12)
         call cpu_convect_scalar_lx12(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case (11)
         call cpu_convect_scalar_lx11(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case (10)
         call cpu_convect_scalar_lx10(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case (9)
         call cpu_convect_scalar_lx9(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case (8)
         call cpu_convect_scalar_lx8(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case (7)
         call cpu_convect_scalar_lx7(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case (6)
         call cpu_convect_scalar_lx6(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case (5)
         call cpu_convect_scalar_lx5(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case (4)
         call cpu_convect_scalar_lx4(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case (3)
         call cpu_convect_scalar_lx3(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case (2)
         call cpu_convect_scalar_lx2(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv)
      case default
         call cpu_convect_scalar_lx(du, u, c, dx, dy, dz, &
              Xh_GLL, coef_GLL, GLL_to_GL, nelv, lx)
      end select
    end associate

  end subroutine opr_cpu_convect_scalar

  subroutine cpu_convect_scalar_lx(du, u, c, dx, dy, dz, Xh_GLL, &
                                coef_GLL, GLL_to_GL, nelv, lx)
    integer, intent(in) :: nelv, lx
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: du(Xh_GLL%lx, Xh_GLL%lx, Xh_GLL%lx, nelv)
    real(kind=rp), intent(in) :: u(lx, lx, lx, nelv)
    real(kind=rp), intent(in) :: c(lx*lx*lx, nelv, 3)
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) :: ur, us, ut
    real(kind=rp) :: ud(lx*lx*lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l, idx, n_GLL

    n_GLL = nelv * Xh_GLL%lxyz

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             ur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                us(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             ut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) &
                + c(i,e,3) * ut(i,1,1)
       end do
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_convect_scalar_lx

  subroutine cpu_convect_scalar_lx18(du, u, c, dx, dy, dz, Xh_GLL, &
                                  coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 18
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: du(Xh_GLL%lx, Xh_GLL%lx, Xh_GLL%lx, nelv)
    real(kind=rp), intent(in) :: u(lx, lx, lx, nelv)
    real(kind=rp), intent(in) :: c(lx*lx*lx, nelv, 3)
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) :: ur, us, ut
    real(kind=rp) :: ud(lx*lx*lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l, idx, n_GLL

    n_GLL = nelv * Xh_GLL%lxyz

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             ur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                us(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             ut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) &
                + c(i,e,3) * ut(i,1,1)
       end do
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_convect_scalar_lx18

  subroutine cpu_convect_scalar_lx17(du, u, c, dx, dy, dz, Xh_GLL, &
                                  coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 17
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: du(Xh_GLL%lx, Xh_GLL%lx, Xh_GLL%lx, nelv)
    real(kind=rp), intent(in) :: u(lx, lx, lx, nelv)
    real(kind=rp), intent(in) :: c(lx*lx*lx, nelv, 3)
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) :: ur, us, ut
    real(kind=rp) :: ud(lx*lx*lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l, idx, n_GLL

    n_GLL = nelv * Xh_GLL%lxyz

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             ur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                us(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             ut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) &
                + c(i,e,3) * ut(i,1,1)
       end do
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_convect_scalar_lx17

  subroutine cpu_convect_scalar_lx16(du, u, c, dx, dy, dz, Xh_GLL, &
                                  coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 16
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: du(Xh_GLL%lx, Xh_GLL%lx, Xh_GLL%lx, nelv)
    real(kind=rp), intent(in) :: u(lx, lx, lx, nelv)
    real(kind=rp), intent(in) :: c(lx*lx*lx, nelv, 3)
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) :: ur, us, ut
    real(kind=rp) :: ud(lx*lx*lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l, idx, n_GLL

    n_GLL = nelv * Xh_GLL%lxyz

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             ur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                us(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             ut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) &
                + c(i,e,3) * ut(i,1,1)
       end do
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_convect_scalar_lx16

  subroutine cpu_convect_scalar_lx15(du, u, c, dx, dy, dz, Xh_GLL, &
                                  coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 15
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: du(Xh_GLL%lx, Xh_GLL%lx, Xh_GLL%lx, nelv)
    real(kind=rp), intent(in) :: u(lx, lx, lx, nelv)
    real(kind=rp), intent(in) :: c(lx*lx*lx, nelv, 3)
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) :: ur, us, ut
    real(kind=rp) :: ud(lx*lx*lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l, idx, n_GLL

    n_GLL = nelv * Xh_GLL%lxyz

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             ur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                us(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             ut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) &
                + c(i,e,3) * ut(i,1,1)
       end do
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_convect_scalar_lx15

  subroutine cpu_convect_scalar_lx14(du, u, c, dx, dy, dz, Xh_GLL, &
                                  coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 14
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: du(Xh_GLL%lx, Xh_GLL%lx, Xh_GLL%lx, nelv)
    real(kind=rp), intent(in) :: u(lx, lx, lx, nelv)
    real(kind=rp), intent(in) :: c(lx*lx*lx, nelv, 3)
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) :: ur, us, ut
    real(kind=rp) :: ud(lx*lx*lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l, idx, n_GLL

    n_GLL = nelv * Xh_GLL%lxyz

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             ur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                us(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             ut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) &
                + c(i,e,3) * ut(i,1,1)
       end do
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_convect_scalar_lx14

  subroutine cpu_convect_scalar_lx13(du, u, c, dx, dy, dz, Xh_GLL, &
                                  coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 13
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: du(Xh_GLL%lx, Xh_GLL%lx, Xh_GLL%lx, nelv)
    real(kind=rp), intent(in) :: u(lx, lx, lx, nelv)
    real(kind=rp), intent(in) :: c(lx*lx*lx, nelv, 3)
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) :: ur, us, ut
    real(kind=rp) :: ud(lx*lx*lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l, idx, n_GLL

    n_GLL = nelv * Xh_GLL%lxyz

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             ur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                us(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             ut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) &
                + c(i,e,3) * ut(i,1,1)
       end do
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_convect_scalar_lx13

  subroutine cpu_convect_scalar_lx12(du, u, c, dx, dy, dz, Xh_GLL, &
                                  coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 12
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: du(Xh_GLL%lx, Xh_GLL%lx, Xh_GLL%lx, nelv)
    real(kind=rp), intent(in) :: u(lx, lx, lx, nelv)
    real(kind=rp), intent(in) :: c(lx*lx*lx, nelv, 3)
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) :: ur, us, ut
    real(kind=rp) :: ud(lx*lx*lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l, idx, n_GLL

    n_GLL = nelv * Xh_GLL%lxyz

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             ur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                us(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             ut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) &
                + c(i,e,3) * ut(i,1,1)
       end do
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_convect_scalar_lx12

  subroutine cpu_convect_scalar_lx11(du, u, c, dx, dy, dz, Xh_GLL, &
                                  coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 11
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: du(Xh_GLL%lx, Xh_GLL%lx, Xh_GLL%lx, nelv)
    real(kind=rp), intent(in) :: u(lx, lx, lx, nelv)
    real(kind=rp), intent(in) :: c(lx*lx*lx, nelv, 3)
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) :: ur, us, ut
    real(kind=rp) :: ud(lx*lx*lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l, idx, n_GLL

    n_GLL = nelv * Xh_GLL%lxyz

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             ur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                us(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             ut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) &
                + c(i,e,3) * ut(i,1,1)
       end do
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_convect_scalar_lx11

  subroutine cpu_convect_scalar_lx10(du, u, c, dx, dy, dz, Xh_GLL, &
                                  coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 10
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: du(Xh_GLL%lx, Xh_GLL%lx, Xh_GLL%lx, nelv)
    real(kind=rp), intent(in) :: u(lx, lx, lx, nelv)
    real(kind=rp), intent(in) :: c(lx*lx*lx, nelv, 3)
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) :: ur, us, ut
    real(kind=rp) :: ud(lx*lx*lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l, idx, n_GLL

    n_GLL = nelv * Xh_GLL%lxyz

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             ur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                us(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             ut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) &
                + c(i,e,3) * ut(i,1,1)
       end do
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_convect_scalar_lx10

  subroutine cpu_convect_scalar_lx9(du, u, c, dx, dy, dz, Xh_GLL, &
                                 coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 9
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: du(Xh_GLL%lx, Xh_GLL%lx, Xh_GLL%lx, nelv)
    real(kind=rp), intent(in) :: u(lx, lx, lx, nelv)
    real(kind=rp), intent(in) :: c(lx*lx*lx, nelv, 3)
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) :: ur, us, ut
    real(kind=rp) :: ud(lx*lx*lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l, idx, n_GLL

    n_GLL = nelv * Xh_GLL%lxyz

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             ur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                us(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             ut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) &
                + c(i,e,3) * ut(i,1,1)
       end do
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_convect_scalar_lx9

  subroutine cpu_convect_scalar_lx8(du, u, c, dx, dy, dz, Xh_GLL, &
                                 coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 8
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: du(Xh_GLL%lx, Xh_GLL%lx, Xh_GLL%lx, nelv)
    real(kind=rp), intent(in) :: u(lx, lx, lx, nelv)
    real(kind=rp), intent(in) :: c(lx*lx*lx, nelv, 3)
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) :: ur, us, ut
    real(kind=rp) :: ud(lx*lx*lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l, idx, n_GLL

    n_GLL = nelv * Xh_GLL%lxyz

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             ur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                us(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             ut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) &
                + c(i,e,3) * ut(i,1,1)
       end do
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_convect_scalar_lx8

  subroutine cpu_convect_scalar_lx7(du, u, c, dx, dy, dz, Xh_GLL, &
                                 coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 7
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: du(Xh_GLL%lx, Xh_GLL%lx, Xh_GLL%lx, nelv)
    real(kind=rp), intent(in) :: u(lx, lx, lx, nelv)
    real(kind=rp), intent(in) :: c(lx*lx*lx, nelv, 3)
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) :: ur, us, ut
    real(kind=rp) :: ud(lx*lx*lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l, idx, n_GLL

    n_GLL = nelv * Xh_GLL%lxyz

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             ur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                us(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             ut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) &
                + c(i,e,3) * ut(i,1,1)
       end do
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_convect_scalar_lx7

  subroutine cpu_convect_scalar_lx6(du, u, c, dx, dy, dz, Xh_GLL, &
                                 coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 6
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: du(Xh_GLL%lx, Xh_GLL%lx, Xh_GLL%lx, nelv)
    real(kind=rp), intent(in) :: u(lx, lx, lx, nelv)
    real(kind=rp), intent(in) :: c(lx*lx*lx, nelv, 3)
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) :: ur, us, ut
    real(kind=rp) :: ud(lx*lx*lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l, idx, n_GLL

    n_GLL = nelv * Xh_GLL%lxyz

    do e = 1, nelv
      do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             ur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                us(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             ut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) &
                + c(i,e,3) * ut(i,1,1)
       end do
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_convect_scalar_lx6

  subroutine cpu_convect_scalar_lx5(du, u, c, dx, dy, dz, Xh_GLL, &
                                 coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 5
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: du(Xh_GLL%lx, Xh_GLL%lx, Xh_GLL%lx, nelv)
    real(kind=rp), intent(in) :: u(lx, lx, lx, nelv)
    real(kind=rp), intent(in) :: c(lx*lx*lx, nelv, 3)
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) :: ur, us, ut
    real(kind=rp) :: ud(lx*lx*lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l, idx, n_GLL

    n_GLL = nelv * Xh_GLL%lxyz

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             ur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                us(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             ut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) &
                + c(i,e,3) * ut(i,1,1)
       end do
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_convect_scalar_lx5

  subroutine cpu_convect_scalar_lx4(du, u, c, dx, dy, dz, Xh_GLL, &
                                 coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 4
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: du(Xh_GLL%lx, Xh_GLL%lx, Xh_GLL%lx, nelv)
    real(kind=rp), intent(in) :: u(lx, lx, lx, nelv)
    real(kind=rp), intent(in) :: c(lx*lx*lx, nelv, 3)
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) :: ur, us, ut
    real(kind=rp) :: ud(lx*lx*lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l, idx, n_GLL

    n_GLL = nelv * Xh_GLL%lxyz

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             ur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                us(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             ut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) &
                + c(i,e,3) * ut(i,1,1)
       end do
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_convect_scalar_lx4

  subroutine cpu_convect_scalar_lx3(du, u, c, dx, dy, dz, Xh_GLL, &
                                 coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 3
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: du(Xh_GLL%lx, Xh_GLL%lx, Xh_GLL%lx, nelv)
    real(kind=rp), intent(in) :: u(lx, lx, lx, nelv)
    real(kind=rp), intent(in) :: c(lx*lx*lx, nelv, 3)
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) :: ur, us, ut
    real(kind=rp) :: ud(lx*lx*lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l, idx, n_GLL

    n_GLL = nelv * Xh_GLL%lxyz

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             ur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                us(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             ut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) &
                + c(i,e,3) * ut(i,1,1)
       end do
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_convect_scalar_lx3

  subroutine cpu_convect_scalar_lx2(du, u, c, dx, dy, dz, Xh_GLL, &
                                 coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 2
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: du(Xh_GLL%lx, Xh_GLL%lx, Xh_GLL%lx, nelv)
    real(kind=rp), intent(in) :: u(lx, lx, lx, nelv)
    real(kind=rp), intent(in) :: c(lx*lx*lx, nelv, 3)
    real(kind=rp), dimension(lx, lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx, lx, lx) :: ur, us, ut
    real(kind=rp) :: ud(lx*lx*lx)
    real(kind=rp) :: tmp
    integer :: e, i, j, k, l, idx, n_GLL

    n_GLL = nelv * Xh_GLL%lxyz

    do e = 1, nelv
       do j = 1, lx * lx
          do i = 1, lx
             tmp = 0.0_rp
             do k = 1, lx
                tmp = tmp + dx(i,k) * u(k,j,1,e)
             end do
             ur(i,j,1) = tmp
          end do
       end do

       do k = 1, lx
          do j = 1, lx
             do i = 1, lx
                tmp = 0.0_rp
                do l = 1, lx
                   tmp = tmp + dy(j,l) * u(i,l,k,e)
                end do
                us(i,j,k) = tmp
             end do
          end do
       end do

       do k = 1, lx
          do i = 1, lx*lx
             tmp = 0.0_rp
             do l = 1, lx
                tmp = tmp + dz(k,l) * u(i,1,l,e)
             end do
             ut(i,1,k) = tmp
          end do
       end do

       do i = 1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) &
                + c(i,e,3) * ut(i,1,1)
       end do
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_convect_scalar_lx2

end submodule cpu_convect_scalar
