! Copyright (c) 2021, The Neko Authors
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
module cpu_conv_fst_3d
  use num_types, only : rp
  use coefs, only : coef_t
  use space, only : space_t
  use interpolation
  use math
  use gather_scatter 
  use mxm_wrapper  
  implicit none

contains

  subroutine cpu_conv_fst_3d_lx(du, u, c, dx, dy, dz, &
        Xh_GLL, coef_GLL, GLL_to_GL, nelv, lx)
    integer, intent(in) :: nelv, lx
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), dimension(Xh_GLL%lx,Xh_GLL%lx,Xh_GLL%lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u
    real(kind=rp), dimension(lx*lx*lx,nelv,3), intent(in) :: c
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: ur, us, ut
    real(kind=rp):: ud(lx*lx*lx)
    integer :: e, i, k, m, idx, n_GLL
    n_GLL = nelv * Xh_GLL%lxyz
    m = lx * lx
    do e = 1, nelv
       call mxm(dx,lx,u(1,1,1,e),lx,ur,m)
       do k=1,lx
          call mxm(u(1,1,k,e),lx,dy,lx,us(1,1,k),lx)
       enddo
       call mxm(u(1,1,1,e),m,dz,lx,ut,lx)
       do i=1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) + c(i,e,3) * ut(i,1,1)
       enddo
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_conv_fst_3d_lx

  subroutine cpu_conv_fst_3d_lx18(du, u, c, dx, dy, dz, &
        Xh_GLL, coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 18
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), dimension(Xh_GLL%lx,Xh_GLL%lx,Xh_GLL%lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u
    real(kind=rp), dimension(lx*lx*lx,nelv,3), intent(in) :: c
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: ur, us, ut
    real(kind=rp):: ud(lx*lx*lx)
    integer :: e, i, k, m, idx, n_GLL
    n_GLL = nelv * Xh_GLL%lxyz
    m = lx * lx
    do e = 1, nelv
       call mxm(dx,lx,u(1,1,1,e),lx,ur,m)
       do k=1,lx
          call mxm(u(1,1,k,e),lx,dy,lx,us(1,1,k),lx)
       enddo
       call mxm(u(1,1,1,e),m,dz,lx,ut,lx)
       do i=1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) + c(i,e,3) * ut(i,1,1)
       enddo
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_conv_fst_3d_lx18

  subroutine cpu_conv_fst_3d_lx17(du, u, c, dx, dy, dz, &
        Xh_GLL, coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 17
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), dimension(Xh_GLL%lx,Xh_GLL%lx,Xh_GLL%lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u
    real(kind=rp), dimension(lx*lx*lx,nelv,3), intent(in) :: c
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: ur, us, ut
    real(kind=rp):: ud(lx*lx*lx)
    integer :: e, i, k, m, idx, n_GLL
    n_GLL = nelv * Xh_GLL%lxyz
    m = lx * lx
    do e = 1, nelv
       call mxm(dx,lx,u(1,1,1,e),lx,ur,m)
       do k=1,lx
          call mxm(u(1,1,k,e),lx,dy,lx,us(1,1,k),lx)
       enddo
       call mxm(u(1,1,1,e),m,dz,lx,ut,lx)
       do i=1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) + c(i,e,3) * ut(i,1,1)
       enddo
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_conv_fst_3d_lx17

  subroutine cpu_conv_fst_3d_lx16(du, u, c, dx, dy, dz, &
        Xh_GLL, coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 16
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), dimension(Xh_GLL%lx,Xh_GLL%lx,Xh_GLL%lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u
    real(kind=rp), dimension(lx*lx*lx,nelv,3), intent(in) :: c
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: ur, us, ut
    real(kind=rp):: ud(lx*lx*lx)
    integer :: e, i, k, m, idx, n_GLL
    n_GLL = nelv * Xh_GLL%lxyz
    m = lx * lx
    do e = 1, nelv
       call mxm(dx,lx,u(1,1,1,e),lx,ur,m)
       do k=1,lx
          call mxm(u(1,1,k,e),lx,dy,lx,us(1,1,k),lx)
       enddo
       call mxm(u(1,1,1,e),m,dz,lx,ut,lx)
       do i=1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) + c(i,e,3) * ut(i,1,1)
       enddo
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_conv_fst_3d_lx16

  subroutine cpu_conv_fst_3d_lx15(du, u, c, dx, dy, dz, &
        Xh_GLL, coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 15
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), dimension(Xh_GLL%lx,Xh_GLL%lx,Xh_GLL%lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u
    real(kind=rp), dimension(lx*lx*lx,nelv,3), intent(in) :: c
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: ur, us, ut
    real(kind=rp):: ud(lx*lx*lx)
    integer :: e, i, k, m, idx, n_GLL
    n_GLL = nelv * Xh_GLL%lxyz
    m = lx * lx
    do e = 1, nelv
       call mxm(dx,lx,u(1,1,1,e),lx,ur,m)
       do k=1,lx
          call mxm(u(1,1,k,e),lx,dy,lx,us(1,1,k),lx)
       enddo
       call mxm(u(1,1,1,e),m,dz,lx,ut,lx)
       do i=1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) + c(i,e,3) * ut(i,1,1)
       enddo
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_conv_fst_3d_lx15

  subroutine cpu_conv_fst_3d_lx14(du, u, c, dx, dy, dz, &
        Xh_GLL, coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 14
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), dimension(Xh_GLL%lx,Xh_GLL%lx,Xh_GLL%lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u
    real(kind=rp), dimension(lx*lx*lx,nelv,3), intent(in) :: c
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: ur, us, ut
    real(kind=rp):: ud(lx*lx*lx)
    integer :: e, i, k, m, idx, n_GLL
    n_GLL = nelv * Xh_GLL%lxyz
    m = lx * lx
    do e = 1, nelv
       call mxm(dx,lx,u(1,1,1,e),lx,ur,m)
       do k=1,lx
          call mxm(u(1,1,k,e),lx,dy,lx,us(1,1,k),lx)
       enddo
       call mxm(u(1,1,1,e),m,dz,lx,ut,lx)
       do i=1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) + c(i,e,3) * ut(i,1,1)
       enddo
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_conv_fst_3d_lx14

  subroutine cpu_conv_fst_3d_lx13(du, u, c, dx, dy, dz, &
        Xh_GLL, coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 13
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), dimension(Xh_GLL%lx,Xh_GLL%lx,Xh_GLL%lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u
    real(kind=rp), dimension(lx*lx*lx,nelv,3), intent(in) :: c
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: ur, us, ut
    real(kind=rp):: ud(lx*lx*lx)
    integer :: e, i, k, m, idx, n_GLL
    n_GLL = nelv * Xh_GLL%lxyz
    m = lx * lx
    do e = 1, nelv
       call mxm(dx,lx,u(1,1,1,e),lx,ur,m)
       do k=1,lx
          call mxm(u(1,1,k,e),lx,dy,lx,us(1,1,k),lx)
       enddo
       call mxm(u(1,1,1,e),m,dz,lx,ut,lx)
       do i=1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) + c(i,e,3) * ut(i,1,1)
       enddo
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_conv_fst_3d_lx13

  subroutine cpu_conv_fst_3d_lx12(du, u, c, dx, dy, dz, &
        Xh_GLL, coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 12
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), dimension(Xh_GLL%lx,Xh_GLL%lx,Xh_GLL%lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u
    real(kind=rp), dimension(lx*lx*lx,nelv,3), intent(in) :: c
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: ur, us, ut
    real(kind=rp):: ud(lx*lx*lx)
    integer :: e, i, k, m, idx, n_GLL
    n_GLL = nelv * Xh_GLL%lxyz
    m = lx * lx
    do e = 1, nelv
       call mxm(dx,lx,u(1,1,1,e),lx,ur,m)
       do k=1,lx
          call mxm(u(1,1,k,e),lx,dy,lx,us(1,1,k),lx)
       enddo
       call mxm(u(1,1,1,e),m,dz,lx,ut,lx)
       do i=1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) + c(i,e,3) * ut(i,1,1)
       enddo
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)
    
  end subroutine cpu_conv_fst_3d_lx12

  subroutine cpu_conv_fst_3d_lx11(du, u, c, dx, dy, dz, &
        Xh_GLL, coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 11
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), dimension(Xh_GLL%lx,Xh_GLL%lx,Xh_GLL%lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u
    real(kind=rp), dimension(lx*lx*lx,nelv,3), intent(in) :: c
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: ur, us, ut
    real(kind=rp):: ud(lx*lx*lx)
    integer :: e, i, k, m, idx, n_GLL
    n_GLL = nelv * Xh_GLL%lxyz
    m = lx * lx
    do e = 1, nelv
       call mxm(dx,lx,u(1,1,1,e),lx,ur,m)
       do k=1,lx
          call mxm(u(1,1,k,e),lx,dy,lx,us(1,1,k),lx)
       enddo
       call mxm(u(1,1,1,e),m,dz,lx,ut,lx)
       do i=1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) + c(i,e,3) * ut(i,1,1)
       enddo
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_conv_fst_3d_lx11

  subroutine cpu_conv_fst_3d_lx10(du, u, c, dx, dy, dz, &
        Xh_GLL, coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 10
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), dimension(Xh_GLL%lx,Xh_GLL%lx,Xh_GLL%lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u
    real(kind=rp), dimension(lx*lx*lx,nelv,3), intent(in) :: c
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: ur, us, ut
    real(kind=rp):: ud(lx*lx*lx)
    integer :: e, i, k, m, idx, n_GLL
    n_GLL = nelv * Xh_GLL%lxyz
    m = lx * lx
    do e = 1, nelv
       call mxm(dx,lx,u(1,1,1,e),lx,ur,m)
       do k=1,lx
          call mxm(u(1,1,k,e),lx,dy,lx,us(1,1,k),lx)
       enddo
       call mxm(u(1,1,1,e),m,dz,lx,ut,lx)
       do i=1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) + c(i,e,3) * ut(i,1,1)
       enddo
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_conv_fst_3d_lx10

  subroutine cpu_conv_fst_3d_lx9(du, u, c, dx, dy, dz, &
        Xh_GLL, coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 9
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), dimension(Xh_GLL%lx,Xh_GLL%lx,Xh_GLL%lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u
    real(kind=rp), dimension(lx*lx*lx,nelv,3), intent(in) :: c
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: ur, us, ut
    real(kind=rp):: ud(lx*lx*lx)
    integer :: e, i, k, m, idx, n_GLL
    n_GLL = nelv * Xh_GLL%lxyz
    m = lx * lx
    do e = 1, nelv
       call mxm(dx,lx,u(1,1,1,e),lx,ur,m)
       do k=1,lx
          call mxm(u(1,1,k,e),lx,dy,lx,us(1,1,k),lx)
       enddo
       call mxm(u(1,1,1,e),m,dz,lx,ut,lx)
       do i=1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) + c(i,e,3) * ut(i,1,1)
       enddo
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_conv_fst_3d_lx9

  subroutine cpu_conv_fst_3d_lx8(du, u, c, dx, dy, dz, &
        Xh_GLL, coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 8
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), dimension(Xh_GLL%lx,Xh_GLL%lx,Xh_GLL%lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u
    real(kind=rp), dimension(lx*lx*lx,nelv,3), intent(in) :: c
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: ur, us, ut
    real(kind=rp):: ud(lx*lx*lx)
    integer :: e, i, k, m, idx, n_GLL
    n_GLL = nelv * Xh_GLL%lxyz
    m = lx * lx
    do e = 1, nelv
       call mxm(dx,lx,u(1,1,1,e),lx,ur,m)
       do k=1,lx
          call mxm(u(1,1,k,e),lx,dy,lx,us(1,1,k),lx)
       enddo
       call mxm(u(1,1,1,e),m,dz,lx,ut,lx)
       do i=1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) + c(i,e,3) * ut(i,1,1)
       enddo
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_conv_fst_3d_lx8

  subroutine cpu_conv_fst_3d_lx7(du, u, c, dx, dy, dz, &
        Xh_GLL, coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 7
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), dimension(Xh_GLL%lx,Xh_GLL%lx,Xh_GLL%lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u
    real(kind=rp), dimension(lx*lx*lx,nelv,3), intent(in) :: c
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: ur, us, ut
    real(kind=rp):: ud(lx*lx*lx)
    integer :: e, i, k, m, idx, n_GLL
    n_GLL = nelv * Xh_GLL%lxyz
    m = lx * lx
    do e = 1, nelv
       call mxm(dx,lx,u(1,1,1,e),lx,ur,m)
       do k=1,lx
          call mxm(u(1,1,k,e),lx,dy,lx,us(1,1,k),lx)
       enddo
       call mxm(u(1,1,1,e),m,dz,lx,ut,lx)
       do i=1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) + c(i,e,3) * ut(i,1,1)
       enddo
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_conv_fst_3d_lx7

  subroutine cpu_conv_fst_3d_lx6(du, u, c, dx, dy, dz, &
        Xh_GLL, coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 6
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), dimension(Xh_GLL%lx,Xh_GLL%lx,Xh_GLL%lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u
    real(kind=rp), dimension(lx*lx*lx,nelv,3), intent(in) :: c
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: ur, us, ut
    real(kind=rp):: ud(lx*lx*lx)
    integer :: e, i, k, m, idx, n_GLL
    n_GLL = nelv * Xh_GLL%lxyz
    m = lx * lx
    do e = 1, nelv
       call mxm(dx,lx,u(1,1,1,e),lx,ur,m)
       do k=1,lx
          call mxm(u(1,1,k,e),lx,dy,lx,us(1,1,k),lx)
       enddo
       call mxm(u(1,1,1,e),m,dz,lx,ut,lx)
       do i=1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) + c(i,e,3) * ut(i,1,1)
       enddo
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_conv_fst_3d_lx6

  subroutine cpu_conv_fst_3d_lx5(du, u, c, dx, dy, dz, &
        Xh_GLL, coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 5
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), dimension(Xh_GLL%lx,Xh_GLL%lx,Xh_GLL%lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u
    real(kind=rp), dimension(lx*lx*lx,nelv,3), intent(in) :: c
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: ur, us, ut
    real(kind=rp):: ud(lx*lx*lx)
    integer :: e, i, k, m, idx, n_GLL
    n_GLL = nelv * Xh_GLL%lxyz
    m = lx * lx
    do e = 1, nelv
       call mxm(dx,lx,u(1,1,1,e),lx,ur,m)
       do k=1,lx
          call mxm(u(1,1,k,e),lx,dy,lx,us(1,1,k),lx)
       enddo
       call mxm(u(1,1,1,e),m,dz,lx,ut,lx)
       do i=1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) + c(i,e,3) * ut(i,1,1)
       enddo
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_conv_fst_3d_lx5

  subroutine cpu_conv_fst_3d_lx4(du, u, c, dx, dy, dz, &
        Xh_GLL, coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 4
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), dimension(Xh_GLL%lx,Xh_GLL%lx,Xh_GLL%lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u
    real(kind=rp), dimension(lx*lx*lx,nelv,3), intent(in) :: c
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: ur, us, ut
    real(kind=rp):: ud(lx*lx*lx)
    integer :: e, i, k, m, idx, n_GLL
    n_GLL = nelv * Xh_GLL%lxyz
    m = lx * lx
    do e = 1, nelv
       call mxm(dx,lx,u(1,1,1,e),lx,ur,m)
       do k=1,lx
          call mxm(u(1,1,k,e),lx,dy,lx,us(1,1,k),lx)
       enddo
       call mxm(u(1,1,1,e),m,dz,lx,ut,lx)
       do i=1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) + c(i,e,3) * ut(i,1,1)
       enddo
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL) 

  end subroutine cpu_conv_fst_3d_lx4

  subroutine cpu_conv_fst_3d_lx3(du, u, c, dx, dy, dz, &
        Xh_GLL, coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 3
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), dimension(Xh_GLL%lx,Xh_GLL%lx,Xh_GLL%lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u
    real(kind=rp), dimension(lx*lx*lx,nelv,3), intent(in) :: c
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: ur, us, ut
    real(kind=rp):: ud(lx*lx*lx)
    integer :: e, i, k, m, idx, n_GLL
    n_GLL = nelv * Xh_GLL%lxyz
    m = lx * lx
    do e = 1, nelv
       call mxm(dx,lx,u(1,1,1,e),lx,ur,m)
       do k=1,lx
          call mxm(u(1,1,k,e),lx,dy,lx,us(1,1,k),lx)
       enddo
       call mxm(u(1,1,1,e),m,dz,lx,ut,lx)
       do i=1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) + c(i,e,3) * ut(i,1,1)
       enddo
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)

  end subroutine cpu_conv_fst_3d_lx3

  subroutine cpu_conv_fst_3d_lx2(du, u, c, dx, dy, dz, &
        Xh_GLL, coef_GLL, GLL_to_GL, nelv)
    integer, parameter :: lx = 2
    integer, intent(in) :: nelv
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), dimension(Xh_GLL%lx,Xh_GLL%lx,Xh_GLL%lx, nelv), intent(inout) :: du
    real(kind=rp), dimension(lx,lx,lx,nelv), intent(in) :: u
    real(kind=rp), dimension(lx*lx*lx,nelv,3), intent(in) :: c
    real(kind=rp), dimension(lx,lx), intent(in) :: dx, dy, dz
    real(kind=rp), dimension(lx,lx,lx) :: ur, us, ut
    real(kind=rp):: ud(lx*lx*lx)
    integer :: e, i, k, m, idx, n_GLL
    n_GLL = nelv * Xh_GLL%lxyz
    m = lx * lx
    do e = 1, nelv
       call mxm(dx,lx,u(1,1,1,e),lx,ur,m)
       do k=1,lx
          call mxm(u(1,1,k,e),lx,dy,lx,us(1,1,k),lx)
       enddo
       call mxm(u(1,1,1,e),m,dz,lx,ut,lx)
       do i=1, lx * lx * lx
          ud(i) = c(i,e,1) * ur(i,1,1) + c(i,e,2) * us(i,1,1) + c(i,e,3) * ut(i,1,1)
       enddo
       idx = (e-1) * Xh_GLL%lxyz+1
       call GLL_to_GL%map(du(idx,1,1,1), ud, 1, Xh_GLL)
    end do
    call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
    call col2(du, coef_GLL%Binv, n_GLL)
    
  end subroutine cpu_conv_fst_3d_lx2

end module cpu_conv_fst_3d
