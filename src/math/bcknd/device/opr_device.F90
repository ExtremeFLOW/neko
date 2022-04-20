! Copyright (c) 2021-2022, The Neko Authors
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
!> Operators accelerator backends
module opr_device
  use gather_scatter  
  use num_types
  use device_math
  use device    
  use space
  use coefs
  use math
  use mesh
  use field
  use mathops
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_HIP
  interface
     subroutine hip_dudxyz(du_d, u_d, dr_d, ds_d, dt_d, &
          dx_d, dy_d, dz_d, jacinv_d, nel, lx) &
          bind(c, name='hip_dudxyz')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: du_d, u_d, dr_d, ds_d, dt_d
       type(c_ptr), value :: dx_d, dy_d, dz_d, jacinv_d
       integer(c_int) :: nel, lx
     end subroutine hip_dudxyz
  end interface

  interface
     subroutine hip_cdtp(dtx_d, x_d, dr_d, ds_d, dt_d, &
          dxt_d, dyt_d, dzt_d, B_d, jac_d, nel, lx) &
          bind(c, name='hip_cdtp')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: dtx_d, x_d, dr_d, ds_d, dt_d
       type(c_ptr), value :: dxt_d, dyt_d, dzt_d, B_d, jac_d
       integer(c_int) :: nel, lx
     end subroutine hip_cdtp
  end interface

  interface
     subroutine hip_conv1(du_d, u_d, vx_d, vy_d, vz_d, &
          dx_d, dy_d, dz_d, drdx_d, dsdx_d, dtdx_d, &
          drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, dtdz_d, &
          jacinv_d, nel, gdim, lx) &
          bind(c, name='hip_conv1')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: du_d, u_d, vx_d, vy_d, vz_d
       type(c_ptr), value :: dx_d, dy_d, dz_d, drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: jacinv_d
       integer(c_int) :: nel, gdim, lx
     end subroutine hip_conv1
  end interface

  interface
     subroutine hip_opgrad(ux_d, uy_d, uz_d, u_d, &
          dx_d, dy_d, dz_d, &
          drdx_d, dsdx_d, dtdx_d, &
          drdy_d, dsdy_d, dtdy_d, &
          drdz_d, dsdz_d, dtdz_d, w3_d, nel, lx) &
          bind(c, name='hip_opgrad')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: ux_d, uy_d, uz_d, u_d
       type(c_ptr), value :: dx_d, dy_d, dz_d
       type(c_ptr), value :: drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d
       type(c_ptr), value :: drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: w3_d
       integer(c_int) :: nel, lx
     end subroutine hip_opgrad
  end interface

#elif HAVE_CUDA
  interface
     subroutine cuda_dudxyz(du_d, u_d, dr_d, ds_d, dt_d, &
          dx_d, dy_d, dz_d, jacinv_d, nel, lx) &
          bind(c, name='cuda_dudxyz')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: du_d, u_d, dr_d, ds_d, dt_d
       type(c_ptr), value :: dx_d, dy_d, dz_d, jacinv_d
       integer(c_int) :: nel, lx
     end subroutine cuda_dudxyz
  end interface

  interface
     subroutine cuda_cdtp(dtx_d, x_d, dr_d, ds_d, dt_d, &
          dxt_d, dyt_d, dzt_d, B_d, jac_d, nel, lx) &
          bind(c, name='cuda_cdtp')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: dtx_d, x_d, dr_d, ds_d, dt_d
       type(c_ptr), value :: dxt_d, dyt_d, dzt_d, B_d, jac_d
       integer(c_int) :: nel, lx
     end subroutine cuda_cdtp
  end interface

  interface
     subroutine cuda_conv1(du_d, u_d, vx_d, vy_d, vz_d, &
          dx_d, dy_d, dz_d, drdx_d, dsdx_d, dtdx_d, &
          drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, dtdz_d, &
          jacinv_d, nel, gdim, lx) &
          bind(c, name='cuda_conv1')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: du_d, u_d, vx_d, vy_d, vz_d
       type(c_ptr), value :: dx_d, dy_d, dz_d, drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: jacinv_d
       integer(c_int) :: nel, gdim, lx
     end subroutine cuda_conv1
  end interface
  
  interface
     subroutine cuda_opgrad(ux_d, uy_d, uz_d, u_d, &
          dx_d, dy_d, dz_d, &
          drdx_d, dsdx_d, dtdx_d, &
          drdy_d, dsdy_d, dtdy_d, &
          drdz_d, dsdz_d, dtdz_d, w3_d, nel, lx) &
          bind(c, name='cuda_opgrad')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: ux_d, uy_d, uz_d, u_d
       type(c_ptr), value :: dx_d, dy_d, dz_d
       type(c_ptr), value :: drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d
       type(c_ptr), value :: drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: w3_d
       integer(c_int) :: nel, lx
     end subroutine cuda_opgrad
  end interface
#elif HAVE_OPENCL
  interface
     subroutine opencl_dudxyz(du_d, u_d, dr_d, ds_d, dt_d, &
          dx_d, dy_d, dz_d, jacinv_d, nel, lx) &
          bind(c, name='opencl_dudxyz')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: du_d, u_d, dr_d, ds_d, dt_d
       type(c_ptr), value :: dx_d, dy_d, dz_d, jacinv_d
       integer(c_int) :: nel, lx
     end subroutine opencl_dudxyz
  end interface

  interface
     subroutine opencl_cdtp(dtx_d, x_d, dr_d, ds_d, dt_d, &
          dxt_d, dyt_d, dzt_d, B_d, jac_d, nel, lx) &
          bind(c, name='opencl_cdtp')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: dtx_d, x_d, dr_d, ds_d, dt_d
       type(c_ptr), value :: dxt_d, dyt_d, dzt_d, B_d, jac_d
       integer(c_int) :: nel, lx
     end subroutine opencl_cdtp
  end interface

  interface
     subroutine opencl_conv1(du_d, u_d, vx_d, vy_d, vz_d, &
          dx_d, dy_d, dz_d, drdx_d, dsdx_d, dtdx_d, &
          drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, dtdz_d, &
          jacinv_d, nel, gdim, lx) &
          bind(c, name='opencl_conv1')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: du_d, u_d, vx_d, vy_d, vz_d
       type(c_ptr), value :: dx_d, dy_d, dz_d, drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: jacinv_d
       integer(c_int) :: nel, gdim, lx
     end subroutine opencl_conv1
  end interface

  interface
     subroutine opencl_opgrad(ux_d, uy_d, uz_d, u_d, &
          dx_d, dy_d, dz_d, &
          drdx_d, dsdx_d, dtdx_d, &
          drdy_d, dsdy_d, dtdy_d, &
          drdz_d, dsdz_d, dtdz_d, w3_d, nel, lx) &
          bind(c, name='opencl_opgrad')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: ux_d, uy_d, uz_d, u_d
       type(c_ptr), value :: dx_d, dy_d, dz_d
       type(c_ptr), value :: drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d
       type(c_ptr), value :: drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: w3_d
       integer(c_int) :: nel, lx
     end subroutine opencl_opgrad
  end interface
#endif  
  
contains

  subroutine opr_device_dudxyz(du, u, dr, ds, dt, coef)
    type(coef_t), intent(in), target :: coef
    real(kind=rp), dimension(coef%Xh%lx,coef%Xh%ly, &
         coef%Xh%lz,coef%msh%nelv), intent(inout) ::  du
    real(kind=rp), dimension(coef%Xh%lx,coef%Xh%ly, &
         coef%Xh%lz,coef%msh%nelv), intent(in) ::  u, dr, ds, dt
    type(c_ptr) :: du_d, u_d, dr_d, ds_d, dt_d

    du_d = device_get_ptr(du)
    u_d = device_get_ptr(u)

    dr_d = device_get_ptr(dr)
    ds_d = device_get_ptr(ds)
    dt_d = device_get_ptr(dt)

    associate(Xh => coef%Xh, msh => coef%msh, dof => coef%dof)    
#ifdef HAVE_HIP
      call hip_dudxyz(du_d, u_d, dr_d, ds_d, dt_d, &
           Xh%dx_d, Xh%dy_d, Xh%dz_d, coef%jacinv_d, &
           msh%nelv, Xh%lx)
#elif HAVE_CUDA
      call cuda_dudxyz(du_d, u_d, dr_d, ds_d, dt_d, &
           Xh%dx_d, Xh%dy_d, Xh%dz_d, coef%jacinv_d, &
           msh%nelv, Xh%lx)
#elif HAVE_OPENCL
      call opencl_dudxyz(du_d, u_d, dr_d, ds_d, dt_d, &
           Xh%dx_d, Xh%dy_d, Xh%dz_d, coef%jacinv_d, &
           msh%nelv, Xh%lx)
#else
      call neko_error('No device backend configured')
#endif
    end associate
  
  end subroutine opr_device_dudxyz

  subroutine opr_device_opgrad(ux, uy, uz, u, coef) 
    type(coef_t), intent(in) :: coef  
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: ux
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: uy
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: uz
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: u
    type(c_ptr) :: ux_d, uy_d, uz_d, u_d

    ux_d = device_get_ptr(ux)
    uy_d = device_get_ptr(uy)
    uz_d = device_get_ptr(uz)

    u_d = device_get_ptr(u)
    
    associate(Xh => coef%Xh, msh => coef%msh)
#ifdef HAVE_HIP
      call hip_opgrad(ux_d, uy_d, uz_d, u_d, &
           Xh%dx_d, Xh%dy_d, Xh%dz_d, &
           coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
           coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
           coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
           Xh%w3_d, msh%nelv, Xh%lx)
#elif HAVE_CUDA
      call cuda_opgrad(ux_d, uy_d, uz_d, u_d, &
           Xh%dx_d, Xh%dy_d, Xh%dz_d, &
           coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
           coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
           coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
           Xh%w3_d, msh%nelv, Xh%lx)
#elif HAVE_OPENCL
      call opencl_opgrad(ux_d, uy_d, uz_d, u_d, &
           Xh%dx_d, Xh%dy_d, Xh%dz_d, &
           coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
           coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
           coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
           Xh%w3_d, msh%nelv, Xh%lx)
#else
      call neko_error('No device backend configured')
#endif
    end associate
    
  end subroutine opr_device_opgrad

  subroutine opr_device_cdtp(dtx, x, dr,ds, dt, coef)
    type(coef_t), intent(in) :: coef
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: dtx
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(inout) :: x
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: dr
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: ds
    real(kind=rp), dimension(coef%Xh%lxyz,coef%msh%nelv), intent(in) :: dt
    type(c_ptr) :: dtx_d, x_d, dr_d, ds_d, dt_d

    dtx_d = device_get_ptr(dtx)
    x_d = device_get_ptr(x)

    dr_d = device_get_ptr(dr)
    ds_d = device_get_ptr(ds)
    dt_d = device_get_ptr(dt)
    
    associate(Xh => coef%Xh, msh => coef%msh, dof => coef%dof)    
#ifdef HAVE_HIP
      call hip_cdtp(dtx_d, x_d, dr_d, ds_d, dt_d, &
           Xh%dxt_d, Xh%dyt_d, Xh%dzt_d, coef%B_d, &
           coef%jac_d, msh%nelv, Xh%lx)
#elif HAVE_CUDA
      call cuda_cdtp(dtx_d, x_d, dr_d, ds_d, dt_d, &
           Xh%dxt_d, Xh%dyt_d, Xh%dzt_d, coef%B_d, &
           coef%jac_d, msh%nelv, Xh%lx)
#elif HAVE_OPENCL
      call opencl_cdtp(dtx_d, x_d, dr_d, ds_d, dt_d, &
           Xh%dxt_d, Xh%dyt_d, Xh%dzt_d, coef%B_d, &
           coef%jac_d, msh%nelv, Xh%lx)
#else
      call neko_error('No device backend configured')
#endif
  end associate

  end subroutine opr_device_cdtp

  subroutine opr_device_conv1(du, u, vx, vy, vz, Xh, coef, nelv, gdim)  
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: nelv, gdim
    real(kind=rp), intent(inout) ::  du(Xh%lxyz,nelv)
    real(kind=rp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  u
    real(kind=rp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  vx
    real(kind=rp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  vy
    real(kind=rp), intent(inout), dimension(Xh%lx,Xh%ly,Xh%lz,nelv) ::  vz
    type(c_ptr) :: du_d, u_d, vx_d, vy_d, vz_d

    du_d = device_get_ptr(du)
    u_d = device_get_ptr(u)

    vx_d = device_get_ptr(vx)
    vy_d = device_get_ptr(vy)
    vz_d = device_get_ptr(vz)
    
    associate(Xh => coef%Xh, msh => coef%msh, dof => coef%dof)    
#ifdef HAVE_HIP
      call hip_conv1(du_d, u_d, vx_d, vy_d, vz_d, &
                     Xh%dx_d, Xh%dy_d, Xh%dz_d, &
                     coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
                     coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
                     coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
                     coef%jacinv_d, msh%nelv, msh%gdim, Xh%lx)
#elif HAVE_CUDA
      call cuda_conv1(du_d, u_d, vx_d, vy_d, vz_d, &
                      Xh%dx_d, Xh%dy_d, Xh%dz_d, &
                      coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
                      coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
                      coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
                      coef%jacinv_d, msh%nelv, msh%gdim, Xh%lx)
#elif HAVE_OPENCL
      call opencl_conv1(du_d, u_d, vx_d, vy_d, vz_d, &
                        Xh%dx_d, Xh%dy_d, Xh%dz_d, &
                        coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
                        coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
                        coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
                        coef%jacinv_d, msh%nelv, msh%gdim, Xh%lx)
#else
      call neko_error('No device backend configured')
#endif
    end associate
    
  end subroutine opr_device_conv1

  subroutine opr_device_curl(w1, w2, w3, u1, u2, u3, work1, work2, c_Xh)
    type(field_t), intent(inout) :: w1
    type(field_t), intent(inout) :: w2
    type(field_t), intent(inout) :: w3
    type(field_t), intent(inout) :: u1
    type(field_t), intent(inout) :: u2
    type(field_t), intent(inout) :: u3
    type(field_t), intent(inout) :: work1
    type(field_t), intent(inout) :: work2
    type(coef_t), intent(in)  :: c_Xh
    integer :: gdim, n, nelv

    n = w1%dof%size()
    gdim = c_Xh%msh%gdim
    nelv = c_Xh%msh%nelv

    !     this%work1=dw/dy ; this%work2=dv/dz
#if defined(HAVE_HIP) || defined(HAVE_CUDA) || defined(HAVE_OPENCL)
#ifdef HAVE_HIP
    call hip_dudxyz(work1%x_d, u3%x_d, &
           c_Xh%drdy_d, c_Xh%dsdy_d, c_Xh%dtdy_d,&
           c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
           c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
#elif HAVE_CUDA
    call cuda_dudxyz(work1%x_d, u3%x_d, &
           c_Xh%drdy_d, c_Xh%dsdy_d, c_Xh%dtdy_d,&
           c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
           c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
#elif HAVE_OPENCL
    call opencl_dudxyz(work1%x_d, u3%x_d, &
           c_Xh%drdy_d, c_Xh%dsdy_d, c_Xh%dtdy_d,&
           c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
           c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
#endif
    if (gdim .eq. 3) then
#ifdef HAVE_HIP
       call hip_dudxyz(work2%x_d, u2%x_d, &
            c_Xh%drdz_d, c_Xh%dsdz_d, c_Xh%dtdz_d,&
            c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
            c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
#elif HAVE_CUDA
       call cuda_dudxyz(work2%x_d, u2%x_d, &
            c_Xh%drdz_d, c_Xh%dsdz_d, c_Xh%dtdz_d,&
            c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
            c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
#elif HAVE_OPENCL
       call opencl_dudxyz(work2%x_d, u2%x_d, &
            c_Xh%drdz_d, c_Xh%dsdz_d, c_Xh%dtdz_d,&
            c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
            c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
#endif
       call device_sub3(w1%x_d, work1%x_d, work2%x_d, n)
    else
       call device_copy(w1%x_d, work1%x_d, n)
    endif
    !     this%work1=du/dz ; this%work2=dw/dx
    if (gdim .eq. 3) then
#ifdef HAVE_HIP
       call hip_dudxyz(work1%x_d, u1%x_d, &
            c_Xh%drdz_d, c_Xh%dsdz_d, c_Xh%dtdz_d,&
            c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
            c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
       call hip_dudxyz(work2%x_d, u3%x_d, &
            c_Xh%drdx_d, c_Xh%dsdx_d, c_Xh%dtdx_d,&
            c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
            c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
#elif HAVE_CUDA
       call cuda_dudxyz(work1%x_d, u1%x_d, &
            c_Xh%drdz_d, c_Xh%dsdz_d, c_Xh%dtdz_d,&
            c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
            c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
       call cuda_dudxyz(work2%x_d, u3%x_d, &
            c_Xh%drdx_d, c_Xh%dsdx_d, c_Xh%dtdx_d,&
            c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
            c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
#elif HAVE_OPENCL
       call opencl_dudxyz(work1%x_d, u1%x_d, &
            c_Xh%drdz_d, c_Xh%dsdz_d, c_Xh%dtdz_d,&
            c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
            c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
       call opencl_dudxyz(work2%x_d, u3%x_d, &
            c_Xh%drdx_d, c_Xh%dsdx_d, c_Xh%dtdx_d,&
            c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
            c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
#endif
       call device_sub3(w2%x_d, work1%x_d, work2%x_d, n)
    else
       call device_rzero (work1%x_d, n)
#ifdef HAVE_HIP
       call hip_dudxyz(work2%x_d, u3%x_d, &
            c_Xh%drdx_d, c_Xh%dsdx_d, c_Xh%dtdx_d,&
            c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
            c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
#elif HAVE_CUDA
       call cuda_dudxyz(work2%x_d, u3%x_d, &
            c_Xh%drdx_d, c_Xh%dsdx_d, c_Xh%dtdx_d,&
            c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
            c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
#elif HAVE_OPENCL
       call opencl_dudxyz(work2%x_d, u3%x_d, &
            c_Xh%drdx_d, c_Xh%dsdx_d, c_Xh%dtdx_d,&
            c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
            c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
#endif
       call device_sub3(w2%x_d, work1%x_d, work2%x_d, n)
    endif
    !     this%work1=dv/dx ; this%work2=du/dy
#ifdef HAVE_HIP
    call hip_dudxyz(work1%x_d, u2%x_d, &
         c_Xh%drdx_d, c_Xh%dsdx_d, c_Xh%dtdx_d,&
         c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
         c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
    call hip_dudxyz(work2%x_d, u1%x_d, &
         c_Xh%drdy_d, c_Xh%dsdy_d, c_Xh%dtdy_d,&
         c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
         c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
#elif HAVE_CUDA
    call cuda_dudxyz(work1%x_d, u2%x_d, &
         c_Xh%drdx_d, c_Xh%dsdx_d, c_Xh%dtdx_d,&
         c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
         c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
    call cuda_dudxyz(work2%x_d, u1%x_d, &
         c_Xh%drdy_d, c_Xh%dsdy_d, c_Xh%dtdy_d,&
         c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
         c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
#elif HAVE_OPENCL
    call opencl_dudxyz(work1%x_d, u2%x_d, &
         c_Xh%drdx_d, c_Xh%dsdx_d, c_Xh%dtdx_d,&
         c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
         c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
    call opencl_dudxyz(work2%x_d, u1%x_d, &
         c_Xh%drdy_d, c_Xh%dsdy_d, c_Xh%dtdy_d,&
         c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
         c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
#endif
    call device_sub3(w3%x_d, work1%x_d, work2%x_d, n)
    !!    BC dependent, Needs to change if cyclic

    !Change to opcolv when there's a device version...
    call device_col2(w1%x_d, c_Xh%B_d, n)
    call device_col2(w2%x_d, c_Xh%B_d, n)
    call device_col2(w3%x_d, c_Xh%B_d, n)
    call gs_op(c_Xh%gs_h, w1, GS_OP_ADD) 
    call gs_op(c_Xh%gs_h, w2, GS_OP_ADD) 
    call gs_op(c_Xh%gs_h, w3, GS_OP_ADD)
    !Change to opcolv when there's a device version...
    call device_col2(w1%x_d, c_Xh%Binv_d, n)
    call device_col2(w2%x_d, c_Xh%Binv_d, n)
    call device_col2(w3%x_d, c_Xh%Binv_d, n)

#else
    call neko_error('No device backend configured')
#endif
        
  end subroutine opr_device_curl



end module opr_device
