! Copyright (c) 2021-2026, The Neko Authors
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
  use gather_scatter, only : GS_OP_ADD
  use num_types, only : rp, c_rp, i8, dp, c_dp, xp, c_xp
  use device, only : device_get_ptr, device_event_sync, device_map, device_unmap
  use space, only : space_t
  use coefs, only : coef_t
  use field, only : field_t
  use vector, only : vector_t
  use scratch_registry, only : neko_scratch_registry
  use utils, only : neko_error
  use interpolation, only : interpolator_t
  use device_math, only : device_sub3, device_rzero, device_copy, device_col2, &
       device_glsum, device_cadd
  use device_mathops, only : device_opcolv
  use, intrinsic :: iso_c_binding
  implicit none
  private

  public :: opr_device_dudxyz, opr_device_opgrad, opr_device_cdtp, &
       opr_device_conv1, opr_device_convect_scalar, opr_device_curl, &
       opr_device_cfl, opr_device_lambda2, opr_device_set_convect_rst, &
       opr_device_rotate_cyc, device_ortho

#ifdef HAVE_HIP
  interface
     subroutine hip_dudxyz(du_d, u_d, dr_d, ds_d, dt_d, &
          dx_d, dy_d, dz_d, jacinv_d, nel, lx) &
          bind(c, name = 'hip_dudxyz')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: du_d, u_d, dr_d, ds_d, dt_d
       type(c_ptr), value :: dx_d, dy_d, dz_d, jacinv_d
       integer(c_int) :: nel, lx
     end subroutine hip_dudxyz
  end interface

  interface
     subroutine hip_cdtp(dtx_d, x_d, dr_d, ds_d, dt_d, &
          dxt_d, dyt_d, dzt_d, w3_d, nel, lx) &
          bind(c, name = 'hip_cdtp')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: dtx_d, x_d, dr_d, ds_d, dt_d
       type(c_ptr), value :: dxt_d, dyt_d, dzt_d, w3_d
       integer(c_int) :: nel, lx
     end subroutine hip_cdtp
  end interface

  interface
     subroutine hip_conv1(du_d, u_d, vx_d, vy_d, vz_d, &
          dx_d, dy_d, dz_d, drdx_d, dsdx_d, dtdx_d, &
          drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, dtdz_d, &
          jacinv_d, nel, gdim, lx) &
          bind(c, name = 'hip_conv1')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: du_d, u_d, vx_d, vy_d, vz_d
       type(c_ptr), value :: dx_d, dy_d, dz_d, drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: jacinv_d
       integer(c_int) :: nel, gdim, lx
     end subroutine hip_conv1
  end interface

  interface
     subroutine hip_convect_scalar(du_d, u_d, cr_d, cs_d, ct_d, &
          dx_d, dy_d, dz_d, nel, lx) bind(c, name = 'hip_convect_scalar')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: du_d, u_d
       type(c_ptr), value :: cr_d, cs_d, ct_d
       type(c_ptr), value :: dx_d, dy_d, dz_d
       integer(c_int) :: nel, lx
     end subroutine hip_convect_scalar
  end interface

  interface
     subroutine hip_opgrad(ux_d, uy_d, uz_d, u_d, &
          dx_d, dy_d, dz_d, &
          drdx_d, dsdx_d, dtdx_d, &
          drdy_d, dsdy_d, dtdy_d, &
          drdz_d, dsdz_d, dtdz_d, w3_d, nel, lx) &
          bind(c, name = 'hip_opgrad')
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

  interface
     subroutine hip_lambda2(lambda2_d, u_d, v_d, w_d, &
          dx_d, dy_d, dz_d, &
          drdx_d, dsdx_d, dtdx_d, &
          drdy_d, dsdy_d, dtdy_d, &
          drdz_d, dsdz_d, dtdz_d, jacinv_d, nel, lx) &
          bind(c, name = 'hip_lambda2')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: lambda2_d, u_d, v_d, w_d
       type(c_ptr), value :: dx_d, dy_d, dz_d
       type(c_ptr), value :: drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d
       type(c_ptr), value :: drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: jacinv_d
       integer(c_int) :: nel, lx
     end subroutine hip_lambda2
  end interface

  interface
     real(c_dp) function hip_cfl(dt, u_d, v_d, w_d, &
          drdx_d, dsdx_d, dtdx_d, drdy_d, dsdy_d, dtdy_d, &
          drdz_d, dsdz_d, dtdz_d, dr_inv_d, ds_inv_d, dt_inv_d, &
          jacinv_d, nel, lx) &
          bind(c, name = 'hip_cfl')
       use, intrinsic :: iso_c_binding
       import c_dp
       type(c_ptr), value :: u_d, v_d, w_d, drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: dr_inv_d, ds_inv_d, dt_inv_d, jacinv_d
       real(c_dp) :: dt
       integer(c_int) :: nel, lx
     end function hip_cfl
  end interface

  interface
     subroutine hip_rotate_cyc(vx_d, vy_d, vz_d, &
          x_d, y_d, z_d, &
          cyc_msk_d, R11_d, R12_d, ncyc, idir) &
          bind(c, name = 'hip_rotate_cyc')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: vx_d, vy_d, vz_d
       type(c_ptr), value :: x_d, y_d, z_d
       type(c_ptr), value :: cyc_msk_d, R11_d, R12_d
       integer(c_int) :: ncyc, idir
     end subroutine hip_rotate_cyc
  end interface

  interface
     subroutine hip_set_convect_rst(cr_d, cs_d, ct_d, cx_d, cy_d, cz_d, &
          drdx_d, dsdx_d, dtdx_d, drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, &
          dtdz_d, w3_d, nel, lx) bind(c, name = 'hip_set_convect_rst')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: cr_d, cs_d, ct_d
       type(c_ptr), value :: cx_d, cy_d, cz_d
       type(c_ptr), value :: drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d
       type(c_ptr), value :: drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: w3_d
       integer(c_int) :: nel, lx
     end subroutine hip_set_convect_rst
  end interface

#elif HAVE_CUDA
  interface
     subroutine cuda_dudxyz(du_d, u_d, dr_d, ds_d, dt_d, &
          dx_d, dy_d, dz_d, jacinv_d, nel, lx) &
          bind(c, name = 'cuda_dudxyz')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: du_d, u_d, dr_d, ds_d, dt_d
       type(c_ptr), value :: dx_d, dy_d, dz_d, jacinv_d
       integer(c_int) :: nel, lx
     end subroutine cuda_dudxyz
  end interface

  interface
     subroutine cuda_cdtp(dtx_d, x_d, dr_d, ds_d, dt_d, &
          dxt_d, dyt_d, dzt_d, w3_d, nel, lx) &
          bind(c, name = 'cuda_cdtp')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: dtx_d, x_d, dr_d, ds_d, dt_d
       type(c_ptr), value :: dxt_d, dyt_d, dzt_d, w3_d
       integer(c_int) :: nel, lx
     end subroutine cuda_cdtp
  end interface

  interface
     subroutine cuda_conv1(du_d, u_d, vx_d, vy_d, vz_d, &
          dx_d, dy_d, dz_d, drdx_d, dsdx_d, dtdx_d, &
          drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, dtdz_d, &
          jacinv_d, nel, gdim, lx) &
          bind(c, name = 'cuda_conv1')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: du_d, u_d, vx_d, vy_d, vz_d
       type(c_ptr), value :: dx_d, dy_d, dz_d, drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: jacinv_d
       integer(c_int) :: nel, gdim, lx
     end subroutine cuda_conv1
  end interface

  interface
     subroutine cuda_convect_scalar(du_d, u_d, cr_d, cs_d, ct_d, &
          dx_d, dy_d, dz_d, nel, lx) bind(c, name = 'cuda_convect_scalar')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: du_d, u_d
       type(c_ptr), value :: cr_d, cs_d, ct_d
       type(c_ptr), value :: dx_d, dy_d, dz_d
       integer(c_int) :: nel, lx
     end subroutine cuda_convect_scalar
  end interface

  interface
     subroutine cuda_opgrad(ux_d, uy_d, uz_d, u_d, &
          dx_d, dy_d, dz_d, &
          drdx_d, dsdx_d, dtdx_d, &
          drdy_d, dsdy_d, dtdy_d, &
          drdz_d, dsdz_d, dtdz_d, w3_d, nel, lx) &
          bind(c, name = 'cuda_opgrad')
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

  interface
     subroutine cuda_lambda2(lambda2_d, u_d, v_d, w_d, &
          dx_d, dy_d, dz_d, &
          drdx_d, dsdx_d, dtdx_d, &
          drdy_d, dsdy_d, dtdy_d, &
          drdz_d, dsdz_d, dtdz_d, jacinv_d, nel, lx) &
          bind(c, name = 'cuda_lambda2')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: lambda2_d, u_d, v_d, w_d
       type(c_ptr), value :: dx_d, dy_d, dz_d
       type(c_ptr), value :: drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d
       type(c_ptr), value :: drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: jacinv_d
       integer(c_int) :: nel, lx
     end subroutine cuda_lambda2
  end interface

  interface
     real(c_dp) function cuda_cfl(dt, u_d, v_d, w_d, &
          drdx_d, dsdx_d, dtdx_d, drdy_d, dsdy_d, dtdy_d, &
          drdz_d, dsdz_d, dtdz_d, dr_inv_d, ds_inv_d, dt_inv_d, &
          jacinv_d, nel, lx) &
          bind(c, name = 'cuda_cfl')
       use, intrinsic :: iso_c_binding
       import c_dp
       type(c_ptr), value :: u_d, v_d, w_d, drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: dr_inv_d, ds_inv_d, dt_inv_d, jacinv_d
       real(c_dp) :: dt
       integer(c_int) :: nel, lx
     end function cuda_cfl
  end interface

  interface
     subroutine cuda_rotate_cyc(vx_d, vy_d, vz_d, &
          x_d, y_d, z_d, &
          cyc_msk_d, R11_d, R12_d, ncyc, idir) &
          bind(c, name = 'cuda_rotate_cyc')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: vx_d, vy_d, vz_d
       type(c_ptr), value :: x_d, y_d, z_d
       type(c_ptr), value :: cyc_msk_d, R11_d, R12_d
       integer(c_int) :: ncyc, idir
     end subroutine cuda_rotate_cyc
  end interface

  interface
     subroutine cuda_set_convect_rst(cr_d, cs_d, ct_d, cx_d, cy_d, cz_d, &
          drdx_d, dsdx_d, dtdx_d, drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, &
          dtdz_d, w3_d, nel, lx) bind(c, name = 'cuda_set_convect_rst')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: cr_d, cs_d, ct_d
       type(c_ptr), value :: cx_d, cy_d, cz_d
       type(c_ptr), value :: drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d
       type(c_ptr), value :: drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: w3_d
       integer(c_int) :: nel, lx
     end subroutine cuda_set_convect_rst
  end interface

#elif HAVE_OPENCL
  interface
     subroutine opencl_dudxyz(du_d, u_d, dr_d, ds_d, dt_d, &
          dx_d, dy_d, dz_d, jacinv_d, nel, lx) &
          bind(c, name = 'opencl_dudxyz')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: du_d, u_d, dr_d, ds_d, dt_d
       type(c_ptr), value :: dx_d, dy_d, dz_d, jacinv_d
       integer(c_int) :: nel, lx
     end subroutine opencl_dudxyz
  end interface

  interface
     subroutine opencl_cdtp(dtx_d, x_d, dr_d, ds_d, dt_d, &
          dxt_d, dyt_d, dzt_d, w3_d, nel, lx) &
          bind(c, name = 'opencl_cdtp')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: dtx_d, x_d, dr_d, ds_d, dt_d
       type(c_ptr), value :: dxt_d, dyt_d, dzt_d, w3_d
       integer(c_int) :: nel, lx
     end subroutine opencl_cdtp
  end interface

  interface
     subroutine opencl_conv1(du_d, u_d, vx_d, vy_d, vz_d, &
          dx_d, dy_d, dz_d, drdx_d, dsdx_d, dtdx_d, &
          drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, dtdz_d, &
          jacinv_d, nel, gdim, lx) &
          bind(c, name = 'opencl_conv1')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: du_d, u_d, vx_d, vy_d, vz_d
       type(c_ptr), value :: dx_d, dy_d, dz_d, drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: jacinv_d
       integer(c_int) :: nel, gdim, lx
     end subroutine opencl_conv1
  end interface

  interface
     subroutine opencl_convect_scalar(du_d, u_d, cr_d, cs_d, ct_d, &
          dx_d, dy_d, dz_d, nel, lx) bind(c, name = 'opencl_convect_scalar')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: du_d, u_d
       type(c_ptr), value :: cr_d, cs_d, ct_d
       type(c_ptr), value :: dx_d, dy_d, dz_d
       integer(c_int) :: nel, lx
     end subroutine opencl_convect_scalar
  end interface

  interface
     subroutine opencl_opgrad(ux_d, uy_d, uz_d, u_d, &
          dx_d, dy_d, dz_d, &
          drdx_d, dsdx_d, dtdx_d, &
          drdy_d, dsdy_d, dtdy_d, &
          drdz_d, dsdz_d, dtdz_d, w3_d, nel, lx) &
          bind(c, name = 'opencl_opgrad')
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

  interface
     real(c_xp) function opencl_cfl(dt, u_d, v_d, w_d, &
          drdx_d, dsdx_d, dtdx_d, drdy_d, dsdy_d, dtdy_d, &
          drdz_d, dsdz_d, dtdz_d, dr_inv_d, ds_inv_d, dt_inv_d, &
          jacinv_d, nel, lx) &
          bind(c, name = 'opencl_cfl')
       use, intrinsic :: iso_c_binding
       import c_xp
       type(c_ptr), value :: u_d, v_d, w_d, drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: dr_inv_d, ds_inv_d, dt_inv_d, jacinv_d
       real(c_xp) :: dt
       integer(c_int) :: nel, lx
     end function opencl_cfl
  end interface

  interface
     subroutine opencl_lambda2(lambda2_d, u_d, v_d, w_d, &
          dx_d, dy_d, dz_d, &
          drdx_d, dsdx_d, dtdx_d, &
          drdy_d, dsdy_d, dtdy_d, &
          drdz_d, dsdz_d, dtdz_d, jacinv_d, nel, lx) &
          bind(c, name = 'opencl_lambda2')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: lambda2_d, u_d, v_d, w_d
       type(c_ptr), value :: dx_d, dy_d, dz_d
       type(c_ptr), value :: drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d
       type(c_ptr), value :: drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: jacinv_d
       integer(c_int) :: nel, lx
     end subroutine opencl_lambda2
  end interface

  interface
     subroutine opencl_rotate_cyc(vx_d, vy_d, vz_d, &
          x_d, y_d, z_d, &
          cyc_msk_d, R11_d, R12_d, ncyc, idir) &
          bind(c, name = 'opencl_rotate_cyc')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: vx_d, vy_d, vz_d
       type(c_ptr), value :: x_d, y_d, z_d
       type(c_ptr), value :: cyc_msk_d, R11_d, R12_d
       integer(c_int) :: ncyc, idir
     end subroutine opencl_rotate_cyc
  end interface

  interface
     subroutine opencl_set_convect_rst(cr_d, cs_d, ct_d, cx_d, cy_d, cz_d, &
          drdx_d, dsdx_d, dtdx_d, drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, &
          dtdz_d, w3_d, nel, lx) bind(c, name = 'opencl_set_convect_rst')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: cr_d, cs_d, ct_d
       type(c_ptr), value :: cx_d, cy_d, cz_d
       type(c_ptr), value :: drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d
       type(c_ptr), value :: drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: w3_d
       integer(c_int) :: nel, lx
     end subroutine opencl_set_convect_rst
  end interface

#elif HAVE_METAL
  interface
     subroutine metal_dudxyz(du_d, u_d, dr_d, ds_d, dt_d, &
          dx_d, dy_d, dz_d, jacinv_d, nel, lx) &
          bind(c, name = 'metal_dudxyz')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: du_d, u_d, dr_d, ds_d, dt_d
       type(c_ptr), value :: dx_d, dy_d, dz_d, jacinv_d
       integer(c_int) :: nel, lx
     end subroutine metal_dudxyz
  end interface

  interface
     subroutine metal_cdtp(dtx_d, x_d, dr_d, ds_d, dt_d, &
          dxt_d, dyt_d, dzt_d, w3_d, nel, lx) &
          bind(c, name = 'metal_cdtp')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: dtx_d, x_d, dr_d, ds_d, dt_d
       type(c_ptr), value :: dxt_d, dyt_d, dzt_d, w3_d
       integer(c_int) :: nel, lx
     end subroutine metal_cdtp
  end interface

  interface
     subroutine metal_conv1(du_d, u_d, vx_d, vy_d, vz_d, &
          dx_d, dy_d, dz_d, drdx_d, dsdx_d, dtdx_d, &
          drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, dtdz_d, &
          jacinv_d, nel, gdim, lx) &
          bind(c, name = 'metal_conv1')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: du_d, u_d, vx_d, vy_d, vz_d
       type(c_ptr), value :: dx_d, dy_d, dz_d, drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: jacinv_d
       integer(c_int) :: nel, gdim, lx
     end subroutine metal_conv1
  end interface

  interface
     subroutine metal_convect_scalar(du_d, u_d, cr_d, cs_d, ct_d, &
          dx_d, dy_d, dz_d, nel, lx) bind(c, name = 'metal_convect_scalar')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: du_d, u_d
       type(c_ptr), value :: cr_d, cs_d, ct_d
       type(c_ptr), value :: dx_d, dy_d, dz_d
       integer(c_int) :: nel, lx
     end subroutine metal_convect_scalar
  end interface

  interface
     subroutine metal_opgrad(ux_d, uy_d, uz_d, u_d, &
          dx_d, dy_d, dz_d, &
          drdx_d, dsdx_d, dtdx_d, &
          drdy_d, dsdy_d, dtdy_d, &
          drdz_d, dsdz_d, dtdz_d, w3_d, nel, lx) &
          bind(c, name = 'metal_opgrad')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: ux_d, uy_d, uz_d, u_d
       type(c_ptr), value :: dx_d, dy_d, dz_d
       type(c_ptr), value :: drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d
       type(c_ptr), value :: drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: w3_d
       integer(c_int) :: nel, lx
     end subroutine metal_opgrad
  end interface

  interface
     real(c_rp) function metal_cfl(dt, u_d, v_d, w_d, &
          drdx_d, dsdx_d, dtdx_d, drdy_d, dsdy_d, dtdy_d, &
          drdz_d, dsdz_d, dtdz_d, dr_inv_d, ds_inv_d, dt_inv_d, &
          jacinv_d, nel, lx) &
          bind(c, name = 'metal_cfl')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: u_d, v_d, w_d, drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: dr_inv_d, ds_inv_d, dt_inv_d, jacinv_d
       real(c_rp) :: dt
       integer(c_int) :: nel, lx
     end function metal_cfl
  end interface

  interface
     subroutine metal_lambda2(lambda2_d, u_d, v_d, w_d, &
          dx_d, dy_d, dz_d, &
          drdx_d, dsdx_d, dtdx_d, &
          drdy_d, dsdy_d, dtdy_d, &
          drdz_d, dsdz_d, dtdz_d, jacinv_d, nel, lx) &
          bind(c, name = 'metal_lambda2')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: lambda2_d, u_d, v_d, w_d
       type(c_ptr), value :: dx_d, dy_d, dz_d
       type(c_ptr), value :: drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d
       type(c_ptr), value :: drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: jacinv_d
       integer(c_int) :: nel, lx
     end subroutine metal_lambda2
  end interface

  interface
     subroutine metal_rotate_cyc(vx_d, vy_d, vz_d, &
          x_d, y_d, z_d, &
          cyc_msk_d, R11_d, R12_d, ncyc, idir) &
          bind(c, name = 'metal_rotate_cyc')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: vx_d, vy_d, vz_d
       type(c_ptr), value :: x_d, y_d, z_d
       type(c_ptr), value :: cyc_msk_d, R11_d, R12_d
       integer(c_int) :: ncyc, idir
     end subroutine metal_rotate_cyc
  end interface

  interface
     subroutine metal_set_convect_rst(cr_d, cs_d, ct_d, cx_d, cy_d, cz_d, &
          drdx_d, dsdx_d, dtdx_d, drdy_d, dsdy_d, dtdy_d, drdz_d, dsdz_d, &
          dtdz_d, w3_d, nel, lx) bind(c, name = 'metal_set_convect_rst')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: cr_d, cs_d, ct_d
       type(c_ptr), value :: cx_d, cy_d, cz_d
       type(c_ptr), value :: drdx_d, dsdx_d, dtdx_d
       type(c_ptr), value :: drdy_d, dsdy_d, dtdy_d
       type(c_ptr), value :: drdz_d, dsdz_d, dtdz_d
       type(c_ptr), value :: w3_d
       integer(c_int) :: nel, lx
     end subroutine metal_set_convect_rst
  end interface

#endif

contains

  subroutine opr_device_dudxyz(du_d, u_d, dr_d, ds_d, dt_d, coef)
    type(coef_t), intent(in), target :: coef
    type(c_ptr), intent(inout) :: du_d
    type(c_ptr), intent(in) :: u_d, dr_d, ds_d, dt_d

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
#elif HAVE_METAL
      call metal_dudxyz(du_d, u_d, dr_d, ds_d, dt_d, &
           Xh%dx_d, Xh%dy_d, Xh%dz_d, coef%jacinv_d, &
           msh%nelv, Xh%lx)
#else
      call neko_error('No device backend configured')
#endif
    end associate

  end subroutine opr_device_dudxyz

  subroutine opr_device_opgrad(ux_d, uy_d, uz_d, u_d, coef)
    type(coef_t), intent(in), target :: coef
    type(c_ptr), intent(inout) :: ux_d, uy_d, uz_d
    type(c_ptr), intent(in) :: u_d

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
#elif HAVE_METAL
      call metal_opgrad(ux_d, uy_d, uz_d, u_d, &
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

  !> Othogonalize with regard to vector (1,1,1,1,1,1...,1)^T.
  !! @param x_d The device vector to orthogonolize.
  !! @param glb_n_points The global number of non-unique gll points in the grid.
  !! @note This is equivalent to subtracting the mean of `x` from each of its
  !! elements.
  subroutine device_ortho(x_d, glb_n_points, n)
    integer, intent(in) :: n
    integer(kind=i8), intent(in) :: glb_n_points
    type(c_ptr), intent(inout) :: x_d
    real(kind=rp) :: c

    c = device_glsum(x_d, n) / glb_n_points
    call device_cadd(x_d, -c, n)

  end subroutine device_ortho

  subroutine opr_device_lambda2(lambda2_d, u_d, v_d, w_d, coef)
    type(coef_t), intent(in) :: coef
    type(c_ptr), intent(inout) :: lambda2_d
    type(c_ptr), intent(in) :: u_d, v_d, w_d
#ifdef HAVE_HIP
    call hip_lambda2(lambda2_d, u_d, v_d, w_d, &
         coef%Xh%dx_d, coef%Xh%dy_d, coef%Xh%dz_d, &
         coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
         coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
         coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
         coef%jacinv_d, coef%msh%nelv, coef%Xh%lx)
#elif HAVE_CUDA
    call cuda_lambda2(lambda2_d, u_d, v_d, w_d, &
         coef%Xh%dx_d, coef%Xh%dy_d, coef%Xh%dz_d, &
         coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
         coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
         coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
         coef%jacinv_d, coef%msh%nelv, coef%Xh%lx)
#elif HAVE_OPENCL
    call opencl_lambda2(lambda2_d, u_d, v_d, w_d, &
         coef%Xh%dx_d, coef%Xh%dy_d, coef%Xh%dz_d, &
         coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
         coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
         coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
         coef%jacinv_d, coef%msh%nelv, coef%Xh%lx)
#elif HAVE_METAL
    call metal_lambda2(lambda2_d, u_d, v_d, w_d, &
         coef%Xh%dx_d, coef%Xh%dy_d, coef%Xh%dz_d, &
         coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
         coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
         coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
         coef%jacinv_d, coef%msh%nelv, coef%Xh%lx)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine opr_device_lambda2

  subroutine opr_device_cdtp(dtx_d, x_d, dr_d, ds_d, dt_d, coef)
    type(coef_t), intent(in), target :: coef
    type(c_ptr), intent(inout) :: dtx_d, x_d
    type(c_ptr), intent(in) :: dr_d, ds_d, dt_d

    associate(Xh => coef%Xh, msh => coef%msh, dof => coef%dof)
#ifdef HAVE_HIP
      call hip_cdtp(dtx_d, x_d, dr_d, ds_d, dt_d, &
           Xh%dxt_d, Xh%dyt_d, Xh%dzt_d, Xh%w3_d, &
           msh%nelv, Xh%lx)
#elif HAVE_CUDA
      call cuda_cdtp(dtx_d, x_d, dr_d, ds_d, dt_d, &
           Xh%dxt_d, Xh%dyt_d, Xh%dzt_d, Xh%w3_d, &
           msh%nelv, Xh%lx)
#elif HAVE_OPENCL
      call opencl_cdtp(dtx_d, x_d, dr_d, ds_d, dt_d, &
           Xh%dxt_d, Xh%dyt_d, Xh%dzt_d, Xh%w3_d, &
           msh%nelv, Xh%lx)
#elif HAVE_METAL
      call metal_cdtp(dtx_d, x_d, dr_d, ds_d, dt_d, &
           Xh%dxt_d, Xh%dyt_d, Xh%dzt_d, Xh%w3_d, &
           msh%nelv, Xh%lx)
#else
      call neko_error('No device backend configured')
#endif
    end associate

  end subroutine opr_device_cdtp

  subroutine opr_device_conv1(du_d, u_d, vx_d, vy_d, vz_d, Xh, coef, nelv, gdim)
    type(space_t), intent(in) :: Xh
    type(coef_t), intent(in), target :: coef
    integer, intent(in) :: nelv, gdim
    type(c_ptr), intent(inout) :: du_d
    type(c_ptr), intent(in) :: u_d, vx_d, vy_d, vz_d

    associate(msh => coef%msh, dof => coef%dof)
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
#elif HAVE_METAL
      call metal_conv1(du_d, u_d, vx_d, vy_d, vz_d, &
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

  subroutine opr_device_convect_scalar(du, u_d, cr_d, cs_d, ct_d, &
       Xh_GLL, Xh_GL, coef_GLL, coef_GL, GLL_to_GL)
    type(space_t), intent(in) :: Xh_GL
    type(space_t), intent(in) :: Xh_GLL
    type(coef_t), intent(in) :: coef_GLL
    type(coef_t), intent(in) :: coef_GL
    type(interpolator_t), intent(inout) :: GLL_to_GL
    real(kind=rp), intent(inout) :: &
         du(Xh_GLL%lx, Xh_GLL%ly, Xh_GLL%lz, coef_GL%msh%nelv)
    type(c_ptr) :: cr_d, cs_d, ct_d, u_d
    type(vector_t), pointer :: ud
    type(c_ptr) :: du_d
    integer :: temp_index
    integer :: n_GL, n_GLL

    n_GL = coef_GL%msh%nelv * Xh_GL%lxyz
    n_GLL = coef_GL%msh%nelv * Xh_GLL%lxyz

    call neko_scratch_registry%request_vector(ud, temp_index, n_GL, .false.)

    du_d = device_get_ptr(du)

    associate(Xh => Xh_GL, nelv => coef_GL%msh%nelv, lx => Xh_GL%lx, &
         ud_d => ud%x_d)
#ifdef HAVE_HIP
      call hip_convect_scalar(ud_d, u_d, cr_d, cs_d, ct_d, &
           Xh%dx_d, Xh%dy_d, Xh%dz_d, nelv, lx)
#elif HAVE_CUDA
      call cuda_convect_scalar(ud_d, u_d, cr_d, cs_d, ct_d, &
           Xh%dx_d, Xh%dy_d, Xh%dz_d, nelv, lx)
#elif HAVE_OPENCL
      call opencl_convect_scalar(ud_d, u_d, cr_d, cs_d, ct_d, &
           Xh%dx_d, Xh%dy_d, Xh%dz_d, nelv, lx)
#elif HAVE_METAL
      call metal_convect_scalar(ud_d, u_d, cr_d, cs_d, ct_d, &
           Xh%dx_d, Xh%dy_d, Xh%dz_d, nelv, lx)
#else
      call neko_error('No device backend configured')
#endif

      call GLL_to_GL%map(du, ud%x, nelv, Xh_GLL)
      call coef_GLL%gs_h%op(du, n_GLL, GS_OP_ADD)
      call device_col2(du_d, coef_GLL%Binv_d, n_GLL)

    end associate

    call neko_scratch_registry%relinquish(temp_index)

  end subroutine opr_device_convect_scalar

  subroutine opr_device_curl(w1, w2, w3, u1, u2, u3, work1, work2, c_Xh, event)
    type(field_t), intent(inout) :: w1
    type(field_t), intent(inout) :: w2
    type(field_t), intent(inout) :: w3
    type(field_t), intent(in) :: u1
    type(field_t), intent(in) :: u2
    type(field_t), intent(in) :: u3
    type(field_t), intent(inout) :: work1
    type(field_t), intent(inout) :: work2
    type(coef_t), intent(in) :: c_Xh
    type(c_ptr), optional, intent(inout) :: event
    integer :: gdim, n, nelv

    n = w1%dof%size()
    gdim = c_Xh%msh%gdim
    nelv = c_Xh%msh%nelv

    !     this%work1=dw/dy ; this%work2=dv/dz
#if defined(HAVE_HIP) || defined(HAVE_CUDA) || defined(HAVE_OPENCL) || defined(HAVE_METAL)
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
#elif HAVE_METAL
    call metal_dudxyz(work1%x_d, u3%x_d, &
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
#elif HAVE_METAL
       call metal_dudxyz(work2%x_d, u2%x_d, &
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
#elif HAVE_METAL
       call metal_dudxyz(work1%x_d, u1%x_d, &
            c_Xh%drdz_d, c_Xh%dsdz_d, c_Xh%dtdz_d,&
            c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
            c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
       call metal_dudxyz(work2%x_d, u3%x_d, &
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
#elif HAVE_METAL
       call metal_dudxyz(work2%x_d, u3%x_d, &
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
#elif HAVE_METAL
    call metal_dudxyz(work1%x_d, u2%x_d, &
         c_Xh%drdx_d, c_Xh%dsdx_d, c_Xh%dtdx_d,&
         c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
         c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
    call metal_dudxyz(work2%x_d, u1%x_d, &
         c_Xh%drdy_d, c_Xh%dsdy_d, c_Xh%dtdy_d,&
         c_Xh%Xh%dx_d, c_Xh%Xh%dy_d, c_Xh%Xh%dz_d, &
         c_Xh%jacinv_d, nelv, c_Xh%Xh%lx)
#endif
    call device_sub3(w3%x_d, work1%x_d, work2%x_d, n)
    !!    BC dependent, Needs to change if cyclic

    call device_opcolv(w1%x_d, w2%x_d, w3%x_d, c_Xh%B_d, gdim, n)

    if(c_Xh%cyclic) call opr_device_rotate_cyc(w1%x_d, w2%x_d, w3%x_d, 1, c_Xh)
    if (present(event)) then
       call c_Xh%gs_h%op(w1%x, w2%x, w3%x, n, GS_OP_ADD, event)
       call device_event_sync(event)
    else
       call c_Xh%gs_h%op(w1%x, w2%x, w3%x, n, GS_OP_ADD)
    end if
    if(c_Xh%cyclic) call opr_device_rotate_cyc(w1%x_d, w2%x_d, w3%x_d, 0, c_Xh)

    call device_opcolv(w1%x_d, w2%x_d, w3%x_d, c_Xh%Binv_d, gdim, n)

#else
    call neko_error('No device backend configured')
#endif

  end subroutine opr_device_curl

  function opr_device_cfl(dt, u_d, v_d, w_d, Xh, coef, nelv, gdim) result(cfl)
    type(space_t) :: Xh
    type(coef_t) :: coef
    integer :: nelv, gdim
    real(kind=dp) :: dt
    type(c_ptr), intent(in) :: u_d, v_d, w_d
    real(kind=dp) :: cfl
    real(kind=rp) :: cfl_rp
    real(kind=xp) :: cfl_xp

#ifdef HAVE_HIP
    cfl = hip_cfl(dt, u_d, v_d, w_d, &
         coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
         coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
         coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
         Xh%dr_inv_d, Xh%ds_inv_d, Xh%dt_inv_d, &
         coef%jacinv_d, nelv, Xh%lx)
#elif HAVE_CUDA
    cfl = cuda_cfl(dt, u_d, v_d, w_d, &
         coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
         coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
         coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
         Xh%dr_inv_d, Xh%ds_inv_d, Xh%dt_inv_d, &
         coef%jacinv_d, nelv, Xh%lx)
#elif HAVE_OPENCL
    cfl_xp = opencl_cfl(real(dt, kind=xp), u_d, v_d, w_d, &
         coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
         coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
         coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
         Xh%dr_inv_d, Xh%ds_inv_d, Xh%dt_inv_d, &
         coef%jacinv_d, nelv, Xh%lx)
    cfl = real(cfl_xp, kind=dp)
#elif HAVE_METAL
    cfl_rp = metal_cfl(real(dt, kind=rp), u_d, v_d, w_d, &
         coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
         coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
         coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
         Xh%dr_inv_d, Xh%ds_inv_d, Xh%dt_inv_d, &
         coef%jacinv_d, nelv, Xh%lx)
    cfl = real(cfl_rp, kind=dp)
#else
    cfl = 0.0_dp
    call neko_error('No device backend configured')
#endif
  end function opr_device_cfl

  subroutine opr_device_rotate_cyc(vx_d, vy_d, vz_d, idir, coef)
    type(c_ptr), intent(inout) :: vx_d, vy_d, vz_d
    integer, intent(in) :: idir
    type(coef_t), intent(in) :: coef
    integer :: ncyc

    ncyc = coef%cyc_msk(0) - 1

    if (ncyc .le. 0) return

#ifdef HAVE_HIP
    call hip_rotate_cyc(vx_d, vy_d, vz_d, &
         coef%dof%x_d, coef%dof%y_d, coef%dof%z_d, &
         coef%cyc_msk_d, coef%R11_d, coef%R12_d, &
         ncyc, idir)
#elif HAVE_CUDA
    call cuda_rotate_cyc(vx_d, vy_d, vz_d, &
         coef%dof%x_d, coef%dof%y_d, coef%dof%z_d, &
         coef%cyc_msk_d, coef%R11_d, coef%R12_d, &
         ncyc, idir)
#elif HAVE_OPENCL
    call opencl_rotate_cyc(vx_d, vy_d, vz_d, &
         coef%dof%x_d, coef%dof%y_d, coef%dof%z_d, &
         coef%cyc_msk_d, coef%R11_d, coef%R12_d, &
         ncyc, idir)
#elif HAVE_METAL
    call metal_rotate_cyc(vx_d, vy_d, vz_d, &
         coef%dof%x_d, coef%dof%y_d, coef%dof%z_d, &
         coef%cyc_msk_d, coef%R11_d, coef%R12_d, &
         ncyc, idir)
#else
    call neko_error('No device backend configured for rotate_cyc')
#endif
  end subroutine opr_device_rotate_cyc

  subroutine opr_device_set_convect_rst(cr_d, cs_d, ct_d, cx_d, cy_d, cz_d, &
       Xh, coef)
    type(space_t), intent(inout) :: Xh
    type(coef_t), intent(inout) :: coef
    type(c_ptr), intent(inout) :: cr_d, cs_d, ct_d, cx_d, cy_d, cz_d

#ifdef HAVE_HIP
    call hip_set_convect_rst(cr_d, cs_d, ct_d, cx_d, cy_d, cz_d, &
         coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
         coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
         coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
         Xh%w3_d, coef%msh%nelv, Xh%lx)
#elif HAVE_CUDA
    call cuda_set_convect_rst(cr_d, cs_d, ct_d, cx_d, cy_d, cz_d, &
         coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
         coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
         coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
         Xh%w3_d, coef%msh%nelv, Xh%lx)
#elif HAVE_OPENCL
    call opencl_set_convect_rst(cr_d, cs_d, ct_d, cx_d, cy_d, cz_d, &
         coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
         coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
         coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
         Xh%w3_d, coef%msh%nelv, Xh%lx)
#elif HAVE_METAL
    call metal_set_convect_rst(cr_d, cs_d, ct_d, cx_d, cy_d, cz_d, &
         coef%drdx_d, coef%dsdx_d, coef%dtdx_d, &
         coef%drdy_d, coef%dsdy_d, coef%dtdy_d, &
         coef%drdz_d, coef%dsdz_d, coef%dtdz_d, &
         Xh%w3_d, coef%msh%nelv, Xh%lx)
#else
    call neko_error('No device backend configured')
#endif

  end subroutine opr_device_set_convect_rst

end module opr_device
