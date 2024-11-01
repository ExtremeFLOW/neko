! Copyright (c) 2022-2024, The Neko Authors
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
module device_local_interpolation
  use num_types, only : rp, c_rp
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int
  implicit none
  private

  public :: device_find_rst_legendre

#ifdef HAVE_HIP
  interface
     subroutine hip_coef_generate_geo(G11, G12, G13, G22, G23, G33, &
          drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
          jacinv, w3, nel, lx, gdim) &
          bind(c, name='hip_coef_generate_geo')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: G11, G12, G13, G22, G23, G33
       type(c_ptr), value :: drdx, drdy, drdz
       type(c_ptr), value :: dsdx, dsdy, dsdz
       type(c_ptr), value :: dtdx, dtdy, dtdz
       type(c_ptr), value :: jacinv, w3
       integer(c_int) :: nel, gdim, lx
     end subroutine hip_coef_generate_geo
  end interface

  interface
     subroutine hip_coef_generate_dxyzdrst(drdx, drdy, drdz, dsdx, dsdy, &
          dsdz, dtdx, dtdy, dtdz, dxdr, dydr, dzdr, dxds, dyds, dzds, dxdt, &
          dydt, dzdt, dx, dy, dz, x, y, z, jacinv, jac, lx, nel) &
          bind(c, name='hip_coef_generate_dxyzdrst')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: drdx, drdy, drdz
       type(c_ptr), value :: dsdx, dsdy, dsdz
       type(c_ptr), value :: dtdx, dtdy, dtdz
       type(c_ptr), value :: dxdr, dydr, dzdr
       type(c_ptr), value :: dxds, dyds, dzds
       type(c_ptr), value :: dxdt, dydt, dzdt
       type(c_ptr), value :: dx, dy, dz, x, y, z
       type(c_ptr), value :: jacinv, jac
       integer(c_int) :: lx, nel
     end subroutine hip_coef_generate_dxyzdrst
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_find_rst_legendre(rst, pt_x, pt_y, pt_z, &
                                       x_hat, y_hat, z_hat, &
                                       resx, resy, resz, &
                                       lx, el_ids, n_pt, tol, conv_pts) &
          bind(c, name='cuda_find_rst_legendre')
       use, intrinsic :: iso_c_binding
       use num_types 
       implicit none
       type(c_ptr), value :: rst
       type(c_ptr), value :: pt_x, pt_y, pt_z
       type(c_ptr), value :: x_hat, y_hat, z_hat
       type(c_ptr), value :: resx, resy, resz
       type(c_ptr), value :: el_ids, conv_pts
       integer(c_int) :: lx, n_pt
       real(c_rp) :: tol
     end subroutine cuda_find_rst_legendre
  end interface
#elif HAVE_OPENCL
  interface
     subroutine opencl_coef_generate_geo(G11, G12, G13, G22, G23, G33, &
          drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
          jacinv, w3, nel, lx, gdim) &
          bind(c, name='opencl_coef_generate_geo')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: G11, G12, G13, G22, G23, G33
       type(c_ptr), value :: drdx, drdy, drdz
       type(c_ptr), value :: dsdx, dsdy, dsdz
       type(c_ptr), value :: dtdx, dtdy, dtdz
       type(c_ptr), value :: jacinv, w3
       integer(c_int) :: nel, gdim, lx
     end subroutine opencl_coef_generate_geo
  end interface

  interface
     subroutine opencl_coef_generate_dxyzdrst(drdx, drdy, drdz, dsdx, dsdy, &
          dsdz, dtdx, dtdy, dtdz, dxdr, dydr, dzdr, dxds, dyds, dzds, dxdt, &
          dydt, dzdt, dx, dy, dz, x, y, z, jacinv, jac, lx, nel) &
          bind(c, name='opencl_coef_generate_dxyzdrst')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: drdx, drdy, drdz
       type(c_ptr), value :: dsdx, dsdy, dsdz
       type(c_ptr), value :: dtdx, dtdy, dtdz
       type(c_ptr), value :: dxdr, dydr, dzdr
       type(c_ptr), value :: dxds, dyds, dzds
       type(c_ptr), value :: dxdt, dydt, dzdt
       type(c_ptr), value :: dx, dy, dz, x, y, z
       type(c_ptr), value :: jacinv, jac
       integer(c_int) :: lx, nel
     end subroutine opencl_coef_generate_dxyzdrst
  end interface
#endif

contains

  subroutine device_find_rst_legendre(rst_d, pt_x_d, pt_y_d, pt_z_d, &
                                      x_hat_d, y_hat_d, z_hat_d, &
                                      resx_d, resy_d, resz_d, &
                                      lx, el_ids_d, n_pts, tol, conv_pts_d)
    type(c_ptr) :: rst_d, pt_x_d, pt_y_d, pt_z_d
    type(c_ptr) :: x_hat_d, y_hat_d, z_hat_d
    type(c_ptr) :: resx_d, resy_d, resz_d
    type(c_ptr) :: el_ids_d, conv_pts_d
    integer :: lx, n_pts
    real(kind=c_rp) :: tol

#ifdef HAVE_HIP
    call hip_coef_generate_geo(G11_d, G12_d, G13_d, G22_d, G23_d, &
         G33_d, drdx_d, drdy_d, drdz_d, dsdx_d, dsdy_d, dsdz_d, &
         dtdx_d, dtdy_d, dtdz_d, jacinv_d, w3_d, nel, lx, gdim)
#elif HAVE_CUDA
     call cuda_find_rst_legendre(rst_d, pt_x_d, pt_y_d, pt_z_d, &
                                       x_hat_d, y_hat_d, z_hat_d, &
                                       resx_d, resy_d, resz_d, &
                                       lx, el_ids_d, n_pts, tol, conv_pts_d)
#elif HAVE_OPENCL
    call opencl_coef_generate_geo(G11_d, G12_d, G13_d, G22_d, G23_d, &
         G33_d, drdx_d, drdy_d, drdz_d, dsdx_d, dsdy_d, dsdz_d, &
         dtdx_d, dtdy_d, dtdz_d, jacinv_d, w3_d, nel, lx, gdim)
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_find_rst_legendre

end module device_local_interpolation
