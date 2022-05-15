! Copyright (c) 2022, The Neko Authors
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
module device_coef
  use num_types
  use utils
  use, intrinsic :: iso_c_binding
  implicit none

#ifdef HAVE_HIP
    interface
     subroutine hip_coef_generate_geo(G11, G12, G13, G22, G23, G33, &
          drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
          jacinv, w3, nel, lx, gdim) &
          bind(c, name='hip_coef_generate_geo')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: G11, G12, G13, G22, G23, G33
       type(c_ptr), value :: drdx, drdy, drdz
       type(c_ptr), value :: dsdx, dsdy, dsdz 
       type(c_ptr), value :: dtdx, dtdy, dtdz
       type(c_ptr), value :: jacinv, w3
       integer(c_int) :: nel, gdim, lx
     end subroutine hip_coef_generate_geo
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_coef_generate_geo(G11, G12, G13, G22, G23, G33, &
          drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
          jacinv, w3, nel, lx, gdim) &
          bind(c, name='cuda_coef_generate_geo')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: G11, G12, G13, G22, G23, G33
       type(c_ptr), value :: drdx, drdy, drdz
       type(c_ptr), value :: dsdx, dsdy, dsdz 
       type(c_ptr), value :: dtdx, dtdy, dtdz
       type(c_ptr), value :: jacinv, w3
       integer(c_int) :: nel, gdim, lx
     end subroutine cuda_coef_generate_geo
  end interface

  interface
     subroutine cuda_coef_generate_dxyzdrst(drdx, drdy, drdz, dsdx, dsdy, &
          dsdz, dtdx, dtdy, dtdz, dxdr, dydr, dzdr, dxds, dyds, dzds, dxdt, &
          dydt, dzdt, dx, dy, dz, x, y, z, jacinv, jac, lx, nel) &
          bind(c, name='cuda_coef_generate_dxyzdrst')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: drdx, drdy, drdz
       type(c_ptr), value :: dsdx, dsdy, dsdz
       type(c_ptr), value :: dtdx, dtdy, dtdz
       type(c_ptr), value :: dxdr, dydr, dzdr
       type(c_ptr), value :: dxds, dyds, dzds
       type(c_ptr), value :: dxdt, dydt, dzdt
       type(c_ptr), value :: dx, dy, dz, x, y, z
       type(c_ptr), value :: jacinv, jac
       integer(c_int) :: lx, nelf
     end subroutine cuda_coef_generate_dxyzdrst
  end interface
#elif HAVE_OPENCL
  interface
     subroutine opencl_coef_generate_geo(G11, G12, G13, G22, G23, G33, &
          drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
          jacinv, w3, nel, lx, gdim) &
          bind(c, name='opencl_coef_generate_geo')
       use, intrinsic :: iso_c_binding
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
       type(c_ptr), value :: drdx, drdy, drdz
       type(c_ptr), value :: dsdx, dsdy, dsdz
       type(c_ptr), value :: dtdx, dtdy, dtdz
       type(c_ptr), value :: dxdr, dydr, dzdr
       type(c_ptr), value :: dxds, dyds, dzds
       type(c_ptr), value :: dxdt, dydt, dzdt
       type(c_ptr), value :: dx, dy, dz, x, y, z
       type(c_ptr), value :: jacinv, jac
       integer(c_int) :: lx, nelf
     end subroutine opencl_coef_generate_dxyzdrst
  end interface
#endif

contains

  subroutine device_coef_generate_geo(G11_d, G12_d, G13_d, G22_d, &
       G23_d, G33_d, drdx_d, drdy_d, drdz_d, dsdx_d, dsdy_d, &
       dsdz_d, dtdx_d, dtdy_d, dtdz_d, jacinv_d, w3_d, nel, lx, gdim)
    type(c_ptr) :: G11_d, G12_d, G13_d, G22_d, G23_d, G33_d
    type(c_ptr) :: drdx_d, drdy_d, drdz_d
    type(c_ptr) :: dsdx_d, dsdy_d, dsdz_d
    type(c_ptr) :: dtdx_d, dtdy_d, dtdz_d
    type(c_ptr) :: jacinv_d, w3_d
    integer :: nel, gdim, lx

#ifdef HAVE_HIP
    call hip_coef_generate_geo(G11_d, G12_d, G13_d, G22_d, G23_d, &
         G33_d, drdx_d, drdy_d, drdz_d, dsdx_d, dsdy_d, dsdz_d, &
         dtdx_d, dtdy_d, dtdz_d, jacinv_d, w3_d, nel, lx, gdim)
#elif HAVE_CUDA
    call cuda_coef_generate_geo(G11_d, G12_d, G13_d, G22_d, G23_d, &
         G33_d, drdx_d, drdy_d, drdz_d, dsdx_d, dsdy_d, dsdz_d, &
         dtdx_d, dtdy_d, dtdz_d, jacinv_d, w3_d, nel, lx, gdim)
#elif HAVE_OPENCL
    call opencl_coef_generate_geo(G11_d, G12_d, G13_d, G22_d, G23_d, &
         G33_d, drdx_d, drdy_d, drdz_d, dsdx_d, dsdy_d, dsdz_d, &
         dtdx_d, dtdy_d, dtdz_d, jacinv_d, w3_d, nel, lx, gdim)
#else
    call neko_error('No device backend configured')
#endif

  end subroutine device_coef_generate_geo

  subroutine device_coef_generate_dxydrst(drdx_d, drdy_d, drdz_d, dsdx_d, dsdy_d,&
       dsdz_d, dtdx_d, dtdy_d, dtdz_d, dxdr_d, dydr_d, dzdr_d, dxds_d, &
       dyds_d, dzds_d, dxdt_d, dydt_d, dzdt_d, dx_d, dy_d, dz_d, &
       x_d, y_d, z_d, jacinv_d, jac_d, lx, nel)
    type(c_ptr) :: drdx_d, drdy_d, drdz_d
    type(c_ptr) :: dsdx_d, dsdy_d, dsdz_d
    type(c_ptr) :: dtdx_d, dtdy_d, dtdz_d
    type(c_ptr) :: dxdr_d, dydr_d, dzdr_d
    type(c_ptr) :: dxds_d, dyds_d, dzds_d
    type(c_ptr) :: dxdt_d, dydt_d, dzdt_d
    type(c_ptr) :: dx_d, dy_d, dz_d, x_d, y_d, z_d
    type(c_ptr) :: jacinv_d, jac_d
    integer :: lx, nel

#ifdef HAVE_HIP
#elif HAVE_CUDA
    call cuda_coef_generate_dxyzdrst(drdx_d, drdy_d, drdz_d, dsdx_d, &
         dsdy_d, dsdz_d, dtdx_d, dtdy_d, dtdz_d, dxdr_d, dydr_d, &
         dzdr_d, dxds_d, dyds_d, dzds_d, dxdt_d, dydt_d, dzdt_d, &
         dx_d, dy_d, dz_d, x_d, y_d, z_d, jacinv_d, jac_d, lx, nel)
#elif HAVE_OPENCL
    call opencl_coef_generate_dxyzdrst(drdx_d, drdy_d, drdz_d, dsdx_d, &
         dsdy_d, dsdz_d, dtdx_d, dtdy_d, dtdz_d, dxdr_d, dydr_d, &
         dzdr_d, dxds_d, dyds_d, dzds_d, dxdt_d, dydt_d, dzdt_d, &
         dx_d, dy_d, dz_d, x_d, y_d, z_d, jacinv_d, jac_d, lx, nel)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_coef_generate_dxydrst

end module device_coef
