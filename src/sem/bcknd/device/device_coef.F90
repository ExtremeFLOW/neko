module device_coef
  use num_types, only : rp, c_rp
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int
  implicit none
  private

  public :: device_coef_generate_geo
  public :: device_coef_generate_dxydrst
  public :: device_coef_generate_mass
  public :: device_coef_generate_area_and_normal

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

  interface
     subroutine hip_coef_generate_mass(B, Binv, jac, w3, lxyz, nel) &
          bind(c, name='hip_coef_generate_mass')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: B, Binv, jac, w3
       integer(c_int) :: lxyz, nel
     end subroutine hip_coef_generate_mass
  end interface

  interface
     subroutine hip_coef_generate_area_and_normal(area, nx, ny, nz, &
          dxdr, dydr, dzdr, dxds, dyds, dzds, dxdt, dydt, dzdt, &
          wx, wy, wz, lx, nel, eps) &
          bind(c, name='hip_coef_generate_area_and_normal')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: area, nx, ny, nz
       type(c_ptr), value :: dxdr, dydr, dzdr
       type(c_ptr), value :: dxds, dyds, dzds
       type(c_ptr), value :: dxdt, dydt, dzdt
       type(c_ptr), value :: wx, wy, wz
       integer(c_int) :: lx, nel
       real(kind=c_rp), value :: eps
     end subroutine hip_coef_generate_area_and_normal
  end interface

#elif HAVE_CUDA
  interface
     subroutine cuda_coef_generate_geo(G11, G12, G13, G22, G23, G33, &
          drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz, &
          jacinv, w3, nel, lx, gdim) &
          bind(c, name='cuda_coef_generate_geo')
       use, intrinsic :: iso_c_binding
       implicit none
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
     end subroutine cuda_coef_generate_dxyzdrst
  end interface

  interface
     subroutine cuda_coef_generate_mass(B, Binv, jac, w3, lxyz, nel) &
          bind(c, name='cuda_coef_generate_mass')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: B, Binv, jac, w3
       integer(c_int) :: lxyz, nel
     end subroutine cuda_coef_generate_mass
  end interface

  interface
     subroutine cuda_coef_generate_area_and_normal(area, nx, ny, nz, &
          dxdr, dydr, dzdr, dxds, dyds, dzds, dxdt, dydt, dzdt, &
          wx, wy, wz, lx, nel, eps) &
          bind(c, name='cuda_coef_generate_area_and_normal')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: area, nx, ny, nz
       type(c_ptr), value :: dxdr, dydr, dzdr
       type(c_ptr), value :: dxds, dyds, dzds
       type(c_ptr), value :: dxdt, dydt, dzdt
       type(c_ptr), value :: wx, wy, wz
       integer(c_int) :: lx, nel
       real(kind=c_rp), value :: eps
     end subroutine cuda_coef_generate_area_and_normal
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

  interface
     subroutine opencl_coef_generate_mass(B, Binv, jac, w3, lxyz, nel) &
          bind(c, name='opencl_coef_generate_mass')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: B, Binv, jac, w3
       integer(c_int) :: lxyz, nel
     end subroutine opencl_coef_generate_mass
  end interface

  interface
     subroutine opencl_coef_generate_area_and_normal(area, nx, ny, nz, &
          dxdr, dydr, dzdr, dxds, dyds, dzds, dxdt, dydt, dzdt, &
          wx, wy, wz, lx, nel, eps) &
          bind(c, name='opencl_coef_generate_area_and_normal')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: area, nx, ny, nz
       type(c_ptr), value :: dxdr, dydr, dzdr
       type(c_ptr), value :: dxds, dyds, dzds
       type(c_ptr), value :: dxdt, dydt, dzdt
       type(c_ptr), value :: wx, wy, wz
       integer(c_int) :: lx, nel
       real(kind=c_rp), value :: eps
     end subroutine opencl_coef_generate_area_and_normal
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
    call hip_coef_generate_dxyzdrst(drdx_d, drdy_d, drdz_d, dsdx_d, &
         dsdy_d, dsdz_d, dtdx_d, dtdy_d, dtdz_d, dxdr_d, dydr_d, &
         dzdr_d, dxds_d, dyds_d, dzds_d, dxdt_d, dydt_d, dzdt_d, &
         dx_d, dy_d, dz_d, x_d, y_d, z_d, jacinv_d, jac_d, lx, nel)
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

  subroutine device_coef_generate_mass(B, Binv, jac, w3, lxyz, nel)
    type(c_ptr) :: B, Binv, jac, w3
    integer :: lxyz, nel

#ifdef HAVE_HIP
    call hip_coef_generate_mass(B, Binv, jac, w3, lxyz, nel)
#elif HAVE_CUDA
    call cuda_coef_generate_mass(B, Binv, jac, w3, lxyz, nel)
#elif HAVE_OPENCL
    call opencl_coef_generate_mass(B, Binv, jac, w3, lxyz, nel)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_coef_generate_mass

  subroutine device_coef_generate_area_and_normal(area_d, nx_d, ny_d, nz_d, &
       dxdr_d, dydr_d, dzdr_d, dxds_d, dyds_d, dzds_d, dxdt_d, dydt_d, dzdt_d, &
       wx_d, wy_d, wz_d, lx, nel, eps)
    type(c_ptr) :: area_d, nx_d, ny_d, nz_d
    type(c_ptr) :: dxdr_d, dydr_d, dzdr_d
    type(c_ptr) :: dxds_d, dyds_d, dzds_d
    type(c_ptr) :: dxdt_d, dydt_d, dzdt_d
    type(c_ptr) :: wx_d, wy_d, wz_d
    integer :: lx, nel
    real(kind=c_rp) :: eps

#ifdef HAVE_HIP
    call hip_coef_generate_area_and_normal(area_d, nx_d, ny_d, nz_d, &
         dxdr_d, dydr_d, dzdr_d, dxds_d, dyds_d, dzds_d, dxdt_d, dydt_d, dzdt_d, &
         wx_d, wy_d, wz_d, lx, nel, eps)
#elif HAVE_CUDA
    call cuda_coef_generate_area_and_normal(area_d, nx_d, ny_d, nz_d, &
         dxdr_d, dydr_d, dzdr_d, dxds_d, dyds_d, dzds_d, dxdt_d, dydt_d, dzdt_d, &
         wx_d, wy_d, wz_d, lx, nel, eps)
#elif HAVE_OPENCL
    call opencl_coef_generate_area_and_normal(area_d, nx_d, ny_d, nz_d, &
         dxdr_d, dydr_d, dzdr_d, dxds_d, dyds_d, dzds_d, dxdt_d, dydt_d, dzdt_d, &
         wx_d, wy_d, wz_d, lx, nel, eps)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine device_coef_generate_area_and_normal

end module device_coef