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
!> Implements the device kernels for the `duprat_t` type.
module duprat_device
  use num_types, only : rp, c_rp
  use, intrinsic :: iso_c_binding, only : c_ptr
  use utils, only : neko_error
  implicit none
  private

#ifdef HAVE_HIP
  interface
     subroutine hip_duprat_compute( &
          u_d, v_d, w_d, n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, dpds_d, &
          tau_x_d, tau_y_d, tau_z_d, alpha_d, n_nodes, kappa, beta, A, tstep) &
          bind(c, name = 'hip_duprat_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d, n_x_d, n_y_d
       type(c_ptr), value :: n_z_d, nu_d, rho_w_d, h_d, dpds_d
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d, alpha_d
       real(c_rp) :: kappa, beta, A
       integer(c_int) :: n_nodes, tstep
     end subroutine hip_duprat_compute
  end interface
  interface
     subroutine hip_duprat_update_dpds( &
          dpx_d, dpy_d, dpz_d, u_d, v_d, w_d, n_x_d, n_y_d, n_z_d, nu_d, &
          rho_w_d, h_d, eps, dpds_d, n_nodes) &
          bind(c, name = 'hip_duprat_update_dpds')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: dpx_d, dpy_d, dpz_d, u_d, v_d
       type(c_ptr), value :: w_d, n_x_d, n_y_d, n_z_d, nu_d
       type(c_ptr), value :: rho_w_d, h_d
       type(c_ptr), value :: dpds_d
       real(c_rp) :: eps
       integer(c_int) :: n_nodes
     end subroutine hip_duprat_update_dpds
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_duprat_compute( &
          u_d, v_d, w_d, n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, dpds_d, &
          tau_x_d, tau_y_d, tau_z_d, alpha_d, n_nodes, kappa, beta, A, tstep) &
          bind(c, name = 'cuda_duprat_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d, n_x_d, n_y_d
       type(c_ptr), value :: n_z_d, nu_d, rho_w_d, h_d, dpds_d
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d, alpha_d
       real(c_rp) :: kappa, beta, A
       integer(c_int) :: n_nodes, tstep
     end subroutine cuda_duprat_compute
  end interface
  interface
     subroutine cuda_duprat_update_dpds( &
          dpx_d, dpy_d, dpz_d, u_d, v_d, w_d, n_x_d, n_y_d, n_z_d, nu_d, &
          rho_w_d, h_d, eps, dpds_d, n_nodes) &
          bind(c, name = 'cuda_duprat_update_dpds')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: dpx_d, dpy_d, dpz_d, u_d, v_d
       type(c_ptr), value :: w_d, n_x_d, n_y_d, n_z_d, nu_d
       type(c_ptr), value :: rho_w_d, h_d
       type(c_ptr), value :: dpds_d
       real(c_rp) :: eps
       integer(c_int) :: n_nodes
     end subroutine cuda_duprat_update_dpds
  end interface
#elif HAVE_OPENCL
  interface
     subroutine opencl_duprat_compute( &
          u_d, v_d, w_d, n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, dpds_d, &
          tau_x_d, tau_y_d, tau_z_d, alpha_d, n_nodes, kappa, beta, A, tstep) &
          bind(c, name = 'opencl_duprat_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d, n_x_d, n_y_d
       type(c_ptr), value :: n_z_d, nu_d, rho_w_d, h_d, dpds_d
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d, alpha_d
       real(c_rp) :: kappa, beta, A
       integer(c_int) :: n_nodes, tstep
     end subroutine opencl_duprat_compute
  end interface
  interface
     subroutine opencl_duprat_update_dpds( &
          dpx_d, dpy_d, dpz_d, u_d, v_d, w_d, n_x_d, n_y_d, n_z_d, nu_d, &
          rho_w_d, h_d, eps, dpds_d, n_nodes) &
          bind(c, name = 'opencl_duprat_update_dpds')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: dpx_d, dpy_d, dpz_d, u_d, v_d
       type(c_ptr), value :: w_d, n_x_d, n_y_d, n_z_d, nu_d
       type(c_ptr), value :: rho_w_d, h_d
       type(c_ptr), value :: dpds_d
       real(c_rp) :: eps
       integer(c_int) :: n_nodes
     end subroutine opencl_duprat_update_dpds
  end interface
#elif HAVE_METAL
  interface
     subroutine metal_duprat_compute( &
          u_d, v_d, w_d, n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, dpds_d, &
          tau_x_d, tau_y_d, tau_z_d, alpha_d, n_nodes, kappa, beta, A, tstep) &
          bind(c, name = 'metal_duprat_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d, n_x_d, n_y_d
       type(c_ptr), value :: n_z_d, nu_d, rho_w_d, h_d, dpds_d
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d, alpha_d
       real(c_rp) :: kappa, beta, A
       integer(c_int) :: n_nodes, tstep
     end subroutine metal_duprat_compute
  end interface
  interface
     subroutine metal_duprat_update_dpds( &
          dpx_d, dpy_d, dpz_d, u_d, v_d, w_d, n_x_d, n_y_d, n_z_d, nu_d, &
          rho_w_d, h_d, eps, dpds_d, n_nodes) &
          bind(c, name = 'metal_duprat_update_dpds')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: dpx_d, dpy_d, dpz_d, u_d, v_d
       type(c_ptr), value :: w_d, n_x_d, n_y_d, n_z_d, nu_d
       type(c_ptr), value :: rho_w_d, h_d
       type(c_ptr), value :: dpds_d
       real(c_rp) :: eps
       integer(c_int) :: n_nodes
     end subroutine metal_duprat_update_dpds
  end interface
#endif
  public :: duprat_compute_device, duprat_update_dpds_device

contains
  !> Compute the wall shear stress on device using the model of Duprat et al.
  !! @param u_d The x component of the sampled velocity.
  !! @param v_d The y component of the sampled velocity.
  !! @param w_d The z component of the sampled velocity.
  !! @param n_x_d The x component of the wall normal.
  !! @param n_y_d The y component of the wall normal.
  !! @param n_z_d The z component of the wall normal.
  !! @param nu_d The kinematic viscosity at the wall.
  !! @param rho_w_d The density at the wall.
  !! @param h_d The wall-normal distance of the sampling point.
  !! @param dpds_d The wall-tangential pressure gradient.
  !! @param tau_x_d The x component of the wall shear stress.
  !! @param tau_y_d The y component of the wall shear stress.
  !! @param tau_z_d The z component of the wall shear stress.
  !! @param alpha_d The ratio u_tau^2 / u_tau_p^2.
  !! @param n_nodes The number of wall nodes.
  !! @param kappa The von Karman coefficient.
  !! @param beta The exponent of the pressure-gradient term.
  !! @param A The damping constant.
  !! @param tstep The current time-step.
  subroutine duprat_compute_device( &
       u_d, v_d, w_d, n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, dpds_d, &
       tau_x_d, tau_y_d, tau_z_d, alpha_d, n_nodes, kappa, beta, A, tstep)
    integer, intent(in) :: n_nodes, tstep
    type(c_ptr), intent(in) :: u_d, v_d, w_d, n_x_d, n_y_d, n_z_d
    type(c_ptr), intent(in) :: nu_d, rho_w_d, h_d, dpds_d
    type(c_ptr), intent(inout) :: tau_x_d, tau_y_d, tau_z_d, alpha_d
    real(kind=rp), intent(in) :: kappa, beta, A

#if HAVE_HIP
    call hip_duprat_compute( &
         u_d, v_d, w_d, n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, dpds_d, &
         tau_x_d, tau_y_d, tau_z_d, alpha_d, n_nodes, kappa, beta, A, tstep)
#elif HAVE_CUDA
    call cuda_duprat_compute( &
         u_d, v_d, w_d, n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, dpds_d, &
         tau_x_d, tau_y_d, tau_z_d, alpha_d, n_nodes, kappa, beta, A, tstep)
#elif HAVE_OPENCL
    call opencl_duprat_compute( &
         u_d, v_d, w_d, n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, dpds_d, &
         tau_x_d, tau_y_d, tau_z_d, alpha_d, n_nodes, kappa, beta, A, tstep)
#elif HAVE_METAL
    call metal_duprat_compute( &
         u_d, v_d, w_d, n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, dpds_d, &
         tau_x_d, tau_y_d, tau_z_d, alpha_d, n_nodes, kappa, beta, A, tstep)
#else
    call neko_error('No device backend configured')
#endif

  end subroutine duprat_compute_device

  !> Update the filtered wall-tangential pressure gradient on device.
  !! @param dpx_d The x component of the sampled pressure gradient.
  !! @param dpy_d The y component of the sampled pressure gradient.
  !! @param dpz_d The z component of the sampled pressure gradient.
  !! @param u_d The x component of the sampled velocity.
  !! @param v_d The y component of the sampled velocity.
  !! @param w_d The z component of the sampled velocity.
  !! @param n_x_d The x component of the wall normal.
  !! @param n_y_d The y component of the wall normal.
  !! @param n_z_d The z component of the wall normal.
  !! @param nu_d The kinematic viscosity at the wall.
  !! @param rho_w_d The density at the wall.
  !! @param h_d The wall-normal distance of the sampling point.
  !! @param eps The filter weight of the current gradient.
  !! @param dpds_d The filtered pressure gradient, updated in place.
  !! @param n_nodes The number of wall nodes.
  subroutine duprat_update_dpds_device( &
       dpx_d, dpy_d, dpz_d, u_d, v_d, w_d, n_x_d, n_y_d, n_z_d, nu_d, &
       rho_w_d, h_d, eps, dpds_d, n_nodes)
    integer, intent(in) :: n_nodes
    type(c_ptr), intent(in) :: dpx_d, dpy_d, dpz_d, u_d, v_d, w_d
    type(c_ptr), intent(in) :: n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d
    type(c_ptr), intent(inout) :: dpds_d
    real(kind=rp), intent(in) :: eps

#if HAVE_HIP
    call hip_duprat_update_dpds( &
         dpx_d, dpy_d, dpz_d, u_d, v_d, w_d, n_x_d, n_y_d, n_z_d, nu_d, &
         rho_w_d, h_d, eps, dpds_d, n_nodes)
#elif HAVE_CUDA
    call cuda_duprat_update_dpds( &
         dpx_d, dpy_d, dpz_d, u_d, v_d, w_d, n_x_d, n_y_d, n_z_d, nu_d, &
         rho_w_d, h_d, eps, dpds_d, n_nodes)
#elif HAVE_OPENCL
    call opencl_duprat_update_dpds( &
         dpx_d, dpy_d, dpz_d, u_d, v_d, w_d, n_x_d, n_y_d, n_z_d, nu_d, &
         rho_w_d, h_d, eps, dpds_d, n_nodes)
#elif HAVE_METAL
    call metal_duprat_update_dpds( &
         dpx_d, dpy_d, dpz_d, u_d, v_d, w_d, n_x_d, n_y_d, n_z_d, nu_d, &
         rho_w_d, h_d, eps, dpds_d, n_nodes)
#else
    call neko_error('No device backend configured')
#endif

  end subroutine duprat_update_dpds_device
end module duprat_device
