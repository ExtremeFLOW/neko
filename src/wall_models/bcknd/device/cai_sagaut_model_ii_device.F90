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
!> Device dispatch for `cai_sagaut_model_ii_t`.
module cai_sagaut_model_ii_device
  use num_types, only : rp, c_rp
  use, intrinsic :: iso_c_binding, only : c_ptr
  use utils, only : neko_error
  implicit none
  private

#ifdef HAVE_HIP
  interface
     subroutine hip_cai_sagaut_model_ii_compute(u_d, v_d, w_d, n_x_d, &
          n_y_d, n_z_d, nu_d, rho_w_d, h_d, tau_x_d, tau_y_d, tau_z_d, &
          n_nodes, kappa, B, p, s) &
          bind(c, name = 'hip_cai_sagaut_model_ii_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d, rho_w_d
       type(c_ptr), value :: n_x_d, n_y_d, n_z_d, h_d, nu_d
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d
       real(c_rp) :: kappa, B, p, s
       integer(c_int) :: n_nodes
     end subroutine hip_cai_sagaut_model_ii_compute
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_cai_sagaut_model_ii_compute(u_d, v_d, w_d, n_x_d, &
          n_y_d, n_z_d, nu_d, rho_w_d, h_d, tau_x_d, tau_y_d, tau_z_d, &
          n_nodes, kappa, B, p, s) &
          bind(c, name = 'cuda_cai_sagaut_model_ii_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d, rho_w_d
       type(c_ptr), value :: n_x_d, n_y_d, n_z_d, h_d, nu_d
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d
       real(c_rp) :: kappa, B, p, s
       integer(c_int) :: n_nodes
     end subroutine cuda_cai_sagaut_model_ii_compute
  end interface
#elif HAVE_OPENCL
  interface
     subroutine opencl_cai_sagaut_model_ii_compute(u_d, v_d, w_d, n_x_d, &
          n_y_d, n_z_d, nu_d, rho_w_d, h_d, tau_x_d, tau_y_d, tau_z_d, &
          n_nodes, kappa, B, p, s) &
          bind(c, name = 'opencl_cai_sagaut_model_ii_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d, rho_w_d
       type(c_ptr), value :: n_x_d, n_y_d, n_z_d, h_d, nu_d
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d
       real(c_rp) :: kappa, B, p, s
       integer(c_int) :: n_nodes
     end subroutine opencl_cai_sagaut_model_ii_compute
  end interface
#endif
  public :: cai_sagaut_model_ii_compute_device

contains
  !> Evaluate the device wall-model kernel for Model-II.
  !! @param u_d The sampled x-velocity field on the device.
  !! @param v_d The sampled y-velocity field on the device.
  !! @param w_d The sampled z-velocity field on the device.
  !! @param u_d Sampled x velocity.
  !! @param v_d Sampled y velocity.
  !! @param w_d Sampled z velocity.
  !! @param n_x_d The x-component of the wall normals.
  !! @param n_y_d The y-component of the wall normals.
  !! @param n_z_d The z-component of the wall normals.
  !! @param nu_d The sampled kinematic viscosity at wall points.
  !! @param rho_w_d The sampled density at wall points.
  !! @param h_d The wall-model sampling distances.
  !! @param tau_x_d The x-component of the wall shear stress.
  !! @param tau_y_d The y-component of the wall shear stress.
  !! @param tau_z_d The z-component of the wall shear stress.
  !! @param n_nodes The number of wall points.
  !! @param lx The number of GLL points per direction.
  !! @param kappa The von Karman coefficient.
  !! @param B The log-law intercept.
  !! @param p The blending exponent.
  !! @param s The blending scale.
  subroutine cai_sagaut_model_ii_compute_device(u_d, v_d, w_d, n_x_d, &
       n_y_d, n_z_d, nu_d, rho_w_d, h_d, tau_x_d, tau_y_d, tau_z_d, &
       n_nodes, kappa, B, p, s)
    integer, intent(in) :: n_nodes
    type(c_ptr), intent(in) :: u_d, v_d, w_d, rho_w_d
    type(c_ptr), intent(in) :: n_x_d, n_y_d, n_z_d, h_d, nu_d
    type(c_ptr), intent(inout) :: tau_x_d, tau_y_d, tau_z_d
    real(kind=rp), intent(in) :: kappa, B, p, s

#if HAVE_HIP
    call hip_cai_sagaut_model_ii_compute(u_d, v_d, w_d, n_x_d, n_y_d, &
         n_z_d, nu_d, rho_w_d, h_d, tau_x_d, tau_y_d, tau_z_d, n_nodes, &
         kappa, B, p, s)
#elif HAVE_CUDA
    call cuda_cai_sagaut_model_ii_compute(u_d, v_d, w_d, n_x_d, n_y_d, &
         n_z_d, nu_d, rho_w_d, h_d, tau_x_d, tau_y_d, tau_z_d, n_nodes, &
         kappa, B, p, s)
#elif HAVE_OPENCL
    call opencl_cai_sagaut_model_ii_compute(u_d, v_d, w_d, n_x_d, n_y_d, &
         n_z_d, nu_d, rho_w_d, h_d, tau_x_d, tau_y_d, tau_z_d, n_nodes, &
         kappa, B, p, s)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine cai_sagaut_model_ii_compute_device
end module cai_sagaut_model_ii_device
