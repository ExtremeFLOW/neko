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
!> Implements the device kernel for the `reichardt_t` type.
module reichardt_device
  use num_types, only : rp, c_rp
  use, intrinsic :: iso_c_binding, only : c_ptr
  use utils, only : neko_error
  implicit none
  private

#ifdef HAVE_HIP
  interface
     subroutine hip_reichardt_compute(u_d, v_d, w_d, &
          n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, &
          tau_x_d, tau_y_d, tau_z_d, n_nodes, &
          kappa, C, B1, B2, tstep) &
          bind(c, name = 'hip_reichardt_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d, rho_w_d
       type(c_ptr), value :: n_x_d, n_y_d, n_z_d, h_d, nu_d
       real(c_rp) :: kappa, C, B1, B2
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d
       integer(c_int) :: n_nodes, tstep
     end subroutine hip_reichardt_compute
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_reichardt_compute(u_d, v_d, w_d, &
          n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, &
          tau_x_d, tau_y_d, tau_z_d, n_nodes, &
          kappa, C, B1, B2, tstep) &
          bind(c, name = 'cuda_reichardt_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d, rho_w_d
       type(c_ptr), value :: n_x_d, n_y_d, n_z_d, h_d, nu_d
       real(c_rp) :: kappa, C, B1, B2
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d
       integer(c_int) :: n_nodes, tstep
     end subroutine cuda_reichardt_compute
  end interface
#elif HAVE_OPENCL
  interface
     subroutine opencl_reichardt_compute(u_d, v_d, w_d, &
          n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, &
          tau_x_d, tau_y_d, tau_z_d, n_nodes, &
          kappa, C, B1, B2, tstep) &
          bind(c, name = 'opencl_reichardt_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d, rho_w_d
       type(c_ptr), value :: n_x_d, n_y_d, n_z_d, h_d, nu_d
       real(c_rp) :: kappa, C, B1, B2
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d
       integer(c_int) :: n_nodes, tstep
     end subroutine opencl_reichardt_compute
  end interface
#elif HAVE_METAL
  interface
     subroutine metal_reichardt_compute(u_d, v_d, w_d, &
          n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, &
          tau_x_d, tau_y_d, tau_z_d, n_nodes, &
          kappa, C, B1, B2, tstep) &
          bind(c, name = 'metal_reichardt_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       use num_types, only : c_rp
       implicit none
       type(c_ptr), value :: u_d, v_d, w_d, rho_w_d
       type(c_ptr), value :: n_x_d, n_y_d, n_z_d, h_d, nu_d
       real(c_rp) :: kappa, C, B1, B2
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d
       integer(c_int) :: n_nodes, tstep
     end subroutine metal_reichardt_compute
  end interface
#endif
  public :: reichardt_compute_device

contains
  !> Compute the wall shear stress on device using Reichardt's law.
  !! @param u_d The x component of the sampled velocity.
  !! @param v_d The y component of the sampled velocity.
  !! @param w_d The z component of the sampled velocity.
  !! @param n_x_d The x component of the wall normal.
  !! @param n_y_d The y component of the wall normal.
  !! @param n_z_d The z component of the wall normal.
  !! @param nu_d The kinematic viscosity at the wall.
  !! @param rho_w_d The density at the wall.
  !! @param h_d The wall-normal distance of the sampling point.
  !! @param tau_x_d The x component of the wall shear stress.
  !! @param tau_y_d The y component of the wall shear stress.
  !! @param tau_z_d The z component of the wall shear stress.
  !! @param n_nodes The number of wall nodes.
  !! @param kappa The von Karman coefficient.
  !! @param C The amplitude of the exponential correction.
  !! @param B1 The damping length scale, in wall units.
  !! @param B2 The decay length scale of the second exponential term, in
  !! wall units.
  !! @param tstep The current time-step.
  subroutine reichardt_compute_device(u_d, v_d, w_d, &
       n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, tau_x_d, tau_y_d, tau_z_d, &
       n_nodes, kappa, C, B1, B2, tstep)
    integer, intent(in) :: n_nodes, tstep
    type(c_ptr), intent(in) :: u_d, v_d, w_d, rho_w_d
    type(c_ptr), intent(in) :: n_x_d, n_y_d, n_z_d, h_d, nu_d
    type(c_ptr), intent(inout) :: tau_x_d, tau_y_d, tau_z_d
    real(kind=rp), intent(in) :: kappa, C, B1, B2

#if HAVE_HIP
    call hip_reichardt_compute(u_d, v_d, w_d, &
         n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, &
         tau_x_d, tau_y_d, tau_z_d, n_nodes, kappa, C, B1, B2, &
         tstep)
#elif HAVE_CUDA
    call cuda_reichardt_compute(u_d, v_d, w_d, &
         n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, &
         tau_x_d, tau_y_d, tau_z_d, n_nodes, kappa, C, B1, B2, &
         tstep)
#elif HAVE_OPENCL
    call opencl_reichardt_compute(u_d, v_d, w_d, &
         n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, &
         tau_x_d, tau_y_d, tau_z_d, n_nodes, kappa, C, B1, B2, &
         tstep)
#elif HAVE_METAL
    call metal_reichardt_compute(u_d, v_d, w_d, &
         n_x_d, n_y_d, n_z_d, nu_d, rho_w_d, h_d, &
         tau_x_d, tau_y_d, tau_z_d, n_nodes, kappa, C, B1, B2, &
         tstep)
#else
    call neko_error('No device backend configured')
#endif

  end subroutine reichardt_compute_device
end module reichardt_device
