! Copyright (c) 2025, The Neko Authors
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
!> Device implementation of compressible flow operations
module compressible_ops_device
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int
  use num_types, only : rp, c_rp
  use field, only : field_t
  use utils, only : neko_error
  implicit none
  private

#ifdef HAVE_HIP
  interface
     subroutine hip_compute_max_wave_speed(max_wave_speed_d, u_d, v_d, w_d, &
          gamma, p_d, rho_d, n) &
          bind(c, name = 'hip_compute_max_wave_speed')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: max_wave_speed_d, u_d, v_d, w_d, p_d, rho_d
       real(c_rp) :: gamma
       integer(c_int) :: n
     end subroutine hip_compute_max_wave_speed
  end interface

  interface
     subroutine hip_compute_entropy(S_d, p_d, rho_d, gamma, n) &
          bind(c, name = 'hip_compute_entropy')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: S_d, p_d, rho_d
       real(c_rp) :: gamma
       integer(c_int) :: n
     end subroutine hip_compute_entropy
  end interface

  interface
     subroutine hip_update_uvw(u_d, v_d, w_d, m_x_d, &
          m_y_d, m_z_d, rho_d, n) &
          bind(c, name = 'hip_update_uvw')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: u_d, v_d, w_d, m_x_d, m_y_d, m_z_d, rho_d
       integer(c_int) :: n
     end subroutine hip_update_uvw
  end interface

  interface
     subroutine hip_update_mxyz_p_ruvw(m_x_d, m_y_d, m_z_d, &
          p_d, ruvw_d, u_d, v_d, w_d, E_d, rho_d, gamma, n) &
          bind(c, name = 'hip_update_mxyz_p_ruvw')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: u_d, v_d, w_d, m_x_d, m_y_d, m_z_d, rho_d
       type(c_ptr), value :: p_d, ruvw_d, E_d
       real(c_rp) :: gamma
       integer(c_int) :: n
     end subroutine hip_update_mxyz_p_ruvw
  end interface

  interface
     subroutine hip_update_e(E_d, p_d, ruvw_d, gamma, n) &
       bind(c, name = 'hip_update_e')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: p_d, E_d, ruvw_d
       real(c_rp) :: gamma
       integer(c_int) :: n
     end subroutine hip_update_e
  end interface

  interface
     subroutine hip_update_temperature(T_d, p_d, rho_d, gamma, n) &
       bind(c, name = 'hip_update_temperature')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: T_d, p_d, rho_d
       real(c_rp) :: gamma
       integer(c_int) :: n
     end subroutine hip_update_temperature
  end interface

  interface
     subroutine hip_ns_flux_prepare(div_flux_d, dissipation_d, h1_d, &
          dudx_d, dudy_d, dudz_d, dvdx_d, dvdy_d, dvdz_d, &
          dwdx_d, dwdy_d, dwdz_d, mu_d, n) &
          bind(c, name = 'hip_ns_flux_prepare')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: div_flux_d, dissipation_d, h1_d
       type(c_ptr), value :: dudx_d, dudy_d, dudz_d
       type(c_ptr), value :: dvdx_d, dvdy_d, dvdz_d
       type(c_ptr), value :: dwdx_d, dwdy_d, dwdz_d, mu_d
       integer(c_int) :: n
     end subroutine hip_ns_flux_prepare
  end interface

  interface
     subroutine hip_ns_flux_finalize(visc_m_x_d, visc_m_y_d, visc_m_z_d, &
          visc_E_d, f_x_d, f_y_d, f_z_d, opgrad_x_d, opgrad_y_d, &
          opgrad_z_d, u_d, v_d, w_d, B_d, dissipation_d, n) &
          bind(c, name = 'hip_ns_flux_finalize')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: visc_m_x_d, visc_m_y_d, visc_m_z_d, visc_E_d
       type(c_ptr), value :: f_x_d, f_y_d, f_z_d
       type(c_ptr), value :: opgrad_x_d, opgrad_y_d, opgrad_z_d
       type(c_ptr), value :: u_d, v_d, w_d, B_d, dissipation_d
       integer(c_int) :: n
     end subroutine hip_ns_flux_finalize
  end interface

  interface
     subroutine hip_ns_flux_temperature(div_flux_d, h1_d, p_d, rho_d, &
          kappa_d, gamma, n) bind(c, name = 'hip_ns_flux_temperature')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: div_flux_d, h1_d, p_d, rho_d, kappa_d
       real(c_rp) :: gamma
       integer(c_int) :: n
     end subroutine hip_ns_flux_temperature
  end interface

#elif HAVE_CUDA
  interface
     subroutine cuda_compute_max_wave_speed(max_wave_speed_d, u_d, v_d, w_d, &
          gamma, p_d, rho_d, n) &
          bind(c, name = 'cuda_compute_max_wave_speed')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: max_wave_speed_d, u_d, v_d, w_d, p_d, rho_d
       real(c_rp) :: gamma
       integer(c_int) :: n
     end subroutine cuda_compute_max_wave_speed
  end interface

  interface
     subroutine cuda_compute_entropy(S_d, p_d, rho_d, gamma, n) &
          bind(c, name = 'cuda_compute_entropy')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: S_d, p_d, rho_d
       real(c_rp) :: gamma
       integer(c_int) :: n
     end subroutine cuda_compute_entropy
  end interface

  interface
     subroutine cuda_update_uvw(u_d, v_d, w_d, m_x_d, &
          m_y_d, m_z_d, rho_d, n) &
          bind(c, name = 'cuda_update_uvw')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: u_d, v_d, w_d, m_x_d, m_y_d, m_z_d, rho_d
       integer(c_int) :: n
     end subroutine cuda_update_uvw
  end interface

  interface
     subroutine cuda_update_mxyz_p_ruvw(m_x_d, m_y_d, m_z_d, &
          p_d, ruvw_d, u_d, v_d, w_d, E_d, rho_d, gamma, n) &
          bind(c, name = 'cuda_update_mxyz_p_ruvw')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: u_d, v_d, w_d, m_x_d, m_y_d, m_z_d, rho_d
       type(c_ptr), value :: p_d, ruvw_d, E_d
       real(c_rp) :: gamma
       integer(c_int) :: n
     end subroutine cuda_update_mxyz_p_ruvw
  end interface

  interface
     subroutine cuda_update_e(E_d, p_d, ruvw_d, gamma, n) &
       bind(c, name = 'cuda_update_e')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: p_d, E_d, ruvw_d
       real(c_rp) :: gamma
       integer(c_int) :: n
     end subroutine cuda_update_e
  end interface

  interface
     subroutine cuda_update_temperature(T_d, p_d, rho_d, gamma, n) &
       bind(c, name = 'cuda_update_temperature')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: T_d, p_d, rho_d
       real(c_rp) :: gamma
       integer(c_int) :: n
     end subroutine cuda_update_temperature
  end interface

  interface
     subroutine cuda_ns_flux_prepare(div_flux_d, dissipation_d, h1_d, &
          dudx_d, dudy_d, dudz_d, dvdx_d, dvdy_d, dvdz_d, &
          dwdx_d, dwdy_d, dwdz_d, mu_d, n) &
          bind(c, name = 'cuda_ns_flux_prepare')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: div_flux_d, dissipation_d, h1_d
       type(c_ptr), value :: dudx_d, dudy_d, dudz_d
       type(c_ptr), value :: dvdx_d, dvdy_d, dvdz_d
       type(c_ptr), value :: dwdx_d, dwdy_d, dwdz_d, mu_d
       integer(c_int) :: n
     end subroutine cuda_ns_flux_prepare
  end interface

  interface
     subroutine cuda_ns_flux_finalize(visc_m_x_d, visc_m_y_d, visc_m_z_d, &
          visc_E_d, f_x_d, f_y_d, f_z_d, opgrad_x_d, opgrad_y_d, &
          opgrad_z_d, u_d, v_d, w_d, B_d, dissipation_d, n) &
          bind(c, name = 'cuda_ns_flux_finalize')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: visc_m_x_d, visc_m_y_d, visc_m_z_d, visc_E_d
       type(c_ptr), value :: f_x_d, f_y_d, f_z_d
       type(c_ptr), value :: opgrad_x_d, opgrad_y_d, opgrad_z_d
       type(c_ptr), value :: u_d, v_d, w_d, B_d, dissipation_d
       integer(c_int) :: n
     end subroutine cuda_ns_flux_finalize
  end interface

  interface
     subroutine cuda_ns_flux_temperature(div_flux_d, h1_d, p_d, rho_d, &
          kappa_d, gamma, n) bind(c, name = 'cuda_ns_flux_temperature')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: div_flux_d, h1_d, p_d, rho_d, kappa_d
       real(c_rp) :: gamma
       integer(c_int) :: n
     end subroutine cuda_ns_flux_temperature
  end interface

#elif HAVE_OPENCL
  interface
     subroutine opencl_compute_max_wave_speed(max_wave_speed_d, u_d, v_d, w_d, &
          gamma, p_d, rho_d, n) &
          bind(c, name = 'opencl_compute_max_wave_speed')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: max_wave_speed_d, u_d, v_d, w_d, p_d, rho_d
       real(c_rp), value :: gamma
       integer(c_int), value :: n
     end subroutine opencl_compute_max_wave_speed
  end interface

  interface
     subroutine opencl_compute_entropy(S_d, p_d, rho_d, gamma, n) &
          bind(c, name = 'opencl_compute_entropy')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: S_d, p_d, rho_d
       real(c_rp), value :: gamma
       integer(c_int), value :: n
     end subroutine opencl_compute_entropy
  end interface

  interface
     subroutine opencl_update_uvw(u_d, v_d, w_d, m_x_d, &
          m_y_d, m_z_d, rho_d, n) &
          bind(c, name = 'opencl_update_uvw')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: u_d, v_d, w_d, m_x_d, m_y_d, m_z_d, rho_d
       integer(c_int), value :: n
     end subroutine opencl_update_uvw
  end interface

  interface
     subroutine opencl_update_mxyz_p_ruvw(m_x_d, m_y_d, m_z_d, &
          p_d, ruvw_d, u_d, v_d, w_d, E_d, rho_d, gamma, n) &
          bind(c, name = 'opencl_update_mxyz_p_ruvw')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: u_d, v_d, w_d, m_x_d, m_y_d, m_z_d, rho_d
       type(c_ptr), value :: p_d, ruvw_d, E_d
       real(c_rp), value :: gamma
       integer(c_int), value :: n
     end subroutine opencl_update_mxyz_p_ruvw
  end interface

  interface
     subroutine opencl_update_e(E_d, p_d, ruvw_d, gamma, n) &
       bind(c, name = 'opencl_update_e')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: p_d, E_d, ruvw_d
       real(c_rp), value :: gamma
       integer(c_int), value :: n
     end subroutine opencl_update_e
  end interface
  interface
     subroutine opencl_update_temperature(T_d, p_d, rho_d, gamma, n) &
       bind(c, name = 'opencl_update_temperature')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: T_d, p_d, rho_d
       real(c_rp), value :: gamma
       integer(c_int), value :: n
     end subroutine opencl_update_temperature
  end interface

  interface
     subroutine opencl_ns_flux_prepare(div_flux_d, dissipation_d, h1_d, &
          dudx_d, dudy_d, dudz_d, dvdx_d, dvdy_d, dvdz_d, &
          dwdx_d, dwdy_d, dwdz_d, mu_d, n) &
          bind(c, name = 'opencl_ns_flux_prepare')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: div_flux_d, dissipation_d, h1_d
       type(c_ptr), value :: dudx_d, dudy_d, dudz_d
       type(c_ptr), value :: dvdx_d, dvdy_d, dvdz_d
       type(c_ptr), value :: dwdx_d, dwdy_d, dwdz_d, mu_d
       integer(c_int), value :: n
     end subroutine opencl_ns_flux_prepare
  end interface

  interface
     subroutine opencl_ns_flux_finalize(visc_m_x_d, visc_m_y_d, visc_m_z_d, &
          visc_E_d, f_x_d, f_y_d, f_z_d, opgrad_x_d, opgrad_y_d, &
          opgrad_z_d, u_d, v_d, w_d, B_d, dissipation_d, n) &
          bind(c, name = 'opencl_ns_flux_finalize')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: visc_m_x_d, visc_m_y_d, visc_m_z_d, visc_E_d
       type(c_ptr), value :: f_x_d, f_y_d, f_z_d
       type(c_ptr), value :: opgrad_x_d, opgrad_y_d, opgrad_z_d
       type(c_ptr), value :: u_d, v_d, w_d, B_d, dissipation_d
       integer(c_int), value :: n
     end subroutine opencl_ns_flux_finalize
  end interface

  interface
     subroutine opencl_ns_flux_temperature(div_flux_d, h1_d, p_d, rho_d, &
          kappa_d, gamma, n) bind(c, name = 'opencl_ns_flux_temperature')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: div_flux_d, h1_d, p_d, rho_d, kappa_d
       real(c_rp), value :: gamma
       integer(c_int), value :: n
     end subroutine opencl_ns_flux_temperature
  end interface
#elif HAVE_METAL
  interface
     subroutine metal_compute_max_wave_speed(max_wave_speed_d, u_d, v_d, w_d, &
          gamma, p_d, rho_d, n) &
          bind(c, name = 'metal_compute_max_wave_speed')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: max_wave_speed_d, u_d, v_d, w_d, p_d, rho_d
       real(c_rp), value :: gamma
       integer(c_int), value :: n
     end subroutine metal_compute_max_wave_speed
  end interface

  interface
     subroutine metal_compute_entropy(S_d, p_d, rho_d, gamma, n) &
          bind(c, name = 'metal_compute_entropy')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: S_d, p_d, rho_d
       real(c_rp), value :: gamma
       integer(c_int), value :: n
     end subroutine metal_compute_entropy
  end interface

  interface
     subroutine metal_update_uvw(u_d, v_d, w_d, m_x_d, &
          m_y_d, m_z_d, rho_d, n) &
          bind(c, name = 'metal_update_uvw')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: u_d, v_d, w_d, m_x_d, m_y_d, m_z_d, rho_d
       integer(c_int), value :: n
     end subroutine metal_update_uvw
  end interface

  interface
     subroutine metal_update_mxyz_p_ruvw(m_x_d, m_y_d, m_z_d, &
          p_d, ruvw_d, u_d, v_d, w_d, E_d, rho_d, gamma, n) &
          bind(c, name = 'metal_update_mxyz_p_ruvw')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: u_d, v_d, w_d, m_x_d, m_y_d, m_z_d, rho_d
       type(c_ptr), value :: p_d, ruvw_d, E_d
       real(c_rp), value :: gamma
       integer(c_int), value :: n
     end subroutine metal_update_mxyz_p_ruvw
  end interface

  interface
     subroutine metal_update_e(E_d, p_d, ruvw_d, gamma, n) &
       bind(c, name = 'metal_update_e')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: p_d, E_d, ruvw_d
       real(c_rp), value :: gamma
       integer(c_int), value :: n
     end subroutine metal_update_e
  end interface
#endif

  public :: compressible_ops_device_compute_max_wave_speed, &
       compressible_ops_device_compute_entropy, &
       compressible_ops_device_update_uvw, &
       compressible_ops_device_update_mxyz_p_ruvw,&
       compressible_ops_device_update_e, &
       compressible_ops_device_update_temperature, &
       compressible_ops_device_ns_flux_prepare, &
       compressible_ops_device_ns_flux_finalize, &
       compressible_ops_device_ns_flux_temperature

contains

  !> Compute maximum wave speed for compressible flows on device
  subroutine compressible_ops_device_compute_max_wave_speed(max_wave_speed, &
       u, v, w, gamma, p, rho, n)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: gamma
    type(field_t), intent(inout) :: max_wave_speed
    type(field_t), intent(in) :: u, v, w, p, rho

#ifdef HAVE_HIP
    call hip_compute_max_wave_speed(max_wave_speed%x_d, u%x_d, v%x_d, &
         w%x_d, gamma, p%x_d, rho%x_d, n)
#elif HAVE_CUDA
    call cuda_compute_max_wave_speed(max_wave_speed%x_d, u%x_d, v%x_d, &
         w%x_d, gamma, p%x_d, rho%x_d, n)
#elif HAVE_OPENCL
    call opencl_compute_max_wave_speed(max_wave_speed%x_d, u%x_d, v%x_d, &
         w%x_d, gamma, p%x_d, rho%x_d, n)
#elif HAVE_METAL
    call metal_compute_max_wave_speed(max_wave_speed%x_d, u%x_d, v%x_d, &
         w%x_d, gamma, p%x_d, rho%x_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine compressible_ops_device_compute_max_wave_speed

  !> Compute entropy field S = 1/(gamma-1) * rho * (log(p) - gamma * log(rho)) on device
  subroutine compressible_ops_device_compute_entropy(S, p, rho, gamma, n)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: gamma
    type(field_t), intent(inout) :: S
    type(field_t), intent(in) :: p, rho

#ifdef HAVE_HIP
    call hip_compute_entropy(S%x_d, p%x_d, rho%x_d, gamma, n)
#elif HAVE_CUDA
    call cuda_compute_entropy(S%x_d, p%x_d, rho%x_d, gamma, n)
#elif HAVE_OPENCL
    call opencl_compute_entropy(S%x_d, p%x_d, rho%x_d, gamma, n)
#elif HAVE_METAL
    call metal_compute_entropy(S%x_d, p%x_d, rho%x_d, gamma, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine compressible_ops_device_compute_entropy

  !> Update u,v,w fields
  subroutine compressible_ops_device_update_uvw(u_d, v_d, w_d, &
       m_x_d, m_y_d, m_z_d, rho_d, n)
    type(c_ptr), intent(inout) :: u_d, v_d, w_d
    type(c_ptr), intent(in) :: m_x_d, m_y_d, m_z_d, rho_d
    integer, intent(in) :: n
    integer :: i

#ifdef HAVE_HIP
    call hip_update_uvw(u_d, v_d, w_d, m_x_d, m_y_d, m_z_d, rho_d, n)
#elif HAVE_CUDA
    call cuda_update_uvw(u_d, v_d, w_d, m_x_d, m_y_d, m_z_d, rho_d, n)
#elif HAVE_OPENCL
    call opencl_update_uvw(u_d, v_d, w_d, m_x_d, m_y_d, m_z_d, rho_d, n)
#elif HAVE_METAL
    call metal_update_uvw(u_d, v_d, w_d, m_x_d, m_y_d, m_z_d, rho_d, n)
#else
    call neko_error('No device backend configured')
#endif

  end subroutine compressible_ops_device_update_uvw

  !> Update m_x, m_y, m_z, p, ruvw, fields
  subroutine compressible_ops_device_update_mxyz_p_ruvw(m_x_d, m_y_d, m_z_d, &
       p_d, ruvw_d, u_d, v_d, w_d, E_d, rho_d, gamma, n)
    integer, intent(in) :: n
    type(c_ptr), intent(inout) :: m_x_d, m_y_d, m_z_d, p_d, ruvw_d
    type(c_ptr), intent(in) :: u_d, v_d, w_d, E_d, rho_d
    real(kind=rp), intent(in) :: gamma


#ifdef HAVE_HIP
    call hip_update_mxyz_p_ruvw(m_x_d, m_y_d, m_z_d, &
         p_d, ruvw_d, u_d, v_d, w_d, E_d, rho_d, gamma, n)
#elif HAVE_CUDA
    call cuda_update_mxyz_p_ruvw(m_x_d, m_y_d, m_z_d, &
         p_d, ruvw_d, u_d, v_d, w_d, E_d, rho_d, gamma, n)
#elif HAVE_OPENCL
    call opencl_update_mxyz_p_ruvw(m_x_d, m_y_d, m_z_d, &
         p_d, ruvw_d, u_d, v_d, w_d, E_d, rho_d, gamma, n)
#elif HAVE_METAL
    call metal_update_mxyz_p_ruvw(m_x_d, m_y_d, m_z_d, &
         p_d, ruvw_d, u_d, v_d, w_d, E_d, rho_d, gamma, n)
#else
    call neko_error('No device backend configured')
#endif

  end subroutine compressible_ops_device_update_mxyz_p_ruvw

  !> Update E field
  subroutine compressible_ops_device_update_e(E_d, p_d, ruvw_d, gamma, n)
    integer, intent(in) :: n
    type(c_ptr), intent(inout) :: E_d, p_d
    ! ruvw = 0.5 * rho * (u^2 + v^2 + w^2)
    type(c_ptr), intent(in) :: ruvw_d
    real(kind=rp), intent(in) :: gamma

#ifdef HAVE_HIP
    call hip_update_e(E_d, p_d, ruvw_d, gamma, n)
#elif HAVE_CUDA
    call cuda_update_e(E_d, p_d, ruvw_d, gamma, n)
#elif HAVE_OPENCL
    call opencl_update_e(E_d, p_d, ruvw_d, gamma, n)
#elif HAVE_METAL
    call metal_update_e(E_d, p_d, ruvw_d, gamma, n)
#else
    call neko_error('No device backend configured')
#endif

  end subroutine compressible_ops_device_update_e

  !> Update temperature field.
  subroutine compressible_ops_device_update_temperature(T_d, p_d, rho_d, &
       gamma, n)
    integer, intent(in) :: n
    type(c_ptr), intent(inout) :: T_d
    type(c_ptr), intent(in) :: p_d, rho_d
    real(kind=rp), intent(in) :: gamma

#ifdef HAVE_HIP
    call hip_update_temperature(T_d, p_d, rho_d, gamma, n)
#elif HAVE_CUDA
    call cuda_update_temperature(T_d, p_d, rho_d, gamma, n)
#elif HAVE_OPENCL
    call opencl_update_temperature(T_d, p_d, rho_d, gamma, n)
#else
    call neko_error('No device backend configured')
#endif

  end subroutine compressible_ops_device_update_temperature

  !> Prepare physical Navier-Stokes flux work arrays.
  subroutine compressible_ops_device_ns_flux_prepare(div_flux_d, &
       dissipation_d, h1_d, dudx_d, dudy_d, dudz_d, dvdx_d, dvdy_d, &
       dvdz_d, dwdx_d, dwdy_d, dwdz_d, mu_d, n)
    integer, intent(in) :: n
    type(c_ptr), intent(inout) :: div_flux_d, dissipation_d, h1_d
    type(c_ptr), intent(in) :: dudx_d, dudy_d, dudz_d
    type(c_ptr), intent(in) :: dvdx_d, dvdy_d, dvdz_d
    type(c_ptr), intent(in) :: dwdx_d, dwdy_d, dwdz_d, mu_d

#ifdef HAVE_HIP
    call hip_ns_flux_prepare(div_flux_d, dissipation_d, h1_d, dudx_d, &
         dudy_d, dudz_d, dvdx_d, dvdy_d, dvdz_d, dwdx_d, dwdy_d, dwdz_d, &
         mu_d, n)
#elif HAVE_CUDA
    call cuda_ns_flux_prepare(div_flux_d, dissipation_d, h1_d, dudx_d, &
         dudy_d, dudz_d, dvdx_d, dvdy_d, dvdz_d, dwdx_d, dwdy_d, dwdz_d, &
         mu_d, n)
#elif HAVE_OPENCL
    call opencl_ns_flux_prepare(div_flux_d, dissipation_d, h1_d, dudx_d, &
         dudy_d, dudz_d, dvdx_d, dvdy_d, dvdz_d, dwdx_d, dwdy_d, dwdz_d, &
         mu_d, n)
#else
    call neko_error('No device backend configured')
#endif

  end subroutine compressible_ops_device_ns_flux_prepare

  !> Finish physical Navier-Stokes flux assembly.
  subroutine compressible_ops_device_ns_flux_finalize(visc_m_x_d, visc_m_y_d, &
       visc_m_z_d, visc_E_d, f_x_d, f_y_d, f_z_d, opgrad_x_d, opgrad_y_d, &
       opgrad_z_d, u_d, v_d, w_d, B_d, dissipation_d, n)
    integer, intent(in) :: n
    type(c_ptr), intent(inout) :: visc_m_x_d, visc_m_y_d, visc_m_z_d
    type(c_ptr), intent(inout) :: visc_E_d, f_x_d, f_y_d, f_z_d
    type(c_ptr), intent(in) :: opgrad_x_d, opgrad_y_d, opgrad_z_d
    type(c_ptr), intent(in) :: u_d, v_d, w_d, B_d, dissipation_d

#ifdef HAVE_HIP
    call hip_ns_flux_finalize(visc_m_x_d, visc_m_y_d, visc_m_z_d, &
         visc_E_d, f_x_d, f_y_d, f_z_d, opgrad_x_d, opgrad_y_d, &
         opgrad_z_d, u_d, v_d, w_d, B_d, dissipation_d, n)
#elif HAVE_CUDA
    call cuda_ns_flux_finalize(visc_m_x_d, visc_m_y_d, visc_m_z_d, &
         visc_E_d, f_x_d, f_y_d, f_z_d, opgrad_x_d, opgrad_y_d, &
         opgrad_z_d, u_d, v_d, w_d, B_d, dissipation_d, n)
#elif HAVE_OPENCL
    call opencl_ns_flux_finalize(visc_m_x_d, visc_m_y_d, visc_m_z_d, &
         visc_E_d, f_x_d, f_y_d, f_z_d, opgrad_x_d, opgrad_y_d, &
         opgrad_z_d, u_d, v_d, w_d, B_d, dissipation_d, n)
#else
    call neko_error('No device backend configured')
#endif

  end subroutine compressible_ops_device_ns_flux_finalize

  !> Prepare temperature and conductivity coefficient for energy flux.
  subroutine compressible_ops_device_ns_flux_temperature(div_flux_d, h1_d, &
       p_d, rho_d, kappa_d, gamma, n)
    integer, intent(in) :: n
    type(c_ptr), intent(inout) :: div_flux_d, h1_d
    type(c_ptr), intent(in) :: p_d, rho_d, kappa_d
    real(kind=rp), intent(in) :: gamma

#ifdef HAVE_HIP
    call hip_ns_flux_temperature(div_flux_d, h1_d, p_d, rho_d, kappa_d, &
         gamma, n)
#elif HAVE_CUDA
    call cuda_ns_flux_temperature(div_flux_d, h1_d, p_d, rho_d, kappa_d, &
         gamma, n)
#elif HAVE_OPENCL
    call opencl_ns_flux_temperature(div_flux_d, h1_d, p_d, rho_d, kappa_d, &
         gamma, n)
#else
    call neko_error('No device backend configured')
#endif

  end subroutine compressible_ops_device_ns_flux_temperature

end module compressible_ops_device
