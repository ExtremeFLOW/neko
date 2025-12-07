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
!> Device backend for entropy viscosity regularization
module entropy_viscosity_device
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int
  use num_types, only: rp, c_rp
  use utils, only: neko_error
  implicit none
  private

#ifdef HAVE_HIP
  interface
     subroutine hip_entropy_visc_compute_residual(entropy_residual_d, &
          S_d, S_lag1_d, S_lag2_d, S_lag3_d, &
          bdf1, bdf2, bdf3, bdf4, dt, n) &
          bind(c, name = 'hip_entropy_visc_compute_residual')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: entropy_residual_d
       type(c_ptr), value :: S_d, S_lag1_d, S_lag2_d, S_lag3_d
       real(c_rp) :: bdf1, bdf2, bdf3, bdf4, dt
       integer(c_int) :: n
     end subroutine hip_entropy_visc_compute_residual
  end interface

  interface
     subroutine hip_entropy_visc_compute_viscosity(reg_coeff_d, &
          entropy_residual_d, h_d, c_entropy, n_S, n) &
          bind(c, name = 'hip_entropy_visc_compute_viscosity')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: reg_coeff_d, entropy_residual_d, h_d
       real(c_rp) :: c_entropy, n_S
       integer(c_int) :: n
     end subroutine hip_entropy_visc_compute_viscosity
  end interface

  interface
     subroutine hip_entropy_visc_apply_element_max(reg_coeff_d, lx, nelv) &
          bind(c, name = 'hip_entropy_visc_apply_element_max')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: reg_coeff_d
       integer(c_int) :: lx, nelv
     end subroutine hip_entropy_visc_apply_element_max
  end interface

  interface
     subroutine hip_entropy_visc_clamp_to_low_order(reg_coeff_d, &
          h_d, max_wave_speed_d, c_max, n) &
          bind(c, name = 'hip_entropy_visc_clamp_to_low_order')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: reg_coeff_d, h_d, max_wave_speed_d
       real(c_rp) :: c_max
       integer(c_int) :: n
     end subroutine hip_entropy_visc_clamp_to_low_order
  end interface

  interface
     subroutine hip_entropy_visc_smooth_divide(reg_coeff_d, &
          temp_field_d, mult_field_d, n) &
          bind(c, name = 'hip_entropy_visc_smooth_divide')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: reg_coeff_d, temp_field_d, mult_field_d
       integer(c_int) :: n
     end subroutine hip_entropy_visc_smooth_divide
  end interface

#elif HAVE_CUDA
  interface
     subroutine cuda_entropy_visc_compute_residual(entropy_residual_d, &
          S_d, S_lag1_d, S_lag2_d, S_lag3_d, &
          bdf1, bdf2, bdf3, bdf4, dt, n) &
          bind(c, name = 'cuda_entropy_visc_compute_residual')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: entropy_residual_d
       type(c_ptr), value :: S_d, S_lag1_d, S_lag2_d, S_lag3_d
       real(c_rp) :: bdf1, bdf2, bdf3, bdf4, dt
       integer(c_int) :: n
     end subroutine cuda_entropy_visc_compute_residual
  end interface

  interface
     subroutine cuda_entropy_visc_compute_viscosity(reg_coeff_d, &
          entropy_residual_d, h_d, c_entropy, n_S, n) &
          bind(c, name = 'cuda_entropy_visc_compute_viscosity')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: reg_coeff_d, entropy_residual_d, h_d
       real(c_rp) :: c_entropy, n_S
       integer(c_int) :: n
     end subroutine cuda_entropy_visc_compute_viscosity
  end interface

  interface
     subroutine cuda_entropy_visc_apply_element_max(reg_coeff_d, lx, nelv) &
          bind(c, name = 'cuda_entropy_visc_apply_element_max')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: reg_coeff_d
       integer(c_int) :: lx, nelv
     end subroutine cuda_entropy_visc_apply_element_max
  end interface

  interface
     subroutine cuda_entropy_visc_clamp_to_low_order(reg_coeff_d, &
          h_d, max_wave_speed_d, c_max, n) &
          bind(c, name = 'cuda_entropy_visc_clamp_to_low_order')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: reg_coeff_d, h_d, max_wave_speed_d
       real(c_rp) :: c_max
       integer(c_int) :: n
     end subroutine cuda_entropy_visc_clamp_to_low_order
  end interface

  interface
     subroutine cuda_entropy_visc_smooth_divide(reg_coeff_d, &
          temp_field_d, mult_field_d, n) &
          bind(c, name = 'cuda_entropy_visc_smooth_divide')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: reg_coeff_d, temp_field_d, mult_field_d
       integer(c_int) :: n
     end subroutine cuda_entropy_visc_smooth_divide
  end interface

#elif HAVE_OPENCL
  interface
     subroutine opencl_entropy_visc_compute_residual(entropy_residual_d, &
          S_d, S_lag1_d, S_lag2_d, S_lag3_d, &
          bdf1, bdf2, bdf3, bdf4, dt, n) &
          bind(c, name = 'opencl_entropy_visc_compute_residual')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: entropy_residual_d
       type(c_ptr), value :: S_d, S_lag1_d, S_lag2_d, S_lag3_d
       real(c_rp), value :: bdf1, bdf2, bdf3, bdf4, dt
       integer(c_int), value :: n
     end subroutine opencl_entropy_visc_compute_residual
  end interface

  interface
     subroutine opencl_entropy_visc_compute_viscosity(reg_coeff_d, &
          entropy_residual_d, h_d, c_entropy, n_S, n) &
          bind(c, name = 'opencl_entropy_visc_compute_viscosity')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: reg_coeff_d, entropy_residual_d, h_d
       real(c_rp), value :: c_entropy, n_S
       integer(c_int), value :: n
     end subroutine opencl_entropy_visc_compute_viscosity
  end interface

  interface
     subroutine opencl_entropy_visc_apply_element_max(reg_coeff_d, lx, nelv) &
          bind(c, name = 'opencl_entropy_visc_apply_element_max')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: reg_coeff_d
       integer(c_int), value :: lx, nelv
     end subroutine opencl_entropy_visc_apply_element_max
  end interface

  interface
     subroutine opencl_entropy_visc_clamp_to_low_order(reg_coeff_d, &
          h_d, max_wave_speed_d, c_max, n) &
          bind(c, name = 'opencl_entropy_visc_clamp_to_low_order')
       use, intrinsic :: iso_c_binding
       import c_rp
       type(c_ptr), value :: reg_coeff_d, h_d, max_wave_speed_d
       real(c_rp), value :: c_max
       integer(c_int), value :: n
     end subroutine opencl_entropy_visc_clamp_to_low_order
  end interface

  interface
     subroutine opencl_entropy_visc_smooth_divide(reg_coeff_d, &
          temp_field_d, mult_field_d, n) &
          bind(c, name = 'opencl_entropy_visc_smooth_divide')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: reg_coeff_d, temp_field_d, mult_field_d
       integer(c_int), value :: n
     end subroutine opencl_entropy_visc_smooth_divide
  end interface
#endif

  public :: entropy_viscosity_compute_residual_device, &
            entropy_viscosity_compute_viscosity_device, &
            entropy_viscosity_apply_element_max_device, &
            entropy_viscosity_clamp_to_low_order_device, &
            entropy_viscosity_smooth_divide_device

contains

  !> Compute entropy residual on device
  subroutine entropy_viscosity_compute_residual_device(entropy_residual_d, &
       S_d, S_lag1_d, S_lag2_d, S_lag3_d, bdf_coeffs, dt, n)
    type(c_ptr), intent(in) :: entropy_residual_d
    type(c_ptr), intent(in) :: S_d, S_lag1_d, S_lag2_d, S_lag3_d
    real(kind=rp), intent(in) :: bdf_coeffs(4)
    real(kind=rp), intent(in) :: dt
    integer, intent(in) :: n

#ifdef HAVE_HIP
    call hip_entropy_visc_compute_residual(entropy_residual_d, &
         S_d, S_lag1_d, S_lag2_d, S_lag3_d, &
         bdf_coeffs(1), bdf_coeffs(2), bdf_coeffs(3), bdf_coeffs(4), dt, n)
#elif HAVE_CUDA
    call cuda_entropy_visc_compute_residual(entropy_residual_d, &
         S_d, S_lag1_d, S_lag2_d, S_lag3_d, &
         bdf_coeffs(1), bdf_coeffs(2), bdf_coeffs(3), bdf_coeffs(4), dt, n)
#elif HAVE_OPENCL
    call opencl_entropy_visc_compute_residual(entropy_residual_d, &
         S_d, S_lag1_d, S_lag2_d, S_lag3_d, &
         bdf_coeffs(1), bdf_coeffs(2), bdf_coeffs(3), bdf_coeffs(4), dt, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine entropy_viscosity_compute_residual_device

  !> Compute viscosity from entropy residual on device
  subroutine entropy_viscosity_compute_viscosity_device(reg_coeff_d, &
       entropy_residual_d, h_d, c_entropy, n_S, n)
    type(c_ptr), intent(in) :: reg_coeff_d, entropy_residual_d, h_d
    real(kind=rp), intent(in) :: c_entropy, n_S
    integer, intent(in) :: n

#ifdef HAVE_HIP
    call hip_entropy_visc_compute_viscosity(reg_coeff_d, &
         entropy_residual_d, h_d, c_entropy, n_S, n)
#elif HAVE_CUDA
    call cuda_entropy_visc_compute_viscosity(reg_coeff_d, &
         entropy_residual_d, h_d, c_entropy, n_S, n)
#elif HAVE_OPENCL
    call opencl_entropy_visc_compute_viscosity(reg_coeff_d, &
         entropy_residual_d, h_d, c_entropy, n_S, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine entropy_viscosity_compute_viscosity_device

  !> Apply element-wise maximum on device
  subroutine entropy_viscosity_apply_element_max_device(reg_coeff_d, lx, nelv)
    type(c_ptr), intent(in) :: reg_coeff_d
    integer, intent(in) :: lx, nelv

#ifdef HAVE_HIP
    call hip_entropy_visc_apply_element_max(reg_coeff_d, lx, nelv)
#elif HAVE_CUDA
    call cuda_entropy_visc_apply_element_max(reg_coeff_d, lx, nelv)
#elif HAVE_OPENCL
    call opencl_entropy_visc_apply_element_max(reg_coeff_d, lx, nelv)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine entropy_viscosity_apply_element_max_device

  !> Clamp regularization coefficient to low-order viscosity on device
  subroutine entropy_viscosity_clamp_to_low_order_device(reg_coeff_d, &
       h_d, max_wave_speed_d, c_max, n)
    type(c_ptr), intent(in) :: reg_coeff_d, h_d, max_wave_speed_d
    real(kind=rp), intent(in) :: c_max
    integer, intent(in) :: n

#ifdef HAVE_HIP
    call hip_entropy_visc_clamp_to_low_order(reg_coeff_d, &
         h_d, max_wave_speed_d, c_max, n)
#elif HAVE_CUDA
    call cuda_entropy_visc_clamp_to_low_order(reg_coeff_d, &
         h_d, max_wave_speed_d, c_max, n)
#elif HAVE_OPENCL
    call opencl_entropy_visc_clamp_to_low_order(reg_coeff_d, &
         h_d, max_wave_speed_d, c_max, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine entropy_viscosity_clamp_to_low_order_device

  !> Divide by multiplicity for smoothing on device
  subroutine entropy_viscosity_smooth_divide_device(reg_coeff_d, &
       temp_field_d, mult_field_d, n)
    type(c_ptr), intent(in) :: reg_coeff_d, temp_field_d, mult_field_d
    integer, intent(in) :: n

#ifdef HAVE_HIP
    call hip_entropy_visc_smooth_divide(reg_coeff_d, &
         temp_field_d, mult_field_d, n)
#elif HAVE_CUDA
    call cuda_entropy_visc_smooth_divide(reg_coeff_d, &
         temp_field_d, mult_field_d, n)
#elif HAVE_OPENCL
    call opencl_entropy_visc_smooth_divide(reg_coeff_d, &
         temp_field_d, mult_field_d, n)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine entropy_viscosity_smooth_divide_device

end module entropy_viscosity_device

