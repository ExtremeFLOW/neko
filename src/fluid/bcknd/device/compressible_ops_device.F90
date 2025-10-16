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
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int
  use num_types, only: rp, c_rp
  use field, only: field_t
  use utils, only: neko_error
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
#endif

  public :: compressible_ops_device_compute_max_wave_speed, &
            compressible_ops_device_compute_entropy

contains

  !> Compute maximum wave speed for compressible flows on device
  subroutine compressible_ops_device_compute_max_wave_speed(max_wave_speed, u, v, w, gamma, p, rho, n)
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: gamma
    type(field_t), intent(inout) :: max_wave_speed
    type(field_t), intent(in) :: u, v, w, p, rho

#ifdef HAVE_HIP
    call hip_compute_max_wave_speed(max_wave_speed%x_d, u%x_d, v%x_d, w%x_d, gamma, p%x_d, rho%x_d, n)
#elif HAVE_CUDA
    call cuda_compute_max_wave_speed(max_wave_speed%x_d, u%x_d, v%x_d, w%x_d, gamma, p%x_d, rho%x_d, n)
#elif HAVE_OPENCL
    call opencl_compute_max_wave_speed(max_wave_speed%x_d, u%x_d, v%x_d, w%x_d, gamma, p%x_d, rho%x_d, n)
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
#else
    call neko_error('No device backend configured')
#endif
  end subroutine compressible_ops_device_compute_entropy

end module compressible_ops_device
