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
!> Implements the device kernel for the `wall_model_t` type.
module wall_model_device
  use num_types, only : rp, c_rp
  use, intrinsic :: iso_c_binding, only : c_ptr
  use utils, only : neko_error
  implicit none
  private

#ifdef HAVE_HIP
  interface
     subroutine hip_wall_model_compute_mag_field(tau_x_d, tau_y_d, tau_z_d, &
          tau_field_d, msk_d, m) &
          bind(c, name = 'hip_wall_model_compute_mag_field')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       implicit none
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d, tau_field_d, msk_d
       integer(c_int) :: m
     end subroutine hip_wall_model_compute_mag_field
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_wall_model_compute_mag_field(tau_x_d, tau_y_d, tau_z_d, &
          tau_field_d, msk_d, m) &
          bind(c, name = 'cuda_wall_model_compute_mag_field')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       implicit none
       type(c_ptr), value :: tau_x_d, tau_y_d, tau_z_d, tau_field_d, msk_d
       integer(c_int) :: m
     end subroutine cuda_wall_model_compute_mag_field
  end interface
#elif HAVE_OPENCL
#endif
  public :: wall_model_compute_mag_field_device

contains
  !> Compute the wall shear stress's magnitude on device.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine wall_model_compute_mag_field_device(tau_x_d, tau_y_d, tau_z_d, &
       tau_field_d, msk_d, m)
    type(c_ptr), intent(in) :: tau_x_d, tau_y_d, tau_z_d
    type(c_ptr), intent(in) :: msk_d
    type(c_ptr), intent(inout) :: tau_field_d
    integer, intent(in) :: m

#if HAVE_HIP
    call hip_wall_model_compute_mag_field(tau_x_d, tau_y_d, tau_z_d, &
         tau_field_d, msk_d, m)
#elif HAVE_CUDA
    call cuda_wall_model_compute_mag_field(tau_x_d, tau_y_d, tau_z_d, &
         tau_field_d, msk_d, m)
#elif HAVE_OPENCL
    call neko_error("OPENCL is not implemented for &
    &wall_model_compute_mag_field")
#else
    call neko_error('No device backend configured')
#endif

  end subroutine wall_model_compute_mag_field_device
end module wall_model_device
