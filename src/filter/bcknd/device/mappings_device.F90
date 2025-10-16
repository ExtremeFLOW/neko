! Copyright (c) 2024, The Neko Authors
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
!> Device implementations of the mapping functions.
module mappings_device
  use num_types, only: rp, c_rp
  use, intrinsic :: iso_c_binding, only: c_ptr
  use utils, only: neko_error

  use cuda_mappings, only: &
       cuda_smooth_step, cuda_step_function, cuda_permeability
  use hip_mappings, only: &
       hip_smooth_step, hip_step_function, hip_permeability
  use opencl_mappings, only: &
       opencl_smooth_step, opencl_step_function, opencl_permeability

  implicit none

contains

  ! ========================================================================== !
  ! Internal functions and subroutines
  ! ========================================================================== !

  !> @brief Apply a smooth step function to a scalar.
  subroutine smooth_step_device(x, edge0, edge1, n)
    type(c_ptr), intent(inout) :: x
    real(kind=rp), intent(in) :: edge0, edge1
    integer, intent(in) :: n

#if HAVE_HIP == 1
    call hip_smooth_step(x, edge0, edge1, n)
#elif HAVE_CUDA == 1
    call cuda_smooth_step(x, edge0, edge1, n)
#elif HAVE_OPENCL == 1
    call opencl_smooth_step(x, edge0, edge1, n)
#else
    call neko_error(&
         "Smooth step function not implemented for the current device.")
#endif

  end subroutine smooth_step_device

  !> @brief Apply a step function to a scalar.
  subroutine step_function_device(x, x_step, value0, value1, n)
    type(c_ptr), intent(inout) :: x
    real(kind=rp), intent(in) :: x_step, value0, value1
    integer, intent(in) :: n

#if HAVE_HIP == 1
    call hip_step_function(x, x_step, value0, value1, n)
#elif HAVE_CUDA == 1
    call cuda_step_function(x, x_step, value0, value1, n)
#elif HAVE_OPENCL == 1
    call opencl_step_function(x, x_step, value0, value1, n)
#else
    call neko_error(&
         "Step function not implemented for the current device.")
#endif

  end subroutine step_function_device

  !> @brief Apply a permeability function to a scalar.
  subroutine permeability_device(x, k_0, k_1, q, n)
    type(c_ptr), intent(inout) :: x
    real(kind=rp), intent(in) :: k_0, k_1, q
    integer, intent(in) :: n

#if HAVE_HIP == 1
    call hip_permeability(x, k_0, k_1, q, n)
#elif HAVE_CUDA == 1
    call cuda_permeability(x, k_0, k_1, q, n)
#elif HAVE_OPENCL == 1
    call opencl_permeability(x, k_0, k_1, q, n)
#else
    call neko_error(&
         "Permeability function not implemented for the current device.")
#endif

  end subroutine permeability_device


end module mappings_device
