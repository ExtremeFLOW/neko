! Copyright (c) 2025-2026, The Neko Authors
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
!> Device kernel wrapper for computing Deardorff SGS quantities.

module device_deardorff_nut
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int
  use num_types, only: rp, c_rp
  use utils, only: neko_error
  use comm, only: NEKO_COMM, pe_size, MPI_REAL_PRECISION
  use mpi_f08, only: MPI_SUM, MPI_IN_PLACE, MPI_Allreduce

  implicit none
  private

#ifdef HAVE_HIP
  interface
     subroutine hip_deardorff_nut_compute(TKE_d, &
          dTdx_d, dTdy_d, dTdz_d, &
          a11_d, a12_d, a13_d, &
          a21_d, a22_d, a23_d, &
          a31_d, a32_d, a33_d, &
          delta_d, nut_d, temperature_alphat, TKE_alphat, TKE_source, &
          c_k, T0, g1, g2, g3, &
          eps, n) bind(C,name="hip_deardorff_nut_compute")
       use, intrinsic :: iso_c_binding, only: c_ptr, c_int
       import c_rp
       type(c_ptr), value :: TKE_d, &
            dTdx_d, dTdy_d, dTdz_d, &
            a11_d, a12_d, a13_d, &
            a21_d, a22_d, a23_d, &
            a31_d, a32_d, a33_d, &
            delta_d, nut_d, temperature_alphat, &
            TKE_alphat, TKE_source
       integer(c_int) :: n
       real(c_rp) :: c_k, T0, g1, g2, g3, eps
     end subroutine hip_deardorff_nut_compute
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_deardorff_nut_compute(TKE_d, &
          dTdx_d, dTdy_d, dTdz_d, &
          a11_d, a12_d, a13_d, &
          a21_d, a22_d, a23_d, &
          a31_d, a32_d, a33_d, &
          delta_d, nut_d, temperature_alphat, TKE_alphat, TKE_source, &
          c_k, T0, g1, g2, g3, &
          eps, n) &
          bind(c, name = 'cuda_deardorff_nut_compute')
       use, intrinsic :: iso_c_binding, only: c_ptr, c_int
       import c_rp
       type(c_ptr), value :: TKE_d, &
            dTdx_d, dTdy_d, dTdz_d, &
            a11_d, a12_d, a13_d, &
            a21_d, a22_d, a23_d, &
            a31_d, a32_d, a33_d, &
            delta_d, nut_d, temperature_alphat, &
            TKE_alphat, TKE_source
       integer(c_int) :: n
       real(c_rp) :: c_k, T0, g1, g2, g3, eps
     end subroutine cuda_deardorff_nut_compute
  end interface
#elif HAVE_OPENCL
#endif

  public :: device_deardorff_nut_compute

contains

  !> Compute Deardorff SGS quantities on the device backend.
  !! @param TKE_d Device pointer to the TKE field.
  !! @param dTdx_d Device pointer to the x-gradient of temperature.
  !! @param dTdy_d Device pointer to the y-gradient of temperature.
  !! @param dTdz_d Device pointer to the z-gradient of temperature.
  !! @param a11_d Device pointer to velocity gradient component a11.
  !! @param a12_d Device pointer to velocity gradient component a12.
  !! @param a13_d Device pointer to velocity gradient component a13.
  !! @param a21_d Device pointer to velocity gradient component a21.
  !! @param a22_d Device pointer to velocity gradient component a22.
  !! @param a23_d Device pointer to velocity gradient component a23.
  !! @param a31_d Device pointer to velocity gradient component a31.
  !! @param a32_d Device pointer to velocity gradient component a32.
  !! @param a33_d Device pointer to velocity gradient component a33.
  !! @param delta_d Device pointer to LES filter width field.
  !! @param nut_d Device pointer to eddy viscosity output field.
  !! @param temperature_alphat Device pointer to temperature eddy diffusivity.
  !! @param TKE_alphat Device pointer to TKE eddy diffusivity.
  !! @param TKE_source Device pointer to TKE source output field.
  !! @param c_k Deardorff model constant.
  !! @param T0 Reference temperature.
  !! @param g Gravity vector.
  !! @param eps Lower bound used to clip TKE.
  !! @param n Number of degrees of freedom.
  subroutine device_deardorff_nut_compute(TKE_d, &
       dTdx_d, dTdy_d, dTdz_d, a11_d, a12_d, a13_d, &
       a21_d, a22_d, a23_d, &
       a31_d, a32_d, a33_d, &
       delta_d, &
       nut_d, temperature_alphat, TKE_alphat, TKE_source, &
       c_k, T0, g, eps, n)
    type(c_ptr) :: TKE_d, dTdx_d, dTdy_d, dTdz_d, &
         a11_d, a12_d, a13_d, &
         a21_d, a22_d, a23_d, &
         a31_d, a32_d, a33_d, &
         delta_d, nut_d, temperature_alphat, TKE_alphat, TKE_source
    integer :: n
    real(kind=rp) :: c_k, T0, g(3), eps
#if HAVE_HIP
    call hip_deardorff_nut_compute(TKE_d, &
         dTdx_d, dTdy_d, dTdz_d, &
         a11_d, a12_d, a13_d, &
         a21_d, a22_d, a23_d, &
         a31_d, a32_d, a33_d, &
         delta_d, nut_d, temperature_alphat, TKE_alphat, TKE_source, &
         c_k, T0, g(1), g(2), g(3), eps, n)
#elif HAVE_CUDA
    call cuda_deardorff_nut_compute(TKE_d, &
         dTdx_d, dTdy_d, dTdz_d, &
         a11_d, a12_d, a13_d, &
         a21_d, a22_d, a23_d, &
         a31_d, a32_d, a33_d, &
         delta_d, nut_d, temperature_alphat, TKE_alphat, TKE_source, &
         c_k, T0, g(1), g(2), g(3), eps, n)
#elif HAVE_OPENCL
    call neko_error('opencl backend is not supported for device_deardorff_nut')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_deardorff_nut_compute


end module device_deardorff_nut
