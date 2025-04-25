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
module device_smagorinsky_nut
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int
  use num_types, only: rp, c_rp
  use utils, only: neko_error
  use comm, only: NEKO_COMM, pe_size, MPI_REAL_PRECISION
  use mpi_f08, only: MPI_SUM, MPI_IN_PLACE, MPI_Allreduce
  
  implicit none
  private

#ifdef HAVE_HIP
  interface
     subroutine hip_smagorinsky_nut_compute(s11_d, s22_d, s33_d, &
                              s12_d, s13_d, s23_d, &
                              delta_d, nut_d, mult_d, c_s, n) &
          bind(c, name = 'hip_smagorinsky_nut_compute')
       use, intrinsic :: iso_c_binding, only: c_ptr, c_int
       import c_rp
       type(c_ptr), value :: s11_d, s22_d, s33_d, &
                             s12_d, s13_d, s23_d, &
                             delta_d, nut_d, mult_d
       integer(c_int) :: n
       real(c_rp) :: c_s
     end subroutine hip_smagorinsky_nut_compute
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_smagorinsky_nut_compute(s11_d, s22_d, s33_d, &
                              s12_d, s13_d, s23_d, &
                              delta_d, nut_d, mult_d, c_s, n) &
          bind(c, name = 'cuda_smagorinsky_nut_compute')
       use, intrinsic :: iso_c_binding, only: c_ptr, c_int
       import c_rp
       type(c_ptr), value :: s11_d, s22_d, s33_d, &
                             s12_d, s13_d, s23_d, &
                             delta_d, nut_d, mult_d
       integer(c_int) :: n
       real(c_rp) :: c_s
     end subroutine cuda_smagorinsky_nut_compute
  end interface
#elif HAVE_OPENCL
#endif

  public :: device_smagorinsky_nut_compute

contains

  !> Compute the eddy viscosity field for the Sigma model indevice
  subroutine device_smagorinsky_nut_compute(s11_d, s22_d, s33_d, &
                              s12_d, s13_d, s23_d, &
                              delta_d, nut_d, mult_d, c_s, n)
    type(c_ptr) :: s11_d, s22_d, s33_d, &
                   s12_d, s13_d, s23_d, &
                   delta_d, nut_d, mult_d
    integer :: n
    real(kind=rp) :: c_s
#if HAVE_HIP
    call hip_smagorinsky_nut_compute(s11_d, s22_d, s33_d, &
                              s12_d, s13_d, s23_d, &
                              delta_d, nut_d, mult_d, c_s, n)
#elif HAVE_CUDA
    call cuda_smagorinsky_nut_compute(s11_d, s22_d, s33_d, &
                              s12_d, s13_d, s23_d, &
                              delta_d, nut_d, mult_d, c_s, n)
#elif HAVE_OPENCL
    call neko_error('opencl backend is not supported for device_smagorinsky_nut')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_smagorinsky_nut_compute

  
end module device_smagorinsky_nut