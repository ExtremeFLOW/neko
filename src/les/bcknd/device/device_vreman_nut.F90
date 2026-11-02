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
module device_vreman_nut
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int
  use num_types, only : rp, c_rp
  use utils, only : neko_error

  implicit none
  private

#ifdef HAVE_HIP
  interface
     subroutine hip_vreman_nut_compute(a11_d, a12_d, a13_d, &
          a21_d, a22_d, a23_d, &
          a31_d, a32_d, a33_d, &
          delta_d, nut_d, mult_d, c, eps, n) &
          bind(c, name = 'hip_vreman_nut_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_rp
       type(c_ptr), value :: a11_d, a12_d, a13_d, &
            a21_d, a22_d, a23_d, &
            a31_d, a32_d, a33_d, &
            delta_d, nut_d, mult_d
       integer(c_int) :: n
       real(c_rp) :: c, eps
     end subroutine hip_vreman_nut_compute

     subroutine hip_vreman_nut_compute_aniso(a11_d, a12_d, a13_d, &
          a21_d, a22_d, a23_d, &
          a31_d, a32_d, a33_d, &
          g11_d, g12_d, g13_d, g22_d, g23_d, g33_d, &
          nut_d, mult_d, c, eps, n) &
          bind(c, name = 'hip_vreman_nut_compute_aniso')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_rp
       type(c_ptr), value :: a11_d, a12_d, a13_d, &
            a21_d, a22_d, a23_d, &
            a31_d, a32_d, a33_d, &
            g11_d, g12_d, g13_d, g22_d, g23_d, g33_d, &
            nut_d, mult_d
       integer(c_int) :: n
       real(c_rp) :: c, eps
     end subroutine hip_vreman_nut_compute_aniso

     subroutine hip_vreman_stability(a11_d, a12_d, a13_d, &
          a21_d, a22_d, a23_d, &
          a31_d, a32_d, a33_d, &
          gdT_d, nut_d, mult_d, &
          g, ri_c, ref_temp, eps, n) &
          bind(c, name = 'hip_vreman_stability')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_rp
       type(c_ptr), value :: a11_d, a12_d, a13_d, &
            a21_d, a22_d, a23_d, &
            a31_d, a32_d, a33_d, &
            gdT_d, nut_d, mult_d
       real(c_rp) :: g(3)
       real(c_rp) :: ri_c, ref_temp, eps
       integer(c_int) :: n
     end subroutine hip_vreman_stability
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_vreman_nut_compute(a11_d, a12_d, a13_d, &
          a21_d, a22_d, a23_d, &
          a31_d, a32_d, a33_d, &
          delta_d, nut_d, mult_d, c, eps, n) &
          bind(c, name = 'cuda_vreman_nut_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_rp
       type(c_ptr), value :: a11_d, a12_d, a13_d, &
            a21_d, a22_d, a23_d, &
            a31_d, a32_d, a33_d, &
            delta_d, nut_d, mult_d
       integer(c_int) :: n
       real(c_rp) :: c, eps
     end subroutine cuda_vreman_nut_compute

     subroutine cuda_vreman_nut_compute_aniso(a11_d, a12_d, a13_d, &
          a21_d, a22_d, a23_d, &
          a31_d, a32_d, a33_d, &
          g11_d, g12_d, g13_d, g22_d, g23_d, g33_d, &
          nut_d, mult_d, c, eps, n) &
          bind(c, name = 'cuda_vreman_nut_compute_aniso')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_rp
       type(c_ptr), value :: a11_d, a12_d, a13_d, &
            a21_d, a22_d, a23_d, &
            a31_d, a32_d, a33_d, &
            g11_d, g12_d, g13_d, g22_d, g23_d, g33_d, &
            nut_d, mult_d
       integer(c_int) :: n
       real(c_rp) :: c, eps
     end subroutine cuda_vreman_nut_compute_aniso

     subroutine cuda_vreman_stability(a11_d, a12_d, a13_d, &
          a21_d, a22_d, a23_d, &
          a31_d, a32_d, a33_d, &
          gdT_d, nut_d, mult_d, &
          g, ri_c, ref_temp, eps, n) &
          bind(c, name = 'cuda_vreman_stability')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_rp
       type(c_ptr), value :: a11_d, a12_d, a13_d, &
            a21_d, a22_d, a23_d, &
            a31_d, a32_d, a33_d, &
            gdT_d, nut_d, mult_d
       real(c_rp) :: g(3)
       real(c_rp) :: ri_c, ref_temp, eps
       integer(c_int) :: n
     end subroutine cuda_vreman_stability
  end interface
#elif HAVE_OPENCL
  interface
     subroutine opencl_vreman_nut_compute(a11_d, a12_d, a13_d, &
          a21_d, a22_d, a23_d, &
          a31_d, a32_d, a33_d, &
          delta_d, nut_d, mult_d, c, eps, n) &
          bind(c, name = 'opencl_vreman_nut_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_rp
       type(c_ptr), value :: a11_d, a12_d, a13_d, &
            a21_d, a22_d, a23_d, &
            a31_d, a32_d, a33_d, &
            delta_d, nut_d, mult_d
       integer(c_int) :: n
       real(c_rp) :: c, eps
     end subroutine opencl_vreman_nut_compute

     subroutine opencl_vreman_nut_compute_aniso(a11_d, a12_d, a13_d, &
          a21_d, a22_d, a23_d, &
          a31_d, a32_d, a33_d, &
          g11_d, g12_d, g13_d, g22_d, g23_d, g33_d, &
          nut_d, mult_d, c, eps, n) &
          bind(c, name = 'opencl_vreman_nut_compute_aniso')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_rp
       type(c_ptr), value :: a11_d, a12_d, a13_d, &
            a21_d, a22_d, a23_d, &
            a31_d, a32_d, a33_d, &
            g11_d, g12_d, g13_d, g22_d, g23_d, g33_d, &
            nut_d, mult_d
       integer(c_int) :: n
       real(c_rp) :: c, eps
     end subroutine opencl_vreman_nut_compute_aniso

     subroutine opencl_vreman_stability(a11_d, a12_d, a13_d, &
          a21_d, a22_d, a23_d, &
          a31_d, a32_d, a33_d, &
          gdT_d, nut_d, mult_d, &
          g, ri_c, ref_temp, eps, n) &
          bind(c, name = 'opencl_vreman_stability')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_rp
       type(c_ptr), value :: a11_d, a12_d, a13_d, &
            a21_d, a22_d, a23_d, &
            a31_d, a32_d, a33_d, &
            gdT_d, nut_d, mult_d
       real(c_rp) :: g(3)
       real(c_rp) :: ri_c, ref_temp, eps
       integer(c_int) :: n
     end subroutine opencl_vreman_stability
  end interface
#elif HAVE_METAL
  interface
     subroutine metal_vreman_nut_compute(a11_d, a12_d, a13_d, &
          a21_d, a22_d, a23_d, &
          a31_d, a32_d, a33_d, &
          delta_d, nut_d, mult_d, c, eps, n) &
          bind(c, name = 'metal_vreman_nut_compute')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_rp
       type(c_ptr), value :: a11_d, a12_d, a13_d, &
            a21_d, a22_d, a23_d, &
            a31_d, a32_d, a33_d, &
            delta_d, nut_d, mult_d
       integer(c_int) :: n
       real(c_rp) :: c, eps
     end subroutine metal_vreman_nut_compute

     subroutine metal_vreman_nut_compute_aniso(a11_d, a12_d, a13_d, &
          a21_d, a22_d, a23_d, &
          a31_d, a32_d, a33_d, &
          g11_d, g12_d, g13_d, g22_d, g23_d, g33_d, &
          nut_d, mult_d, c, eps, n) &
          bind(c, name = 'metal_vreman_nut_compute_aniso')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_rp
       type(c_ptr), value :: a11_d, a12_d, a13_d, &
            a21_d, a22_d, a23_d, &
            a31_d, a32_d, a33_d, &
            g11_d, g12_d, g13_d, g22_d, g23_d, g33_d, &
            nut_d, mult_d
       integer(c_int) :: n
       real(c_rp) :: c, eps
     end subroutine metal_vreman_nut_compute_aniso

     subroutine metal_vreman_stability(a11_d, a12_d, a13_d, &
          a21_d, a22_d, a23_d, &
          a31_d, a32_d, a33_d, &
          gdT_d, nut_d, mult_d, &
          g, ri_c, ref_temp, eps, n) &
          bind(c, name = 'metal_vreman_stability')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_rp
       type(c_ptr), value :: a11_d, a12_d, a13_d, &
            a21_d, a22_d, a23_d, &
            a31_d, a32_d, a33_d, &
            gdT_d, nut_d, mult_d
       real(c_rp) :: g(3)
       real(c_rp) :: ri_c, ref_temp, eps
       integer(c_int) :: n
     end subroutine metal_vreman_stability
  end interface
#endif

  public :: device_vreman_nut_compute
  public :: device_vreman_nut_compute_aniso
  public :: device_vreman_stability

contains

  !> Compute the isotropic Vreman eddy viscosity on the device.
  subroutine device_vreman_nut_compute(a11_d, a12_d, a13_d, &
       a21_d, a22_d, a23_d, &
       a31_d, a32_d, a33_d, &
       delta_d, nut_d, mult_d, c, eps, n)
    type(c_ptr) :: a11_d, a12_d, a13_d, &
         a21_d, a22_d, a23_d, &
         a31_d, a32_d, a33_d, &
         delta_d, nut_d, mult_d
    integer :: n
    real(kind=rp) :: c, eps
#if HAVE_HIP
    call hip_vreman_nut_compute(a11_d, a12_d, a13_d, &
         a21_d, a22_d, a23_d, &
         a31_d, a32_d, a33_d, &
         delta_d, nut_d, mult_d, c, eps, n)
#elif HAVE_CUDA
    call cuda_vreman_nut_compute(a11_d, a12_d, a13_d, &
         a21_d, a22_d, a23_d, &
         a31_d, a32_d, a33_d, &
         delta_d, nut_d, mult_d, c, eps, n)
#elif HAVE_OPENCL
    call opencl_vreman_nut_compute(a11_d, a12_d, a13_d, &
         a21_d, a22_d, a23_d, &
         a31_d, a32_d, a33_d, &
         delta_d, nut_d, mult_d, c, eps, n)
#elif HAVE_METAL
    call metal_vreman_nut_compute(a11_d, a12_d, a13_d, &
         a21_d, a22_d, a23_d, &
         a31_d, a32_d, a33_d, &
         delta_d, nut_d, mult_d, c, eps, n)
#else
    call neko_error('no device backend configured (vreman_nut_compute)')
#endif
  end subroutine device_vreman_nut_compute

  !> Compute the anisotropic Vreman eddy viscosity on the device.
  !! @param g11_d ... g33_d The filter-width tensor G (symmetric).
  subroutine device_vreman_nut_compute_aniso(a11_d, a12_d, a13_d, &
       a21_d, a22_d, a23_d, &
       a31_d, a32_d, a33_d, &
       g11_d, g12_d, g13_d, g22_d, g23_d, g33_d, &
       nut_d, mult_d, c, eps, n)
    type(c_ptr) :: a11_d, a12_d, a13_d, &
         a21_d, a22_d, a23_d, &
         a31_d, a32_d, a33_d, &
         g11_d, g12_d, g13_d, g22_d, g23_d, g33_d, &
         nut_d, mult_d
    integer :: n
    real(kind=rp) :: c, eps
#if HAVE_HIP
    call hip_vreman_nut_compute_aniso(a11_d, a12_d, a13_d, &
         a21_d, a22_d, a23_d, &
         a31_d, a32_d, a33_d, &
         g11_d, g12_d, g13_d, g22_d, g23_d, g33_d, &
         nut_d, mult_d, c, eps, n)
#elif HAVE_CUDA
    call cuda_vreman_nut_compute_aniso(a11_d, a12_d, a13_d, &
         a21_d, a22_d, a23_d, &
         a31_d, a32_d, a33_d, &
         g11_d, g12_d, g13_d, g22_d, g23_d, g33_d, &
         nut_d, mult_d, c, eps, n)
#elif HAVE_OPENCL
    call opencl_vreman_nut_compute_aniso(a11_d, a12_d, a13_d, &
         a21_d, a22_d, a23_d, &
         a31_d, a32_d, a33_d, &
         g11_d, g12_d, g13_d, g22_d, g23_d, g33_d, &
         nut_d, mult_d, c, eps, n)
#elif HAVE_METAL
    call metal_vreman_nut_compute_aniso(a11_d, a12_d, a13_d, &
         a21_d, a22_d, a23_d, &
         a31_d, a32_d, a33_d, &
         g11_d, g12_d, g13_d, g22_d, g23_d, g33_d, &
         nut_d, mult_d, c, eps, n)
#else
    call neko_error('no device backend configured (vreman_nut_compute_aniso)')
#endif
  end subroutine device_vreman_nut_compute_aniso

  !> Apply the Richardson-number correction of the Vreman model on the
  !! device (Lilly 1962; Moeng & Sullivan 2015).
  !! @param gdT_d The interface-averaged g . grad(T).
  !! @note g points upwards (reversed from the case-file input at init).
  subroutine device_vreman_stability(a11_d, a12_d, a13_d, &
       a21_d, a22_d, a23_d, &
       a31_d, a32_d, a33_d, &
       gdT_d, nut_d, mult_d, &
       g, ri_c, ref_temp, eps, n)
    type(c_ptr) :: a11_d, a12_d, a13_d, &
         a21_d, a22_d, a23_d, &
         a31_d, a32_d, a33_d, &
         gdT_d, nut_d, mult_d
    real(kind=rp) :: g(3), ri_c, ref_temp, eps
    integer :: n
#if HAVE_HIP
    call hip_vreman_stability(a11_d, a12_d, a13_d, &
         a21_d, a22_d, a23_d, &
         a31_d, a32_d, a33_d, &
         gdT_d, nut_d, mult_d, &
         g, ri_c, ref_temp, eps, n)
#elif HAVE_CUDA
    call cuda_vreman_stability(a11_d, a12_d, a13_d, &
         a21_d, a22_d, a23_d, &
         a31_d, a32_d, a33_d, &
         gdT_d, nut_d, mult_d, &
         g, ri_c, ref_temp, eps, n)
#elif HAVE_OPENCL
    call opencl_vreman_stability(a11_d, a12_d, a13_d, &
         a21_d, a22_d, a23_d, &
         a31_d, a32_d, a33_d, &
         gdT_d, nut_d, mult_d, &
         g, ri_c, ref_temp, eps, n)
#elif HAVE_METAL
    call metal_vreman_stability(a11_d, a12_d, a13_d, &
         a21_d, a22_d, a23_d, &
         a31_d, a32_d, a33_d, &
         gdT_d, nut_d, mult_d, &
         g, ri_c, ref_temp, eps, n)
#else
    call neko_error('no device backend configured (vreman_stability)')
#endif
  end subroutine device_vreman_stability

end module device_vreman_nut
