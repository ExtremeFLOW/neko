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
module device_dynamic_smagorinsky_nut
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int
  use num_types, only: rp, c_rp
  use utils, only: neko_error
  use comm, only: NEKO_COMM, pe_size, MPI_REAL_PRECISION
  use mpi_f08, only: MPI_SUM, MPI_IN_PLACE, MPI_Allreduce

  implicit none
  private

#ifdef HAVE_HIP
  interface
     subroutine hip_s_abs_compute(s_abs_d, s11_d, s22_d, s33_d, &
          s12_d, s13_d, s23_d, &
          n) &
          bind(c, name = 'hip_s_abs_compute')
       use, intrinsic :: iso_c_binding, only: c_ptr, c_int
       import c_rp
       type(c_ptr), value :: s_abs_d, s11_d, s22_d, s33_d, &
            s12_d, s13_d, s23_d
       integer(c_int) :: n
     end subroutine hip_s_abs_compute
  end interface
  interface
     subroutine hip_lij_compute_part1(l11_d, l22_d, l33_d, &
          l12_d, l13_d, l23_d, &
          u_d, v_d, w_d, &
          fu_d, fv_d, fw_d, &
          fuu_d, fvv_d, fww_d, &
          fuv_d, fuw_d, fvw_d, n) &
          bind(c, name = 'hip_lij_compute_part1')
       use, intrinsic :: iso_c_binding, only: c_ptr, c_int
       import c_rp
       type(c_ptr), value :: l11_d, l22_d, l33_d, l12_d, l13_d, l23_d, &
            u_d, v_d, w_d, fu_d, fv_d, fw_d, &
            fuu_d, fvv_d, fww_d, fuv_d, fuw_d, fvw_d
       integer(c_int) :: n
     end subroutine hip_lij_compute_part1
  end interface
  interface
     subroutine hip_lij_compute_part2(l11_d, l22_d, l33_d, &
          l12_d, l13_d, l23_d, &
          fuu_d, fvv_d, fww_d, &
          fuv_d, fuw_d, fvw_d, n) &
          bind(c, name = 'hip_lij_compute_part2')
       use, intrinsic :: iso_c_binding, only: c_ptr, c_int
       import c_rp
       type(c_ptr), value :: l11_d, l22_d, l33_d, l12_d, l13_d, l23_d, &
            fuu_d, fvv_d, fww_d, fuv_d, fuw_d, fvw_d
       integer(c_int) :: n
     end subroutine hip_lij_compute_part2
  end interface
  interface
     subroutine hip_mij_compute_part1(m11_d, m22_d, m33_d, &
          m12_d, m13_d, m23_d, &
          s_abs_d, s11_d, s22_d, s33_d, &
          s12_d, s13_d, s23_d, &
          fs_abs_d, fs11_d, fs22_d, fs33_d, &
          fs12_d, fs13_d, fs23_d, &
          fsabss11_d, fsabss22_d, fsabss33_d, &
          fsabss12_d, fsabss13_d, fsabss23_d, &
          delta_ratio2, n) &
          bind(c, name = 'hip_mij_compute_part1')
       use, intrinsic :: iso_c_binding, only: c_ptr, c_int
       import c_rp
       type(c_ptr), value :: m11_d, m22_d, m33_d, &
            m12_d, m13_d, m23_d, &
            s_abs_d, s11_d, s22_d, s33_d, &
            s12_d, s13_d, s23_d, &
            fs_abs_d, fs11_d, fs22_d, fs33_d, &
            fs12_d, fs13_d, fs23_d, &
            fsabss11_d, fsabss22_d, fsabss33_d, &
            fsabss12_d, fsabss13_d, fsabss23_d
       real(c_rp) :: delta_ratio2
       integer(c_int) :: n
     end subroutine hip_mij_compute_part1
  end interface
  interface
     subroutine hip_mij_nut_compute_part2(m11_d, m22_d, m33_d, &
          m12_d, m13_d, m23_d, &
          l11_d, l22_d, l33_d, &
          l12_d, l13_d, l23_d, &
          fsabss11_d, fsabss22_d, fsabss33_d, &
          fsabss12_d, fsabss13_d, fsabss23_d, &
          num_d, den_d, c_dyn_d, delta_d, &
          s_abs_d, nut_d, alpha, n) &
          bind(c, name = 'hip_mij_nut_compute_part2')
       use, intrinsic :: iso_c_binding, only: c_ptr, c_int
       import c_rp
       type(c_ptr), value :: m11_d, m22_d, m33_d, &
            m12_d, m13_d, m23_d, &
            l11_d, l22_d, l33_d, &
            l12_d, l13_d, l23_d, &
            fsabss11_d, fsabss22_d, fsabss33_d, &
            fsabss12_d, fsabss13_d, fsabss23_d, &
            num_d, den_d, c_dyn_d, delta_d, s_abs_d, nut_d
       real(c_rp) :: alpha
       integer(c_int) :: n
     end subroutine hip_mij_nut_compute_part2
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_s_abs_compute(s_abs_d, s11_d, s22_d, s33_d, &
          s12_d, s13_d, s23_d, &
          n) &
          bind(c, name = 'cuda_s_abs_compute')
       use, intrinsic :: iso_c_binding, only: c_ptr, c_int
       import c_rp
       type(c_ptr), value :: s_abs_d, s11_d, s22_d, s33_d, &
            s12_d, s13_d, s23_d
       integer(c_int) :: n
     end subroutine cuda_s_abs_compute
  end interface
  interface
     subroutine cuda_lij_compute_part1(l11_d, l22_d, l33_d, &
          l12_d, l13_d, l23_d, &
          u_d, v_d, w_d, &
          fu_d, fv_d, fw_d, &
          fuu_d, fvv_d, fww_d, &
          fuv_d, fuw_d, fvw_d, n) &
          bind(c, name = 'cuda_lij_compute_part1')
       use, intrinsic :: iso_c_binding, only: c_ptr, c_int
       import c_rp
       type(c_ptr), value :: l11_d, l22_d, l33_d, l12_d, l13_d, l23_d, &
            u_d, v_d, w_d, fu_d, fv_d, fw_d, &
            fuu_d, fvv_d, fww_d, fuv_d, fuw_d, fvw_d
       integer(c_int) :: n
     end subroutine cuda_lij_compute_part1
  end interface
  interface
     subroutine cuda_lij_compute_part2(l11_d, l22_d, l33_d, &
          l12_d, l13_d, l23_d, &
          fuu_d, fvv_d, fww_d, &
          fuv_d, fuw_d, fvw_d, n) &
          bind(c, name = 'cuda_lij_compute_part2')
       use, intrinsic :: iso_c_binding, only: c_ptr, c_int
       import c_rp
       type(c_ptr), value :: l11_d, l22_d, l33_d, l12_d, l13_d, l23_d, &
            fuu_d, fvv_d, fww_d, fuv_d, fuw_d, fvw_d
       integer(c_int) :: n
     end subroutine cuda_lij_compute_part2
  end interface
  interface
     subroutine cuda_mij_compute_part1(m11_d, m22_d, m33_d, &
          m12_d, m13_d, m23_d, &
          s_abs_d, s11_d, s22_d, s33_d, &
          s12_d, s13_d, s23_d, &
          fs_abs_d, fs11_d, fs22_d, fs33_d, &
          fs12_d, fs13_d, fs23_d, &
          fsabss11_d, fsabss22_d, fsabss33_d, &
          fsabss12_d, fsabss13_d, fsabss23_d, &
          delta_ratio2, n) &
          bind(c, name = 'cuda_mij_compute_part1')
       use, intrinsic :: iso_c_binding, only: c_ptr, c_int
       import c_rp
       type(c_ptr), value :: m11_d, m22_d, m33_d, &
            m12_d, m13_d, m23_d, &
            s_abs_d, s11_d, s22_d, s33_d, &
            s12_d, s13_d, s23_d, &
            fs_abs_d, fs11_d, fs22_d, fs33_d, &
            fs12_d, fs13_d, fs23_d, &
            fsabss11_d, fsabss22_d, fsabss33_d, &
            fsabss12_d, fsabss13_d, fsabss23_d
       real(c_rp) :: delta_ratio2
       integer(c_int) :: n
     end subroutine cuda_mij_compute_part1
  end interface
  interface
     subroutine cuda_mij_nut_compute_part2(m11_d, m22_d, m33_d, &
          m12_d, m13_d, m23_d, &
          l11_d, l22_d, l33_d, &
          l12_d, l13_d, l23_d, &
          fsabss11_d, fsabss22_d, fsabss33_d, &
          fsabss12_d, fsabss13_d, fsabss23_d, &
          num_d, den_d, c_dyn_d, delta_d, &
          s_abs_d, nut_d, alpha, n) &
          bind(c, name = 'cuda_mij_nut_compute_part2')
       use, intrinsic :: iso_c_binding, only: c_ptr, c_int
       import c_rp
       type(c_ptr), value :: m11_d, m22_d, m33_d, &
            m12_d, m13_d, m23_d, &
            l11_d, l22_d, l33_d, &
            l12_d, l13_d, l23_d, &
            fsabss11_d, fsabss22_d, fsabss33_d, &
            fsabss12_d, fsabss13_d, fsabss23_d, &
            num_d, den_d, c_dyn_d, delta_d, s_abs_d, nut_d
       real(c_rp) :: alpha
       integer(c_int) :: n
     end subroutine cuda_mij_nut_compute_part2
  end interface
#elif HAVE_OPENCL
#endif

  public :: device_s_abs_compute, device_lij_compute_part1, &
       device_lij_compute_part2, device_mij_compute_part1, &
       device_mij_nut_compute_part2
contains

  !> Compute the s_abs field for the Sigma model indevice
  subroutine device_s_abs_compute(s_abs_d, s11_d, s22_d, s33_d, &
       s12_d, s13_d, s23_d, &
       n)
    type(c_ptr) :: s_abs_d, s11_d, s22_d, s33_d, &
         s12_d, s13_d, s23_d
    integer :: n
#if HAVE_HIP
    call hip_s_abs_compute(s_abs_d, s11_d, s22_d, s33_d, &
         s12_d, s13_d, s23_d, &
         n)
#elif HAVE_CUDA
    call cuda_s_abs_compute(s_abs_d, s11_d, s22_d, s33_d, &
         s12_d, s13_d, s23_d, &
         n)
#elif HAVE_OPENCL
    call neko_error('opencl backend is not supported for device_s_abs_compute')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_s_abs_compute

  !> part 1 of the computing of the lij field
  subroutine device_lij_compute_part1(l11_d, l22_d, l33_d, &
       l12_d, l13_d, l23_d, &
       u_d, v_d, w_d, &
       fu_d, fv_d, fw_d, &
       fuu_d, fvv_d, fww_d, &
       fuv_d, fuw_d, fvw_d, n)
    type(c_ptr) :: l11_d, l22_d, l33_d, l12_d, l13_d, l23_d, &
         u_d, v_d, w_d, fu_d, fv_d, fw_d, &
         fuu_d, fvv_d, fww_d, fuv_d, fuw_d, fvw_d
    integer :: n
#if HAVE_HIP
    call hip_lij_compute_part1(l11_d, l22_d, l33_d, &
         l12_d, l13_d, l23_d, &
         u_d, v_d, w_d, &
         fu_d, fv_d, fw_d, &
         fuu_d, fvv_d, fww_d, &
         fuv_d, fuw_d, fvw_d, n)
#elif HAVE_CUDA
    call cuda_lij_compute_part1(l11_d, l22_d, l33_d, &
         l12_d, l13_d, l23_d, &
         u_d, v_d, w_d, &
         fu_d, fv_d, fw_d, &
         fuu_d, fvv_d, fww_d, &
         fuv_d, fuw_d, fvw_d, n)
#elif HAVE_OPENCL
    call neko_error('opencl backend is not supported &
    &for device_lij_compute_part1')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_lij_compute_part1

!> part 2 of the computing of the lij field
  subroutine device_lij_compute_part2(l11_d, l22_d, l33_d, &
       l12_d, l13_d, l23_d, &
       fuu_d, fvv_d, fww_d, &
       fuv_d, fuw_d, fvw_d, n)
    type(c_ptr) :: l11_d, l22_d, l33_d, l12_d, l13_d, l23_d, &
         fuu_d, fvv_d, fww_d, fuv_d, fuw_d, fvw_d
    integer :: n
#if HAVE_HIP
    call hip_lij_compute_part2(l11_d, l22_d, l33_d, &
         l12_d, l13_d, l23_d, &
         fuu_d, fvv_d, fww_d, &
         fuv_d, fuw_d, fvw_d, n)
#elif HAVE_CUDA
    call cuda_lij_compute_part2(l11_d, l22_d, l33_d, &
         l12_d, l13_d, l23_d, &
         fuu_d, fvv_d, fww_d, &
         fuv_d, fuw_d, fvw_d, n)
#elif HAVE_OPENCL
    call neko_error('opencl backend is not supported &
    &for device_lij_compute_part2')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_lij_compute_part2

  !> part 1 of the computing of the mij field
  subroutine device_mij_compute_part1(m11_d, m22_d, m33_d, &
       m12_d, m13_d, m23_d, &
       s_abs_d, s11_d, s22_d, s33_d, &
       s12_d, s13_d, s23_d, &
       fs_abs_d, fs11_d, fs22_d, fs33_d, &
       fs12_d, fs13_d, fs23_d, &
       fsabss11_d, fsabss22_d, fsabss33_d, &
       fsabss12_d, fsabss13_d, fsabss23_d, &
       delta_ratio2, n)
    type(c_ptr) :: m11_d, m22_d, m33_d, &
         m12_d, m13_d, m23_d, &
         s_abs_d, s11_d, s22_d, s33_d, &
         s12_d, s13_d, s23_d, &
         fs_abs_d, fs11_d, fs22_d, fs33_d, &
         fs12_d, fs13_d, fs23_d, &
         fsabss11_d, fsabss22_d, fsabss33_d, &
         fsabss12_d, fsabss13_d, fsabss23_d
    real(kind=rp) :: delta_ratio2
    integer :: n

#if HAVE_HIP
    call hip_mij_compute_part1(m11_d, m22_d, m33_d, &
         m12_d, m13_d, m23_d, &
         s_abs_d, s11_d, s22_d, s33_d, &
         s12_d, s13_d, s23_d, &
         fs_abs_d, fs11_d, fs22_d, fs33_d, &
         fs12_d, fs13_d, fs23_d, &
         fsabss11_d, fsabss22_d, fsabss33_d, &
         fsabss12_d, fsabss13_d, fsabss23_d, &
         delta_ratio2, n)
#elif HAVE_CUDA
    call cuda_mij_compute_part1(m11_d, m22_d, m33_d, &
         m12_d, m13_d, m23_d, &
         s_abs_d, s11_d, s22_d, s33_d, &
         s12_d, s13_d, s23_d, &
         fs_abs_d, fs11_d, fs22_d, fs33_d, &
         fs12_d, fs13_d, fs23_d, &
         fsabss11_d, fsabss22_d, fsabss33_d, &
         fsabss12_d, fsabss13_d, fsabss23_d, &
         delta_ratio2, n)
#elif HAVE_OPENCL
    call neko_error('opencl backend is not supported for &
    &device_mij_compute_part1')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_mij_compute_part1

  !> part 1 of the computing of the mij field
  subroutine device_mij_nut_compute_part2(m11_d, m22_d, m33_d, &
       m12_d, m13_d, m23_d, &
       l11_d, l22_d, l33_d, &
       l12_d, l13_d, l23_d, &
       fsabss11_d, fsabss22_d, fsabss33_d, &
       fsabss12_d, fsabss13_d, fsabss23_d, &
       num_d, den_d, c_dyn_d, delta_d, &
       s_abs_d, nut_d, alpha, n)
    type(c_ptr) :: m11_d, m22_d, m33_d, &
         m12_d, m13_d, m23_d, &
         l11_d, l22_d, l33_d, &
         l12_d, l13_d, l23_d, &
         fsabss11_d, fsabss22_d, fsabss33_d, &
         fsabss12_d, fsabss13_d, fsabss23_d, &
         num_d, den_d, c_dyn_d, delta_d, s_abs_d, nut_d
    real(kind=rp) :: alpha
    integer :: n

#if HAVE_HIP
    call hip_mij_nut_compute_part2(m11_d, m22_d, m33_d, &
         m12_d, m13_d, m23_d, &
         l11_d, l22_d, l33_d, &
         l12_d, l13_d, l23_d, &
         fsabss11_d, fsabss22_d, fsabss33_d, &
         fsabss12_d, fsabss13_d, fsabss23_d, &
         num_d, den_d, c_dyn_d, delta_d, &
         s_abs_d, nut_d, alpha, n)
#elif HAVE_CUDA
    call cuda_mij_nut_compute_part2(m11_d, m22_d, m33_d, &
         m12_d, m13_d, m23_d, &
         l11_d, l22_d, l33_d, &
         l12_d, l13_d, l23_d, &
         fsabss11_d, fsabss22_d, fsabss33_d, &
         fsabss12_d, fsabss13_d, fsabss23_d, &
         num_d, den_d, c_dyn_d, delta_d, &
         s_abs_d, nut_d, alpha, n)
#elif HAVE_OPENCL
    call neko_error('opencl backend is not supported for &
    &device_mij_nut_compute_part2')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_mij_nut_compute_part2


end module device_dynamic_smagorinsky_nut
