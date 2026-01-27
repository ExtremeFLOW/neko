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
!> Implements device kernels for use with TreeAMG smoothers
module device_tree_amg_smoother
  use num_types, only : rp, c_rp
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_bool
  use device, only : glb_cmd_queue
  use utils, only : neko_error
  implicit none
  private

#ifdef HAVE_HIP
  interface
     subroutine hip_amg_cheby_solve_part1(r_d, f_d, w_d, x_d, d_d, &
          inv_thet, n, zero_initial, strm) &
          bind(c, name='hip_amg_cheby_solve_part1')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_bool
       import c_rp
       type(c_ptr), value :: r_d, f_d, w_d, x_d, d_d, strm
       real(c_rp) :: inv_thet
       logical(c_bool) :: zero_initial
     end subroutine hip_amg_cheby_solve_part1
  end interface

  interface
     subroutine hip_amg_cheby_solve_part2(r_d, w_d, d_d, x_d, &
          tmp1, tmp2, n, strm) bind(c, name='hip_amg_cheby_solve_part2')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_rp
       type(c_ptr), value :: r_d, w_d, d_d, x_d, strm
       real(c_rp) :: tmp1, tmp2
       integer(c_int) :: n
     end subroutine hip_amg_cheby_solve_part2
  end interface
#elif HAVE_CUDA
  interface
     subroutine cuda_amg_cheby_solve_part1(r_d, f_d, w_d, x_d, d_d, &
          inv_thet, n, zero_initial, strm) &
          bind(c, name='cuda_amg_cheby_solve_part1')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int, c_bool
       import c_rp
       type(c_ptr), value :: r_d, f_d, w_d, x_d, d_d, strm
       real(c_rp) :: inv_thet
       logical(c_bool) :: zero_initial
     end subroutine cuda_amg_cheby_solve_part1
  end interface

  interface
     subroutine cuda_amg_cheby_solve_part2(r_d, w_d, d_d, x_d, &
          tmp1, tmp2, n, strm) bind(c, name='cuda_amg_cheby_solve_part2')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_rp
       type(c_ptr), value :: r_d, w_d, d_d, x_d, strm
       real(c_rp) :: tmp1, tmp2
       integer(c_int) :: n
     end subroutine cuda_amg_cheby_solve_part2
  end interface
#elif HAVE_OPENCL
#endif

  public :: amg_device_cheby_solve_part1, amg_device_cheby_solve_part2

contains

  subroutine amg_device_cheby_solve_part1(r_d, f_d, w_d, x_d, d_d, &
       inv_thet, n, zero_initial)
    type(c_ptr), intent(inout) :: r_d, f_d, w_d, d_d, x_d
    real(kind=rp), intent(in) :: inv_thet
    integer, intent(in) :: n
    logical, intent(in) :: zero_initial
    logical(c_bool) :: zinit

    zinit = zero_initial

#ifdef HAVE_HIP
    call hip_amg_cheby_solve_part1(r_d, f_d, w_d, x_d, d_d, &
         inv_thet, n, zinit, glb_cmd_queue)
#elif HAVE_CUDA
    call cuda_amg_cheby_solve_part1(r_d, f_d, w_d, x_d, d_d, &
         inv_thet, n, zinit, glb_cmd_queue)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine amg_device_cheby_solve_part1

  subroutine amg_device_cheby_solve_part2(r_d, w_d, d_d, x_d, tmp1, tmp2, n)
    type(c_ptr), intent(inout) :: r_d, w_d, d_d, x_d
    real(kind=rp), intent(in) :: tmp1, tmp2
    integer, intent(in) :: n
#ifdef HAVE_HIP
    call hip_amg_cheby_solve_part2(r_d, w_d, d_d, x_d, &
         tmp1, tmp2, n, glb_cmd_queue)
#elif HAVE_CUDA
    call cuda_amg_cheby_solve_part2(r_d, w_d, d_d, x_d, &
         tmp1, tmp2, n, glb_cmd_queue)
#else
    call neko_error('No device backend configured')
#endif
  end subroutine amg_device_cheby_solve_part2

end module device_tree_amg_smoother
