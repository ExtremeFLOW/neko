! Copyright (c) 2026, The Neko Authors
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
!> Device dispatch for LPT periodic boundary-condition wrapping.
module lpt_periodic_bc_device
  use num_types, only : rp, c_rp
  use vector, only : vector_t
  use utils, only : neko_error
  use device, only : glb_cmd_queue
  use, intrinsic :: iso_c_binding, only : c_int, c_ptr, C_NULL_PTR
  implicit none
  private

#ifdef HAVE_HIP
  interface
     !> HIP kernel entry point for translational periodic wrapping.
     subroutine hip_lpt_periodic_bc_wrap_translational(x, y, z, n, &
          n_periodic_dirs, periodic_dir_x1, periodic_dir_y1, periodic_dir_z1, &
          periodic_dir_x2, periodic_dir_y2, periodic_dir_z2, periodic_dir_x3, &
          periodic_dir_y3, periodic_dir_z3, periodic_min1, periodic_min2, &
          periodic_min3, periodic_max1, periodic_max2, periodic_max3, &
          periodic_shift_x1, periodic_shift_y1, periodic_shift_z1, &
          periodic_shift_x2, periodic_shift_y2, periodic_shift_z2, &
          periodic_shift_x3, periodic_shift_y3, periodic_shift_z3, &
          periodic_len1, periodic_len2, periodic_len3, strm) &
          bind(c, name = 'hip_lpt_periodic_bc_wrap_translational')
       use, intrinsic :: iso_c_binding
       use num_types, only : c_rp
       implicit none
       integer(c_int) :: n, n_periodic_dirs
       real(c_rp) :: periodic_dir_x1, periodic_dir_y1, periodic_dir_z1
       real(c_rp) :: periodic_dir_x2, periodic_dir_y2, periodic_dir_z2
       real(c_rp) :: periodic_dir_x3, periodic_dir_y3, periodic_dir_z3
       real(c_rp) :: periodic_min1, periodic_min2, periodic_min3
       real(c_rp) :: periodic_max1, periodic_max2, periodic_max3
       real(c_rp) :: periodic_shift_x1, periodic_shift_y1, periodic_shift_z1
       real(c_rp) :: periodic_shift_x2, periodic_shift_y2, periodic_shift_z2
       real(c_rp) :: periodic_shift_x3, periodic_shift_y3, periodic_shift_z3
       real(c_rp) :: periodic_len1, periodic_len2, periodic_len3
       type(c_ptr), value :: x, y, z, strm
     end subroutine hip_lpt_periodic_bc_wrap_translational

     !> HIP kernel entry point for rotational periodic wrapping.
     subroutine hip_lpt_periodic_bc_wrap_rotational(x, y, z, n, theta_min, &
          theta_max, theta_len, u, v, w, u_lag, v_lag, w_lag, u_laglag, &
          v_laglag, w_laglag, acc_xlag, acc_ylag, acc_zlag, acc_xlaglag, &
          acc_ylaglag, acc_zlaglag, strm) &
          bind(c, name = 'hip_lpt_periodic_bc_wrap_rotational')
       use, intrinsic :: iso_c_binding
       use num_types, only : c_rp
       implicit none
       integer(c_int) :: n
       real(c_rp) :: theta_min, theta_max, theta_len
       type(c_ptr), value :: x, y, z, u, v, w, u_lag, v_lag, w_lag
       type(c_ptr), value :: u_laglag, v_laglag, w_laglag
       type(c_ptr), value :: acc_xlag, acc_ylag, acc_zlag
       type(c_ptr), value :: acc_xlaglag, acc_ylaglag, acc_zlaglag, strm
     end subroutine hip_lpt_periodic_bc_wrap_rotational
  end interface
#elif HAVE_CUDA
  interface
     !> CUDA kernel entry point for translational periodic wrapping.
     subroutine cuda_lpt_periodic_bc_wrap_translational(x, y, z, n, &
          n_periodic_dirs, periodic_dir_x1, periodic_dir_y1, periodic_dir_z1, &
          periodic_dir_x2, periodic_dir_y2, periodic_dir_z2, periodic_dir_x3, &
          periodic_dir_y3, periodic_dir_z3, periodic_min1, periodic_min2, &
          periodic_min3, periodic_max1, periodic_max2, periodic_max3, &
          periodic_shift_x1, periodic_shift_y1, periodic_shift_z1, &
          periodic_shift_x2, periodic_shift_y2, periodic_shift_z2, &
          periodic_shift_x3, periodic_shift_y3, periodic_shift_z3, &
          periodic_len1, periodic_len2, periodic_len3, strm) &
          bind(c, name = 'cuda_lpt_periodic_bc_wrap_translational')
       use, intrinsic :: iso_c_binding
       use num_types, only : c_rp
       implicit none
       integer(c_int) :: n, n_periodic_dirs
       real(c_rp) :: periodic_dir_x1, periodic_dir_y1, periodic_dir_z1
       real(c_rp) :: periodic_dir_x2, periodic_dir_y2, periodic_dir_z2
       real(c_rp) :: periodic_dir_x3, periodic_dir_y3, periodic_dir_z3
       real(c_rp) :: periodic_min1, periodic_min2, periodic_min3
       real(c_rp) :: periodic_max1, periodic_max2, periodic_max3
       real(c_rp) :: periodic_shift_x1, periodic_shift_y1, periodic_shift_z1
       real(c_rp) :: periodic_shift_x2, periodic_shift_y2, periodic_shift_z2
       real(c_rp) :: periodic_shift_x3, periodic_shift_y3, periodic_shift_z3
       real(c_rp) :: periodic_len1, periodic_len2, periodic_len3
       type(c_ptr), value :: x, y, z, strm
     end subroutine cuda_lpt_periodic_bc_wrap_translational

     !> CUDA kernel entry point for rotational periodic wrapping.
     subroutine cuda_lpt_periodic_bc_wrap_rotational(x, y, z, n, theta_min, &
          theta_max, theta_len, u, v, w, u_lag, v_lag, w_lag, u_laglag, &
          v_laglag, w_laglag, acc_xlag, acc_ylag, acc_zlag, acc_xlaglag, &
          acc_ylaglag, acc_zlaglag, strm) &
          bind(c, name = 'cuda_lpt_periodic_bc_wrap_rotational')
       use, intrinsic :: iso_c_binding
       use num_types, only : c_rp
       implicit none
       integer(c_int) :: n
       real(c_rp) :: theta_min, theta_max, theta_len
       type(c_ptr), value :: x, y, z, u, v, w, u_lag, v_lag, w_lag
       type(c_ptr), value :: u_laglag, v_laglag, w_laglag
       type(c_ptr), value :: acc_xlag, acc_ylag, acc_zlag
       type(c_ptr), value :: acc_xlaglag, acc_ylaglag, acc_zlaglag, strm
     end subroutine cuda_lpt_periodic_bc_wrap_rotational
  end interface
#endif

  public :: lpt_periodic_bc_wrap_translational_device
  public :: lpt_periodic_bc_wrap_rotational_device

contains

  !> Launch device kernel for translational periodic wrapping.
  !! @param x Particle x coordinates on device.
  !! @param y Particle y coordinates on device.
  !! @param z Particle z coordinates on device.
  !! @param n Number of local particles.
  !! @param n_periodic_dirs Number of translational periodic directions.
  subroutine lpt_periodic_bc_wrap_translational_device(x, y, z, n, &
       n_periodic_dirs, periodic_dir_x1, periodic_dir_y1, periodic_dir_z1, &
       periodic_dir_x2, periodic_dir_y2, periodic_dir_z2, periodic_dir_x3, &
       periodic_dir_y3, periodic_dir_z3, periodic_min1, periodic_min2, &
       periodic_min3, periodic_max1, periodic_max2, periodic_max3, &
       periodic_shift_x1, periodic_shift_y1, periodic_shift_z1, &
       periodic_shift_x2, periodic_shift_y2, periodic_shift_z2, &
       periodic_shift_x3, periodic_shift_y3, periodic_shift_z3, &
       periodic_len1, periodic_len2, periodic_len3, strm)
    type(vector_t), intent(inout) :: x, y, z
    integer, intent(in) :: n
    integer, intent(in) :: n_periodic_dirs
    real(kind=rp), intent(in) :: periodic_dir_x1, periodic_dir_y1
    real(kind=rp), intent(in) :: periodic_dir_z1, periodic_dir_x2
    real(kind=rp), intent(in) :: periodic_dir_y2, periodic_dir_z2
    real(kind=rp), intent(in) :: periodic_dir_x3, periodic_dir_y3
    real(kind=rp), intent(in) :: periodic_dir_z3
    real(kind=rp), intent(in) :: periodic_min1, periodic_min2, periodic_min3
    real(kind=rp), intent(in) :: periodic_max1, periodic_max2, periodic_max3
    real(kind=rp), intent(in) :: periodic_shift_x1, periodic_shift_y1
    real(kind=rp), intent(in) :: periodic_shift_z1, periodic_shift_x2
    real(kind=rp), intent(in) :: periodic_shift_y2, periodic_shift_z2
    real(kind=rp), intent(in) :: periodic_shift_x3, periodic_shift_y3
    real(kind=rp), intent(in) :: periodic_shift_z3
    real(kind=rp), intent(in) :: periodic_len1, periodic_len2, periodic_len3
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_

    if (n .lt. 1 .or. n_periodic_dirs .lt. 1) return
    strm_ = glb_cmd_queue
    if (present(strm)) strm_ = strm

#ifdef HAVE_HIP
    call hip_lpt_periodic_bc_wrap_translational(x%x_d, y%x_d, z%x_d, n, &
         n_periodic_dirs, periodic_dir_x1, periodic_dir_y1, periodic_dir_z1, &
         periodic_dir_x2, periodic_dir_y2, periodic_dir_z2, periodic_dir_x3, &
         periodic_dir_y3, periodic_dir_z3, periodic_min1, periodic_min2, &
         periodic_min3, periodic_max1, periodic_max2, periodic_max3, &
         periodic_shift_x1, periodic_shift_y1, periodic_shift_z1, &
         periodic_shift_x2, periodic_shift_y2, periodic_shift_z2, &
         periodic_shift_x3, periodic_shift_y3, periodic_shift_z3, &
         periodic_len1, periodic_len2, periodic_len3, strm_)
#elif HAVE_CUDA
    call cuda_lpt_periodic_bc_wrap_translational(x%x_d, y%x_d, z%x_d, n, &
         n_periodic_dirs, periodic_dir_x1, periodic_dir_y1, periodic_dir_z1, &
         periodic_dir_x2, periodic_dir_y2, periodic_dir_z2, periodic_dir_x3, &
         periodic_dir_y3, periodic_dir_z3, periodic_min1, periodic_min2, &
         periodic_min3, periodic_max1, periodic_max2, periodic_max3, &
         periodic_shift_x1, periodic_shift_y1, periodic_shift_z1, &
         periodic_shift_x2, periodic_shift_y2, periodic_shift_z2, &
         periodic_shift_x3, periodic_shift_y3, periodic_shift_z3, &
         periodic_len1, periodic_len2, periodic_len3, strm_)
#else
    call neko_error('LPT periodic BC device wrapping requires CUDA or HIP')
#endif
  end subroutine lpt_periodic_bc_wrap_translational_device

  !> Launch device kernel for rotational periodic wrapping.
  !! @param x Particle x coordinates on device.
  !! @param y Particle y coordinates on device.
  !! @param z Particle z coordinates on device.
  !! @param n Number of local particles.
  !! @param theta_min Minimum sector angle.
  !! @param theta_max Maximum sector angle.
  !! @param theta_len Angular sector length.
  subroutine lpt_periodic_bc_wrap_rotational_device(x, y, z, n, theta_min, &
       theta_max, theta_len, u, v, w, u_lag, v_lag, w_lag, u_laglag, &
       v_laglag, w_laglag, acc_xlag, acc_ylag, acc_zlag, acc_xlaglag, &
       acc_ylaglag, acc_zlaglag, strm)
    type(vector_t), intent(inout) :: x, y, z
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: theta_min
    real(kind=rp), intent(in) :: theta_max
    real(kind=rp), intent(in) :: theta_len
    type(vector_t), intent(inout), optional :: u, v, w
    type(vector_t), intent(inout), optional :: u_lag, v_lag, w_lag
    type(vector_t), intent(inout), optional :: u_laglag, v_laglag, w_laglag
    type(vector_t), intent(inout), optional :: acc_xlag, acc_ylag, acc_zlag
    type(vector_t), intent(inout), optional :: acc_xlaglag
    type(vector_t), intent(inout), optional :: acc_ylaglag
    type(vector_t), intent(inout), optional :: acc_zlaglag
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_
    type(c_ptr) :: u_d, v_d, w_d
    type(c_ptr) :: u_lag_d, v_lag_d, w_lag_d
    type(c_ptr) :: u_laglag_d, v_laglag_d, w_laglag_d
    type(c_ptr) :: acc_xlag_d, acc_ylag_d, acc_zlag_d
    type(c_ptr) :: acc_xlaglag_d, acc_ylaglag_d, acc_zlaglag_d

    if (n .lt. 1) return
    strm_ = glb_cmd_queue
    if (present(strm)) strm_ = strm

    u_d = C_NULL_PTR
    v_d = C_NULL_PTR
    w_d = C_NULL_PTR
    u_lag_d = C_NULL_PTR
    v_lag_d = C_NULL_PTR
    w_lag_d = C_NULL_PTR
    u_laglag_d = C_NULL_PTR
    v_laglag_d = C_NULL_PTR
    w_laglag_d = C_NULL_PTR
    acc_xlag_d = C_NULL_PTR
    acc_ylag_d = C_NULL_PTR
    acc_zlag_d = C_NULL_PTR
    acc_xlaglag_d = C_NULL_PTR
    acc_ylaglag_d = C_NULL_PTR
    acc_zlaglag_d = C_NULL_PTR

    if (present(u)) u_d = u%x_d
    if (present(v)) v_d = v%x_d
    if (present(w)) w_d = w%x_d
    if (present(u_lag)) u_lag_d = u_lag%x_d
    if (present(v_lag)) v_lag_d = v_lag%x_d
    if (present(w_lag)) w_lag_d = w_lag%x_d
    if (present(u_laglag)) u_laglag_d = u_laglag%x_d
    if (present(v_laglag)) v_laglag_d = v_laglag%x_d
    if (present(w_laglag)) w_laglag_d = w_laglag%x_d
    if (present(acc_xlag)) acc_xlag_d = acc_xlag%x_d
    if (present(acc_ylag)) acc_ylag_d = acc_ylag%x_d
    if (present(acc_zlag)) acc_zlag_d = acc_zlag%x_d
    if (present(acc_xlaglag)) acc_xlaglag_d = acc_xlaglag%x_d
    if (present(acc_ylaglag)) acc_ylaglag_d = acc_ylaglag%x_d
    if (present(acc_zlaglag)) acc_zlaglag_d = acc_zlaglag%x_d

#ifdef HAVE_HIP
    call hip_lpt_periodic_bc_wrap_rotational(x%x_d, y%x_d, z%x_d, n, &
         theta_min, theta_max, theta_len, u_d, v_d, w_d, u_lag_d, v_lag_d, &
         w_lag_d, u_laglag_d, v_laglag_d, w_laglag_d, acc_xlag_d, &
         acc_ylag_d, acc_zlag_d, acc_xlaglag_d, acc_ylaglag_d, &
         acc_zlaglag_d, strm_)
#elif HAVE_CUDA
    call cuda_lpt_periodic_bc_wrap_rotational(x%x_d, y%x_d, z%x_d, n, &
         theta_min, theta_max, theta_len, u_d, v_d, w_d, u_lag_d, v_lag_d, &
         w_lag_d, u_laglag_d, v_laglag_d, w_laglag_d, acc_xlag_d, &
         acc_ylag_d, acc_zlag_d, acc_xlaglag_d, acc_ylaglag_d, &
         acc_zlaglag_d, strm_)
#else
    call neko_error('LPT periodic BC device wrapping requires CUDA or HIP')
#endif
  end subroutine lpt_periodic_bc_wrap_rotational_device

end module lpt_periodic_bc_device
