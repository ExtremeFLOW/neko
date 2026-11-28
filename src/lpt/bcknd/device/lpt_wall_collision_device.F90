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
!> Device dispatch for the LPT elastic wall collisions with the high-order
!! wall geometry. CUDA and HIP use one kernel; OpenCL and Metal (at most 31
!! buffers per Metal function) use a collision-search kernel that stores
!! the contact normals and a second kernel that reflects the velocity and
!! its history.
module lpt_wall_collision_device
  use vector, only : vector_t
  use utils, only : neko_error
  use device, only : glb_cmd_queue
  use num_types, only : rp, c_rp
  use lpt_wall_geometry, only : lpt_wall_geometry_t
  use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_associated
  implicit none
  private

#ifdef HAVE_HIP
  interface
     !> HIP kernel entry point for the high-order wall collisions.
     subroutine hip_lpt_wall_collide(xf, fr, fi, zg, el_ptr, el_facets, &
          el_radius, lx, nelv, el_list, rst, x_old, y_old, z_old, x, y, z, &
          d, u, v, w, u_lag, v_lag, w_lag, u_laglag, v_laglag, w_laglag, &
          acc_xlag, acc_ylag, acc_zlag, acc_xlaglag, acc_ylaglag, &
          acc_zlaglag, u_old, v_old, w_old, acc_x, acc_y, acc_z, lag_len, &
          e_n, n, stats, strm) bind(c, name = 'hip_lpt_wall_collide')
       use, intrinsic :: iso_c_binding
       import :: c_rp
       implicit none
       type(c_ptr), value :: xf, fr, fi, zg, el_ptr, el_facets, el_radius
       type(c_ptr), value :: el_list, rst
       type(c_ptr), value :: x_old, y_old, z_old, x, y, z, d, u, v, w
       type(c_ptr), value :: u_lag, v_lag, w_lag, u_laglag, v_laglag
       type(c_ptr), value :: w_laglag, acc_xlag, acc_ylag, acc_zlag
       type(c_ptr), value :: acc_xlaglag, acc_ylaglag, acc_zlaglag
       type(c_ptr), value :: u_old, v_old, w_old, acc_x, acc_y, acc_z
       type(c_ptr), value :: stats, strm
       integer(c_int) :: lx, nelv, lag_len, n
       real(kind=c_rp) :: e_n
     end subroutine hip_lpt_wall_collide
  end interface
#elif HAVE_CUDA
  interface
     !> CUDA kernel entry point for the high-order wall collisions.
     subroutine cuda_lpt_wall_collide(xf, fr, fi, zg, el_ptr, el_facets, &
          el_radius, lx, nelv, el_list, rst, x_old, y_old, z_old, x, y, z, &
          d, u, v, w, u_lag, v_lag, w_lag, u_laglag, v_laglag, w_laglag, &
          acc_xlag, acc_ylag, acc_zlag, acc_xlaglag, acc_ylaglag, &
          acc_zlaglag, u_old, v_old, w_old, acc_x, acc_y, acc_z, lag_len, &
          e_n, n, stats, strm) bind(c, name = 'cuda_lpt_wall_collide')
       use, intrinsic :: iso_c_binding
       import :: c_rp
       implicit none
       type(c_ptr), value :: xf, fr, fi, zg, el_ptr, el_facets, el_radius
       type(c_ptr), value :: el_list, rst
       type(c_ptr), value :: x_old, y_old, z_old, x, y, z, d, u, v, w
       type(c_ptr), value :: u_lag, v_lag, w_lag, u_laglag, v_laglag
       type(c_ptr), value :: w_laglag, acc_xlag, acc_ylag, acc_zlag
       type(c_ptr), value :: acc_xlaglag, acc_ylaglag, acc_zlaglag
       type(c_ptr), value :: u_old, v_old, w_old, acc_x, acc_y, acc_z
       type(c_ptr), value :: stats, strm
       integer(c_int) :: lx, nelv, lag_len, n
       real(kind=c_rp) :: e_n
     end subroutine cuda_lpt_wall_collide
  end interface
#elif HAVE_OPENCL
  interface
     !> OPENCL entry point of the collision search.
     subroutine opencl_lpt_wall_collide(xf, fr, fi, zg, el_ptr, el_facets, &
          el_radius, lx, nelv, el_list, rst, x_old, y_old, z_old, x, y, z, &
          d, hits, normals, stats, n, e_n, has_rst, strm) &
          bind(c, name = 'opencl_lpt_wall_collide')
       use, intrinsic :: iso_c_binding
       import :: c_rp
       implicit none
       type(c_ptr), value :: xf, fr, fi, zg, el_ptr, el_facets, el_radius
       type(c_ptr), value :: el_list, rst, x_old, y_old, z_old, x, y, z, d
       type(c_ptr), value :: hits, normals, stats, strm
       integer(c_int) :: lx, nelv, n, has_rst
       real(kind=c_rp) :: e_n
     end subroutine opencl_lpt_wall_collide

     !> OPENCL entry point of the reflection of velocity and history.
     subroutine opencl_lpt_wall_reflect(hits, normals, u, v, w, u_old, v_old, &
          w_old, acc_x, acc_y, acc_z, u_lag, v_lag, w_lag, acc_xlag, &
          acc_ylag, acc_zlag, u_laglag, v_laglag, w_laglag, acc_xlaglag, &
          acc_ylaglag, acc_zlaglag, n, lag_len, e_n, strm) &
          bind(c, name = 'opencl_lpt_wall_reflect')
       use, intrinsic :: iso_c_binding
       import :: c_rp
       implicit none
       type(c_ptr), value :: hits, normals, u, v, w, u_old, v_old, w_old
       type(c_ptr), value :: acc_x, acc_y, acc_z, u_lag, v_lag, w_lag
       type(c_ptr), value :: acc_xlag, acc_ylag, acc_zlag, u_laglag
       type(c_ptr), value :: v_laglag, w_laglag, acc_xlaglag, acc_ylaglag
       type(c_ptr), value :: acc_zlaglag, strm
       integer(c_int) :: n, lag_len
       real(kind=c_rp) :: e_n
     end subroutine opencl_lpt_wall_reflect
  end interface
#elif HAVE_METAL
  interface
     !> METAL entry point of the collision search.
     subroutine metal_lpt_wall_collide(xf, fr, fi, zg, el_ptr, el_facets, &
          el_radius, lx, nelv, el_list, rst, x_old, y_old, z_old, x, y, z, &
          d, hits, normals, stats, n, e_n, has_rst, strm) &
          bind(c, name = 'metal_lpt_wall_collide')
       use, intrinsic :: iso_c_binding
       import :: c_rp
       implicit none
       type(c_ptr), value :: xf, fr, fi, zg, el_ptr, el_facets, el_radius
       type(c_ptr), value :: el_list, rst, x_old, y_old, z_old, x, y, z, d
       type(c_ptr), value :: hits, normals, stats, strm
       integer(c_int) :: lx, nelv, n, has_rst
       real(kind=c_rp) :: e_n
     end subroutine metal_lpt_wall_collide

     !> METAL entry point of the reflection of velocity and history.
     subroutine metal_lpt_wall_reflect(hits, normals, u, v, w, u_old, v_old, &
          w_old, acc_x, acc_y, acc_z, u_lag, v_lag, w_lag, acc_xlag, &
          acc_ylag, acc_zlag, u_laglag, v_laglag, w_laglag, acc_xlaglag, &
          acc_ylaglag, acc_zlaglag, n, lag_len, e_n, strm) &
          bind(c, name = 'metal_lpt_wall_reflect')
       use, intrinsic :: iso_c_binding
       import :: c_rp
       implicit none
       type(c_ptr), value :: hits, normals, u, v, w, u_old, v_old, w_old
       type(c_ptr), value :: acc_x, acc_y, acc_z, u_lag, v_lag, w_lag
       type(c_ptr), value :: acc_xlag, acc_ylag, acc_zlag, u_laglag
       type(c_ptr), value :: v_laglag, w_laglag, acc_xlaglag, acc_ylaglag
       type(c_ptr), value :: acc_zlaglag, strm
       integer(c_int) :: n, lag_len
       real(kind=c_rp) :: e_n
     end subroutine metal_lpt_wall_reflect
  end interface
#endif

  public :: lpt_wall_collide_device

contains

  !> Detect and resolve wall collisions of all local particles on the
  !! device. The geometry must have been copied with
  !! lpt_wall_geometry_t%to_device.
  !! @param g High-order wall geometry with device copies.
  !! @param el_list_d Device list of local owner elements (0-based).
  !! @param rst_d Device reference coordinates (3, n) of the start
  !! positions in those elements, used as initial guesses (may be null).
  !! @param e_n Normal coefficient of restitution.
  !! @param lag_len Number of history levels in use.
  !! @param n Number of local particles.
  subroutine lpt_wall_collide_device(g, el_list_d, rst_d, x_old, y_old, &
       z_old, x, y, z, d, u, v, w, u_lag, v_lag, w_lag, u_laglag, v_laglag, &
       w_laglag, acc_xlag, acc_ylag, acc_zlag, acc_xlaglag, acc_ylaglag, &
       acc_zlaglag, u_old, v_old, w_old, acc_x, acc_y, acc_z, e_n, lag_len, &
       n, strm)
    type(lpt_wall_geometry_t), intent(inout) :: g
    type(c_ptr), intent(in) :: el_list_d, rst_d
    type(vector_t), intent(in) :: x_old, y_old, z_old
    type(vector_t), intent(inout) :: x, y, z
    type(vector_t), intent(in) :: d
    type(vector_t), intent(inout) :: u, v, w
    type(vector_t), intent(inout) :: u_lag, v_lag, w_lag
    type(vector_t), intent(inout) :: u_laglag, v_laglag, w_laglag
    type(vector_t), intent(inout) :: acc_xlag, acc_ylag, acc_zlag
    type(vector_t), intent(inout) :: acc_xlaglag, acc_ylaglag
    type(vector_t), intent(inout) :: acc_zlaglag
    type(vector_t), intent(inout) :: u_old, v_old, w_old
    type(vector_t), intent(inout) :: acc_x, acc_y, acc_z
    real(kind=rp), intent(in) :: e_n
    integer, intent(in) :: lag_len
    integer, intent(in) :: n
    type(c_ptr), optional :: strm
    type(c_ptr) :: strm_
    integer(c_int) :: lx_, nelv_, lag_len_, n_, has_rst_
    real(kind=c_rp) :: e_n_

    if (n .lt. 1) return
    strm_ = glb_cmd_queue
    if (present(strm)) strm_ = strm
    lx_ = int(g%lx, c_int)
    nelv_ = int(size(g%el_radius), c_int)
    lag_len_ = int(lag_len, c_int)
    n_ = int(n, c_int)
    e_n_ = real(e_n, c_rp)
    has_rst_ = merge(1_c_int, 0_c_int, c_associated(rst_d))

#ifdef HAVE_HIP
    call hip_lpt_wall_collide(g%xf_d, g%fr_d, g%fi_d, g%zg_d, &
         g%el_ptr_d, g%el_facets_d, g%el_radius_d, lx_, nelv_, el_list_d, &
         rst_d, x_old%x_d, y_old%x_d, z_old%x_d, x%x_d, y%x_d, z%x_d, &
         d%x_d, u%x_d, v%x_d, w%x_d, u_lag%x_d, v_lag%x_d, w_lag%x_d, &
         u_laglag%x_d, v_laglag%x_d, w_laglag%x_d, acc_xlag%x_d, &
         acc_ylag%x_d, acc_zlag%x_d, acc_xlaglag%x_d, acc_ylaglag%x_d, &
         acc_zlaglag%x_d, u_old%x_d, v_old%x_d, w_old%x_d, acc_x%x_d, &
         acc_y%x_d, acc_z%x_d, lag_len_, e_n_, n_, g%stats_d, strm_)
#elif HAVE_CUDA
    call cuda_lpt_wall_collide(g%xf_d, g%fr_d, g%fi_d, g%zg_d, &
         g%el_ptr_d, g%el_facets_d, g%el_radius_d, lx_, nelv_, el_list_d, &
         rst_d, x_old%x_d, y_old%x_d, z_old%x_d, x%x_d, y%x_d, z%x_d, &
         d%x_d, u%x_d, v%x_d, w%x_d, u_lag%x_d, v_lag%x_d, w_lag%x_d, &
         u_laglag%x_d, v_laglag%x_d, w_laglag%x_d, acc_xlag%x_d, &
         acc_ylag%x_d, acc_zlag%x_d, acc_xlaglag%x_d, acc_ylaglag%x_d, &
         acc_zlaglag%x_d, u_old%x_d, v_old%x_d, w_old%x_d, acc_x%x_d, &
         acc_y%x_d, acc_z%x_d, lag_len_, e_n_, n_, g%stats_d, strm_)
#elif HAVE_OPENCL
    call g%reserve_scratch(n)
    call opencl_lpt_wall_collide(g%xf_d, g%fr_d, g%fi_d, g%zg_d, &
         g%el_ptr_d, g%el_facets_d, g%el_radius_d, lx_, nelv_, el_list_d, &
         rst_d, x_old%x_d, y_old%x_d, z_old%x_d, x%x_d, y%x_d, z%x_d, &
         d%x_d, g%hits_d, g%normals_d, g%stats_d, n_, e_n_, has_rst_, &
         strm_)
    call opencl_lpt_wall_reflect(g%hits_d, g%normals_d, u%x_d, v%x_d, &
         w%x_d, u_old%x_d, v_old%x_d, w_old%x_d, acc_x%x_d, acc_y%x_d, &
         acc_z%x_d, u_lag%x_d, v_lag%x_d, w_lag%x_d, acc_xlag%x_d, &
         acc_ylag%x_d, acc_zlag%x_d, u_laglag%x_d, v_laglag%x_d, &
         w_laglag%x_d, acc_xlaglag%x_d, acc_ylaglag%x_d, acc_zlaglag%x_d, &
         n_, lag_len_, e_n_, strm_)
#elif HAVE_METAL
    call g%reserve_scratch(n)
    call metal_lpt_wall_collide(g%xf_d, g%fr_d, g%fi_d, g%zg_d, &
         g%el_ptr_d, g%el_facets_d, g%el_radius_d, lx_, nelv_, el_list_d, &
         rst_d, x_old%x_d, y_old%x_d, z_old%x_d, x%x_d, y%x_d, z%x_d, &
         d%x_d, g%hits_d, g%normals_d, g%stats_d, n_, e_n_, has_rst_, &
         strm_)
    call metal_lpt_wall_reflect(g%hits_d, g%normals_d, u%x_d, v%x_d, &
         w%x_d, u_old%x_d, v_old%x_d, w_old%x_d, acc_x%x_d, acc_y%x_d, &
         acc_z%x_d, u_lag%x_d, v_lag%x_d, w_lag%x_d, acc_xlag%x_d, &
         acc_ylag%x_d, acc_zlag%x_d, u_laglag%x_d, v_laglag%x_d, &
         w_laglag%x_d, acc_xlaglag%x_d, acc_ylaglag%x_d, acc_zlaglag%x_d, &
         n_, lag_len_, e_n_, strm_)
#else
    call neko_error('No device backend configured')
#endif
    call g%gather_stats()
  end subroutine lpt_wall_collide_device

end module lpt_wall_collision_device
