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
module device_cartesian_el_finder
  use num_types, only : xp, c_xp
  use utils, only : neko_error
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int
  implicit none
  private

  public :: device_cartesian_el_finder_count
  public :: device_cartesian_el_finder_fill

#ifdef HAVE_CUDA
  interface
     subroutine cuda_cartesian_el_finder_count(points_d, el_map_offset_d, &
          point_box_d, n_el_cands_d, min_x, min_y, min_z, &
          x_res, y_res, z_res, n_boxes, n_points) &
          bind(c, name='cuda_cartesian_el_finder_count')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_xp
       type(c_ptr), value :: points_d, el_map_offset_d
       type(c_ptr), value :: point_box_d, n_el_cands_d
       real(c_xp), intent(in) :: min_x, min_y, min_z
       real(c_xp), intent(in) :: x_res, y_res, z_res
       integer(c_int), intent(in) :: n_boxes, n_points
     end subroutine cuda_cartesian_el_finder_count

     subroutine cuda_cartesian_el_finder_fill(point_box_d, &
          candidate_offsets_d, el_map_offset_d, el_map_data_d, &
          candidate_array_d, n_points) &
          bind(c, name='cuda_cartesian_el_finder_fill')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       type(c_ptr), value :: point_box_d, candidate_offsets_d
       type(c_ptr), value :: el_map_offset_d, el_map_data_d
       type(c_ptr), value :: candidate_array_d
       integer(c_int), intent(in) :: n_points
     end subroutine cuda_cartesian_el_finder_fill
  end interface
#elif HAVE_HIP
  interface
     subroutine hip_cartesian_el_finder_count(points_d, el_map_offset_d, &
          point_box_d, n_el_cands_d, min_x, min_y, min_z, &
          x_res, y_res, z_res, n_boxes, n_points) &
          bind(c, name='hip_cartesian_el_finder_count')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_xp
       type(c_ptr), value :: points_d, el_map_offset_d
       type(c_ptr), value :: point_box_d, n_el_cands_d
       real(c_xp), intent(in) :: min_x, min_y, min_z
       real(c_xp), intent(in) :: x_res, y_res, z_res
       integer(c_int), intent(in) :: n_boxes, n_points
     end subroutine hip_cartesian_el_finder_count

     subroutine hip_cartesian_el_finder_fill(point_box_d, &
          candidate_offsets_d, el_map_offset_d, el_map_data_d, &
          candidate_array_d, n_points) &
          bind(c, name='hip_cartesian_el_finder_fill')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       type(c_ptr), value :: point_box_d, candidate_offsets_d
       type(c_ptr), value :: el_map_offset_d, el_map_data_d
       type(c_ptr), value :: candidate_array_d
       integer(c_int), intent(in) :: n_points
     end subroutine hip_cartesian_el_finder_fill
  end interface
#elif HAVE_OPENCL
  interface
     subroutine opencl_cartesian_el_finder_count(points_d, el_map_offset_d, &
          point_box_d, n_el_cands_d, min_x, min_y, min_z, &
          x_res, y_res, z_res, n_boxes, n_points, xp_bytes) &
          bind(c, name='opencl_cartesian_el_finder_count')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_xp
       type(c_ptr), value :: points_d, el_map_offset_d
       type(c_ptr), value :: point_box_d, n_el_cands_d
       real(c_xp), intent(in) :: min_x, min_y, min_z
       real(c_xp), intent(in) :: x_res, y_res, z_res
       integer(c_int), intent(in) :: n_boxes, n_points, xp_bytes
     end subroutine opencl_cartesian_el_finder_count

     subroutine opencl_cartesian_el_finder_fill(point_box_d, &
          candidate_offsets_d, el_map_offset_d, el_map_data_d, &
          candidate_array_d, n_points) &
          bind(c, name='opencl_cartesian_el_finder_fill')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       type(c_ptr), value :: point_box_d, candidate_offsets_d
       type(c_ptr), value :: el_map_offset_d, el_map_data_d
       type(c_ptr), value :: candidate_array_d
       integer(c_int), intent(in) :: n_points
     end subroutine opencl_cartesian_el_finder_fill
  end interface
#elif HAVE_METAL
  interface
     subroutine metal_cartesian_el_finder_count(points_d, el_map_offset_d, &
          point_box_d, n_el_cands_d, min_x, min_y, min_z, &
          x_res, y_res, z_res, n_boxes, n_points, xp_bytes) &
          bind(c, name='metal_cartesian_el_finder_count')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       import c_xp
       type(c_ptr), value :: points_d, el_map_offset_d
       type(c_ptr), value :: point_box_d, n_el_cands_d
       real(c_xp), intent(in) :: min_x, min_y, min_z
       real(c_xp), intent(in) :: x_res, y_res, z_res
       integer(c_int), intent(in) :: n_boxes, n_points, xp_bytes
     end subroutine metal_cartesian_el_finder_count

     subroutine metal_cartesian_el_finder_fill(point_box_d, &
          candidate_offsets_d, el_map_offset_d, el_map_data_d, &
          candidate_array_d, n_points) &
          bind(c, name='metal_cartesian_el_finder_fill')
       use, intrinsic :: iso_c_binding, only : c_ptr, c_int
       type(c_ptr), value :: point_box_d, candidate_offsets_d
       type(c_ptr), value :: el_map_offset_d, el_map_data_d
       type(c_ptr), value :: candidate_array_d
       integer(c_int), intent(in) :: n_points
     end subroutine metal_cartesian_el_finder_fill
  end interface
#endif

contains

  subroutine device_cartesian_el_finder_count(points_d, el_map_offset_d, &
       point_box_d, n_el_cands_d, min_x, min_y, min_z, &
       x_res, y_res, z_res, n_boxes, n_points)
    type(c_ptr), intent(in) :: points_d, el_map_offset_d
    type(c_ptr), intent(in) :: point_box_d, n_el_cands_d
    real(kind=xp), intent(in) :: min_x, min_y, min_z
    real(kind=xp), intent(in) :: x_res, y_res, z_res
    integer, intent(in) :: n_boxes, n_points

#ifdef HAVE_CUDA
    call cuda_cartesian_el_finder_count(points_d, el_map_offset_d, &
         point_box_d, n_el_cands_d, min_x, min_y, min_z, &
         x_res, y_res, z_res, n_boxes, n_points)
#elif HAVE_HIP
    call hip_cartesian_el_finder_count(points_d, el_map_offset_d, &
         point_box_d, n_el_cands_d, min_x, min_y, min_z, &
         x_res, y_res, z_res, n_boxes, n_points)
#elif HAVE_OPENCL
    call opencl_cartesian_el_finder_count(points_d, el_map_offset_d, &
         point_box_d, n_el_cands_d, min_x, min_y, min_z, &
         x_res, y_res, z_res, n_boxes, n_points, storage_size(min_x)/8)
#elif HAVE_METAL
    call metal_cartesian_el_finder_count(points_d, el_map_offset_d, &
         point_box_d, n_el_cands_d, min_x, min_y, min_z, &
         x_res, y_res, z_res, n_boxes, n_points, storage_size(min_x)/8)
#else
    call neko_error('No Cartesian element finder device backend configured')
#endif
  end subroutine device_cartesian_el_finder_count

  subroutine device_cartesian_el_finder_fill(point_box_d, &
       candidate_offsets_d, el_map_offset_d, el_map_data_d, &
       candidate_array_d, n_points)
    type(c_ptr), intent(in) :: point_box_d, candidate_offsets_d
    type(c_ptr), intent(in) :: el_map_offset_d, el_map_data_d
    type(c_ptr), intent(in) :: candidate_array_d
    integer, intent(in) :: n_points

#ifdef HAVE_CUDA
    call cuda_cartesian_el_finder_fill(point_box_d, candidate_offsets_d, &
         el_map_offset_d, el_map_data_d, candidate_array_d, n_points)
#elif HAVE_HIP
    call hip_cartesian_el_finder_fill(point_box_d, candidate_offsets_d, &
         el_map_offset_d, el_map_data_d, candidate_array_d, n_points)
#elif HAVE_OPENCL
    call opencl_cartesian_el_finder_fill(point_box_d, candidate_offsets_d, &
         el_map_offset_d, el_map_data_d, candidate_array_d, n_points)
#elif HAVE_METAL
    call metal_cartesian_el_finder_fill(point_box_d, candidate_offsets_d, &
         el_map_offset_d, el_map_data_d, candidate_array_d, n_points)
#else
    call neko_error('No Cartesian element finder device backend configured')
#endif
  end subroutine device_cartesian_el_finder_fill

end module device_cartesian_el_finder
