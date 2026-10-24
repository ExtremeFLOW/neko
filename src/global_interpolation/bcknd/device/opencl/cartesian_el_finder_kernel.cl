/*
 Copyright (c) 2026, The Neko Authors
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

   * Redistributions in binary form must reproduce the above
     copyright notice, this list of conditions and the following
     disclaimer in the documentation and/or other materials provided
     with the distribution.

   * Neither the name of the authors nor the names of its
     contributors may be used to endorse or promote products derived
     from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
*/

__kernel void cartesian_el_finder_count_kernel(
    __global const POINT_T *points, __global const int *el_map_offset,
    __global int *point_box, __global int *n_el_cands,
    const XP_T min_x, const XP_T min_y, const XP_T min_z,
    const XP_T x_res, const XP_T y_res, const XP_T z_res,
    const int n_boxes, const int n_points) {
  const int p = get_global_id(0);
  if (p >= n_points) return;

  int x_id = (int)(((XP_T) points[3*p] - min_x) / x_res);
  int y_id = (int)(((XP_T) points[3*p + 1] - min_y) / y_res);
  int z_id = (int)(((XP_T) points[3*p + 2] - min_z) / z_res);

  if (x_id == -1) x_id = 0;
  if (x_id == n_boxes) x_id = n_boxes - 1;
  if (y_id == -1) y_id = 0;
  if (y_id == n_boxes) y_id = n_boxes - 1;
  if (z_id == -1) z_id = 0;
  if (z_id == n_boxes) z_id = n_boxes - 1;

  if (x_id < 0 || x_id >= n_boxes ||
      y_id < 0 || y_id >= n_boxes ||
      z_id < 0 || z_id >= n_boxes) {
    point_box[p] = -1;
    n_el_cands[p] = 0;
    return;
  }

  const int box = x_id + n_boxes * (y_id + n_boxes * z_id);
  point_box[p] = box;
  n_el_cands[p] = el_map_offset[box + 1] - el_map_offset[box];
}

__kernel void cartesian_el_finder_fill_kernel(
    __global const int *point_box, __global const int *candidate_offsets,
    __global const int *el_map_offset, __global const int *el_map_data,
    __global int *candidate_array, const int n_points) {
  const int p = get_global_id(0);
  if (p >= n_points) return;

  const int box = point_box[p];
  if (box < 0) return;

  const int src_begin = el_map_offset[box];
  const int src_end = el_map_offset[box + 1];
  const int dst_begin = candidate_offsets[p];
  for (int j = 0; j < src_end - src_begin; ++j) {
    candidate_array[dst_begin + j] = el_map_data[src_begin + j] - 1;
  }
}
