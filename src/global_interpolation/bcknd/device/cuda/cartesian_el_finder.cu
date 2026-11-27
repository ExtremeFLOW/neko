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

#include <device/device_config.h>
#include <device/cuda/check.h>
#include <global_interpolation/bcknd/device/cuda/cartesian_el_finder.h>

extern "C" void cuda_cartesian_el_finder_count(
    const void *points, const void *el_map_offset,
    void *point_box, void *n_el_cands,
    const real_xp *min_x, const real_xp *min_y, const real_xp *min_z,
    const real_xp *x_res, const real_xp *y_res, const real_xp *z_res,
    const int *n_boxes, const int *n_points) {
  const int threads = 256;
  const int blocks = (*n_points + threads - 1) / threads;
  const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;
  cartesian_el_finder_count_kernel<<<blocks, threads, 0, stream>>>(
      (const real *) points, (const int *) el_map_offset,
      (int *) point_box, (int *) n_el_cands,
      *min_x, *min_y, *min_z, *x_res, *y_res, *z_res,
      *n_boxes, *n_points);
  CUDA_CHECK(cudaGetLastError());
}

extern "C" void cuda_cartesian_el_finder_fill(
    const void *point_box, const void *candidate_offsets,
    const void *el_map_offset, const void *el_map_data,
    void *candidate_array, const int *n_points) {
  const int threads = 256;
  const int blocks = (*n_points + threads - 1) / threads;
  const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;
  cartesian_el_finder_fill_kernel<<<blocks, threads, 0, stream>>>(
      (const int *) point_box, (const int *) candidate_offsets,
      (const int *) el_map_offset, (const int *) el_map_data,
      (int *) candidate_array, *n_points);
  CUDA_CHECK(cudaGetLastError());
}
