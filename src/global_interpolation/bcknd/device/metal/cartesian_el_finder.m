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

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#include <stdlib.h>
#include <device/metal/kernel_utils.h>

static float metal_xp_to_float(const void *value, int bytes) {
  if (bytes == sizeof(float)) return *(const float *)value;
  if (bytes == sizeof(double)) return (float)*(const double *)value;
  abort();
}

void metal_cartesian_el_finder_count(
    void *points, void *el_map_offset, void *point_box, void *n_el_cands,
    const void *min_x, const void *min_y, const void *min_z,
    const void *x_res, const void *y_res, const void *z_res,
    const int *n_boxes, const int *n_points, const int *xp_bytes) {
  if (*n_points <= 0) return;

  const float min_x_r = metal_xp_to_float(min_x, *xp_bytes);
  const float min_y_r = metal_xp_to_float(min_y, *xp_bytes);
  const float min_z_r = metal_xp_to_float(min_z, *xp_bytes);
  const float x_res_r = metal_xp_to_float(x_res, *xp_bytes);
  const float y_res_r = metal_xp_to_float(y_res, *xp_bytes);
  const float z_res_r = metal_xp_to_float(z_res, *xp_bytes);
  const int n_boxes_r = *n_boxes;
  const int n_points_r = *n_points;

  neko_metal_dispatch_1d(
    neko_metal_pipeline(@"cartesian_el_finder_count_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>)points offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>)el_map_offset offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>)point_box offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>)n_el_cands offset:0 atIndex:3];
      [enc setBytes:&min_x_r length:sizeof(float) atIndex:4];
      [enc setBytes:&min_y_r length:sizeof(float) atIndex:5];
      [enc setBytes:&min_z_r length:sizeof(float) atIndex:6];
      [enc setBytes:&x_res_r length:sizeof(float) atIndex:7];
      [enc setBytes:&y_res_r length:sizeof(float) atIndex:8];
      [enc setBytes:&z_res_r length:sizeof(float) atIndex:9];
      [enc setBytes:&n_boxes_r length:sizeof(int) atIndex:10];
      [enc setBytes:&n_points_r length:sizeof(int) atIndex:11];
    }, (NSUInteger)n_points_r);
}

void metal_cartesian_el_finder_fill(
    void *point_box, void *candidate_offsets, void *el_map_offset,
    void *el_map_data, void *candidate_array, const int *n_points) {
  if (*n_points <= 0) return;
  const int n_points_r = *n_points;

  neko_metal_dispatch_1d(
    neko_metal_pipeline(@"cartesian_el_finder_fill_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>)point_box offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>)candidate_offsets offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>)el_map_offset offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>)el_map_data offset:0 atIndex:3];
      [enc setBuffer:(__bridge id<MTLBuffer>)candidate_array offset:0 atIndex:4];
      [enc setBytes:&n_points_r length:sizeof(int) atIndex:5];
    }, (NSUInteger)n_points_r);
}
