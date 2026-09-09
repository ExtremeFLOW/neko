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

/**
 * Metal host-side dispatch for elastic particle-wall collisions.
 *
 * Metal allows at most 31 buffer arguments per function, so the lagged
 * state is reflected by a separate kernel. Hit detection depends on the
 * pre-update position, so the lag dispatches must precede the main one,
 * which overwrites it.
 *
 * @note Apple GPUs do not support FP64. This backend operates in FP32.
 */

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include <device/device_config.h>
#include <device/metal/kernel_utils.h>

/** Reflect one lagged velocity and acceleration triple */
static void lpt_wall_collision_reflect_lag(
    id<MTLCommandQueue> queue,
    void *wall_facet_mask, void *el_list,
    void *x_old, void *y_old, void *z_old,
    void *x, void *y, void *z, void *d,
    void *u_lag, void *v_lag, void *w_lag,
    void *acc_xlag, void *acc_ylag, void *acc_zlag,
    void *dm_x, void *dm_y, void *dm_z,
    void *nx, void *ny, void *nz,
    int n, int gdim, int nelv, int lx, int ly, int lz) {

  neko_metal_dispatch_1d_queue(
    queue, neko_metal_pipeline(@"lpt_wall_collision_reflect_lag_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) wall_facet_mask
              offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) el_list offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) x_old offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>) y_old offset:0 atIndex:3];
      [enc setBuffer:(__bridge id<MTLBuffer>) z_old offset:0 atIndex:4];
      [enc setBuffer:(__bridge id<MTLBuffer>) x offset:0 atIndex:5];
      [enc setBuffer:(__bridge id<MTLBuffer>) y offset:0 atIndex:6];
      [enc setBuffer:(__bridge id<MTLBuffer>) z offset:0 atIndex:7];
      [enc setBuffer:(__bridge id<MTLBuffer>) d offset:0 atIndex:8];
      [enc setBuffer:(__bridge id<MTLBuffer>) u_lag offset:0 atIndex:9];
      [enc setBuffer:(__bridge id<MTLBuffer>) v_lag offset:0 atIndex:10];
      [enc setBuffer:(__bridge id<MTLBuffer>) w_lag offset:0 atIndex:11];
      [enc setBuffer:(__bridge id<MTLBuffer>) acc_xlag offset:0 atIndex:12];
      [enc setBuffer:(__bridge id<MTLBuffer>) acc_ylag offset:0 atIndex:13];
      [enc setBuffer:(__bridge id<MTLBuffer>) acc_zlag offset:0 atIndex:14];
      [enc setBuffer:(__bridge id<MTLBuffer>) dm_x offset:0 atIndex:15];
      [enc setBuffer:(__bridge id<MTLBuffer>) dm_y offset:0 atIndex:16];
      [enc setBuffer:(__bridge id<MTLBuffer>) dm_z offset:0 atIndex:17];
      [enc setBuffer:(__bridge id<MTLBuffer>) nx offset:0 atIndex:18];
      [enc setBuffer:(__bridge id<MTLBuffer>) ny offset:0 atIndex:19];
      [enc setBuffer:(__bridge id<MTLBuffer>) nz offset:0 atIndex:20];
      [enc setBytes:&n length:sizeof(int) atIndex:21];
      [enc setBytes:&gdim length:sizeof(int) atIndex:22];
      [enc setBytes:&nelv length:sizeof(int) atIndex:23];
      [enc setBytes:&lx length:sizeof(int) atIndex:24];
      [enc setBytes:&ly length:sizeof(int) atIndex:25];
      [enc setBytes:&lz length:sizeof(int) atIndex:26];
    }, (NSUInteger) n);
}

/** Fortran wrapper for the Metal elastic particle-wall collision kernels */
void metal_lpt_handle_elastic_wall_collisions(
    void *wall_facet_mask, void *el_list,
    void *x_old, void *y_old, void *z_old,
    void *x, void *y, void *z,
    void *d, void *u, void *v, void *w,
    void *u_lag, void *v_lag, void *w_lag,
    void *u_laglag, void *v_laglag, void *w_laglag,
    void *acc_xlag, void *acc_ylag, void *acc_zlag,
    void *acc_xlaglag, void *acc_ylaglag, void *acc_zlaglag,
    void *u_old, void *v_old, void *w_old,
    void *acc_x, void *acc_y, void *acc_z,
    void *dm_x, void *dm_y, void *dm_z,
    void *nx, void *ny, void *nz,
    int *n, int *gdim, int *nelv,
    int *lx, int *ly, int *lz, int *lag_len, void *strm) {
  if (*n < 1) return;

  const int n_r = *n;
  const int gdim_r = *gdim;
  const int nelv_r = *nelv;
  const int lx_r = *lx;
  const int ly_r = *ly;
  const int lz_r = *lz;

  id<MTLCommandQueue> queue = (__bridge id<MTLCommandQueue>) strm;

  /* The lagged state must be reflected before the position is updated */
  if (*lag_len >= 2) {
    lpt_wall_collision_reflect_lag(queue, wall_facet_mask, el_list,
                                   x_old, y_old, z_old, x, y, z, d,
                                   u_laglag, v_laglag, w_laglag,
                                   acc_xlaglag, acc_ylaglag, acc_zlaglag,
                                   dm_x, dm_y, dm_z, nx, ny, nz,
                                   n_r, gdim_r, nelv_r, lx_r, ly_r, lz_r);
  }

  if (*lag_len >= 1) {
    lpt_wall_collision_reflect_lag(queue, wall_facet_mask, el_list,
                                   x_old, y_old, z_old, x, y, z, d,
                                   u_lag, v_lag, w_lag,
                                   acc_xlag, acc_ylag, acc_zlag,
                                   dm_x, dm_y, dm_z, nx, ny, nz,
                                   n_r, gdim_r, nelv_r, lx_r, ly_r, lz_r);
  }

  neko_metal_dispatch_1d_queue(
    queue,
    neko_metal_pipeline(@"lpt_handle_elastic_wall_collisions_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) wall_facet_mask
              offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) el_list offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) x_old offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>) y_old offset:0 atIndex:3];
      [enc setBuffer:(__bridge id<MTLBuffer>) z_old offset:0 atIndex:4];
      [enc setBuffer:(__bridge id<MTLBuffer>) x offset:0 atIndex:5];
      [enc setBuffer:(__bridge id<MTLBuffer>) y offset:0 atIndex:6];
      [enc setBuffer:(__bridge id<MTLBuffer>) z offset:0 atIndex:7];
      [enc setBuffer:(__bridge id<MTLBuffer>) d offset:0 atIndex:8];
      [enc setBuffer:(__bridge id<MTLBuffer>) u offset:0 atIndex:9];
      [enc setBuffer:(__bridge id<MTLBuffer>) v offset:0 atIndex:10];
      [enc setBuffer:(__bridge id<MTLBuffer>) w offset:0 atIndex:11];
      [enc setBuffer:(__bridge id<MTLBuffer>) u_old offset:0 atIndex:12];
      [enc setBuffer:(__bridge id<MTLBuffer>) v_old offset:0 atIndex:13];
      [enc setBuffer:(__bridge id<MTLBuffer>) w_old offset:0 atIndex:14];
      [enc setBuffer:(__bridge id<MTLBuffer>) acc_x offset:0 atIndex:15];
      [enc setBuffer:(__bridge id<MTLBuffer>) acc_y offset:0 atIndex:16];
      [enc setBuffer:(__bridge id<MTLBuffer>) acc_z offset:0 atIndex:17];
      [enc setBuffer:(__bridge id<MTLBuffer>) dm_x offset:0 atIndex:18];
      [enc setBuffer:(__bridge id<MTLBuffer>) dm_y offset:0 atIndex:19];
      [enc setBuffer:(__bridge id<MTLBuffer>) dm_z offset:0 atIndex:20];
      [enc setBuffer:(__bridge id<MTLBuffer>) nx offset:0 atIndex:21];
      [enc setBuffer:(__bridge id<MTLBuffer>) ny offset:0 atIndex:22];
      [enc setBuffer:(__bridge id<MTLBuffer>) nz offset:0 atIndex:23];
      [enc setBytes:&n_r length:sizeof(int) atIndex:24];
      [enc setBytes:&gdim_r length:sizeof(int) atIndex:25];
      [enc setBytes:&nelv_r length:sizeof(int) atIndex:26];
      [enc setBytes:&lx_r length:sizeof(int) atIndex:27];
      [enc setBytes:&ly_r length:sizeof(int) atIndex:28];
      [enc setBytes:&lz_r length:sizeof(int) atIndex:29];
    }, (NSUInteger) n_r);
}
