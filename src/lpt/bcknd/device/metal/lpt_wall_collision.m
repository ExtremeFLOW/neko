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
 * Metal host-side dispatch for elastic collisions of particles with the
 * high-order wall facets.
 *
 * Metal allows at most 31 buffer arguments per function, so the work is
 * split into the collision search, which corrects the positions and stores
 * the contact normals, and the reflection of the velocity and its history.
 *
 * @note Apple GPUs do not support FP64. This backend operates in FP32.
 */

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include <device/device_config.h>
#include <device/metal/kernel_utils.h>

/** Fortran wrapper for the Metal collision search */
void metal_lpt_wall_collide(void *xf, void *fr, void *fi, void *zg,
                            void *el_ptr, void *el_facets, void *el_radius,
                            int *lx, int *nelv, void *el_list, void *rst,
                            void *x_old, void *y_old, void *z_old,
                            void *x, void *y, void *z, void *d,
                            void *hits, void *normals, void *stats,
                            int *n, float *e_n, int *has_rst, void *strm) {
  if (*n < 1) return;

  const int lx_r = *lx;
  const int nelv_r = *nelv;
  const int n_r = *n;
  const float e_n_r = *e_n;
  const int has_rst_r = *has_rst;
  /* rst is only read when has_rst is set, but needs a valid buffer */
  if (rst == NULL) rst = x;

  id<MTLCommandQueue> queue = (__bridge id<MTLCommandQueue>) strm;

  neko_metal_dispatch_1d_queue(
    queue, neko_metal_pipeline(@"lpt_wall_collide_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) xf offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) fr offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) fi offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>) zg offset:0 atIndex:3];
      [enc setBuffer:(__bridge id<MTLBuffer>) el_ptr offset:0 atIndex:4];
      [enc setBuffer:(__bridge id<MTLBuffer>) el_facets offset:0 atIndex:5];
      [enc setBuffer:(__bridge id<MTLBuffer>) el_radius offset:0 atIndex:6];
      [enc setBuffer:(__bridge id<MTLBuffer>) el_list offset:0 atIndex:7];
      [enc setBuffer:(__bridge id<MTLBuffer>) rst offset:0 atIndex:8];
      [enc setBuffer:(__bridge id<MTLBuffer>) x_old offset:0 atIndex:9];
      [enc setBuffer:(__bridge id<MTLBuffer>) y_old offset:0 atIndex:10];
      [enc setBuffer:(__bridge id<MTLBuffer>) z_old offset:0 atIndex:11];
      [enc setBuffer:(__bridge id<MTLBuffer>) x offset:0 atIndex:12];
      [enc setBuffer:(__bridge id<MTLBuffer>) y offset:0 atIndex:13];
      [enc setBuffer:(__bridge id<MTLBuffer>) z offset:0 atIndex:14];
      [enc setBuffer:(__bridge id<MTLBuffer>) d offset:0 atIndex:15];
      [enc setBuffer:(__bridge id<MTLBuffer>) hits offset:0 atIndex:16];
      [enc setBuffer:(__bridge id<MTLBuffer>) normals offset:0 atIndex:17];
      [enc setBuffer:(__bridge id<MTLBuffer>) stats offset:0 atIndex:18];
      [enc setBytes:&lx_r length:sizeof(int) atIndex:19];
      [enc setBytes:&nelv_r length:sizeof(int) atIndex:20];
      [enc setBytes:&n_r length:sizeof(int) atIndex:21];
      [enc setBytes:&e_n_r length:sizeof(float) atIndex:22];
      [enc setBytes:&has_rst_r length:sizeof(int) atIndex:23];
    }, (NSUInteger) n_r);
}

/** Fortran wrapper for the Metal reflection of the velocity and history */
void metal_lpt_wall_reflect(void *hits, void *normals,
                            void *u, void *v, void *w,
                            void *u_old, void *v_old, void *w_old,
                            void *acc_x, void *acc_y, void *acc_z,
                            void *u_lag, void *v_lag, void *w_lag,
                            void *acc_xlag, void *acc_ylag, void *acc_zlag,
                            void *u_laglag, void *v_laglag, void *w_laglag,
                            void *acc_xlaglag, void *acc_ylaglag,
                            void *acc_zlaglag, int *n, int *lag_len,
                            float *e_n, void *strm) {
  if (*n < 1) return;

  const int n_r = *n;
  const int lag_len_r = *lag_len;
  const float e_n_r = *e_n;
  void *lag[12] = { u_lag, v_lag, w_lag, acc_xlag, acc_ylag, acc_zlag,
                    u_laglag, v_laglag, w_laglag, acc_xlaglag, acc_ylaglag,
                    acc_zlaglag };
  /* The lag arrays are only read when lag_len allows it, but every slot
     still needs a valid buffer */
  for (int i = 0; i < 12; i++) {
    if (lag[i] == NULL) lag[i] = u;
  }
  u_lag = lag[0]; v_lag = lag[1]; w_lag = lag[2];
  acc_xlag = lag[3]; acc_ylag = lag[4]; acc_zlag = lag[5];
  u_laglag = lag[6]; v_laglag = lag[7]; w_laglag = lag[8];
  acc_xlaglag = lag[9]; acc_ylaglag = lag[10]; acc_zlaglag = lag[11];

  id<MTLCommandQueue> queue = (__bridge id<MTLCommandQueue>) strm;

  neko_metal_dispatch_1d_queue(
    queue, neko_metal_pipeline(@"lpt_wall_reflect_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) hits offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) normals offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) u offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>) v offset:0 atIndex:3];
      [enc setBuffer:(__bridge id<MTLBuffer>) w offset:0 atIndex:4];
      [enc setBuffer:(__bridge id<MTLBuffer>) u_old offset:0 atIndex:5];
      [enc setBuffer:(__bridge id<MTLBuffer>) v_old offset:0 atIndex:6];
      [enc setBuffer:(__bridge id<MTLBuffer>) w_old offset:0 atIndex:7];
      [enc setBuffer:(__bridge id<MTLBuffer>) acc_x offset:0 atIndex:8];
      [enc setBuffer:(__bridge id<MTLBuffer>) acc_y offset:0 atIndex:9];
      [enc setBuffer:(__bridge id<MTLBuffer>) acc_z offset:0 atIndex:10];
      [enc setBuffer:(__bridge id<MTLBuffer>) u_lag offset:0 atIndex:11];
      [enc setBuffer:(__bridge id<MTLBuffer>) v_lag offset:0 atIndex:12];
      [enc setBuffer:(__bridge id<MTLBuffer>) w_lag offset:0 atIndex:13];
      [enc setBuffer:(__bridge id<MTLBuffer>) acc_xlag offset:0 atIndex:14];
      [enc setBuffer:(__bridge id<MTLBuffer>) acc_ylag offset:0 atIndex:15];
      [enc setBuffer:(__bridge id<MTLBuffer>) acc_zlag offset:0 atIndex:16];
      [enc setBuffer:(__bridge id<MTLBuffer>) u_laglag offset:0 atIndex:17];
      [enc setBuffer:(__bridge id<MTLBuffer>) v_laglag offset:0 atIndex:18];
      [enc setBuffer:(__bridge id<MTLBuffer>) w_laglag offset:0 atIndex:19];
      [enc setBuffer:(__bridge id<MTLBuffer>) acc_xlaglag offset:0 atIndex:20];
      [enc setBuffer:(__bridge id<MTLBuffer>) acc_ylaglag offset:0 atIndex:21];
      [enc setBuffer:(__bridge id<MTLBuffer>) acc_zlaglag offset:0 atIndex:22];
      [enc setBytes:&n_r length:sizeof(int) atIndex:23];
      [enc setBytes:&lag_len_r length:sizeof(int) atIndex:24];
      [enc setBytes:&e_n_r length:sizeof(float) atIndex:25];
    }, (NSUInteger) n_r);
}
