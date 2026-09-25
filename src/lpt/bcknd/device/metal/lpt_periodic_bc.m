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
 * Metal host-side dispatch for LPT periodic boundary wrapping.
 *
 * @note Apple GPUs do not support FP64. This backend operates in FP32.
 */

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include <device/device_config.h>
#include <device/metal/kernel_utils.h>
#include <lpt/bcknd/device/lpt_periodic_params.h>

/** Fortran wrapper for the Metal translational periodic wrapping kernel */
void metal_lpt_periodic_bc_wrap_translational(
    void *x, void *y, void *z, int *n, int *n_periodic_dirs,
    real *periodic_dir_x1, real *periodic_dir_y1, real *periodic_dir_z1,
    real *periodic_dir_x2, real *periodic_dir_y2, real *periodic_dir_z2,
    real *periodic_dir_x3, real *periodic_dir_y3, real *periodic_dir_z3,
    real *periodic_min1, real *periodic_min2, real *periodic_min3,
    real *periodic_max1, real *periodic_max2, real *periodic_max3,
    real *periodic_shift_x1, real *periodic_shift_y1,
    real *periodic_shift_z1, real *periodic_shift_x2,
    real *periodic_shift_y2, real *periodic_shift_z2,
    real *periodic_shift_x3, real *periodic_shift_y3,
    real *periodic_shift_z3, real *periodic_len1, real *periodic_len2,
    real *periodic_len3, void *strm) {
  if (*n < 1) return;

  const int n_r = *n;
  const int n_dirs = *n_periodic_dirs;

  lpt_periodic_params prm;
  lpt_periodic_params_pack(&prm,
      periodic_dir_x1, periodic_dir_y1, periodic_dir_z1,
      periodic_dir_x2, periodic_dir_y2, periodic_dir_z2,
      periodic_dir_x3, periodic_dir_y3, periodic_dir_z3,
      periodic_min1, periodic_min2, periodic_min3,
      periodic_max1, periodic_max2, periodic_max3,
      periodic_shift_x1, periodic_shift_y1, periodic_shift_z1,
      periodic_shift_x2, periodic_shift_y2, periodic_shift_z2,
      periodic_shift_x3, periodic_shift_y3, periodic_shift_z3,
      periodic_len1, periodic_len2, periodic_len3);

  neko_metal_dispatch_1d_queue(
    (__bridge id<MTLCommandQueue>) strm,
    neko_metal_pipeline(@"lpt_periodic_bc_wrap_translational_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) x offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) y offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) z offset:0 atIndex:2];
      [enc setBytes:&n_r length:sizeof(int) atIndex:3];
      [enc setBytes:&n_dirs length:sizeof(int) atIndex:4];
      [enc setBytes:&prm length:sizeof(lpt_periodic_params) atIndex:5];
    }, (NSUInteger) n_r);
}

/**
 * Fortran wrapper for the Metal rotational periodic wrapping kernel
 * @note Absent optional arrays arrive as NULL; they are replaced by a
 * bound placeholder and masked out so the kernel never touches them.
 */
void metal_lpt_periodic_bc_wrap_rotational(
    void *x, void *y, void *z, int *n,
    real *theta_min, real *theta_max, real *theta_len,
    void *u, void *v, void *w,
    void *u_lag, void *v_lag, void *w_lag,
    void *u_laglag, void *v_laglag, void *w_laglag,
    void *acc_xlag, void *acc_ylag, void *acc_zlag,
    void *acc_xlaglag, void *acc_ylaglag, void *acc_zlaglag,
    void *strm) {
  if (*n < 1) return;

  const int n_r = *n;
  const real theta_min_r = *theta_min;
  const real theta_max_r = *theta_max;
  const real theta_len_r = *theta_len;

  int present_mask = 0;
  if (u != NULL && v != NULL) present_mask |= 1;
  if (u_lag != NULL && v_lag != NULL) present_mask |= 2;
  if (u_laglag != NULL && v_laglag != NULL) present_mask |= 4;
  if (acc_xlag != NULL && acc_ylag != NULL) present_mask |= 8;
  if (acc_xlaglag != NULL && acc_ylaglag != NULL) present_mask |= 16;

  /* Unbound slots still need a valid buffer */
  void *pairs[10] = { u, v, u_lag, v_lag, u_laglag, v_laglag,
                      acc_xlag, acc_ylag, acc_xlaglag, acc_ylaglag };
  for (int i = 0; i < 10; i++) {
    if (pairs[i] == NULL) pairs[i] = x;
  }
  /* blocks cannot capture a C array, capture a pointer to it instead */
  void * const *pair = pairs;

  neko_metal_dispatch_1d_queue(
    (__bridge id<MTLCommandQueue>) strm,
    neko_metal_pipeline(@"lpt_periodic_bc_wrap_rotational_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) x offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) y offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) z offset:0 atIndex:2];
      [enc setBytes:&n_r length:sizeof(int) atIndex:3];
      [enc setBytes:&theta_min_r length:sizeof(real) atIndex:4];
      [enc setBytes:&theta_max_r length:sizeof(real) atIndex:5];
      [enc setBytes:&theta_len_r length:sizeof(real) atIndex:6];
      for (int i = 0; i < 10; i++) {
        [enc setBuffer:(__bridge id<MTLBuffer>) pair[i]
                offset:0 atIndex:(NSUInteger) (7 + i)];
      }
      [enc setBytes:&present_mask length:sizeof(int) atIndex:17];
    }, (NSUInteger) n_r);
}
