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
 * Metal host-side dispatch for the Deardorff eddy viscosity.
 *
 * @note Apple GPUs do not support FP64. This backend operates in FP32.
 */

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include <device/device_config.h>
#include <device/metal/kernel_utils.h>

/** Fortran wrapper for the Metal Deardorff eddy viscosity kernel */
void metal_deardorff_nut_compute(void *TKE,
                                 void *dTdx, void *dTdy, void *dTdz,
                                 void *a11, void *a12, void *a13,
                                 void *a21, void *a22, void *a23,
                                 void *a31, void *a32, void *a33,
                                 void *delta, void *nut,
                                 void *temperature_alphat,
                                 void *TKE_alphat, void *TKE_source,
                                 real *c_k, real *T0,
                                 real *g1, real *g2, real *g3, real *eps,
                                 int *n) {
  if (*n < 1) return;

  const real c_k_r = *c_k;
  const real T0_r = *T0;
  const real g1_r = *g1;
  const real g2_r = *g2;
  const real g3_r = *g3;
  const real eps_r = *eps;
  const int n_r = *n;

  neko_metal_dispatch_1d(neko_metal_pipeline(@"deardorff_nut_compute_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) TKE offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) dTdx offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) dTdy offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>) dTdz offset:0 atIndex:3];
      [enc setBuffer:(__bridge id<MTLBuffer>) a11 offset:0 atIndex:4];
      [enc setBuffer:(__bridge id<MTLBuffer>) a12 offset:0 atIndex:5];
      [enc setBuffer:(__bridge id<MTLBuffer>) a13 offset:0 atIndex:6];
      [enc setBuffer:(__bridge id<MTLBuffer>) a21 offset:0 atIndex:7];
      [enc setBuffer:(__bridge id<MTLBuffer>) a22 offset:0 atIndex:8];
      [enc setBuffer:(__bridge id<MTLBuffer>) a23 offset:0 atIndex:9];
      [enc setBuffer:(__bridge id<MTLBuffer>) a31 offset:0 atIndex:10];
      [enc setBuffer:(__bridge id<MTLBuffer>) a32 offset:0 atIndex:11];
      [enc setBuffer:(__bridge id<MTLBuffer>) a33 offset:0 atIndex:12];
      [enc setBuffer:(__bridge id<MTLBuffer>) delta offset:0 atIndex:13];
      [enc setBuffer:(__bridge id<MTLBuffer>) nut offset:0 atIndex:14];
      [enc setBuffer:(__bridge id<MTLBuffer>) temperature_alphat
              offset:0 atIndex:15];
      [enc setBuffer:(__bridge id<MTLBuffer>) TKE_alphat offset:0 atIndex:16];
      [enc setBuffer:(__bridge id<MTLBuffer>) TKE_source offset:0 atIndex:17];
      [enc setBytes:&c_k_r length:sizeof(real) atIndex:18];
      [enc setBytes:&T0_r length:sizeof(real) atIndex:19];
      [enc setBytes:&g1_r length:sizeof(real) atIndex:20];
      [enc setBytes:&g2_r length:sizeof(real) atIndex:21];
      [enc setBytes:&g3_r length:sizeof(real) atIndex:22];
      [enc setBytes:&eps_r length:sizeof(real) atIndex:23];
      [enc setBytes:&n_r length:sizeof(int) atIndex:24];
    }, (NSUInteger) n_r);
}
