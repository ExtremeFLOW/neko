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
 * Metal host-side dispatch for the Vreman eddy viscosity.
 *
 * @note Apple GPUs do not support FP64. This backend operates in FP32.
 */

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include <math.h>

#include <device/device_config.h>
#include <device/metal/kernel_utils.h>

/** Fortran wrapper for the Metal Vreman eddy viscosity kernel */
void metal_vreman_nut_compute(void *a11, void *a12, void *a13,
                              void *a21, void *a22, void *a23,
                              void *a31, void *a32, void *a33,
                              void *delta, void *nut, void *mult,
                              real *c, real *eps, int *n) {
  if (*n < 1) return;

  const real c_r = *c;
  const real eps_r = *eps;
  const int n_r = *n;

  neko_metal_dispatch_1d(neko_metal_pipeline(@"vreman_nut_compute_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) a11 offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) a12 offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) a13 offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>) a21 offset:0 atIndex:3];
      [enc setBuffer:(__bridge id<MTLBuffer>) a22 offset:0 atIndex:4];
      [enc setBuffer:(__bridge id<MTLBuffer>) a23 offset:0 atIndex:5];
      [enc setBuffer:(__bridge id<MTLBuffer>) a31 offset:0 atIndex:6];
      [enc setBuffer:(__bridge id<MTLBuffer>) a32 offset:0 atIndex:7];
      [enc setBuffer:(__bridge id<MTLBuffer>) a33 offset:0 atIndex:8];
      [enc setBuffer:(__bridge id<MTLBuffer>) delta offset:0 atIndex:9];
      [enc setBuffer:(__bridge id<MTLBuffer>) nut offset:0 atIndex:10];
      [enc setBuffer:(__bridge id<MTLBuffer>) mult offset:0 atIndex:11];
      [enc setBytes:&c_r length:sizeof(real) atIndex:12];
      [enc setBytes:&eps_r length:sizeof(real) atIndex:13];
      [enc setBytes:&n_r length:sizeof(int) atIndex:14];
    }, (NSUInteger) n_r);
}

/**
 * Fortran wrapper for the Metal buoyancy-corrected Vreman kernel
 * @note @a g is the (unnormalised) gravity vector; the shear direction
 * uses the normalised vector, the buoyancy term the original one.
 */
void metal_vreman_nut_compute_buoy(void *a11, void *a12, void *a13,
                                   void *a21, void *a22, void *a23,
                                   void *a31, void *a32, void *a33,
                                   void *delta, void *nut, void *mult,
                                   real *c, real *eps, int *n,
                                   void *dTdx, void *dTdy, void *dTdz,
                                   real *g, real *ri_c, real *ref_temp) {
  if (*n < 1) return;

  const real c_r = *c;
  const real eps_r = *eps;
  const int n_r = *n;
  const real g1 = g[0];
  const real g2 = g[1];
  const real g3 = g[2];
  const real gmag = sqrt(g1 * g1 + g2 * g2 + g3 * g3);
  const real n1 = g1 / gmag;
  const real n2 = g2 / gmag;
  const real n3 = g3 / gmag;
  const real ri_c_r = *ri_c;
  const real ref_temp_r = *ref_temp;

  neko_metal_dispatch_1d(neko_metal_pipeline(@"vreman_nut_compute_buoy_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) a11 offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) a12 offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) a13 offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>) a21 offset:0 atIndex:3];
      [enc setBuffer:(__bridge id<MTLBuffer>) a22 offset:0 atIndex:4];
      [enc setBuffer:(__bridge id<MTLBuffer>) a23 offset:0 atIndex:5];
      [enc setBuffer:(__bridge id<MTLBuffer>) a31 offset:0 atIndex:6];
      [enc setBuffer:(__bridge id<MTLBuffer>) a32 offset:0 atIndex:7];
      [enc setBuffer:(__bridge id<MTLBuffer>) a33 offset:0 atIndex:8];
      [enc setBuffer:(__bridge id<MTLBuffer>) delta offset:0 atIndex:9];
      [enc setBuffer:(__bridge id<MTLBuffer>) nut offset:0 atIndex:10];
      [enc setBuffer:(__bridge id<MTLBuffer>) mult offset:0 atIndex:11];
      [enc setBytes:&c_r length:sizeof(real) atIndex:12];
      [enc setBytes:&eps_r length:sizeof(real) atIndex:13];
      [enc setBytes:&n_r length:sizeof(int) atIndex:14];
      [enc setBuffer:(__bridge id<MTLBuffer>) dTdx offset:0 atIndex:15];
      [enc setBuffer:(__bridge id<MTLBuffer>) dTdy offset:0 atIndex:16];
      [enc setBuffer:(__bridge id<MTLBuffer>) dTdz offset:0 atIndex:17];
      [enc setBytes:&n1 length:sizeof(real) atIndex:18];
      [enc setBytes:&n2 length:sizeof(real) atIndex:19];
      [enc setBytes:&n3 length:sizeof(real) atIndex:20];
      [enc setBytes:&g1 length:sizeof(real) atIndex:21];
      [enc setBytes:&g2 length:sizeof(real) atIndex:22];
      [enc setBytes:&g3 length:sizeof(real) atIndex:23];
      [enc setBytes:&ri_c_r length:sizeof(real) atIndex:24];
      [enc setBytes:&ref_temp_r length:sizeof(real) atIndex:25];
    }, (NSUInteger) n_r);
}
