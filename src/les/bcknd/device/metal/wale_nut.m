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
 * Metal host-side dispatch for the wale eddy viscosity.
 *
 * @note Apple GPUs do not support FP64. This backend operates in FP32.
 */

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include <device/device_config.h>
#include <device/metal/kernel_utils.h>

/** Fortran wrapper for the Metal wale eddy viscosity kernel */
void metal_wale_nut_compute(void *g11, void *g12, void *g13,
                            void *g21, void *g22, void *g23,
                            void *g31, void *g32, void *g33,
                            void *delta, void *nut, void *mult,
                            real *c, real *eps, int *n) {
  if (*n < 1) return;

  const real c_r = *c;
  const real eps_r = *eps;
  const int n_r = *n;

  neko_metal_dispatch_1d(neko_metal_pipeline(@"wale_nut_compute_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) g11 offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) g12 offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) g13 offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>) g21 offset:0 atIndex:3];
      [enc setBuffer:(__bridge id<MTLBuffer>) g22 offset:0 atIndex:4];
      [enc setBuffer:(__bridge id<MTLBuffer>) g23 offset:0 atIndex:5];
      [enc setBuffer:(__bridge id<MTLBuffer>) g31 offset:0 atIndex:6];
      [enc setBuffer:(__bridge id<MTLBuffer>) g32 offset:0 atIndex:7];
      [enc setBuffer:(__bridge id<MTLBuffer>) g33 offset:0 atIndex:8];
      [enc setBuffer:(__bridge id<MTLBuffer>) delta offset:0 atIndex:9];
      [enc setBuffer:(__bridge id<MTLBuffer>) nut offset:0 atIndex:10];
      [enc setBuffer:(__bridge id<MTLBuffer>) mult offset:0 atIndex:11];
      [enc setBytes:&c_r length:sizeof(real) atIndex:12];
      [enc setBytes:&eps_r length:sizeof(real) atIndex:13];
      [enc setBytes:&n_r length:sizeof(int) atIndex:14];
    }, (NSUInteger) n_r);
}
