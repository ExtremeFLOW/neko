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
 * Metal host-side dispatch for the Smagorinsky eddy viscosity.
 *
 * @note Apple GPUs do not support FP64. This backend operates in FP32.
 */

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include <device/device_config.h>
#include <device/metal/kernel_utils.h>

/** Fortran wrapper for the Metal Smagorinsky eddy viscosity kernel */
void metal_smagorinsky_nut_compute(void *s11, void *s22, void *s33,
                                   void *s12, void *s13, void *s23,
                                   void *delta, void *nut, void *mult,
                                   real *c_s, int *n) {
  if (*n < 1) return;

  const real c_s_r = *c_s;
  const int n_r = *n;

  neko_metal_dispatch_1d(neko_metal_pipeline(@"smagorinsky_nut_compute_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) s11 offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) s22 offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) s33 offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>) s12 offset:0 atIndex:3];
      [enc setBuffer:(__bridge id<MTLBuffer>) s13 offset:0 atIndex:4];
      [enc setBuffer:(__bridge id<MTLBuffer>) s23 offset:0 atIndex:5];
      [enc setBuffer:(__bridge id<MTLBuffer>) delta offset:0 atIndex:6];
      [enc setBuffer:(__bridge id<MTLBuffer>) nut offset:0 atIndex:7];
      [enc setBuffer:(__bridge id<MTLBuffer>) mult offset:0 atIndex:8];
      [enc setBytes:&c_s_r length:sizeof(real) atIndex:9];
      [enc setBytes:&n_r length:sizeof(int) atIndex:10];
    }, (NSUInteger) n_r);
}
