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
 * Metal host-side dispatch for the gradient jump penalty source term.
 *
 * @note Apple GPUs do not support FP64. This backend operates in FP32.
 */

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include <device/device_config.h>
#include <device/metal/kernel_utils.h>

#define GJP_NTHREADS 256

/** Fortran wrapper for the Metal facet value gather kernel */
void metal_pick_facet_value_hex(void *b, void *a, int *nx, int *nel) {
  if (*nel < 1) return;

  const int nx_r = *nx;

  neko_metal_dispatch_groups(
    neko_metal_pipeline(@"pick_facet_value_hex_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) b offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) a offset:0 atIndex:1];
      [enc setBytes:&nx_r length:sizeof(int) atIndex:2];
    }, (NSUInteger) *nel, GJP_NTHREADS);
}

/** Fortran wrapper for the Metal gradient jump penalty finalize kernel */
void metal_gradient_jump_penalty_finalize(void *penalty_d,
                                          void *penalty_facet_d,
                                          void *dphidxi_d,
                                          int *nx, int *nel) {
  if (*nel < 1) return;

  const int nx_r = *nx;

  neko_metal_dispatch_groups(
    neko_metal_pipeline(@"gradient_jump_penalty_finalize_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) penalty_d offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) penalty_facet_d
              offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) dphidxi_d offset:0 atIndex:2];
      [enc setBytes:&nx_r length:sizeof(int) atIndex:3];
    }, (NSUInteger) *nel, GJP_NTHREADS);
}
