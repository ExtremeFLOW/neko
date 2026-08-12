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
 * Metal host-side dispatch for CFL computation.
 *
 * @note Apple GPUs do not support FP64. This backend operates in FP32.
 */

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <device/device_config.h>
#include <device/metal/check.h>

/* Defined in device/metal/metal.m */
extern id<MTLDevice> neko_metal_device(void);
extern id<MTLLibrary> neko_metal_library(void);

/* Cached pipeline states, indexed by LX */
static id<MTLComputePipelineState> pso_cfl[13] = { nil };

/**
 * Create a compute pipeline state for the named kernel.
 */
static id<MTLComputePipelineState>
get_cfl_pipeline(const char *name) {
  id<MTLDevice> device = neko_metal_device();
  id<MTLLibrary> lib = neko_metal_library();
  NSString *nsName = [NSString stringWithUTF8String:name];
  id<MTLFunction> func = [lib newFunctionWithName:nsName];
  if (func == nil) {
    fprintf(stderr, "Metal: kernel '%s' not found in metallib\n", name);
    exit(EXIT_FAILURE);
  }
  NSError *error = nil;
  id<MTLComputePipelineState> pso =
    [device newComputePipelineStateWithFunction:func error:&error];
  METAL_CHECK(error);
  return pso;
}

/**
 * Compute the maximum CFL number on the Metal GPU.
 */
real metal_cfl(real *dt, void *u, void *v, void *w,
               void *drdx, void *dsdx, void *dtdx,
               void *drdy, void *dsdy, void *dtdy,
               void *drdz, void *dsdz, void *dtdz,
               void *dr_inv, void *ds_inv, void *dt_inv,
               void *jacinv, int *nel, int *lx) {

  if (*lx < 2 || *lx > 12) {
    fprintf(stderr, "%s: size not supported: %d\n", __FILE__, *lx);
    exit(1);
  }

  if (pso_cfl[*lx] == nil) {
    char name[64];
    snprintf(name, sizeof(name), "cfl_kernel_lx%d", *lx);
    pso_cfl[*lx] = get_cfl_pipeline(name);
  }

  id<MTLDevice> device = neko_metal_device();
  id<MTLBuffer> cfl_d =
    [device newBufferWithLength:((NSUInteger)(*nel) * sizeof(real))
                        options:MTLResourceStorageModeShared];

  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_cfl[*lx]];

  [enc setBytes:dt length:sizeof(real) atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)u      offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)v      offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)w      offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)drdx   offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)dsdx   offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)dtdx   offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)drdy   offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)dsdy   offset:0 atIndex:8];
  [enc setBuffer:(__bridge id<MTLBuffer>)dtdy   offset:0 atIndex:9];
  [enc setBuffer:(__bridge id<MTLBuffer>)drdz   offset:0 atIndex:10];
  [enc setBuffer:(__bridge id<MTLBuffer>)dsdz   offset:0 atIndex:11];
  [enc setBuffer:(__bridge id<MTLBuffer>)dtdz   offset:0 atIndex:12];
  [enc setBuffer:(__bridge id<MTLBuffer>)dr_inv offset:0 atIndex:13];
  [enc setBuffer:(__bridge id<MTLBuffer>)ds_inv offset:0 atIndex:14];
  [enc setBuffer:(__bridge id<MTLBuffer>)dt_inv offset:0 atIndex:15];
  [enc setBuffer:(__bridge id<MTLBuffer>)jacinv offset:0 atIndex:16];
  [enc setBuffer:cfl_d offset:0 atIndex:17];

  NSUInteger nthrds = 256;
  MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);
  MTLSize numGroups = MTLSizeMake((NSUInteger)(*nel), 1, 1);

  [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];

  const real *cfl = (const real *)[cfl_d contents];
  real cfl_max = 0.0f;
  for (int i = 0; i < (*nel); i++) {
    cfl_max = fmaxf(cfl_max, cfl[i]);
  }

  return cfl_max;
}
