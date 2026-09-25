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
 * Metal host-side dispatch for the design-field mapping functions
 * (smooth step, step, inverse permeability) used by the Brinkman term.
 *
 * @note Apple GPUs do not support FP64. This backend operates in FP32.
 */

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#include <stdio.h>
#include <stdlib.h>
#include <device/device_config.h>
#include <device/metal/check.h>

/* Defined in device/metal/metal.m */
extern id<MTLDevice> neko_metal_device(void);
extern id<MTLLibrary> neko_metal_library(void);
extern void *glb_cmd_queue;

/**
 * Create a compute pipeline state for the named kernel.
 */
static id<MTLComputePipelineState>
get_mappings_pipeline(const char *name) {
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
 * Dispatch one thread per element and wait for completion.
 */
static void mappings_dispatch(id<MTLComputePipelineState> pso,
                              void (^encode)(id<MTLComputeCommandEncoder>),
                              int n) {
  if (n <= 0)
    return;
  @autoreleasepool {
    id<MTLCommandQueue> queue = (__bridge id<MTLCommandQueue>)glb_cmd_queue;
    id<MTLCommandBuffer> cb = [queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
    [enc setComputePipelineState:pso];
    encode(enc);
    const NSUInteger nthrds = 1024;
    MTLSize threads = MTLSizeMake((NSUInteger)n, 1, 1);
    NSUInteger tpg = (NSUInteger)n < nthrds ? (NSUInteger)n : nthrds;
    MTLSize tg = MTLSizeMake(tpg, 1, 1);
    [enc dispatchThreads:threads threadsPerThreadgroup:tg];
    [enc endEncoding];
    [cb commit];
    [cb waitUntilCompleted];
  }
}

/* Cached pipeline states */
static id<MTLComputePipelineState> pso_smooth_step = nil;
static id<MTLComputePipelineState> pso_step = nil;
static id<MTLComputePipelineState> pso_permeability = nil;

/**
 * Apply a second order smooth step function to a scalar field.
 */
void metal_smooth_step(void *x, real *edge0, real *edge1, int *n) {
  if (pso_smooth_step == nil)
    pso_smooth_step = get_mappings_pipeline("smooth_step_kernel");

  float fe0 = *edge0, fe1 = *edge1;
  int fn = *n;
  mappings_dispatch(pso_smooth_step,
      ^(id<MTLComputeCommandEncoder> enc) {
        [enc setBuffer:(__bridge id<MTLBuffer>)x offset:0 atIndex:0];
        [enc setBytes:&fe0 length:sizeof(float) atIndex:1];
        [enc setBytes:&fe1 length:sizeof(float) atIndex:2];
        [enc setBytes:&fn  length:sizeof(int)   atIndex:3];
      }, *n);
}

/**
 * Apply a step function to a scalar field.
 */
void metal_step_function(void *x, real *edge, real *left, real *right, int *n) {
  if (pso_step == nil)
    pso_step = get_mappings_pipeline("step_kernel");

  float fe = *edge, fl = *left, fr = *right;
  int fn = *n;
  mappings_dispatch(pso_step,
      ^(id<MTLComputeCommandEncoder> enc) {
        [enc setBuffer:(__bridge id<MTLBuffer>)x offset:0 atIndex:0];
        [enc setBytes:&fe length:sizeof(float) atIndex:1];
        [enc setBytes:&fl length:sizeof(float) atIndex:2];
        [enc setBytes:&fr length:sizeof(float) atIndex:3];
        [enc setBytes:&fn length:sizeof(int)   atIndex:4];
      }, *n);
}

/**
 * Apply the inverse permeability mapping to a scalar field.
 */
void metal_permeability(void *x, real *k_0, real *k_1, real *q, int *n) {
  if (pso_permeability == nil)
    pso_permeability = get_mappings_pipeline("permeability_kernel");

  float fk0 = *k_0, fk1 = *k_1, fq = *q;
  int fn = *n;
  mappings_dispatch(pso_permeability,
      ^(id<MTLComputeCommandEncoder> enc) {
        [enc setBuffer:(__bridge id<MTLBuffer>)x offset:0 atIndex:0];
        [enc setBytes:&fk0 length:sizeof(float) atIndex:1];
        [enc setBytes:&fk1 length:sizeof(float) atIndex:2];
        [enc setBytes:&fq  length:sizeof(float) atIndex:3];
        [enc setBytes:&fn  length:sizeof(int)   atIndex:4];
      }, *n);
}
