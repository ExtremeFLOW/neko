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
 * Metal host-side dispatch for the Chebyshev iteration.
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

/* Cached pipeline states */
static id<MTLComputePipelineState> pso_part1 = nil;
static id<MTLComputePipelineState> pso_part2 = nil;

/**
 * Create a compute pipeline state for the named kernel.
 */
static id<MTLComputePipelineState>
get_cheby_pipeline(const char *name) {
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
 * Resolve the command queue for a (possibly null) stream handle.
 */
static id<MTLCommandQueue> cheby_queue(void *stream) {
  if (stream != NULL)
    return (__bridge id<MTLCommandQueue>)stream;
  return (__bridge id<MTLCommandQueue>)glb_cmd_queue;
}

/**
 * Chebyshev iteration, part 1.
 */
void metal_cheby_part1(void *d, void *x, real *inv_tha, int *n, void *stream) {

  if (*n <= 0)
    return;

  if (pso_part1 == nil)
    pso_part1 = get_cheby_pipeline("cheby_part1_kernel");

  id<MTLCommandQueue> queue = cheby_queue(stream);
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_part1];

  [enc setBuffer:(__bridge id<MTLBuffer>)d offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)x offset:0 atIndex:1];
  [enc setBytes:inv_tha length:sizeof(real) atIndex:2];
  [enc setBytes:n length:sizeof(int) atIndex:3];

  NSUInteger nthrds = 256;
  MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);
  MTLSize numGroups = MTLSizeMake(((NSUInteger)(*n) + nthrds - 1) / nthrds, 1, 1);

  [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}

/**
 * Chebyshev iteration, part 2.
 */
void metal_cheby_part2(void *d, void *w, void *x,
                       real *tmp1, real *tmp2, int *n, void *stream) {

  if (*n <= 0)
    return;

  if (pso_part2 == nil)
    pso_part2 = get_cheby_pipeline("cheby_part2_kernel");

  id<MTLCommandQueue> queue = cheby_queue(stream);
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_part2];

  [enc setBuffer:(__bridge id<MTLBuffer>)d offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)w offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)x offset:0 atIndex:2];
  [enc setBytes:tmp1 length:sizeof(real) atIndex:3];
  [enc setBytes:tmp2 length:sizeof(real) atIndex:4];
  [enc setBytes:n length:sizeof(int) atIndex:5];

  NSUInteger nthrds = 256;
  MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);
  MTLSize numGroups = MTLSizeMake(((NSUInteger)(*n) + nthrds - 1) / nthrds, 1, 1);

  [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}
