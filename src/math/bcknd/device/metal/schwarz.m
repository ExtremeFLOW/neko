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
 * Metal host-side dispatch for the additive Schwarz method.
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
static id<MTLComputePipelineState> pso_extrude[17] = { nil };
static id<MTLComputePipelineState> pso_toext3d = nil;
static id<MTLComputePipelineState> pso_toreg3d = nil;

/**
 * Create a compute pipeline state for the named kernel.
 */
static id<MTLComputePipelineState>
get_schwarz_pipeline(const char *name) {
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
static id<MTLCommandQueue> schwarz_queue(void *stream) {
  if (stream != NULL)
    return (__bridge id<MTLCommandQueue>)stream;
  return (__bridge id<MTLCommandQueue>)glb_cmd_queue;
}

/**
 * Sum the element shells (schwarz_extrude) on the Metal GPU.
 */
void metal_schwarz_extrude(void *arr1, int *l1, real *f1,
                           void *arr2, int *l2, real *f2,
                           int *nx, int *nelv, void *stream) {

  if (*nx < 2 || *nx > 16) {
    fprintf(stderr, "%s: size not supported: %d\n", __FILE__, *nx);
    exit(1);
  }
  if (*nelv <= 0)
    return;

  if (pso_extrude[*nx] == nil) {
    char name[64];
    snprintf(name, sizeof(name), "schwarz_extrude_kernel_nx%d", *nx);
    pso_extrude[*nx] = get_schwarz_pipeline(name);
  }

  int fl1 = *l1, fl2 = *l2;
  float ff1 = *f1, ff2 = *f2;

  id<MTLCommandQueue> queue = schwarz_queue(stream);
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_extrude[*nx]];

  [enc setBuffer:(__bridge id<MTLBuffer>)arr1 offset:0 atIndex:0];
  [enc setBytes:&fl1 length:sizeof(int) atIndex:1];
  [enc setBytes:&ff1 length:sizeof(float) atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)arr2 offset:0 atIndex:3];
  [enc setBytes:&fl2 length:sizeof(int) atIndex:4];
  [enc setBytes:&ff2 length:sizeof(float) atIndex:5];

  NSUInteger nthrds = (NSUInteger)((*nx - 2) * (*nx - 2));
  MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);
  MTLSize numGroups = MTLSizeMake((NSUInteger)(*nelv), 1, 1);

  [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}

/**
 * Copy the regular array into the extended array (schwarz_toext3d).
 */
void metal_schwarz_toext3d(void *a, void *b, int *nx, int *nelv, void *stream) {

  if (*nelv <= 0)
    return;

  if (pso_toext3d == nil)
    pso_toext3d = get_schwarz_pipeline("schwarz_toext3d_kernel");

  int fnx = *nx;

  id<MTLCommandQueue> queue = schwarz_queue(stream);
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_toext3d];

  [enc setBuffer:(__bridge id<MTLBuffer>)a offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)b offset:0 atIndex:1];
  [enc setBytes:&fnx length:sizeof(int) atIndex:2];

  MTLSize groupSize = MTLSizeMake(256, 1, 1);
  MTLSize numGroups = MTLSizeMake((NSUInteger)(*nelv), 1, 1);

  [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}

/**
 * Copy the extended array back to the regular array (schwarz_toreg3d).
 */
void metal_schwarz_toreg3d(void *b, void *a, int *nx, int *nelv, void *stream) {

  if (*nelv <= 0)
    return;

  if (pso_toreg3d == nil)
    pso_toreg3d = get_schwarz_pipeline("schwarz_toreg3d_kernel");

  int fnx = *nx;

  id<MTLCommandQueue> queue = schwarz_queue(stream);
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_toreg3d];

  [enc setBuffer:(__bridge id<MTLBuffer>)b offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)a offset:0 atIndex:1];
  [enc setBytes:&fnx length:sizeof(int) atIndex:2];

  MTLSize groupSize = MTLSizeMake(256, 1, 1);
  MTLSize numGroups = MTLSizeMake((NSUInteger)(*nelv), 1, 1);

  [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}
