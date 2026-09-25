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
 * Metal host-side dispatch for tensor product application (tnsr3d).
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

/* Cached pipeline states, indexed by max(nu, nv) */
static id<MTLComputePipelineState> pso_tnsr3d[17]    = { nil };
static id<MTLComputePipelineState> pso_tnsr3d_el[17] = { nil };

/**
 * Create a compute pipeline state for the named kernel.
 */
static id<MTLComputePipelineState>
get_tensor_pipeline(const char *name) {
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
 * Apply the tensor product v = A x Bt x Ct u on the Metal GPU.
 */
void metal_tnsr3d(void *v, int *nv, void *u, int *nu,
                  void *A, void *Bt, void *Ct, int *nel) {

  const int n = MAX(*nu, *nv);
  if (n < 2 || n > 16) {
    fprintf(stderr, "%s: size not supported: %d\n", __FILE__, n);
    exit(1);
  }

  if (pso_tnsr3d[n] == nil) {
    char name[64];
    snprintf(name, sizeof(name), "tnsr3d_kernel_n%d", n);
    pso_tnsr3d[n] = get_tensor_pipeline(name);
  }

  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_tnsr3d[n]];

  [enc setBuffer:(__bridge id<MTLBuffer>)v  offset:0 atIndex:0];
  [enc setBytes:nv length:sizeof(int) atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)u  offset:0 atIndex:2];
  [enc setBytes:nu length:sizeof(int) atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)A  offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)Bt offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)Ct offset:0 atIndex:6];

  NSUInteger nthrds = 256;
  MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);
  MTLSize numGroups = MTLSizeMake((NSUInteger)(*nel), 1, 1);

  [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}

/**
 * Apply per-point tensor products for a list of elements on the Metal GPU.
 */
void metal_tnsr3d_el_list(void *v, int *nv, void *u, int *nu,
                          void *A, void *Bt, void *Ct,
                          void *elements, int *n_points) {

  const int n = MAX(*nu, *nv);
  if (n < 2 || n > 16) {
    fprintf(stderr, "%s: size not supported: %d\n", __FILE__, n);
    exit(1);
  }

  if (*n_points == 0)
    return;

  if (pso_tnsr3d_el[n] == nil) {
    char name[64];
    snprintf(name, sizeof(name), "tnsr3d_el_kernel_n%d", n);
    pso_tnsr3d_el[n] = get_tensor_pipeline(name);
  }

  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_tnsr3d_el[n]];

  [enc setBuffer:(__bridge id<MTLBuffer>)v        offset:0 atIndex:0];
  [enc setBytes:nv length:sizeof(int) atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)u        offset:0 atIndex:2];
  [enc setBytes:nu length:sizeof(int) atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)A        offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)Bt       offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)Ct       offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)elements offset:0 atIndex:7];

  NSUInteger nthrds = 256;
  MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);
  MTLSize numGroups = MTLSizeMake((NSUInteger)(*n_points), 1, 1);

  [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}
