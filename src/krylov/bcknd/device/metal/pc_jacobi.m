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
 * Metal host-side dispatch for the Jacobi preconditioner update.
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

/* Cached pipeline states, indexed by LX */
static id<MTLComputePipelineState> pso_jacobi[17] = { nil };

/**
 * Create a compute pipeline state for the named kernel.
 */
static id<MTLComputePipelineState>
get_jacobi_pipeline(const char *name) {
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
 * Update the diagonal of the stiffness matrix on the Metal GPU.
 */
void metal_jacobi_update(void *d_d, void *dxt_d, void *dyt_d, void *dzt_d,
                         void *G11_d, void *G22_d, void *G33_d,
                         void *G12_d, void *G13_d, void *G23_d,
                         int *nelv, int *lx) {

  int LX = *lx;
  if (LX < 2 || LX > 16) {
    fprintf(stderr, "%s: size not supported: %d\n", __FILE__, LX);
    exit(1);
  }

  if (pso_jacobi[LX] == nil) {
    char name[64];
    snprintf(name, sizeof(name), "jacobi_kernel_lx%d", LX);
    pso_jacobi[LX] = get_jacobi_pipeline(name);
  }

  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_jacobi[LX]];

  [enc setBuffer:(__bridge id<MTLBuffer>)d_d   offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)dxt_d offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)dyt_d offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)dzt_d offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)G11_d offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)G22_d offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)G33_d offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)G12_d offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)G13_d offset:0 atIndex:8];
  [enc setBuffer:(__bridge id<MTLBuffer>)G23_d offset:0 atIndex:9];
  [enc setBytes:nelv length:sizeof(int) atIndex:10];

  const NSUInteger n =
    (NSUInteger)(*nelv) * (NSUInteger)(LX*LX*LX);
  NSUInteger nthrds = 1024;
  MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);
  MTLSize numGroups = MTLSizeMake((n + nthrds - 1) / nthrds, 1, 1);

  [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}
