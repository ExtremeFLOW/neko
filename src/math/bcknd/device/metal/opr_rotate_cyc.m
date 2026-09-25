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
 * Metal host-side dispatch for the cyclic boundary velocity rotation.
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

/* Cached pipeline state */
static id<MTLComputePipelineState> pso_rotate_cyc = nil;

/**
 * Create a compute pipeline state for the named kernel.
 */
static id<MTLComputePipelineState>
get_rotate_cyc_pipeline(const char *name) {
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
 * Rotate the velocity field at the cyclic mask points on the Metal GPU.
 */
void metal_rotate_cyc(void *vx, void *vy, void *vz,
                      void *x, void *y, void *z,
                      void *cyc_msk, void *R11, void *R12,
                      int *ncyc, int *idir) {

  if (*ncyc <= 0)
    return;

  if (pso_rotate_cyc == nil)
    pso_rotate_cyc = get_rotate_cyc_pipeline("rotate_cyc_kernel");

  int fncyc = *ncyc, fidir = *idir;

  @autoreleasepool {
    id<MTLCommandQueue> queue = (__bridge id<MTLCommandQueue>)glb_cmd_queue;
    id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

    [enc setComputePipelineState:pso_rotate_cyc];

    [enc setBuffer:(__bridge id<MTLBuffer>)vx      offset:0 atIndex:0];
    [enc setBuffer:(__bridge id<MTLBuffer>)vy      offset:0 atIndex:1];
    [enc setBuffer:(__bridge id<MTLBuffer>)vz      offset:0 atIndex:2];
    [enc setBuffer:(__bridge id<MTLBuffer>)x       offset:0 atIndex:3];
    [enc setBuffer:(__bridge id<MTLBuffer>)y       offset:0 atIndex:4];
    [enc setBuffer:(__bridge id<MTLBuffer>)z       offset:0 atIndex:5];
    [enc setBuffer:(__bridge id<MTLBuffer>)cyc_msk offset:0 atIndex:6];
    [enc setBuffer:(__bridge id<MTLBuffer>)R11     offset:0 atIndex:7];
    [enc setBuffer:(__bridge id<MTLBuffer>)R12     offset:0 atIndex:8];
    [enc setBytes:&fncyc length:sizeof(int) atIndex:9];
    [enc setBytes:&fidir length:sizeof(int) atIndex:10];

    const NSUInteger nthrds = 1024;
    MTLSize threads = MTLSizeMake((NSUInteger)(*ncyc), 1, 1);
    NSUInteger tpg = (NSUInteger)(*ncyc) < nthrds ? (NSUInteger)(*ncyc) : nthrds;
    MTLSize tg = MTLSizeMake(tpg, 1, 1);

    [enc dispatchThreads:threads threadsPerThreadgroup:tg];
    [enc endEncoding];
    [cmdBuf commit];
    [cmdBuf waitUntilCompleted];
  }
}
