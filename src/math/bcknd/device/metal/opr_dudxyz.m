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
 * Metal host-side dispatch for derivative computation (dudxyz).
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
static id<MTLComputePipelineState> pso_dudxyz[17] = { nil };

/**
 * Create a compute pipeline state for the named kernel.
 */
static id<MTLComputePipelineState>
get_dudxyz_pipeline(const char *name) {
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
 * Compute the derivative du = (dr * du/dr + ds * du/ds + dt * du/dt) * jacinv
 * on the Metal GPU.
 */
void metal_dudxyz(void *du, void *u,
                  void *dr, void *ds, void *dt,
                  void *dx, void *dy, void *dz,
                  void *jacinv, int *nel, int *lx) {

  if (*lx < 2 || *lx > 16) {
    fprintf(stderr, "%s: size not supported: %d\n", __FILE__, *lx);
    exit(1);
  }

  if (pso_dudxyz[*lx] == nil) {
    char name[64];
    snprintf(name, sizeof(name), "dudxyz_kernel_lx%d", *lx);
    pso_dudxyz[*lx] = get_dudxyz_pipeline(name);
  }

  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_dudxyz[*lx]];

  [enc setBuffer:(__bridge id<MTLBuffer>)du     offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)u      offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)dr     offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)ds     offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)dt     offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)dx     offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)dy     offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)dz     offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)jacinv offset:0 atIndex:8];

  MTLSize groupSize = MTLSizeMake((NSUInteger)(*lx), (NSUInteger)(*lx), 1);
  MTLSize numGroups = MTLSizeMake((NSUInteger)(*nel), 1, 1);

  [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}
