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
 * Metal host-side dispatch for the stress formulation Ax kernels.
 *
 * Metal buffers are passed as opaque void* (id<MTLBuffer> cast to void*).
 * The Metal command queue is stored in the global glb_cmd_queue.
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

/** Cached pipeline states for each LX (indices 2..16) */
static id<MTLComputePipelineState> pso_stress_full[17] = { nil };
static id<MTLComputePipelineState> pso_stress_part2 = nil;

/**
 * Get (or create) a compute pipeline state for a given kernel name.
 */
static id<MTLComputePipelineState>
get_pipeline(const char *name) {
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
 * Fortran wrapper for device Metal Ax helm full (stress formulation)
 */
void metal_ax_helm_stress_vector(void *au, void *av, void *aw,
                                 void *u, void *v, void *w,
                                 void *dx, void *dy, void *dz,
                                 void *dxt, void *dyt, void *dzt,
                                 void *h1,
                                 void *drdx, void *drdy, void *drdz,
                                 void *dsdx, void *dsdy, void *dsdz,
                                 void *dtdx, void *dtdy, void *dtdz,
                                 void *jacinv, void *w3, int *nelv, int *lx) {

  if (*lx < 2 || *lx > 16) {
    fprintf(stderr, "%s: size not supported: %d\n", __FILE__, *lx);
    exit(1);
  }

  if (pso_stress_full[*lx] == nil) {
    char name[64];
    snprintf(name, sizeof(name), "ax_helm_stress_kernel_full_lx%d", *lx);
    pso_stress_full[*lx] = get_pipeline(name);
  }

  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_stress_full[*lx]];

  /* dxt, dyt, dzt unused — the kernel uses dx/dy/dz transposed
     access patterns directly */
  [enc setBuffer:(__bridge id<MTLBuffer>)au     offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)av     offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)aw     offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)u      offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)v      offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)w      offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)dx     offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)dy     offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)dz     offset:0 atIndex:8];
  [enc setBuffer:(__bridge id<MTLBuffer>)h1     offset:0 atIndex:9];
  [enc setBuffer:(__bridge id<MTLBuffer>)drdx   offset:0 atIndex:10];
  [enc setBuffer:(__bridge id<MTLBuffer>)drdy   offset:0 atIndex:11];
  [enc setBuffer:(__bridge id<MTLBuffer>)drdz   offset:0 atIndex:12];
  [enc setBuffer:(__bridge id<MTLBuffer>)dsdx   offset:0 atIndex:13];
  [enc setBuffer:(__bridge id<MTLBuffer>)dsdy   offset:0 atIndex:14];
  [enc setBuffer:(__bridge id<MTLBuffer>)dsdz   offset:0 atIndex:15];
  [enc setBuffer:(__bridge id<MTLBuffer>)dtdx   offset:0 atIndex:16];
  [enc setBuffer:(__bridge id<MTLBuffer>)dtdy   offset:0 atIndex:17];
  [enc setBuffer:(__bridge id<MTLBuffer>)dtdz   offset:0 atIndex:18];
  [enc setBuffer:(__bridge id<MTLBuffer>)jacinv offset:0 atIndex:19];
  [enc setBuffer:(__bridge id<MTLBuffer>)w3     offset:0 atIndex:20];

  MTLSize threadsPerGroup = MTLSizeMake(*lx, *lx, 1);
  MTLSize numGroups = MTLSizeMake(*nelv, 1, 1);

  [enc dispatchThreadgroups:numGroups
      threadsPerThreadgroup:threadsPerGroup];
  [enc endEncoding];
  [cmdBuf commit];
}

/**
 * Fortran wrapper for device Metal Ax helm full (stress formulation) part2
 */
void metal_ax_helm_stress_vector_part2(void *au, void *av, void *aw,
                                       void *u, void *v, void *w,
                                       void *h2, void *B, int *n) {

  if (pso_stress_part2 == nil) {
    pso_stress_part2 = get_pipeline("ax_helm_stress_kernel_vector_part2");
  }

  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_stress_part2];

  [enc setBuffer:(__bridge id<MTLBuffer>)au offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)av offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)aw offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)u  offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)v  offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)w  offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)h2 offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)B  offset:0 atIndex:7];
  [enc setBytes:n length:sizeof(int) atIndex:8];

  NSUInteger threadGroupSize = MIN(
      (NSUInteger)pso_stress_part2.maxTotalThreadsPerThreadgroup,
      (NSUInteger)1024);
  MTLSize threadsPerGroup = MTLSizeMake(threadGroupSize, 1, 1);
  MTLSize numGroups = MTLSizeMake(
      ((NSUInteger)*n + threadGroupSize - 1) / threadGroupSize, 1, 1);

  [enc dispatchThreadgroups:numGroups
      threadsPerThreadgroup:threadsPerGroup];
  [enc endEncoding];
  [cmdBuf commit];
}
