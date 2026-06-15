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
 * Metal host-side dispatch for Ax Helm kernels.
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
#include <string.h>
#include <device/device_config.h>
#include <device/metal/check.h>
#include <common/neko_log.h>

/* Defined in device/metal/metal.m */
extern id<MTLDevice> neko_metal_device(void);
extern id<MTLLibrary> neko_metal_library(void);

/** Cached pipeline states for each LX (indices 2..16) */
static id<MTLComputePipelineState> pso_kstep[17] = { nil };
static id<MTLComputePipelineState> pso_vector_kstep[17] = { nil };
static id<MTLComputePipelineState> pso_vector_part2 = nil;
/**
 * Return the shared Metal shader library (embedded in the binary).
 */
static id<MTLLibrary> get_metal_library(void) {
  return neko_metal_library();
}

/**
 * Get (or create) a compute pipeline state for a given kernel name.
 */
static id<MTLComputePipelineState>
get_pipeline(const char *name) {
  id<MTLDevice> device = neko_metal_device();
  id<MTLLibrary> lib = get_metal_library();
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
 * Dispatch a scalar ax_helm kernel for a given LX.
 */
static void dispatch_ax_helm(int lx, int nelv,
                             void *w_buf, void *u_buf,
                             void *dx_buf, void *dy_buf, void *dz_buf,
                             void *h1_buf,
                             void *g11_buf, void *g22_buf, void *g33_buf,
                             void *g12_buf, void *g13_buf, void *g23_buf) {

  if (pso_kstep[lx] == nil) {
    char name[64];
    if (lx == 4 || lx == 8 || lx == 16) {
      snprintf(name, sizeof(name),
               "ax_helm_kernel_kstep_padded_lx%d", lx);
    } else {
      snprintf(name, sizeof(name),
               "ax_helm_kernel_kstep_lx%d", lx);
    }
    pso_kstep[lx] = get_pipeline(name);
  }

  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_kstep[lx]];

  [enc setBuffer:(__bridge id<MTLBuffer>)w_buf   offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)u_buf   offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)dx_buf  offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)dy_buf  offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)dz_buf  offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)h1_buf  offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)g11_buf offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)g22_buf offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)g33_buf offset:0 atIndex:8];
  [enc setBuffer:(__bridge id<MTLBuffer>)g12_buf offset:0 atIndex:9];
  [enc setBuffer:(__bridge id<MTLBuffer>)g13_buf offset:0 atIndex:10];
  [enc setBuffer:(__bridge id<MTLBuffer>)g23_buf offset:0 atIndex:11];

  MTLSize threadsPerGroup = MTLSizeMake(lx, lx, 1);
  MTLSize numGroups = MTLSizeMake(nelv, 1, 1);

  [enc dispatchThreadgroups:numGroups
      threadsPerThreadgroup:threadsPerGroup];
  [enc endEncoding];
  [cmdBuf commit];
}

/**
 * Dispatch a vector ax_helm kernel for a given LX.
 */
static void dispatch_ax_helm_vector(int lx, int nelv,
                                    void *au_buf, void *av_buf, void *aw_buf,
                                    void *u_buf, void *v_buf, void *w_buf,
                                    void *dx_buf, void *dy_buf, void *dz_buf,
                                    void *h1_buf,
                                    void *g11_buf, void *g22_buf,
                                    void *g33_buf, void *g12_buf,
                                    void *g13_buf, void *g23_buf) {

  if (pso_vector_kstep[lx] == nil) {
    char name[64];
    if (lx == 4 || lx == 8 || lx == 16) {
      snprintf(name, sizeof(name),
               "ax_helm_kernel_vector_kstep_padded_lx%d", lx);
    } else {
      snprintf(name, sizeof(name),
               "ax_helm_kernel_vector_kstep_lx%d", lx);
    }
    pso_vector_kstep[lx] = get_pipeline(name);
  }

  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_vector_kstep[lx]];

  [enc setBuffer:(__bridge id<MTLBuffer>)au_buf  offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)av_buf  offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)aw_buf  offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)u_buf   offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)v_buf   offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)w_buf   offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)dx_buf  offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)dy_buf  offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)dz_buf  offset:0 atIndex:8];
  [enc setBuffer:(__bridge id<MTLBuffer>)h1_buf  offset:0 atIndex:9];
  [enc setBuffer:(__bridge id<MTLBuffer>)g11_buf offset:0 atIndex:10];
  [enc setBuffer:(__bridge id<MTLBuffer>)g22_buf offset:0 atIndex:11];
  [enc setBuffer:(__bridge id<MTLBuffer>)g33_buf offset:0 atIndex:12];
  [enc setBuffer:(__bridge id<MTLBuffer>)g12_buf offset:0 atIndex:13];
  [enc setBuffer:(__bridge id<MTLBuffer>)g13_buf offset:0 atIndex:14];
  [enc setBuffer:(__bridge id<MTLBuffer>)g23_buf offset:0 atIndex:15];

  MTLSize threadsPerGroup = MTLSizeMake(lx, lx, 1);
  MTLSize numGroups = MTLSizeMake(nelv, 1, 1);

  [enc dispatchThreadgroups:numGroups
      threadsPerThreadgroup:threadsPerGroup];
  [enc endEncoding];
  [cmdBuf commit];
}

  /**
   * Fortran wrapper for device Metal Ax helm
   */
  void metal_ax_helm(void *w, void *u, void *dx, void *dy, void *dz,
                     void *dxt, void *dyt, void *dzt, void *h1,
                     void *g11, void *g22, void *g33, void *g12,
                     void *g13, void *g23, int *nelv, int *lx) {

    if (*lx < 2 || *lx > 16) {
      fprintf(stderr, "%s: size not supported: %d\n", __FILE__, *lx);
      exit(1);
    }

    /*
     * dxt, dyt, dzt unused — the kstep kernels use dx/dy/dz transposed
     * access patterns directly
     */
    dispatch_ax_helm(*lx, *nelv, w, u, dx, dy, dz, h1,
                     g11, g22, g33, g12, g13, g23);
  }

  /**
   * Fortran wrapper for device Metal Ax helm vector version
   */
  void metal_ax_helm_vector(void *au, void *av, void *aw,
                            void *u, void *v, void *w,
                            void *dx, void *dy, void *dz,
                            void *dxt, void *dyt, void *dzt,
                            void *h1, void *g11, void *g22,
                            void *g33, void *g12, void *g13,
                            void *g23, int *nelv, int *lx) {

    if (*lx < 2 || *lx > 16) {
      fprintf(stderr, "%s: size not supported: %d\n", __FILE__, *lx);
      exit(1);
    }

    dispatch_ax_helm_vector(*lx, *nelv, au, av, aw, u, v, w,
                            dx, dy, dz, h1,
                            g11, g22, g33, g12, g13, g23);
  }

  /**
   * Fortran wrapper for device Metal Ax helm vector part2
   */
  void metal_ax_helm_vector_part2(void *au, void *av, void *aw,
                                  void *u, void *v, void *w,
                                  void *h2, void *B, int *n) {

    if (pso_vector_part2 == nil) {
      pso_vector_part2 = get_pipeline("ax_helm_kernel_vector_part2");
    }

    id<MTLCommandQueue> queue =
      (__bridge id<MTLCommandQueue>)glb_cmd_queue;
    id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

    [enc setComputePipelineState:pso_vector_part2];

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
        (NSUInteger)pso_vector_part2.maxTotalThreadsPerThreadgroup,
        (NSUInteger)1024);
    MTLSize threadsPerGroup = MTLSizeMake(threadGroupSize, 1, 1);
    MTLSize numGroups = MTLSizeMake(
        ((NSUInteger)*n + threadGroupSize - 1) / threadGroupSize, 1, 1);

    [enc dispatchThreadgroups:numGroups
        threadsPerThreadgroup:threadsPerGroup];
    [enc endEncoding];
    [cmdBuf commit];
  }
