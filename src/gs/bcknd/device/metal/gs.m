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
 * Metal host-side dispatch for gather-scatter kernels.
 *
 * Metal buffers are passed as opaque void* (id<MTLBuffer> cast to void*).
 * The command queue is stored in the global glb_cmd_queue.
 *
 * @note Apple GPUs do not support FP64. This backend operates in FP32.
 */

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#include <stdio.h>
#include <stdlib.h>
#include <device/device_config.h>
#include <device/metal/check.h>

#define GS_OP_ADD  1
#define GS_OP_MUL  2
#define GS_OP_MIN  3
#define GS_OP_MAX  4

/* Defined in device/metal/metal.m */
extern id<MTLDevice> neko_metal_device(void);
extern id<MTLLibrary> neko_metal_library(void);

/* Cached pipeline states */
static id<MTLComputePipelineState> pso_gather_add    = nil;
static id<MTLComputePipelineState> pso_gather_mul    = nil;
static id<MTLComputePipelineState> pso_gather_min    = nil;
static id<MTLComputePipelineState> pso_gather_max    = nil;
static id<MTLComputePipelineState> pso_scatter       = nil;
static id<MTLComputePipelineState> pso_pack          = nil;
static id<MTLComputePipelineState> pso_unpack_add    = nil;
static id<MTLComputePipelineState> pso_unpack_min    = nil;
static id<MTLComputePipelineState> pso_unpack_max    = nil;
/**
 * Return the shared Metal shader library (embedded in the binary).
 */
static id<MTLLibrary> get_gs_metal_library(void) {
  return neko_metal_library();
}

/**
 * Get (or create) a compute pipeline state for a given kernel name.
 */
static id<MTLComputePipelineState>
get_gs_pipeline(const char *name) {
  id<MTLDevice> device = neko_metal_device();
  id<MTLLibrary> lib = get_gs_metal_library();
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
 * Lazily initialize all gather/scatter pipeline states.
 */
static void ensure_pipelines(void) {
  if (pso_gather_add == nil) {
    pso_gather_add = get_gs_pipeline("gather_kernel_add");
    pso_gather_mul = get_gs_pipeline("gather_kernel_mul");
    pso_gather_min = get_gs_pipeline("gather_kernel_min");
    pso_gather_max = get_gs_pipeline("gather_kernel_max");
    pso_scatter    = get_gs_pipeline("scatter_kernel");
    pso_pack       = get_gs_pipeline("gs_pack_kernel");
    pso_unpack_add = get_gs_pipeline("gs_unpack_add_kernel");
    pso_unpack_min = get_gs_pipeline("gs_unpack_min_kernel");
    pso_unpack_max = get_gs_pipeline("gs_unpack_max_kernel");
  }
}

/*
 * Helper: encode gather kernel arguments
 */

static void encode_gather_args(id<MTLComputeCommandEncoder> enc,
                               void *v, int m, int o,
                               void *dg, void *u, int n,
                               void *gd, int nb,
                               void *b, void *bo) {
  [enc setBuffer:(__bridge id<MTLBuffer>)v  offset:0 atIndex:0];
  [enc setBytes:&m length:sizeof(int) atIndex:1];
  [enc setBytes:&o length:sizeof(int) atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)dg offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)u  offset:0 atIndex:4];
  [enc setBytes:&n length:sizeof(int) atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)gd offset:0 atIndex:6];
  [enc setBytes:&nb length:sizeof(int) atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)b  offset:0 atIndex:8];
  [enc setBuffer:(__bridge id<MTLBuffer>)bo offset:0 atIndex:9];
}

/*
 * Fortran-callable wrappers
 */

/**
 * Fortran wrapper for device gather kernels
 */
void metal_gather_kernel(void *v, int *m, int *o, void *dg,
                         void *u, int *n, void *gd, int *nb,
                         void *b, void *bo, int *op,
                         void *stream) {

  if ((*m) == 0) return;

  ensure_pipelines();

  id<MTLComputePipelineState> pso;
  switch (*op) {
  case GS_OP_ADD: pso = pso_gather_add; break;
  case GS_OP_MUL: pso = pso_gather_mul; break;
  case GS_OP_MIN: pso = pso_gather_min; break;
  case GS_OP_MAX: pso = pso_gather_max; break;
  default:
    fprintf(stderr, "%s: unknown gs op %d\n", __FILE__, *op);
    abort();
  }

  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)(stream ? stream : glb_cmd_queue);
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso];
  encode_gather_args(enc, v, *m, *o, dg, u, *n, gd, *nb, b, bo);

  NSUInteger nthrds = 1024;
  NSUInteger total = (NSUInteger)(*m);
  MTLSize gridSize = MTLSizeMake(((total + nthrds - 1) / nthrds) * nthrds,
                                 1, 1);
  MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);

  [enc dispatchThreads:gridSize threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}

/**
 * Fortran wrapper for device scatter kernel
 */
void metal_scatter_kernel(void *v, int *m, void *dg,
                          void *u, int *n, void *gd,
                          int *nb, void *b, void *bo,
                          void *stream) {

  if ((*m) == 0) return;

  ensure_pipelines();

  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)(stream ? stream : glb_cmd_queue);
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_scatter];

  [enc setBuffer:(__bridge id<MTLBuffer>)v  offset:0 atIndex:0];
  [enc setBytes:m length:sizeof(int) atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)dg offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)u  offset:0 atIndex:3];
  [enc setBytes:n length:sizeof(int) atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)gd offset:0 atIndex:5];
  [enc setBytes:nb length:sizeof(int) atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)b  offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)bo offset:0 atIndex:8];

  NSUInteger nthrds = 1024;
  NSUInteger total = (NSUInteger)(*m);
  MTLSize gridSize = MTLSizeMake(((total + nthrds - 1) / nthrds) * nthrds,
                                 1, 1);
  MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);

  [enc dispatchThreads:gridSize threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}

/**
 * Pack send buffer on device
 */
void metal_gs_pack(void *u_d, void *buf_d, void *dof_d,
                   int offset, int n, void *stream) {

  if (n == 0) return;

  ensure_pipelines();

  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)(stream ? stream : glb_cmd_queue);
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_pack];

  [enc setBuffer:(__bridge id<MTLBuffer>)u_d   offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)buf_d offset:(offset * sizeof(float))
        atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)dof_d offset:(offset * sizeof(int))
        atIndex:2];
  [enc setBytes:&n length:sizeof(int) atIndex:3];

  NSUInteger nthrds = 1024;
  MTLSize gridSize = MTLSizeMake((((NSUInteger)n + nthrds - 1) / nthrds) * nthrds,
                                 1, 1);
  MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);

  [enc dispatchThreads:gridSize threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}

/**
 * Unpack receive buffer on device
 */
void metal_gs_unpack(void *u_d, int op, void *buf_d, void *dof_d,
                     int offset, int n, void *stream) {

  if (n == 0) return;

  ensure_pipelines();

  id<MTLComputePipelineState> pso;
  switch (op) {
  case GS_OP_ADD: pso = pso_unpack_add; break;
  case GS_OP_MIN: pso = pso_unpack_min; break;
  case GS_OP_MAX: pso = pso_unpack_max; break;
  default:
    fprintf(stderr, "%s: unknown gs op %d\n", __FILE__, op);
    abort();
  }

  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)(stream ? stream : glb_cmd_queue);
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso];

  [enc setBuffer:(__bridge id<MTLBuffer>)u_d   offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)buf_d offset:(offset * sizeof(float))
        atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)dof_d offset:(offset * sizeof(int))
        atIndex:2];
  [enc setBytes:&n length:sizeof(int) atIndex:3];

  NSUInteger nthrds = 1024;
  MTLSize gridSize = MTLSizeMake((((NSUInteger)n + nthrds - 1) / nthrds) * nthrds,
                                 1, 1);
  MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);

  [enc dispatchThreads:gridSize threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}
