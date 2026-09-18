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
 * Metal host-side dispatch for the fused BiCGStab kernels.
 *
 * The reductions return the process local sums; the Fortran side combines
 * them across ranks. Partial sums are accumulated on the host rather than
 * with a second reduction kernel, which keeps every routine here to a
 * single command buffer.
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

#define BICGSTAB_NTHREADS 1024
#define BICGSTAB_ELEM_NTHREADS 256

/* Cached pipeline states */
static id<MTLComputePipelineState> pso_update_p = nil;
static id<MTLComputePipelineState> pso_product_and_norm = nil;
static id<MTLComputePipelineState> pso_part1 = nil;
static id<MTLComputePipelineState> pso_part2 = nil;

/* Partial sums, two values per threadgroup */
static id<MTLBuffer> bicgstab_buf = nil;
static int bicgstab_buf_cap = 0;

/**
 * Create a compute pipeline state for the named kernel.
 */
static id<MTLComputePipelineState>
get_bicgstab_pipeline(const char *name) {
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
static id<MTLCommandQueue> bicgstab_queue(void *stream) {
  if (stream != NULL)
    return (__bridge id<MTLCommandQueue>)stream;
  return (__bridge id<MTLCommandQueue>)glb_cmd_queue;
}

/**
 * A shared buffer holding at least @a nfloats partial sums.
 */
static id<MTLBuffer> get_bicgstab_buf(int nfloats) {
  if (nfloats > bicgstab_buf_cap) {
    int cap = nfloats > 2048 ? nfloats : 2048;
    bicgstab_buf = [neko_metal_device()
        newBufferWithLength:cap * sizeof(float)
                    options:MTLResourceStorageModeShared];
    bicgstab_buf_cap = cap;
  }
  return bicgstab_buf;
}

/**
 * Sum @a count partial sums, accumulating in double precision.
 */
static double sum_partials(const float *partials, int count) {
  double sum = 0.0;
  for (int i = 0; i < count; i++)
    sum += (double) partials[i];
  return sum;
}

/**
 * Search direction update \f$ p = r + \beta (p - \omega v) \f$.
 */
void metal_bicgstab_update_p(void *p, void *r, void *v, real *beta,
                             real *omega, int *n, void *stream) {

  if (*n <= 0)
    return;

  if (pso_update_p == nil)
    pso_update_p = get_bicgstab_pipeline("bicgstab_update_p_kernel");

  @autoreleasepool {
    const float fbeta = (float) *beta;
    const float fomega = (float) *omega;

    id<MTLCommandQueue> queue = bicgstab_queue(stream);
    id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

    [enc setComputePipelineState:pso_update_p];
    [enc setBuffer:(__bridge id<MTLBuffer>)p offset:0 atIndex:0];
    [enc setBuffer:(__bridge id<MTLBuffer>)r offset:0 atIndex:1];
    [enc setBuffer:(__bridge id<MTLBuffer>)v offset:0 atIndex:2];
    [enc setBytes:&fbeta length:sizeof(float) atIndex:3];
    [enc setBytes:&fomega length:sizeof(float) atIndex:4];
    [enc setBytes:n length:sizeof(int) atIndex:5];

    const NSUInteger nthrds = BICGSTAB_ELEM_NTHREADS;
    MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);
    MTLSize numGroups =
      MTLSizeMake(((NSUInteger)(*n) + nthrds - 1) / nthrds, 1, 1);

    [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
    [enc endEncoding];
    [cmdBuf commit];
    [cmdBuf waitUntilCompleted];
  }
}

/**
 * Weighted inner product and squared norm, \f$ (a^T M b, b^T M b) \f$.
 * The two process local sums are returned in @a res.
 */
void metal_bicgstab_product_and_norm(void *a, void *b, void *mult,
                                     real *res, int *n, void *stream) {

  res[0] = 0.0;
  res[1] = 0.0;
  if (*n <= 0)
    return;

  if (pso_product_and_norm == nil)
    pso_product_and_norm =
      get_bicgstab_pipeline("bicgstab_product_and_norm_kernel");

  @autoreleasepool {
    const NSUInteger nn = (NSUInteger)(*n);
    const NSUInteger tg_size = nn < BICGSTAB_NTHREADS ? nn
                                                      : BICGSTAB_NTHREADS;
    const NSUInteger nblocks = (nn + tg_size - 1) / tg_size;
    const int nb = (int) nblocks;

    id<MTLBuffer> buf = get_bicgstab_buf(2 * nb);

    id<MTLCommandQueue> queue = bicgstab_queue(stream);
    id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

    [enc setComputePipelineState:pso_product_and_norm];
    [enc setBuffer:(__bridge id<MTLBuffer>)a offset:0 atIndex:0];
    [enc setBuffer:(__bridge id<MTLBuffer>)b offset:0 atIndex:1];
    [enc setBuffer:(__bridge id<MTLBuffer>)mult offset:0 atIndex:2];
    [enc setBuffer:buf offset:0 atIndex:3];
    [enc setBytes:n length:sizeof(int) atIndex:4];
    [enc setBytes:&nb length:sizeof(int) atIndex:5];

    [enc dispatchThreadgroups:MTLSizeMake(nblocks, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
    [enc endEncoding];
    [cmdBuf commit];
    [cmdBuf waitUntilCompleted];

    const float *partials = (const float *) [buf contents];
    res[0] = (real) sum_partials(partials, nb);
    res[1] = (real) sum_partials(partials + nb, nb);
  }
}

/**
 * BiCGStab part 1: \f$ s = r - \alpha v \f$, returning the process local
 * \f$ s^T M s \f$.
 */
real metal_bicgstab_part1(void *s, void *r, void *v, void *mult,
                          real *alpha, int *n, void *stream) {

  if (*n <= 0)
    return 0.0;

  if (pso_part1 == nil)
    pso_part1 = get_bicgstab_pipeline("bicgstab_part1_kernel");

  real res = 0.0;
  @autoreleasepool {
    const float falpha = (float) *alpha;
    const NSUInteger nn = (NSUInteger)(*n);
    const NSUInteger tg_size = nn < BICGSTAB_NTHREADS ? nn
                                                      : BICGSTAB_NTHREADS;
    const NSUInteger nblocks = (nn + tg_size - 1) / tg_size;
    const int nb = (int) nblocks;

    id<MTLBuffer> buf = get_bicgstab_buf(nb);

    id<MTLCommandQueue> queue = bicgstab_queue(stream);
    id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

    [enc setComputePipelineState:pso_part1];
    [enc setBuffer:(__bridge id<MTLBuffer>)s offset:0 atIndex:0];
    [enc setBuffer:(__bridge id<MTLBuffer>)r offset:0 atIndex:1];
    [enc setBuffer:(__bridge id<MTLBuffer>)v offset:0 atIndex:2];
    [enc setBuffer:(__bridge id<MTLBuffer>)mult offset:0 atIndex:3];
    [enc setBuffer:buf offset:0 atIndex:4];
    [enc setBytes:&falpha length:sizeof(float) atIndex:5];
    [enc setBytes:n length:sizeof(int) atIndex:6];

    [enc dispatchThreadgroups:MTLSizeMake(nblocks, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
    [enc endEncoding];
    [cmdBuf commit];
    [cmdBuf waitUntilCompleted];

    res = (real) sum_partials((const float *) [buf contents], nb);
  }
  return res;
}

/**
 * BiCGStab part 2: \f$ x = x + \alpha \hat{p} + \omega \hat{s} \f$ and
 * \f$ r = s - \omega t \f$. The process local \f$ r^T M r \f$ and
 * \f$ f^T M r \f$ are returned in @a res.
 */
void metal_bicgstab_part2(void *x, void *r, void *p_hat, void *s_hat,
                          void *s, void *t, void *f, void *mult,
                          real *alpha, real *omega, real *res, int *n,
                          void *stream) {

  res[0] = 0.0;
  res[1] = 0.0;
  if (*n <= 0)
    return;

  if (pso_part2 == nil)
    pso_part2 = get_bicgstab_pipeline("bicgstab_part2_kernel");

  @autoreleasepool {
    const float falpha = (float) *alpha;
    const float fomega = (float) *omega;
    const NSUInteger nn = (NSUInteger)(*n);
    const NSUInteger tg_size = nn < BICGSTAB_NTHREADS ? nn
                                                      : BICGSTAB_NTHREADS;
    const NSUInteger nblocks = (nn + tg_size - 1) / tg_size;
    const int nb = (int) nblocks;

    id<MTLBuffer> buf = get_bicgstab_buf(2 * nb);

    id<MTLCommandQueue> queue = bicgstab_queue(stream);
    id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

    [enc setComputePipelineState:pso_part2];
    [enc setBuffer:(__bridge id<MTLBuffer>)x offset:0 atIndex:0];
    [enc setBuffer:(__bridge id<MTLBuffer>)r offset:0 atIndex:1];
    [enc setBuffer:(__bridge id<MTLBuffer>)p_hat offset:0 atIndex:2];
    [enc setBuffer:(__bridge id<MTLBuffer>)s_hat offset:0 atIndex:3];
    [enc setBuffer:(__bridge id<MTLBuffer>)s offset:0 atIndex:4];
    [enc setBuffer:(__bridge id<MTLBuffer>)t offset:0 atIndex:5];
    [enc setBuffer:(__bridge id<MTLBuffer>)f offset:0 atIndex:6];
    [enc setBuffer:(__bridge id<MTLBuffer>)mult offset:0 atIndex:7];
    [enc setBuffer:buf offset:0 atIndex:8];
    [enc setBytes:&falpha length:sizeof(float) atIndex:9];
    [enc setBytes:&fomega length:sizeof(float) atIndex:10];
    [enc setBytes:n length:sizeof(int) atIndex:11];
    [enc setBytes:&nb length:sizeof(int) atIndex:12];

    [enc dispatchThreadgroups:MTLSizeMake(nblocks, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
    [enc endEncoding];
    [cmdBuf commit];
    [cmdBuf waitUntilCompleted];

    const float *partials = (const float *) [buf contents];
    res[0] = (real) sum_partials(partials, nb);
    res[1] = (real) sum_partials(partials + nb, nb);
  }
}
