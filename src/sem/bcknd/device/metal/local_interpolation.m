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
 * Metal host-side dispatch for local interpolation (find rst via Legendre).
 *
 * @note Apple GPUs do not support FP64. This backend operates in FP32.
 */

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#include <stdio.h>
#include <stdlib.h>
#include <device/device_config.h>
#include <device/metal/check.h>
#include <device/metal/kernel_utils.h>

void metal_local_interpolation_compute_weights(
    void *rst, void *zg, void *wr, void *ws, void *wt,
    const int *lx, const int *n) {
  if (*n <= 0) return;
  if (*lx < 1 || *lx > 16) {
    fprintf(stderr, "Unsupported interpolation order: %d\n", *lx);
    exit(1);
  }
  const int lx_r = *lx, n_r = *n;
  neko_metal_dispatch_1d(
    neko_metal_pipeline(@"local_interpolation_compute_weights_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>)rst offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>)zg offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>)wr offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>)ws offset:0 atIndex:3];
      [enc setBuffer:(__bridge id<MTLBuffer>)wt offset:0 atIndex:4];
      [enc setBytes:&lx_r length:sizeof(int) atIndex:5];
      [enc setBytes:&n_r length:sizeof(int) atIndex:6];
    }, (NSUInteger)n_r);
}

void metal_local_interpolation_compute_weights_3arrays(
    void *r, void *s, void *t, void *zg, void *wr, void *ws, void *wt,
    const int *lx, const int *n) {
  if (*n <= 0) return;
  if (*lx < 1 || *lx > 16) {
    fprintf(stderr, "Unsupported interpolation order: %d\n", *lx);
    exit(1);
  }
  const int lx_r = *lx, n_r = *n;
  neko_metal_dispatch_1d(
    neko_metal_pipeline(@"local_interpolation_compute_weights_3arrays_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>)r offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>)s offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>)t offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>)zg offset:0 atIndex:3];
      [enc setBuffer:(__bridge id<MTLBuffer>)wr offset:0 atIndex:4];
      [enc setBuffer:(__bridge id<MTLBuffer>)ws offset:0 atIndex:5];
      [enc setBuffer:(__bridge id<MTLBuffer>)wt offset:0 atIndex:6];
      [enc setBytes:&lx_r length:sizeof(int) atIndex:7];
      [enc setBytes:&n_r length:sizeof(int) atIndex:8];
    }, (NSUInteger)n_r);
}

/* Defined in device/metal/metal.m */
extern id<MTLDevice> neko_metal_device(void);
extern id<MTLLibrary> neko_metal_library(void);

/** Cached pipeline states for each LX (indices 2..16) */
static id<MTLComputePipelineState> pso_find_rst[17] = { nil };

/**
 * Create a compute pipeline state for the named kernel.
 */
static id<MTLComputePipelineState>
get_local_interpolation_pipeline(const char *name) {
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
 * Fortran wrapper for local interpolation
 */
void metal_find_rst_legendre(void *rst,
                             void *pt_x, void *pt_y, void *pt_z,
                             void *x_hat, void *y_hat, void *z_hat,
                             void *resx, void *resy, void *resz,
                             int *lx, void *el_ids, int *n_pt, real *tol,
                             void *conv_pts) {

  if (*n_pt <= 0)
    return;

  if (*lx < 2 || *lx > 16) {
    fprintf(stderr, "%s: size not supported: %d\n", __FILE__, *lx);
    exit(1);
  }

  if (pso_find_rst[*lx] == nil) {
    char name[64];
    snprintf(name, sizeof(name), "find_rst_legendre_kernel_lx%d", *lx);
    pso_find_rst[*lx] = get_local_interpolation_pipeline(name);
  }

  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_find_rst[*lx]];

  [enc setBuffer:(__bridge id<MTLBuffer>)rst   offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)pt_x  offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)pt_y  offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)pt_z  offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)x_hat offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)y_hat offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)z_hat offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)resx  offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)resy  offset:0 atIndex:8];
  [enc setBuffer:(__bridge id<MTLBuffer>)resz  offset:0 atIndex:9];
  [enc setBuffer:(__bridge id<MTLBuffer>)el_ids offset:0 atIndex:10];
  [enc setBytes:n_pt length:sizeof(int) atIndex:11];
  [enc setBytes:tol length:sizeof(real) atIndex:12];
  [enc setBuffer:(__bridge id<MTLBuffer>)conv_pts offset:0 atIndex:13];

  MTLSize threadsPerGroup = MTLSizeMake(1, 128, 1);
  MTLSize numGroups = MTLSizeMake(*n_pt, 1, 1);

  [enc dispatchThreadgroups:numGroups
      threadsPerThreadgroup:threadsPerGroup];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}
