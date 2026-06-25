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
 * Metal host-side dispatch for the rhs makers (sumab/ext/bdf/oifs).
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
static id<MTLComputePipelineState> pso_sumab = nil;
static id<MTLComputePipelineState> pso_makeext = nil;
static id<MTLComputePipelineState> pso_scalar_makeext = nil;
static id<MTLComputePipelineState> pso_makebdf = nil;
static id<MTLComputePipelineState> pso_scalar_makebdf = nil;
static id<MTLComputePipelineState> pso_makeoifs = nil;
static id<MTLComputePipelineState> pso_scalar_makeoifs = nil;

/**
 * Create a compute pipeline state for the named kernel.
 */
static id<MTLComputePipelineState>
get_rhs_maker_pipeline(const char *name) {
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
 * Create a compute command encoder on the global command queue.
 */
static id<MTLComputeCommandEncoder>
rhs_maker_encoder(id<MTLCommandBuffer> *cmdBuf) {
  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  *cmdBuf = [queue commandBuffer];
  return [*cmdBuf computeCommandEncoder];
}

/**
 * Dispatch one thread per point and wait for completion.
 */
static void rhs_maker_dispatch(id<MTLCommandBuffer> cmdBuf,
                               id<MTLComputeCommandEncoder> enc, int n) {
  const NSUInteger nthrds = 256;
  MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);
  MTLSize numGroups = MTLSizeMake(((NSUInteger)n + nthrds - 1) / nthrds, 1, 1);

  [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}

/**
 * Adams-Bashforth extrapolation of the velocity field.
 */
void rhs_maker_sumab_metal(void *u, void *v, void *w,
                           void *uu, void *vv, void *ww,
                           void *ulag1, void *ulag2, void *vlag1, void *vlag2,
                           void *wlag1, void *wlag2,
                           real *ab1, real *ab2, real *ab3, int *nab, int *n) {

  if (*n <= 0)
    return;

  if (pso_sumab == nil)
    pso_sumab = get_rhs_maker_pipeline("sumab_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = rhs_maker_encoder(&cmdBuf);
  [enc setComputePipelineState:pso_sumab];

  [enc setBuffer:(__bridge id<MTLBuffer>)u     offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)v     offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)w     offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)uu    offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)vv    offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)ww    offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)ulag1 offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)ulag2 offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)vlag1 offset:0 atIndex:8];
  [enc setBuffer:(__bridge id<MTLBuffer>)vlag2 offset:0 atIndex:9];
  [enc setBuffer:(__bridge id<MTLBuffer>)wlag1 offset:0 atIndex:10];
  [enc setBuffer:(__bridge id<MTLBuffer>)wlag2 offset:0 atIndex:11];
  [enc setBytes:ab1 length:sizeof(real) atIndex:12];
  [enc setBytes:ab2 length:sizeof(real) atIndex:13];
  [enc setBytes:ab3 length:sizeof(real) atIndex:14];
  [enc setBytes:nab length:sizeof(int) atIndex:15];
  [enc setBytes:n length:sizeof(int) atIndex:16];

  rhs_maker_dispatch(cmdBuf, enc, *n);
}

/**
 * Extrapolate the velocity source terms.
 */
void rhs_maker_ext_metal(void *abx1, void *aby1, void *abz1,
                         void *abx2, void *aby2, void *abz2,
                         void *bfx, void *bfy, void *bfz,
                         real *rho, real *ab1, real *ab2, real *ab3, int *n) {

  if (*n <= 0)
    return;

  if (pso_makeext == nil)
    pso_makeext = get_rhs_maker_pipeline("makeext_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = rhs_maker_encoder(&cmdBuf);
  [enc setComputePipelineState:pso_makeext];

  [enc setBuffer:(__bridge id<MTLBuffer>)abx1 offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)aby1 offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)abz1 offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)abx2 offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)aby2 offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)abz2 offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)bfx  offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)bfy  offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)bfz  offset:0 atIndex:8];
  [enc setBytes:rho length:sizeof(real) atIndex:9];
  [enc setBytes:ab1 length:sizeof(real) atIndex:10];
  [enc setBytes:ab2 length:sizeof(real) atIndex:11];
  [enc setBytes:ab3 length:sizeof(real) atIndex:12];
  [enc setBytes:n length:sizeof(int) atIndex:13];

  rhs_maker_dispatch(cmdBuf, enc, *n);
}

/**
 * Extrapolate the scalar source term.
 */
void scalar_rhs_maker_ext_metal(void *fs_lag, void *fs_laglag, void *fs,
                                real *rho, real *ext1, real *ext2, real *ext3,
                                int *n) {

  if (*n <= 0)
    return;

  if (pso_scalar_makeext == nil)
    pso_scalar_makeext = get_rhs_maker_pipeline("scalar_makeext_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = rhs_maker_encoder(&cmdBuf);
  [enc setComputePipelineState:pso_scalar_makeext];

  [enc setBuffer:(__bridge id<MTLBuffer>)fs_lag    offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)fs_laglag offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)fs        offset:0 atIndex:2];
  [enc setBytes:rho length:sizeof(real) atIndex:3];
  [enc setBytes:ext1 length:sizeof(real) atIndex:4];
  [enc setBytes:ext2 length:sizeof(real) atIndex:5];
  [enc setBytes:ext3 length:sizeof(real) atIndex:6];
  [enc setBytes:n length:sizeof(int) atIndex:7];

  rhs_maker_dispatch(cmdBuf, enc, *n);
}

/**
 * BDF contribution to the velocity source terms.
 */
void rhs_maker_bdf_metal(void *ulag1, void *ulag2, void *vlag1, void *vlag2,
                         void *wlag1, void *wlag2,
                         void *bfx, void *bfy, void *bfz,
                         void *u, void *v, void *w, void *B,
                         real *rho, real *dt, real *bd2, real *bd3, real *bd4,
                         int *nbd, int *n) {

  if (*n <= 0)
    return;

  if (pso_makebdf == nil)
    pso_makebdf = get_rhs_maker_pipeline("makebdf_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = rhs_maker_encoder(&cmdBuf);
  [enc setComputePipelineState:pso_makebdf];

  [enc setBuffer:(__bridge id<MTLBuffer>)ulag1 offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)ulag2 offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)vlag1 offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)vlag2 offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)wlag1 offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)wlag2 offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)bfx   offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)bfy   offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)bfz   offset:0 atIndex:8];
  [enc setBuffer:(__bridge id<MTLBuffer>)u     offset:0 atIndex:9];
  [enc setBuffer:(__bridge id<MTLBuffer>)v     offset:0 atIndex:10];
  [enc setBuffer:(__bridge id<MTLBuffer>)w     offset:0 atIndex:11];
  [enc setBuffer:(__bridge id<MTLBuffer>)B     offset:0 atIndex:12];
  [enc setBytes:rho length:sizeof(real) atIndex:13];
  [enc setBytes:dt length:sizeof(real) atIndex:14];
  [enc setBytes:bd2 length:sizeof(real) atIndex:15];
  [enc setBytes:bd3 length:sizeof(real) atIndex:16];
  [enc setBytes:bd4 length:sizeof(real) atIndex:17];
  [enc setBytes:nbd length:sizeof(int) atIndex:18];
  [enc setBytes:n length:sizeof(int) atIndex:19];

  rhs_maker_dispatch(cmdBuf, enc, *n);
}

/**
 * BDF contribution to the scalar source term.
 */
void scalar_rhs_maker_bdf_metal(void *s_lag, void *s_laglag, void *fs,
                                void *s, void *B,
                                real *rho, real *dt, real *bd2, real *bd3,
                                real *bd4, int *nbd, int *n) {

  if (*n <= 0)
    return;

  if (pso_scalar_makebdf == nil)
    pso_scalar_makebdf = get_rhs_maker_pipeline("scalar_makebdf_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = rhs_maker_encoder(&cmdBuf);
  [enc setComputePipelineState:pso_scalar_makebdf];

  [enc setBuffer:(__bridge id<MTLBuffer>)s_lag    offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)s_laglag offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)fs       offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)s        offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)B        offset:0 atIndex:4];
  [enc setBytes:rho length:sizeof(real) atIndex:5];
  [enc setBytes:dt length:sizeof(real) atIndex:6];
  [enc setBytes:bd2 length:sizeof(real) atIndex:7];
  [enc setBytes:bd3 length:sizeof(real) atIndex:8];
  [enc setBytes:bd4 length:sizeof(real) atIndex:9];
  [enc setBytes:nbd length:sizeof(int) atIndex:10];
  [enc setBytes:n length:sizeof(int) atIndex:11];

  rhs_maker_dispatch(cmdBuf, enc, *n);
}

/**
 * OIFS contribution to the velocity source terms.
 */
void rhs_maker_oifs_metal(void *phi_x, void *phi_y, void *phi_z,
                          void *bf_x, void *bf_y, void *bf_z,
                          real *rho, real *dt, int *n) {

  if (*n <= 0)
    return;

  if (pso_makeoifs == nil)
    pso_makeoifs = get_rhs_maker_pipeline("makeoifs_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = rhs_maker_encoder(&cmdBuf);
  [enc setComputePipelineState:pso_makeoifs];

  [enc setBuffer:(__bridge id<MTLBuffer>)phi_x offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)phi_y offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)phi_z offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)bf_x  offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)bf_y  offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)bf_z  offset:0 atIndex:5];
  [enc setBytes:rho length:sizeof(real) atIndex:6];
  [enc setBytes:dt length:sizeof(real) atIndex:7];
  [enc setBytes:n length:sizeof(int) atIndex:8];

  rhs_maker_dispatch(cmdBuf, enc, *n);
}

/**
 * OIFS contribution to the scalar source term.
 */
void scalar_rhs_maker_oifs_metal(void *phi_s, void *bf_s,
                                 real *rho, real *dt, int *n) {

  if (*n <= 0)
    return;

  if (pso_scalar_makeoifs == nil)
    pso_scalar_makeoifs = get_rhs_maker_pipeline("scalar_makeoifs_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = rhs_maker_encoder(&cmdBuf);
  [enc setComputePipelineState:pso_scalar_makeoifs];

  [enc setBuffer:(__bridge id<MTLBuffer>)phi_s offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)bf_s  offset:0 atIndex:1];
  [enc setBytes:rho length:sizeof(real) atIndex:2];
  [enc setBytes:dt length:sizeof(real) atIndex:3];
  [enc setBytes:n length:sizeof(int) atIndex:4];

  rhs_maker_dispatch(cmdBuf, enc, *n);
}
