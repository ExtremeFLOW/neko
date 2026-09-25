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
 * Metal host-side dispatch for the incompressible Pn/Pn residuals.
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
static id<MTLComputePipelineState> pso_prs_part1 = nil;
static id<MTLComputePipelineState> pso_prs_part2 = nil;
static id<MTLComputePipelineState> pso_prs_part3 = nil;
static id<MTLComputePipelineState> pso_vel_update = nil;

/**
 * Create a compute pipeline state for the named kernel.
 */
static id<MTLComputePipelineState>
get_pnpn_res_pipeline(const char *name) {
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
pnpn_res_encoder(id<MTLCommandBuffer> *cmdBuf) {
  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  *cmdBuf = [queue commandBuffer];
  return [*cmdBuf computeCommandEncoder];
}

/**
 * Dispatch one thread per point and wait for completion.
 */
static void pnpn_res_dispatch(id<MTLCommandBuffer> cmdBuf,
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
 * Pressure residual, part 1.
 */
void pnpn_prs_res_part1_metal(void *ta1, void *ta2, void *ta3,
                              void *wa1, void *wa2, void *wa3,
                              void *f_u, void *f_v, void *f_w,
                              void *B, void *h1, real *mu, real *rho, int *n) {

  if (*n <= 0)
    return;

  if (pso_prs_part1 == nil)
    pso_prs_part1 = get_pnpn_res_pipeline("prs_res_part1_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = pnpn_res_encoder(&cmdBuf);
  [enc setComputePipelineState:pso_prs_part1];

  [enc setBuffer:(__bridge id<MTLBuffer>)ta1 offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)ta2 offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)ta3 offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)wa1 offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)wa2 offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)wa3 offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)f_u offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)f_v offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)f_w offset:0 atIndex:8];
  [enc setBuffer:(__bridge id<MTLBuffer>)B   offset:0 atIndex:9];
  [enc setBuffer:(__bridge id<MTLBuffer>)h1  offset:0 atIndex:10];
  [enc setBytes:mu length:sizeof(real) atIndex:11];
  [enc setBytes:rho length:sizeof(real) atIndex:12];
  [enc setBytes:n length:sizeof(int) atIndex:13];

  pnpn_res_dispatch(cmdBuf, enc, *n);
}

/**
 * Pressure residual, part 2.
 */
void pnpn_prs_res_part2_metal(void *p_res, void *wa1, void *wa2, void *wa3,
                              int *n) {

  if (*n <= 0)
    return;

  if (pso_prs_part2 == nil)
    pso_prs_part2 = get_pnpn_res_pipeline("prs_res_part2_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = pnpn_res_encoder(&cmdBuf);
  [enc setComputePipelineState:pso_prs_part2];

  [enc setBuffer:(__bridge id<MTLBuffer>)p_res offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)wa1   offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)wa2   offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)wa3   offset:0 atIndex:3];
  [enc setBytes:n length:sizeof(int) atIndex:4];

  pnpn_res_dispatch(cmdBuf, enc, *n);
}

/**
 * Pressure residual, part 3.
 */
void pnpn_prs_res_part3_metal(void *p_res, void *ta1, void *ta2, void *ta3,
                              real *dtbd, int *n) {

  if (*n <= 0)
    return;

  if (pso_prs_part3 == nil)
    pso_prs_part3 = get_pnpn_res_pipeline("prs_res_part3_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = pnpn_res_encoder(&cmdBuf);
  [enc setComputePipelineState:pso_prs_part3];

  [enc setBuffer:(__bridge id<MTLBuffer>)p_res offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)ta1   offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)ta2   offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)ta3   offset:0 atIndex:3];
  [enc setBytes:dtbd length:sizeof(real) atIndex:4];
  [enc setBytes:n length:sizeof(int) atIndex:5];

  pnpn_res_dispatch(cmdBuf, enc, *n);
}

/**
 * Velocity residual update.
 */
void pnpn_vel_res_update_metal(void *u_res, void *v_res, void *w_res,
                               void *ta1, void *ta2, void *ta3,
                               void *f_u, void *f_v, void *f_w, int *n) {

  if (*n <= 0)
    return;

  if (pso_vel_update == nil)
    pso_vel_update = get_pnpn_res_pipeline("vel_res_update_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = pnpn_res_encoder(&cmdBuf);
  [enc setComputePipelineState:pso_vel_update];

  [enc setBuffer:(__bridge id<MTLBuffer>)u_res offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)v_res offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)w_res offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)ta1   offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)ta2   offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)ta3   offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)f_u   offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)f_v   offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)f_w   offset:0 atIndex:8];
  [enc setBytes:n length:sizeof(int) atIndex:9];

  pnpn_res_dispatch(cmdBuf, enc, *n);
}
