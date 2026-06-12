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
 * Metal host-side dispatch for entropy viscosity regularization.
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
static id<MTLComputePipelineState> pso_ev_residual = nil;
static id<MTLComputePipelineState> pso_ev_viscosity = nil;
static id<MTLComputePipelineState> pso_ev_element_max = nil;
static id<MTLComputePipelineState> pso_ev_clamp = nil;
static id<MTLComputePipelineState> pso_ev_smooth = nil;

/**
 * Create a compute pipeline state for the named kernel.
 */
static id<MTLComputePipelineState>
get_entropy_visc_pipeline(const char *name) {
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
entropy_visc_encoder(id<MTLCommandBuffer> *cmdBuf) {
  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  *cmdBuf = [queue commandBuffer];
  return [*cmdBuf computeCommandEncoder];
}

/**
 * Dispatch one thread per point and wait for completion.
 */
static void entropy_visc_dispatch(id<MTLCommandBuffer> cmdBuf,
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
 * Compute the entropy residual (BDF time derivative) on the Metal GPU.
 */
void metal_entropy_visc_compute_residual(void *entropy_residual,
                                         void *S, void *S_lag1,
                                         void *S_lag2, void *S_lag3,
                                         real bdf1, real bdf2,
                                         real bdf3, real bdf4,
                                         real dt, int n) {

  if (n <= 0)
    return;

  if (pso_ev_residual == nil)
    pso_ev_residual =
      get_entropy_visc_pipeline("entropy_visc_compute_residual_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = entropy_visc_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_ev_residual];

  [enc setBuffer:(__bridge id<MTLBuffer>)entropy_residual offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)S      offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)S_lag1 offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)S_lag2 offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)S_lag3 offset:0 atIndex:4];
  [enc setBytes:&bdf1 length:sizeof(real) atIndex:5];
  [enc setBytes:&bdf2 length:sizeof(real) atIndex:6];
  [enc setBytes:&bdf3 length:sizeof(real) atIndex:7];
  [enc setBytes:&bdf4 length:sizeof(real) atIndex:8];
  [enc setBytes:&dt length:sizeof(real) atIndex:9];
  [enc setBytes:&n length:sizeof(int) atIndex:10];

  entropy_visc_dispatch(cmdBuf, enc, n);
}

/**
 * Compute viscosity from the entropy residual on the Metal GPU.
 */
void metal_entropy_visc_compute_viscosity(void *reg_coeff,
                                          void *entropy_residual,
                                          void *h, real c_avisc_entropy,
                                          real n_S, int n) {

  if (n <= 0)
    return;

  if (pso_ev_viscosity == nil)
    pso_ev_viscosity =
      get_entropy_visc_pipeline("entropy_visc_compute_viscosity_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = entropy_visc_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_ev_viscosity];

  [enc setBuffer:(__bridge id<MTLBuffer>)reg_coeff        offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)entropy_residual offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)h                offset:0 atIndex:2];
  [enc setBytes:&c_avisc_entropy length:sizeof(real) atIndex:3];
  [enc setBytes:&n_S length:sizeof(real) atIndex:4];
  [enc setBytes:&n length:sizeof(int) atIndex:5];

  entropy_visc_dispatch(cmdBuf, enc, n);
}

/**
 * Apply the element-wise maximum on the Metal GPU.
 */
void metal_entropy_visc_apply_element_max(void *reg_coeff, int lx, int nelv) {

  if (nelv <= 0)
    return;

  if (pso_ev_element_max == nil)
    pso_ev_element_max =
      get_entropy_visc_pipeline("entropy_visc_apply_element_max_kernel");

  const int lx3 = lx * lx * lx;

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = entropy_visc_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_ev_element_max];

  [enc setBuffer:(__bridge id<MTLBuffer>)reg_coeff offset:0 atIndex:0];
  [enc setBytes:&lx3 length:sizeof(int) atIndex:1];

  MTLSize groupSize = MTLSizeMake(256, 1, 1);
  MTLSize numGroups = MTLSizeMake((NSUInteger)nelv, 1, 1);

  [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}

/**
 * Clamp to the low-order viscosity on the Metal GPU.
 */
void metal_entropy_visc_clamp_to_low_order(void *reg_coeff, void *h,
                                           void *max_wave_speed,
                                           real c_avisc_low, int n) {

  if (n <= 0)
    return;

  if (pso_ev_clamp == nil)
    pso_ev_clamp =
      get_entropy_visc_pipeline("entropy_visc_clamp_to_low_order_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = entropy_visc_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_ev_clamp];

  [enc setBuffer:(__bridge id<MTLBuffer>)reg_coeff      offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)h              offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)max_wave_speed offset:0 atIndex:2];
  [enc setBytes:&c_avisc_low length:sizeof(real) atIndex:3];
  [enc setBytes:&n length:sizeof(int) atIndex:4];

  entropy_visc_dispatch(cmdBuf, enc, n);
}

/**
 * Divide by the multiplicity (smoothing) on the Metal GPU.
 */
void metal_entropy_visc_smooth_divide(void *reg_coeff, void *temp_field,
                                      void *mult_field, int n) {

  if (n <= 0)
    return;

  if (pso_ev_smooth == nil)
    pso_ev_smooth =
      get_entropy_visc_pipeline("entropy_visc_smooth_divide_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = entropy_visc_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_ev_smooth];

  [enc setBuffer:(__bridge id<MTLBuffer>)reg_coeff  offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)temp_field offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)mult_field offset:0 atIndex:2];
  [enc setBytes:&n length:sizeof(int) atIndex:3];

  entropy_visc_dispatch(cmdBuf, enc, n);
}
