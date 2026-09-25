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
 * Metal host-side dispatch for compressible flow update operations.
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
static id<MTLComputePipelineState> pso_update_uvw = nil;
static id<MTLComputePipelineState> pso_update_mxyz_p_ruvw = nil;
static id<MTLComputePipelineState> pso_update_e = nil;
static id<MTLComputePipelineState> pso_update_temperature = nil;
static id<MTLComputePipelineState> pso_ns_flux_prepare = nil;
static id<MTLComputePipelineState> pso_ns_flux_finalize = nil;
static id<MTLComputePipelineState> pso_ns_flux_temperature = nil;

/**
 * Create a compute pipeline state for the named kernel.
 */
static id<MTLComputePipelineState>
get_compressible_ops_pipeline(const char *name) {
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
compressible_ops_encoder(id<MTLCommandBuffer> *cmdBuf) {
  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  *cmdBuf = [queue commandBuffer];
  return [*cmdBuf computeCommandEncoder];
}

/**
 * Dispatch one thread per point and wait for completion.
 */
static void compressible_ops_dispatch(id<MTLCommandBuffer> cmdBuf,
                                      id<MTLComputeCommandEncoder> enc,
                                      int n) {
  const NSUInteger nthrds = 256;
  MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);
  MTLSize numGroups = MTLSizeMake(((NSUInteger)n + nthrds - 1) / nthrds, 1, 1);

  [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}

/**
 * Update the velocity components u = m / rho on the Metal GPU.
 */
void metal_update_uvw(void *u, void *v, void *w,
                      void *m_x, void *m_y, void *m_z,
                      void *rho, int n) {

  if (n <= 0)
    return;

  if (pso_update_uvw == nil)
    pso_update_uvw = get_compressible_ops_pipeline("update_uvw_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = compressible_ops_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_update_uvw];

  [enc setBuffer:(__bridge id<MTLBuffer>)u   offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)v   offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)w   offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)m_x offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)m_y offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)m_z offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)rho offset:0 atIndex:6];
  [enc setBytes:&n length:sizeof(int) atIndex:7];

  compressible_ops_dispatch(cmdBuf, enc, n);
}

/**
 * Update momentum, pressure and kinetic energy on the Metal GPU.
 */
void metal_update_mxyz_p_ruvw(void *m_x, void *m_y, void *m_z,
                              void *p, void *ruvw,
                              void *u, void *v, void *w,
                              void *E, void *rho, real gamma, int n) {

  if (n <= 0)
    return;

  if (pso_update_mxyz_p_ruvw == nil)
    pso_update_mxyz_p_ruvw =
      get_compressible_ops_pipeline("update_mxyz_p_ruvw_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = compressible_ops_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_update_mxyz_p_ruvw];

  [enc setBuffer:(__bridge id<MTLBuffer>)m_x  offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)m_y  offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)m_z  offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)p    offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)ruvw offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)u    offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)v    offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)w    offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)E    offset:0 atIndex:8];
  [enc setBuffer:(__bridge id<MTLBuffer>)rho  offset:0 atIndex:9];
  [enc setBytes:&gamma length:sizeof(real) atIndex:10];
  [enc setBytes:&n length:sizeof(int) atIndex:11];

  compressible_ops_dispatch(cmdBuf, enc, n);
}

/**
 * Update the total energy on the Metal GPU.
 */
void metal_update_e(void *E, void *p, void *ruvw, real gamma, int n) {

  if (n <= 0)
    return;

  if (pso_update_e == nil)
    pso_update_e = get_compressible_ops_pipeline("update_e_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = compressible_ops_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_update_e];

  [enc setBuffer:(__bridge id<MTLBuffer>)E    offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)p    offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)ruvw offset:0 atIndex:2];
  [enc setBytes:&gamma length:sizeof(real) atIndex:3];
  [enc setBytes:&n length:sizeof(int) atIndex:4];

  compressible_ops_dispatch(cmdBuf, enc, n);
}

/**
 * Update the temperature T = p / (rho * (gamma - 1)) on the Metal GPU.
 */
void metal_update_temperature(void *T, void *p, void *rho,
                              real gamma, int n) {

  if (n <= 0)
    return;

  if (pso_update_temperature == nil)
    pso_update_temperature =
      get_compressible_ops_pipeline("update_temperature_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = compressible_ops_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_update_temperature];

  [enc setBuffer:(__bridge id<MTLBuffer>)T   offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)p   offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)rho offset:0 atIndex:2];
  [enc setBytes:&gamma length:sizeof(real) atIndex:3];
  [enc setBytes:&n length:sizeof(int) atIndex:4];

  compressible_ops_dispatch(cmdBuf, enc, n);
}

/**
 * Prepare physical Navier-Stokes flux work arrays on the Metal GPU.
 */
void metal_ns_flux_prepare(void *div_flux, void *dissipation, void *h1,
                           void *dudx, void *dudy, void *dudz,
                           void *dvdx, void *dvdy, void *dvdz,
                           void *dwdx, void *dwdy, void *dwdz,
                           void *mu, int n) {

  if (n <= 0)
    return;

  if (pso_ns_flux_prepare == nil)
    pso_ns_flux_prepare =
      get_compressible_ops_pipeline("ns_flux_prepare_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = compressible_ops_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_ns_flux_prepare];

  [enc setBuffer:(__bridge id<MTLBuffer>)div_flux    offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)dissipation offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)h1          offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)dudx offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)dudy offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)dudz offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)dvdx offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)dvdy offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)dvdz offset:0 atIndex:8];
  [enc setBuffer:(__bridge id<MTLBuffer>)dwdx offset:0 atIndex:9];
  [enc setBuffer:(__bridge id<MTLBuffer>)dwdy offset:0 atIndex:10];
  [enc setBuffer:(__bridge id<MTLBuffer>)dwdz offset:0 atIndex:11];
  [enc setBuffer:(__bridge id<MTLBuffer>)mu   offset:0 atIndex:12];
  [enc setBytes:&n length:sizeof(int) atIndex:13];

  compressible_ops_dispatch(cmdBuf, enc, n);
}

/**
 * Finish physical Navier-Stokes flux assembly on the Metal GPU.
 */
void metal_ns_flux_finalize(void *visc_m_x, void *visc_m_y, void *visc_m_z,
                            void *visc_E, void *f_x, void *f_y, void *f_z,
                            void *opgrad_x, void *opgrad_y, void *opgrad_z,
                            void *u, void *v, void *w, void *B,
                            void *dissipation, int n) {

  if (n <= 0)
    return;

  if (pso_ns_flux_finalize == nil)
    pso_ns_flux_finalize =
      get_compressible_ops_pipeline("ns_flux_finalize_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = compressible_ops_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_ns_flux_finalize];

  [enc setBuffer:(__bridge id<MTLBuffer>)visc_m_x offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)visc_m_y offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)visc_m_z offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)visc_E   offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)f_x offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)f_y offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)f_z offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)opgrad_x offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)opgrad_y offset:0 atIndex:8];
  [enc setBuffer:(__bridge id<MTLBuffer>)opgrad_z offset:0 atIndex:9];
  [enc setBuffer:(__bridge id<MTLBuffer>)u offset:0 atIndex:10];
  [enc setBuffer:(__bridge id<MTLBuffer>)v offset:0 atIndex:11];
  [enc setBuffer:(__bridge id<MTLBuffer>)w offset:0 atIndex:12];
  [enc setBuffer:(__bridge id<MTLBuffer>)B offset:0 atIndex:13];
  [enc setBuffer:(__bridge id<MTLBuffer>)dissipation offset:0 atIndex:14];
  [enc setBytes:&n length:sizeof(int) atIndex:15];

  compressible_ops_dispatch(cmdBuf, enc, n);
}

/**
 * Prepare temperature and conductivity for the energy flux on the Metal GPU.
 */
void metal_ns_flux_temperature(void *div_flux, void *h1, void *p, void *rho,
                               void *kappa, real gamma, int n) {

  if (n <= 0)
    return;

  if (pso_ns_flux_temperature == nil)
    pso_ns_flux_temperature =
      get_compressible_ops_pipeline("ns_flux_temperature_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = compressible_ops_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_ns_flux_temperature];

  [enc setBuffer:(__bridge id<MTLBuffer>)div_flux offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)h1       offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)p     offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)rho   offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)kappa offset:0 atIndex:4];
  [enc setBytes:&gamma length:sizeof(real) atIndex:5];
  [enc setBytes:&n length:sizeof(int) atIndex:6];

  compressible_ops_dispatch(cmdBuf, enc, n);
}
