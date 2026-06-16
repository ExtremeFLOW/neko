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
 * Metal host-side dispatch for the compressible Euler residual.
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
static id<MTLComputePipelineState> pso_visc = nil;
static id<MTLComputePipelineState> pso_mx_flux = nil;
static id<MTLComputePipelineState> pso_my_flux = nil;
static id<MTLComputePipelineState> pso_mz_flux = nil;
static id<MTLComputePipelineState> pso_E_flux = nil;
static id<MTLComputePipelineState> pso_coef_mult = nil;
static id<MTLComputePipelineState> pso_rk_sum = nil;

/**
 * Create a compute pipeline state for the named kernel.
 */
static id<MTLComputePipelineState>
get_euler_res_pipeline(const char *name) {
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
euler_res_encoder(id<MTLCommandBuffer> *cmdBuf) {
  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  *cmdBuf = [queue commandBuffer];
  return [*cmdBuf computeCommandEncoder];
}

/**
 * Dispatch one thread per point and wait for completion.
 */
static void euler_res_dispatch(id<MTLCommandBuffer> cmdBuf,
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
 * Viscous part of the Euler residual on the Metal GPU.
 */
void euler_res_part_visc_metal(void *rhs, void *Binv, void *lap_sol,
                               void *effective_visc, int *n) {

  if (*n <= 0)
    return;

  if (pso_visc == nil)
    pso_visc = get_euler_res_pipeline("euler_res_part_visc_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = euler_res_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_visc];

  [enc setBuffer:(__bridge id<MTLBuffer>)rhs            offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)Binv           offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)lap_sol        offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)effective_visc offset:0 atIndex:3];
  [enc setBytes:n length:sizeof(int) atIndex:4];

  euler_res_dispatch(cmdBuf, enc, *n);
}

/**
 * Shared launcher for the three momentum flux kernels.
 */
static void euler_res_flux(id<MTLComputePipelineState> pso,
                           void *f_x, void *f_y, void *f_z,
                           void *m_x, void *m_y, void *m_z,
                           void *rho_field, void *p, int *n) {

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = euler_res_encoder(&cmdBuf);

  [enc setComputePipelineState:pso];

  [enc setBuffer:(__bridge id<MTLBuffer>)f_x       offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)f_y       offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)f_z       offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)m_x       offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)m_y       offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)m_z       offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)rho_field offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)p         offset:0 atIndex:7];
  [enc setBytes:n length:sizeof(int) atIndex:8];

  euler_res_dispatch(cmdBuf, enc, *n);
}

/**
 * x-momentum flux on the Metal GPU.
 */
void euler_res_part_mx_flux_metal(void *f_x, void *f_y, void *f_z,
                                  void *m_x, void *m_y, void *m_z,
                                  void *rho_field, void *p, int *n) {

  if (*n <= 0)
    return;

  if (pso_mx_flux == nil)
    pso_mx_flux = get_euler_res_pipeline("euler_res_part_mx_flux_kernel");

  euler_res_flux(pso_mx_flux, f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p, n);
}

/**
 * y-momentum flux on the Metal GPU.
 */
void euler_res_part_my_flux_metal(void *f_x, void *f_y, void *f_z,
                                  void *m_x, void *m_y, void *m_z,
                                  void *rho_field, void *p, int *n) {

  if (*n <= 0)
    return;

  if (pso_my_flux == nil)
    pso_my_flux = get_euler_res_pipeline("euler_res_part_my_flux_kernel");

  euler_res_flux(pso_my_flux, f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p, n);
}

/**
 * z-momentum flux on the Metal GPU.
 */
void euler_res_part_mz_flux_metal(void *f_x, void *f_y, void *f_z,
                                  void *m_x, void *m_y, void *m_z,
                                  void *rho_field, void *p, int *n) {

  if (*n <= 0)
    return;

  if (pso_mz_flux == nil)
    pso_mz_flux = get_euler_res_pipeline("euler_res_part_mz_flux_kernel");

  euler_res_flux(pso_mz_flux, f_x, f_y, f_z, m_x, m_y, m_z, rho_field, p, n);
}

/**
 * Energy flux on the Metal GPU.
 */
void euler_res_part_E_flux_metal(void *f_x, void *f_y, void *f_z,
                                 void *m_x, void *m_y, void *m_z,
                                 void *rho_field, void *p, void *E, int *n) {

  if (*n <= 0)
    return;

  if (pso_E_flux == nil)
    pso_E_flux = get_euler_res_pipeline("euler_res_part_E_flux_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = euler_res_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_E_flux];

  [enc setBuffer:(__bridge id<MTLBuffer>)f_x       offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)f_y       offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)f_z       offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)m_x       offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)m_y       offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)m_z       offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)rho_field offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)p         offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)E         offset:0 atIndex:8];
  [enc setBytes:n length:sizeof(int) atIndex:9];

  euler_res_dispatch(cmdBuf, enc, *n);
}

/**
 * Multiply the residual with the inverse multiplicity on the Metal GPU.
 */
void euler_res_part_coef_mult_metal(void *rhs_rho, void *rhs_m_x,
                                    void *rhs_m_y, void *rhs_m_z,
                                    void *rhs_E, void *mult, int *n) {

  if (*n <= 0)
    return;

  if (pso_coef_mult == nil)
    pso_coef_mult = get_euler_res_pipeline("euler_res_part_coef_mult_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = euler_res_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_coef_mult];

  [enc setBuffer:(__bridge id<MTLBuffer>)rhs_rho offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)rhs_m_x offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)rhs_m_y offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)rhs_m_z offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)rhs_E   offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)mult    offset:0 atIndex:5];
  [enc setBytes:n length:sizeof(int) atIndex:6];

  euler_res_dispatch(cmdBuf, enc, *n);
}

/**
 * Runge-Kutta stage summation on the Metal GPU.
 */
void euler_res_part_rk_sum_metal(void *rho, void *m_x, void *m_y, void *m_z,
                                 void *E, void *k_rho_i, void *k_m_x_i,
                                 void *k_m_y_i, void *k_m_z_i, void *k_E_i,
                                 real *dt, real *c, int *n) {

  if (*n <= 0)
    return;

  if (pso_rk_sum == nil)
    pso_rk_sum = get_euler_res_pipeline("euler_res_part_rk_sum_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = euler_res_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_rk_sum];

  [enc setBuffer:(__bridge id<MTLBuffer>)rho     offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)m_x     offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)m_y     offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)m_z     offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)E       offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)k_rho_i offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)k_m_x_i offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)k_m_y_i offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)k_m_z_i offset:0 atIndex:8];
  [enc setBuffer:(__bridge id<MTLBuffer>)k_E_i   offset:0 atIndex:9];
  [enc setBytes:dt length:sizeof(real) atIndex:10];
  [enc setBytes:c length:sizeof(real) atIndex:11];
  [enc setBytes:n length:sizeof(int) atIndex:12];

  euler_res_dispatch(cmdBuf, enc, *n);
}
