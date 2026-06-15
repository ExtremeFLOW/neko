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
 * Metal host-side dispatch for SEM coefficient generation kernels.
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
static id<MTLComputePipelineState> pso_geo[17]  = { nil };  /* index by LX */
static id<MTLComputePipelineState> pso_dxyz[17] = { nil };
static id<MTLComputePipelineState> pso_drst     = nil;
/**
 * Return the shared Metal shader library (embedded in the binary).
 */
static id<MTLLibrary> get_coef_metal_library(void) {
  return neko_metal_library();
}

/**
 * Create a compute pipeline state for the named kernel.
 */
static id<MTLComputePipelineState>
get_coef_pipeline(const char *name) {
  id<MTLDevice> device = neko_metal_device();
  id<MTLLibrary> lib = get_coef_metal_library();
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

/*
 * Fortran-callable wrappers
 */

/**
 * Generate geometric factors G11..G33 on the Metal GPU.
 */
void metal_coef_generate_geo(void *G11, void *G12, void *G13,
                             void *G22, void *G23, void *G33,
                             void *drdx, void *drdy, void *drdz,
                             void *dsdx, void *dsdy, void *dsdz,
                             void *dtdx, void *dtdy, void *dtdz,
                             void *jacinv, void *w3,
                             int *nel, int *lx, int *gdim) {

  int LX = *lx;
  if (LX < 2 || LX > 16) {
    fprintf(stderr, "%s: size not supported: %d\n", __FILE__, LX);
    exit(1);
  }

  if (pso_geo[LX] == nil) {
    char name[64];
    snprintf(name, sizeof(name), "coef_generate_geo_kernel_lx%d", LX);
    pso_geo[LX] = get_coef_pipeline(name);
  }

  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_geo[LX]];

  [enc setBuffer:(__bridge id<MTLBuffer>)G11    offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)G12    offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)G13    offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)G22    offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)G23    offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)G33    offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)drdx   offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)drdy   offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)drdz   offset:0 atIndex:8];
  [enc setBuffer:(__bridge id<MTLBuffer>)dsdx   offset:0 atIndex:9];
  [enc setBuffer:(__bridge id<MTLBuffer>)dsdy   offset:0 atIndex:10];
  [enc setBuffer:(__bridge id<MTLBuffer>)dsdz   offset:0 atIndex:11];
  [enc setBuffer:(__bridge id<MTLBuffer>)dtdx   offset:0 atIndex:12];
  [enc setBuffer:(__bridge id<MTLBuffer>)dtdy   offset:0 atIndex:13];
  [enc setBuffer:(__bridge id<MTLBuffer>)dtdz   offset:0 atIndex:14];
  [enc setBuffer:(__bridge id<MTLBuffer>)jacinv offset:0 atIndex:15];
  [enc setBuffer:(__bridge id<MTLBuffer>)w3     offset:0 atIndex:16];
  [enc setBytes:gdim length:sizeof(int) atIndex:17];

  NSUInteger nthrds = 1024;
  MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);
  MTLSize numGroups = MTLSizeMake((NSUInteger)(*nel), 1, 1);

  [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}

/**
 * Generate dxyz/drst mappings and Jacobian on the Metal GPU.
 */
void metal_coef_generate_dxyzdrst(void *drdx, void *drdy, void *drdz,
                                  void *dsdx, void *dsdy, void *dsdz,
                                  void *dtdx, void *dtdy, void *dtdz,
                                  void *dxdr, void *dydr, void *dzdr,
                                  void *dxds, void *dyds, void *dzds,
                                  void *dxdt, void *dydt, void *dzdt,
                                  void *dx, void *dy, void *dz,
                                  void *x, void *y, void *z,
                                  void *jacinv, void *jac,
                                  int *lx, int *nel) {

  int LX = *lx;
  int n = (*nel) * LX * LX * LX;

  if (LX < 2 || LX > 16) {
    fprintf(stderr, "%s: size not supported: %d\n", __FILE__, LX);
    exit(1);
  }

  /* --- dxyz kernel (per-element, LX-specific) --- */
  if (pso_dxyz[LX] == nil) {
    char name[64];
    snprintf(name, sizeof(name), "coef_generate_dxyz_kernel_lx%d", LX);
    pso_dxyz[LX] = get_coef_pipeline(name);
  }

  {
    id<MTLCommandQueue> queue =
      (__bridge id<MTLCommandQueue>)glb_cmd_queue;
    id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

    [enc setComputePipelineState:pso_dxyz[LX]];

    [enc setBuffer:(__bridge id<MTLBuffer>)dxdr offset:0 atIndex:0];
    [enc setBuffer:(__bridge id<MTLBuffer>)dydr offset:0 atIndex:1];
    [enc setBuffer:(__bridge id<MTLBuffer>)dzdr offset:0 atIndex:2];
    [enc setBuffer:(__bridge id<MTLBuffer>)dxds offset:0 atIndex:3];
    [enc setBuffer:(__bridge id<MTLBuffer>)dyds offset:0 atIndex:4];
    [enc setBuffer:(__bridge id<MTLBuffer>)dzds offset:0 atIndex:5];
    [enc setBuffer:(__bridge id<MTLBuffer>)dxdt offset:0 atIndex:6];
    [enc setBuffer:(__bridge id<MTLBuffer>)dydt offset:0 atIndex:7];
    [enc setBuffer:(__bridge id<MTLBuffer>)dzdt offset:0 atIndex:8];
    [enc setBuffer:(__bridge id<MTLBuffer>)dx   offset:0 atIndex:9];
    [enc setBuffer:(__bridge id<MTLBuffer>)dy   offset:0 atIndex:10];
    [enc setBuffer:(__bridge id<MTLBuffer>)dz   offset:0 atIndex:11];
    [enc setBuffer:(__bridge id<MTLBuffer>)x    offset:0 atIndex:12];
    [enc setBuffer:(__bridge id<MTLBuffer>)y    offset:0 atIndex:13];
    [enc setBuffer:(__bridge id<MTLBuffer>)z    offset:0 atIndex:14];

    NSUInteger nthrds = 1024;
    MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);
    MTLSize numGroups = MTLSizeMake((NSUInteger)(*nel), 1, 1);

    [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
    [enc endEncoding];
    [cmdBuf commit];
    [cmdBuf waitUntilCompleted];
  }

  /* --- drst kernel (flat, not LX-specific) --- */
  if (pso_drst == nil) {
    pso_drst = get_coef_pipeline("coef_generate_drst_kernel");
  }

  {
    id<MTLCommandQueue> queue =
      (__bridge id<MTLCommandQueue>)glb_cmd_queue;
    id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

    [enc setComputePipelineState:pso_drst];

    [enc setBuffer:(__bridge id<MTLBuffer>)jac    offset:0 atIndex:0];
    [enc setBuffer:(__bridge id<MTLBuffer>)jacinv offset:0 atIndex:1];
    [enc setBuffer:(__bridge id<MTLBuffer>)drdx   offset:0 atIndex:2];
    [enc setBuffer:(__bridge id<MTLBuffer>)drdy   offset:0 atIndex:3];
    [enc setBuffer:(__bridge id<MTLBuffer>)drdz   offset:0 atIndex:4];
    [enc setBuffer:(__bridge id<MTLBuffer>)dsdx   offset:0 atIndex:5];
    [enc setBuffer:(__bridge id<MTLBuffer>)dsdy   offset:0 atIndex:6];
    [enc setBuffer:(__bridge id<MTLBuffer>)dsdz   offset:0 atIndex:7];
    [enc setBuffer:(__bridge id<MTLBuffer>)dtdx   offset:0 atIndex:8];
    [enc setBuffer:(__bridge id<MTLBuffer>)dtdy   offset:0 atIndex:9];
    [enc setBuffer:(__bridge id<MTLBuffer>)dtdz   offset:0 atIndex:10];
    [enc setBuffer:(__bridge id<MTLBuffer>)dxdr   offset:0 atIndex:11];
    [enc setBuffer:(__bridge id<MTLBuffer>)dydr   offset:0 atIndex:12];
    [enc setBuffer:(__bridge id<MTLBuffer>)dzdr   offset:0 atIndex:13];
    [enc setBuffer:(__bridge id<MTLBuffer>)dxds   offset:0 atIndex:14];
    [enc setBuffer:(__bridge id<MTLBuffer>)dyds   offset:0 atIndex:15];
    [enc setBuffer:(__bridge id<MTLBuffer>)dzds   offset:0 atIndex:16];
    [enc setBuffer:(__bridge id<MTLBuffer>)dxdt   offset:0 atIndex:17];
    [enc setBuffer:(__bridge id<MTLBuffer>)dydt   offset:0 atIndex:18];
    [enc setBuffer:(__bridge id<MTLBuffer>)dzdt   offset:0 atIndex:19];
    [enc setBytes:&n length:sizeof(int) atIndex:20];

    NSUInteger nthrds = 1024;
    NSUInteger total = (NSUInteger)n;
    MTLSize gridSize = MTLSizeMake(((total + nthrds - 1) / nthrds) * nthrds,
                                   1, 1);
    MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);

    [enc dispatchThreads:gridSize threadsPerThreadgroup:groupSize];
    [enc endEncoding];
    [cmdBuf commit];
    [cmdBuf waitUntilCompleted];
  }
}

/* Cached pipeline states for mass and area/normal kernels */
static id<MTLComputePipelineState> pso_mass     = nil;
static id<MTLComputePipelineState> pso_area[17] = { nil };

/**
 * Generate the mass matrix (diagonal) on the Metal GPU.
 */
void metal_coef_generate_mass(void *B, void *Binv, void *jac, void *w3,
                              int *lxyz, int *nel) {

  if (pso_mass == nil) {
    pso_mass = get_coef_pipeline("coef_generate_mass_kernel");
  }

  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_mass];

  [enc setBuffer:(__bridge id<MTLBuffer>)B    offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)Binv offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)jac  offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)w3   offset:0 atIndex:3];
  [enc setBytes:lxyz length:sizeof(int) atIndex:4];
  [enc setBytes:nel length:sizeof(int) atIndex:5];

  const NSUInteger n = (NSUInteger)(*lxyz) * (NSUInteger)(*nel);
  NSUInteger nthrds = 1024;
  MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);
  MTLSize numGroups = MTLSizeMake((n + nthrds - 1) / nthrds, 1, 1);

  [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}

/**
 * Generate facet areas and unit normals on the Metal GPU.
 */
void metal_coef_generate_area_and_normal(void *area, void *nx, void *ny,
                                         void *nz, void *dxdr, void *dydr,
                                         void *dzdr, void *dxds, void *dyds,
                                         void *dzds, void *dxdt, void *dydt,
                                         void *dzdt, void *wx, void *wy,
                                         void *wz, int *lx, int *nel,
                                         real eps) {

  int LX = *lx;
  if (LX < 2 || LX > 16) {
    fprintf(stderr, "%s: size not supported: %d\n", __FILE__, LX);
    exit(1);
  }

  if (pso_area[LX] == nil) {
    char name[64];
    snprintf(name, sizeof(name),
             "coef_generate_area_and_normal_kernel_lx%d", LX);
    pso_area[LX] = get_coef_pipeline(name);
  }

  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_area[LX]];

  [enc setBuffer:(__bridge id<MTLBuffer>)area offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)nx   offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)ny   offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)nz   offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)dxdr offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)dydr offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)dzdr offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)dxds offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)dyds offset:0 atIndex:8];
  [enc setBuffer:(__bridge id<MTLBuffer>)dzds offset:0 atIndex:9];
  [enc setBuffer:(__bridge id<MTLBuffer>)dxdt offset:0 atIndex:10];
  [enc setBuffer:(__bridge id<MTLBuffer>)dydt offset:0 atIndex:11];
  [enc setBuffer:(__bridge id<MTLBuffer>)dzdt offset:0 atIndex:12];
  [enc setBuffer:(__bridge id<MTLBuffer>)wx   offset:0 atIndex:13];
  [enc setBuffer:(__bridge id<MTLBuffer>)wy   offset:0 atIndex:14];
  [enc setBuffer:(__bridge id<MTLBuffer>)wz   offset:0 atIndex:15];
  [enc setBytes:&eps length:sizeof(real) atIndex:16];

  NSUInteger nthrds = 1024;
  MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);
  MTLSize numGroups = MTLSizeMake((NSUInteger)(*nel), 1, 1);

  [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}
