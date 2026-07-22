/*
 Copyright (c) 2024-2026, The Neko Authors
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
 * Metal host-side dispatch for the atomic-free IDW gather kernel.
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

/* Cached pipeline state */
static id<MTLComputePipelineState> pso_idw_gather = nil;

/**
 * Fortran wrapper for the atomic-free IDW gather kernel.
 */
void metal_idw_gather_one_sided(void *fu, void *fv, void *fw,
                                void *fu_ib, void *fv_ib, void *fw_ib,
                                void *fum_ib, void *fvm_ib, void *fwm_ib,
                                void *x, void *y, void *z, void *ds,
                                void *pmsk, void *w, void *wm,
                                void *lpx, void *lpy, void *lpz,
                                void *active_el, void *el_off, void *el_lag,
                                int *n_active, int *lx3, real *dt,
                                real *rmax, real *pwr, real *eps, real *wtol) {

  const int total = (*n_active) * (*lx3);
  if (total < 1)
    return;

  if (pso_idw_gather == nil) {
    id<MTLDevice> device = neko_metal_device();
    id<MTLLibrary> lib = neko_metal_library();
    id<MTLFunction> func = [lib newFunctionWithName:@"idw_gather_one_sided"];
    if (func == nil) {
      fprintf(stderr, "Metal: kernel 'idw_gather_one_sided' "
              "not found in metallib\n");
      exit(EXIT_FAILURE);
    }
    NSError *error = nil;
    pso_idw_gather =
      [device newComputePipelineStateWithFunction:func error:&error];
    METAL_CHECK(error);
  }

  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  id<MTLCommandBuffer> cmdBuf = [queue commandBuffer];
  id<MTLComputeCommandEncoder> enc = [cmdBuf computeCommandEncoder];

  [enc setComputePipelineState:pso_idw_gather];

  [enc setBuffer:(__bridge id<MTLBuffer>)fu     offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)fv     offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)fw     offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)fu_ib  offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)fv_ib  offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)fw_ib  offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)fum_ib offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)fvm_ib offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)fwm_ib offset:0 atIndex:8];
  [enc setBuffer:(__bridge id<MTLBuffer>)x      offset:0 atIndex:9];
  [enc setBuffer:(__bridge id<MTLBuffer>)y      offset:0 atIndex:10];
  [enc setBuffer:(__bridge id<MTLBuffer>)z      offset:0 atIndex:11];
  [enc setBuffer:(__bridge id<MTLBuffer>)ds     offset:0 atIndex:12];
  [enc setBuffer:(__bridge id<MTLBuffer>)pmsk   offset:0 atIndex:13];
  [enc setBuffer:(__bridge id<MTLBuffer>)w      offset:0 atIndex:14];
  [enc setBuffer:(__bridge id<MTLBuffer>)wm     offset:0 atIndex:15];
  [enc setBuffer:(__bridge id<MTLBuffer>)lpx    offset:0 atIndex:16];
  [enc setBuffer:(__bridge id<MTLBuffer>)lpy    offset:0 atIndex:17];
  [enc setBuffer:(__bridge id<MTLBuffer>)lpz    offset:0 atIndex:18];
  [enc setBuffer:(__bridge id<MTLBuffer>)active_el offset:0 atIndex:19];
  [enc setBuffer:(__bridge id<MTLBuffer>)el_off    offset:0 atIndex:20];
  [enc setBuffer:(__bridge id<MTLBuffer>)el_lag    offset:0 atIndex:21];
  [enc setBytes:n_active length:sizeof(int)  atIndex:22];
  [enc setBytes:lx3      length:sizeof(int)  atIndex:23];
  [enc setBytes:dt       length:sizeof(real) atIndex:24];
  [enc setBytes:rmax     length:sizeof(real) atIndex:25];
  [enc setBytes:pwr      length:sizeof(real) atIndex:26];
  [enc setBytes:eps      length:sizeof(real) atIndex:27];
  [enc setBytes:wtol     length:sizeof(real) atIndex:28];

  NSUInteger nthrds = 256;
  MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);
  MTLSize numGroups =
    MTLSizeMake(((NSUInteger)total + nthrds - 1) / nthrds, 1, 1);

  [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
}
