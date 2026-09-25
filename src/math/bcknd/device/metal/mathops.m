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
 * Metal host-side dispatch for field operations (mathops).
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
static id<MTLComputePipelineState> pso_opchsign = nil;
static id<MTLComputePipelineState> pso_opcolv = nil;
static id<MTLComputePipelineState> pso_opcolv3c = nil;
static id<MTLComputePipelineState> pso_opadd2cm = nil;
static id<MTLComputePipelineState> pso_opadd2col = nil;

/**
 * Create a compute pipeline state for the named kernel.
 */
static id<MTLComputePipelineState>
get_mathops_pipeline(const char *name) {
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
mathops_encoder(id<MTLCommandBuffer> *cmdBuf) {
  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  *cmdBuf = [queue commandBuffer];
  return [*cmdBuf computeCommandEncoder];
}

/**
 * Dispatch one thread per point and wait for completion.
 */
static void mathops_dispatch(id<MTLCommandBuffer> cmdBuf,
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
 * \f$ a = -a \f$ for vector fields on the Metal GPU.
 */
void metal_opchsign(void *a1, void *a2, void *a3, int *gdim, int *n) {

  if (*n <= 0)
    return;

  if (pso_opchsign == nil)
    pso_opchsign = get_mathops_pipeline("opchsign_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = mathops_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_opchsign];

  [enc setBuffer:(__bridge id<MTLBuffer>)a1 offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)a2 offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)a3 offset:0 atIndex:2];
  [enc setBytes:gdim length:sizeof(int) atIndex:3];
  [enc setBytes:n length:sizeof(int) atIndex:4];

  mathops_dispatch(cmdBuf, enc, *n);
}

/**
 * \f$ a = a \cdot c \f$ for vector fields on the Metal GPU.
 */
void metal_opcolv(void *a1, void *a2, void *a3, void *c, int *gdim, int *n) {

  if (*n <= 0)
    return;

  if (pso_opcolv == nil)
    pso_opcolv = get_mathops_pipeline("opcolv_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = mathops_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_opcolv];

  [enc setBuffer:(__bridge id<MTLBuffer>)a1 offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)a2 offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)a3 offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)c  offset:0 atIndex:3];
  [enc setBytes:gdim length:sizeof(int) atIndex:4];
  [enc setBytes:n length:sizeof(int) atIndex:5];

  mathops_dispatch(cmdBuf, enc, *n);
}

/**
 * \f$ a = b \cdot c \cdot d \f$ for vector fields on the Metal GPU.
 */
void metal_opcolv3c(void *a1, void *a2, void *a3,
                    void *b1, void *b2, void *b3,
                    void *c, real *d, int *gdim, int *n) {

  if (*n <= 0)
    return;

  if (pso_opcolv3c == nil)
    pso_opcolv3c = get_mathops_pipeline("opcolv3c_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = mathops_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_opcolv3c];

  [enc setBuffer:(__bridge id<MTLBuffer>)a1 offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)a2 offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)a3 offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)b1 offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)b2 offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)b3 offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)c  offset:0 atIndex:6];
  [enc setBytes:d length:sizeof(real) atIndex:7];
  [enc setBytes:gdim length:sizeof(int) atIndex:8];
  [enc setBytes:n length:sizeof(int) atIndex:9];

  mathops_dispatch(cmdBuf, enc, *n);
}

/**
 * \f$ a = a + b \cdot c \f$ (scalar c) for vector fields on the Metal GPU.
 */
void metal_opadd2cm(void *a1, void *a2, void *a3,
                    void *b1, void *b2, void *b3,
                    real *c, int *gdim, int *n) {

  if (*n <= 0)
    return;

  if (pso_opadd2cm == nil)
    pso_opadd2cm = get_mathops_pipeline("opadd2cm_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = mathops_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_opadd2cm];

  [enc setBuffer:(__bridge id<MTLBuffer>)a1 offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)a2 offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)a3 offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)b1 offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)b2 offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)b3 offset:0 atIndex:5];
  [enc setBytes:c length:sizeof(real) atIndex:6];
  [enc setBytes:gdim length:sizeof(int) atIndex:7];
  [enc setBytes:n length:sizeof(int) atIndex:8];

  mathops_dispatch(cmdBuf, enc, *n);
}

/**
 * \f$ a = a + b \cdot c \f$ (vector c) for vector fields on the Metal GPU.
 */
void metal_opadd2col(void *a1, void *a2, void *a3,
                     void *b1, void *b2, void *b3,
                     void *c, int *gdim, int *n) {

  if (*n <= 0)
    return;

  if (pso_opadd2col == nil)
    pso_opadd2col = get_mathops_pipeline("opadd2col_kernel");

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = mathops_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_opadd2col];

  [enc setBuffer:(__bridge id<MTLBuffer>)a1 offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)a2 offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)a3 offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)b1 offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)b2 offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)b3 offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)c  offset:0 atIndex:6];
  [enc setBytes:gdim length:sizeof(int) atIndex:7];
  [enc setBytes:n length:sizeof(int) atIndex:8];

  mathops_dispatch(cmdBuf, enc, *n);
}
