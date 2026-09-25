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

#ifndef __DEVICE_METAL_KERNEL_UTILS_H__
#define __DEVICE_METAL_KERNEL_UTILS_H__

#ifdef __APPLE__

#import <Metal/Metal.h>

extern id<MTLDevice>  neko_metal_device(void);
extern id<MTLLibrary> neko_metal_library(void);
extern void *glb_cmd_queue;

/* ------------------------------------------------------------------ */
/*  Pipeline state cache (one per translation unit)                    */
/* ------------------------------------------------------------------ */

static NSMutableDictionary<NSString *, id<MTLComputePipelineState>>
  *neko_metal_pipelines = nil;

/**
 * Look up (and cache) the compute pipeline state of a named kernel.
 */
static inline id<MTLComputePipelineState>
neko_metal_pipeline(NSString *name) {
    if (!neko_metal_pipelines) {
        neko_metal_pipelines = [NSMutableDictionary new];
    }
    id<MTLComputePipelineState> pso = neko_metal_pipelines[name];
    if (!pso) {
        id<MTLDevice> dev = neko_metal_device();
        NSError *err = nil;
        id<MTLLibrary> lib = neko_metal_library();
        id<MTLFunction> fn = [lib newFunctionWithName:name];
        if (!fn) {
            NSLog(@"Metal: kernel '%@' not found in library", name);
            abort();
        }
        pso = [dev newComputePipelineStateWithFunction:fn error:&err];
        if (err) {
            NSLog(@"Metal: pipeline error for '%@': %@", name, err);
            abort();
        }
        neko_metal_pipelines[name] = pso;
    }
    return pso;
}

/* ------------------------------------------------------------------ */
/*  Dispatch helper                                                    */
/* ------------------------------------------------------------------ */

#define NEKO_METAL_NTHREADS 1024

/**
 * Encode and run a one-dimensional kernel over @a n threads on @a queue,
 * waiting for completion.
 */
static inline void
neko_metal_dispatch_1d_queue(id<MTLCommandQueue> queue,
                             id<MTLComputePipelineState> pso,
                             void (^encode)(id<MTLComputeCommandEncoder>),
                             NSUInteger n) {
    if (n == 0) return;
    @autoreleasepool {
        id<MTLCommandBuffer> cb = [queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
        [enc setComputePipelineState:pso];
        encode(enc);
        MTLSize threads = MTLSizeMake(n, 1, 1);
        NSUInteger tpg = n < NEKO_METAL_NTHREADS ? n : NEKO_METAL_NTHREADS;
        MTLSize tg = MTLSizeMake(tpg, 1, 1);
        [enc dispatchThreads:threads threadsPerThreadgroup:tg];
        [enc endEncoding];
        [cb commit];
        [cb waitUntilCompleted];
    }
}

/**
 * Encode and run a kernel as @a ngroups threadgroups of @a tpg threads
 * on the global command queue, waiting for completion. Use this for
 * kernels that treat a threadgroup as one element.
 */
static inline void
neko_metal_dispatch_groups(id<MTLComputePipelineState> pso,
                           void (^encode)(id<MTLComputeCommandEncoder>),
                           NSUInteger ngroups, NSUInteger tpg) {
    if (ngroups == 0) return;
    @autoreleasepool {
        id<MTLCommandQueue> queue =
          (__bridge id<MTLCommandQueue>) glb_cmd_queue;
        id<MTLCommandBuffer> cb = [queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
        [enc setComputePipelineState:pso];
        encode(enc);
        [enc dispatchThreadgroups:MTLSizeMake(ngroups, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(tpg, 1, 1)];
        [enc endEncoding];
        [cb commit];
        [cb waitUntilCompleted];
    }
}

/**
 * Encode and run a one-dimensional kernel over @a n threads on the
 * global command queue, waiting for completion.
 */
static inline void
neko_metal_dispatch_1d(id<MTLComputePipelineState> pso,
                       void (^encode)(id<MTLComputeCommandEncoder>),
                       NSUInteger n) {
    neko_metal_dispatch_1d_queue(
      (__bridge id<MTLCommandQueue>) glb_cmd_queue, pso, encode, n);
}

#endif /* __APPLE__ */
#endif /* __DEVICE_METAL_KERNEL_UTILS_H__ */
