/*
 Copyright (c) 2025, The Neko Authors
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

#ifndef __BC_METAL_UTILS_H__
#define __BC_METAL_UTILS_H__

#ifdef __APPLE__

#import <Metal/Metal.h>

extern id<MTLDevice>  neko_metal_device(void);
extern id<MTLLibrary> neko_metal_library(void);
extern void *glb_cmd_queue;

/* ------------------------------------------------------------------ */
/*  Pipeline state cache (shared across all BC files)                  */
/* ------------------------------------------------------------------ */

static NSMutableDictionary<NSString *, id<MTLComputePipelineState>> *bc_pipelines = nil;

static inline id<MTLComputePipelineState> bc_get_pipeline(NSString *name) {
    if (!bc_pipelines) {
        bc_pipelines = [NSMutableDictionary new];
    }
    id<MTLComputePipelineState> pso = bc_pipelines[name];
    if (!pso) {
        id<MTLDevice> dev = neko_metal_device();
        NSError *err = nil;
        id<MTLLibrary> lib = neko_metal_library();
        id<MTLFunction> fn = [lib newFunctionWithName:name];
        if (!fn) {
            NSLog(@"Metal BC: kernel '%@' not found in library", name);
            abort();
        }
        pso = [dev newComputePipelineStateWithFunction:fn error:&err];
        if (err) {
            NSLog(@"Metal BC: pipeline error for '%@': %@", name, err);
            abort();
        }
        bc_pipelines[name] = pso;
    }
    return pso;
}

/* ------------------------------------------------------------------ */
/*  Dispatch helper                                                    */
/* ------------------------------------------------------------------ */

#define METAL_BC_NTHREADS 1024

static inline void bc_dispatch_1d(id<MTLCommandQueue> queue,
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
        NSUInteger tpg = n < METAL_BC_NTHREADS ? n : METAL_BC_NTHREADS;
        MTLSize tg = MTLSizeMake(tpg, 1, 1);
        [enc dispatchThreads:threads threadsPerThreadgroup:tg];
        [enc endEncoding];
        [cb commit];
        [cb waitUntilCompleted];
    }
}

#endif /* __APPLE__ */
#endif /* __BC_METAL_UTILS_H__ */
