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

#import <Metal/Metal.h>
#include <device/device_config.h>

extern id<MTLDevice> neko_metal_device(void);
extern id<MTLLibrary> neko_metal_library(void);
extern void *glb_cmd_queue;

// ============================================================================
// Pipeline state management
// ============================================================================

static NSMutableDictionary<NSString *, id<MTLComputePipelineState>> *pipelines = nil;

static id<MTLComputePipelineState> get_pipeline(NSString *name) {
    if (!pipelines) {
        pipelines = [NSMutableDictionary new];
    }
    id<MTLComputePipelineState> pso = pipelines[name];
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
        pipelines[name] = pso;
    }
    return pso;
}

// ============================================================================
// Reduction buffer (shared memory, readable by CPU after GPU completes)
// ============================================================================

static id<MTLBuffer> reduction_buf = nil;
static int reduction_buf_cap = 0;

static id<MTLBuffer> get_reduction_buf(int nblocks) {
    if (nblocks > reduction_buf_cap) {
        int cap = nblocks > 1024 ? nblocks : 1024;
        reduction_buf = [neko_metal_device()
            newBufferWithLength:cap * sizeof(float)
                       options:MTLResourceStorageModeShared];
        reduction_buf_cap = cap;
    }
    return reduction_buf;
}

// ============================================================================
// Helper: dispatch a simple element-wise kernel
// ============================================================================

#define METAL_NTHREADS 1024

static inline void dispatch_simple(id<MTLCommandQueue> queue,
                                   id<MTLComputePipelineState> pso,
                                   void (^encode)(id<MTLComputeCommandEncoder>),
                                   NSUInteger n) {
    @autoreleasepool {
        id<MTLCommandBuffer> cb = [queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
        [enc setComputePipelineState:pso];
        encode(enc);
        MTLSize threads = MTLSizeMake(n, 1, 1);
        NSUInteger tpg = n < METAL_NTHREADS ? n : METAL_NTHREADS;
        MTLSize tg = MTLSizeMake(tpg, 1, 1);
        [enc dispatchThreads:threads threadsPerThreadgroup:tg];
        [enc endEncoding];
        [cb commit];
        [cb waitUntilCompleted];
    }
}

// ============================================================================
// Copy / fill / zero
// ============================================================================

void metal_copy(void *a_ptr, void *b_ptr, int *n, void *strm) {
    @autoreleasepool {
        id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
        id<MTLCommandBuffer> cb = [q commandBuffer];
        id<MTLBlitCommandEncoder> blit = [cb blitCommandEncoder];
        NSUInteger sz = (NSUInteger)(*n) * sizeof(float);
        [blit copyFromBuffer:(__bridge id<MTLBuffer>)(b_ptr) sourceOffset:0
                    toBuffer:(__bridge id<MTLBuffer>)(a_ptr) destinationOffset:0
                        size:sz];
        [blit endEncoding];
        [cb commit];
        [cb waitUntilCompleted];
    }
}

void metal_rzero(void *a_ptr, int *n, void *strm) {
    @autoreleasepool {
        id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
        id<MTLCommandBuffer> cb = [q commandBuffer];
        id<MTLBlitCommandEncoder> blit = [cb blitCommandEncoder];
        NSUInteger sz = (NSUInteger)(*n) * sizeof(float);
        [blit fillBuffer:(__bridge id<MTLBuffer>)(a_ptr)
                   range:NSMakeRange(0, sz)
                   value:0];
        [blit endEncoding];
        [cb commit];
        [cb waitUntilCompleted];
    }
}

void metal_cfill(void *a_ptr, real *c, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    float fc = (float)*c;
    dispatch_simple(q, get_pipeline(@"cfill_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBytes:&fc length:sizeof(float) atIndex:1];
            [enc setBytes:n length:sizeof(int) atIndex:2];
        }, (NSUInteger)*n);
}

void metal_cfill_mask(void *a_ptr, real *c, int *n, void *mask_ptr,
                      int *n_mask, void *strm) {
    if (*n_mask < 1) return;
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    float fc = (float)*c;
    dispatch_simple(q, get_pipeline(@"cfill_mask_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBytes:&fc length:sizeof(float) atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(mask_ptr) offset:0 atIndex:2];
            [enc setBytes:n_mask length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n_mask);
}

// ============================================================================
// Masked copy variants
// ============================================================================

void metal_masked_copy_0(void *a_ptr, void *b_ptr, void *mask_ptr,
                       int *n, int *n_mask, void *strm) {
    if (*n_mask < 1) return;
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"masked_copy_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(mask_ptr) offset:0 atIndex:2];
            [enc setBytes:n_mask length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n_mask);
}

void metal_masked_gather_copy(void *a_ptr, void *b_ptr, void *mask_ptr,
                              int *n, int *n_mask, void *strm) {
    if (*n_mask < 1) return;
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"masked_gather_copy_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(mask_ptr) offset:0 atIndex:2];
            [enc setBytes:n_mask length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n_mask);
}

void metal_masked_gather_copy_aligned(void *a_ptr, void *b_ptr, void *mask_ptr,
                                      int *n, int *n_mask, void *strm) {
    if (*n_mask < 1) return;
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"masked_gather_copy_aligned_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(mask_ptr) offset:0 atIndex:2];
            [enc setBytes:n_mask length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n_mask);
}

void metal_masked_scatter_copy(void *a_ptr, void *b_ptr, void *mask_ptr,
                               int *n, int *n_mask, void *strm) {
    if (*n_mask < 1) return;
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"masked_scatter_copy_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(mask_ptr) offset:0 atIndex:2];
            [enc setBytes:n_mask length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n_mask);
}

void metal_masked_atomic_reduction(void *a_ptr, void *b_ptr, void *mask_ptr,
                                   int *n, int *m, void *strm) {
    if (*m < 1) return;
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"masked_atomic_reduction_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(mask_ptr) offset:0 atIndex:2];
            [enc setBytes:m length:sizeof(int) atIndex:3];
        }, (NSUInteger)*m);
}

void metal_masked_scatter_copy_aligned(void *a_ptr, void *b_ptr, void *mask_ptr,
                                       int *n, int *n_mask, void *strm) {
    if (*n_mask < 1) return;
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"masked_scatter_copy_aligned_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(mask_ptr) offset:0 atIndex:2];
            [enc setBytes:n_mask length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n_mask);
}

void metal_face_masked_gather_copy(void *a_ptr, void *b_ptr, void *mask_ptr,
                                   void *facet_ptr, int *n1, int *n2,
                                   int *lx, int *ly, int *lz, int *n_mask,
                                   void *strm) {
    if (*n_mask < 1) return;
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"face_masked_gather_copy_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(mask_ptr) offset:0 atIndex:2];
            [enc setBuffer:(__bridge id<MTLBuffer>)(facet_ptr) offset:0 atIndex:3];
            [enc setBytes:n1 length:sizeof(int) atIndex:4];
            [enc setBytes:n2 length:sizeof(int) atIndex:5];
            [enc setBytes:lx length:sizeof(int) atIndex:6];
            [enc setBytes:ly length:sizeof(int) atIndex:7];
            [enc setBytes:lz length:sizeof(int) atIndex:8];
            [enc setBytes:n_mask length:sizeof(int) atIndex:9];
        }, (NSUInteger)*n_mask);
}

void metal_cwrap(void *a_ptr, real *min_val, real *max_val, int *n,
                 void *strm) {
    if (*n < 1) return;
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    float fmin = (float)*min_val, fmax = (float)*max_val;
    dispatch_simple(q, get_pipeline(@"cwrap_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBytes:&fmin length:sizeof(float) atIndex:1];
            [enc setBytes:&fmax length:sizeof(float) atIndex:2];
            [enc setBytes:n length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n);
}

// ============================================================================
// Scalar-vector operations
// ============================================================================

void metal_cmult(void *a_ptr, real *c, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    float fc = (float)*c;
    dispatch_simple(q, get_pipeline(@"cmult_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBytes:&fc length:sizeof(float) atIndex:1];
            [enc setBytes:n length:sizeof(int) atIndex:2];
        }, (NSUInteger)*n);
}

void metal_cmult2(void *a_ptr, void *b_ptr, real *c, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    float fc = (float)*c;
    dispatch_simple(q, get_pipeline(@"cmult2_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBytes:&fc length:sizeof(float) atIndex:2];
            [enc setBytes:n length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n);
}

void metal_cdiv(void *a_ptr, real *c, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    float fc = (float)*c;
    dispatch_simple(q, get_pipeline(@"cdiv_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBytes:&fc length:sizeof(float) atIndex:1];
            [enc setBytes:n length:sizeof(int) atIndex:2];
        }, (NSUInteger)*n);
}

void metal_cdiv2(void *a_ptr, void *b_ptr, real *c, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    float fc = (float)*c;
    dispatch_simple(q, get_pipeline(@"cdiv2_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBytes:&fc length:sizeof(float) atIndex:2];
            [enc setBytes:n length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n);
}

void metal_radd(void *a_ptr, real *c, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    float fc = (float)*c;
    dispatch_simple(q, get_pipeline(@"cadd_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBytes:&fc length:sizeof(float) atIndex:1];
            [enc setBytes:n length:sizeof(int) atIndex:2];
        }, (NSUInteger)*n);
}

void metal_cadd2(void *a_ptr, void *b_ptr, real *c, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    float fc = (float)*c;
    dispatch_simple(q, get_pipeline(@"cadd2_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBytes:&fc length:sizeof(float) atIndex:2];
            [enc setBytes:n length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n);
}

// ============================================================================
// Vector-vector add operations
// ============================================================================

void metal_add2(void *a_ptr, void *b_ptr, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"add2_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBytes:n length:sizeof(int) atIndex:2];
        }, (NSUInteger)*n);
}

void metal_add3(void *a_ptr, void *b_ptr, void *c_ptr, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"add3_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(c_ptr) offset:0 atIndex:2];
            [enc setBytes:n length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n);
}

void metal_add4(void *a_ptr, void *b_ptr, void *c_ptr, void *d_ptr,
                int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"add4_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(c_ptr) offset:0 atIndex:2];
            [enc setBuffer:(__bridge id<MTLBuffer>)(d_ptr) offset:0 atIndex:3];
            [enc setBytes:n length:sizeof(int) atIndex:4];
        }, (NSUInteger)*n);
}

void metal_add2s1(void *a_ptr, void *b_ptr, real *c1, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    float fc1 = (float)*c1;
    dispatch_simple(q, get_pipeline(@"add2s1_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBytes:&fc1 length:sizeof(float) atIndex:2];
            [enc setBytes:n length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n);
}

void metal_add2s2(void *a_ptr, void *b_ptr, real *c1, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    float fc1 = (float)*c1;
    dispatch_simple(q, get_pipeline(@"add2s2_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBytes:&fc1 length:sizeof(float) atIndex:2];
            [enc setBytes:n length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n);
}

void metal_add2s2_many(void *y_ptr, void *x_d_d_ptr, void *a_ptr,
                       int *j, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    int p_cur = *j;
    NSUInteger nn = (NSUInteger)*n;

    /*
     * x_d_d_ptr is a Metal buffer whose contents are p_cur device pointers.
     * a_ptr is a Metal buffer whose contents are p_cur float alpha values.
     * Because Metal shared memory is CPU-visible, we read them on the host
     * and dispatch one kernel per column.
     */
    id<MTLBuffer> ptrs_buf = (__bridge id<MTLBuffer>)(x_d_d_ptr);
    void **ptrs = (void **)[ptrs_buf contents];
    id<MTLBuffer> alpha_buf = (__bridge id<MTLBuffer>)(a_ptr);
    float *alpha = (float *)[alpha_buf contents];

    for (int col = 0; col < p_cur; col++) {
        float a_col = alpha[col];
        void *p_col = ptrs[col];
        dispatch_simple(q, get_pipeline(@"add2s2_many_kernel"),
            ^(id<MTLComputeCommandEncoder> enc) {
                [enc setBuffer:(__bridge id<MTLBuffer>)(y_ptr)
                        offset:0 atIndex:0];
                [enc setBuffer:(__bridge id<MTLBuffer>)(p_col)
                        offset:0 atIndex:1];
                [enc setBytes:&a_col length:sizeof(float) atIndex:2];
                [enc setBytes:n length:sizeof(int) atIndex:3];
            }, nn);
    }
}

void metal_addsqr2s2(void *a_ptr, void *b_ptr, real *c1, int *n,
                     void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    float fc1 = (float)*c1;
    dispatch_simple(q, get_pipeline(@"addsqr2s2_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBytes:&fc1 length:sizeof(float) atIndex:2];
            [enc setBytes:n length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n);
}

void metal_add3s2(void *a_ptr, void *b_ptr, void *c_ptr,
                  real *c1, real *c2, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    float fc1 = (float)*c1, fc2 = (float)*c2;
    dispatch_simple(q, get_pipeline(@"add3s2_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(c_ptr) offset:0 atIndex:2];
            [enc setBytes:&fc1 length:sizeof(float) atIndex:3];
            [enc setBytes:&fc2 length:sizeof(float) atIndex:4];
            [enc setBytes:n length:sizeof(int) atIndex:5];
        }, (NSUInteger)*n);
}

void metal_add4s3(void *a_ptr, void *b_ptr, void *c_ptr, void *d_ptr,
                  real *c1, real *c2, real *c3, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    float fc1 = (float)*c1, fc2 = (float)*c2, fc3 = (float)*c3;
    dispatch_simple(q, get_pipeline(@"add4s3_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(c_ptr) offset:0 atIndex:2];
            [enc setBuffer:(__bridge id<MTLBuffer>)(d_ptr) offset:0 atIndex:3];
            [enc setBytes:&fc1 length:sizeof(float) atIndex:4];
            [enc setBytes:&fc2 length:sizeof(float) atIndex:5];
            [enc setBytes:&fc3 length:sizeof(float) atIndex:6];
            [enc setBytes:n length:sizeof(int) atIndex:7];
        }, (NSUInteger)*n);
}

void metal_add5s4(void *a_ptr, void *b_ptr, void *c_ptr, void *d_ptr,
                  void *e_ptr, real *c1, real *c2, real *c3,
                  real *c4, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    float fc1 = (float)*c1, fc2 = (float)*c2;
    float fc3 = (float)*c3, fc4 = (float)*c4;
    dispatch_simple(q, get_pipeline(@"add5s4_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(c_ptr) offset:0 atIndex:2];
            [enc setBuffer:(__bridge id<MTLBuffer>)(d_ptr) offset:0 atIndex:3];
            [enc setBuffer:(__bridge id<MTLBuffer>)(e_ptr) offset:0 atIndex:4];
            [enc setBytes:&fc1 length:sizeof(float) atIndex:5];
            [enc setBytes:&fc2 length:sizeof(float) atIndex:6];
            [enc setBytes:&fc3 length:sizeof(float) atIndex:7];
            [enc setBytes:&fc4 length:sizeof(float) atIndex:8];
            [enc setBytes:n length:sizeof(int) atIndex:9];
        }, (NSUInteger)*n);
}

// ============================================================================
// Inverse / column operations
// ============================================================================

void metal_invcol1(void *a_ptr, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"invcol1_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBytes:n length:sizeof(int) atIndex:1];
        }, (NSUInteger)*n);
}

void metal_invcol2(void *a_ptr, void *b_ptr, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"invcol2_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBytes:n length:sizeof(int) atIndex:2];
        }, (NSUInteger)*n);
}

void metal_invcol3(void *a_ptr, void *b_ptr, void *c_ptr, int *n,
                   void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"invcol3_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(c_ptr) offset:0 atIndex:2];
            [enc setBytes:n length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n);
}

void metal_col2(void *a_ptr, void *b_ptr, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"col2_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBytes:n length:sizeof(int) atIndex:2];
        }, (NSUInteger)*n);
}

void metal_col3(void *a_ptr, void *b_ptr, void *c_ptr, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"col3_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(c_ptr) offset:0 atIndex:2];
            [enc setBytes:n length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n);
}

void metal_subcol3(void *a_ptr, void *b_ptr, void *c_ptr, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"subcol3_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(c_ptr) offset:0 atIndex:2];
            [enc setBytes:n length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n);
}

void metal_sub2(void *a_ptr, void *b_ptr, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"sub2_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBytes:n length:sizeof(int) atIndex:2];
        }, (NSUInteger)*n);
}

void metal_sub3(void *a_ptr, void *b_ptr, void *c_ptr, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"sub3_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(c_ptr) offset:0 atIndex:2];
            [enc setBytes:n length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n);
}

void metal_addcol3(void *a_ptr, void *b_ptr, void *c_ptr, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"addcol3_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(c_ptr) offset:0 atIndex:2];
            [enc setBytes:n length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n);
}

void metal_addcol4(void *a_ptr, void *b_ptr, void *c_ptr, void *d_ptr,
                   int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"addcol4_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(c_ptr) offset:0 atIndex:2];
            [enc setBuffer:(__bridge id<MTLBuffer>)(d_ptr) offset:0 atIndex:3];
            [enc setBytes:n length:sizeof(int) atIndex:4];
        }, (NSUInteger)*n);
}

void metal_addcol3s2(void *a_ptr, void *b_ptr, void *c_ptr, real *s,
                     int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    float fs = (float)*s;
    dispatch_simple(q, get_pipeline(@"addcol3s2_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(c_ptr) offset:0 atIndex:2];
            [enc setBytes:&fs length:sizeof(float) atIndex:3];
            [enc setBytes:n length:sizeof(int) atIndex:4];
        }, (NSUInteger)*n);
}

// ============================================================================
// Vector dot / cross
// ============================================================================

void metal_vdot3(void *dot_ptr, void *u1_ptr, void *u2_ptr, void *u3_ptr,
                 void *v1_ptr, void *v2_ptr, void *v3_ptr, int *n,
                 void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"vdot3_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(dot_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(u1_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(u2_ptr) offset:0 atIndex:2];
            [enc setBuffer:(__bridge id<MTLBuffer>)(u3_ptr) offset:0 atIndex:3];
            [enc setBuffer:(__bridge id<MTLBuffer>)(v1_ptr) offset:0 atIndex:4];
            [enc setBuffer:(__bridge id<MTLBuffer>)(v2_ptr) offset:0 atIndex:5];
            [enc setBuffer:(__bridge id<MTLBuffer>)(v3_ptr) offset:0 atIndex:6];
            [enc setBytes:n length:sizeof(int) atIndex:7];
        }, (NSUInteger)*n);
}

void metal_vcross(void *u1_ptr, void *u2_ptr, void *u3_ptr,
                  void *v1_ptr, void *v2_ptr, void *v3_ptr,
                  void *w1_ptr, void *w2_ptr, void *w3_ptr,
                  int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"vcross_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(u1_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(u2_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(u3_ptr) offset:0 atIndex:2];
            [enc setBuffer:(__bridge id<MTLBuffer>)(v1_ptr) offset:0 atIndex:3];
            [enc setBuffer:(__bridge id<MTLBuffer>)(v2_ptr) offset:0 atIndex:4];
            [enc setBuffer:(__bridge id<MTLBuffer>)(v3_ptr) offset:0 atIndex:5];
            [enc setBuffer:(__bridge id<MTLBuffer>)(w1_ptr) offset:0 atIndex:6];
            [enc setBuffer:(__bridge id<MTLBuffer>)(w2_ptr) offset:0 atIndex:7];
            [enc setBuffer:(__bridge id<MTLBuffer>)(w3_ptr) offset:0 atIndex:8];
            [enc setBytes:n length:sizeof(int) atIndex:9];
        }, (NSUInteger)*n);
}

void metal_absval(void *a_ptr, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"absval_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBytes:n length:sizeof(int) atIndex:1];
        }, (NSUInteger)*n);
}

void metal_iadd(void *a_ptr, int *c, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"iadd_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBytes:c length:sizeof(int) atIndex:1];
            [enc setBytes:n length:sizeof(int) atIndex:2];
        }, (NSUInteger)*n);
}

// ============================================================================
// Pointwise max/min
// ============================================================================

void metal_pwmax_vec2(void *a_ptr, void *b_ptr, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"pwmax_vec2_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBytes:n length:sizeof(int) atIndex:2];
        }, (NSUInteger)*n);
}

void metal_pwmax_vec3(void *a_ptr, void *b_ptr, void *c_ptr, int *n,
                      void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"pwmax_vec3_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(c_ptr) offset:0 atIndex:2];
            [enc setBytes:n length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n);
}

void metal_pwmax_sca2(void *a_ptr, real *c, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    float fc = (float)*c;
    dispatch_simple(q, get_pipeline(@"pwmax_sca2_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBytes:&fc length:sizeof(float) atIndex:1];
            [enc setBytes:n length:sizeof(int) atIndex:2];
        }, (NSUInteger)*n);
}

void metal_pwmax_sca3(void *a_ptr, void *b_ptr, real *c, int *n,
                      void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    float fc = (float)*c;
    dispatch_simple(q, get_pipeline(@"pwmax_sca3_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBytes:&fc length:sizeof(float) atIndex:2];
            [enc setBytes:n length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n);
}

void metal_pwmin_vec2(void *a_ptr, void *b_ptr, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"pwmin_vec2_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBytes:n length:sizeof(int) atIndex:2];
        }, (NSUInteger)*n);
}

void metal_pwmin_vec3(void *a_ptr, void *b_ptr, void *c_ptr, int *n,
                      void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    dispatch_simple(q, get_pipeline(@"pwmin_vec3_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(c_ptr) offset:0 atIndex:2];
            [enc setBytes:n length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n);
}

void metal_pwmin_sca2(void *a_ptr, real *c, int *n, void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    float fc = (float)*c;
    dispatch_simple(q, get_pipeline(@"pwmin_sca2_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBytes:&fc length:sizeof(float) atIndex:1];
            [enc setBytes:n length:sizeof(int) atIndex:2];
        }, (NSUInteger)*n);
}

void metal_pwmin_sca3(void *a_ptr, void *b_ptr, real *c, int *n,
                      void *strm) {
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    float fc = (float)*c;
    dispatch_simple(q, get_pipeline(@"pwmin_sca3_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:1];
            [enc setBytes:&fc length:sizeof(float) atIndex:2];
            [enc setBytes:n length:sizeof(int) atIndex:3];
        }, (NSUInteger)*n);
}

// ============================================================================
// Reduction operations
// ============================================================================

static float metal_reduce_sum(void *a_ptr, void *b_ptr, void *c_ptr,
                               int n, void *strm, NSString *kernel_name) {
    @autoreleasepool {
        id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
        NSUInteger nn = (NSUInteger)n;
        NSUInteger tg_size = nn < METAL_NTHREADS ? nn : METAL_NTHREADS;
        NSUInteger nblocks = (nn + tg_size - 1) / tg_size;

        id<MTLBuffer> buf = get_reduction_buf((int)nblocks);

        // Phase 1: per-threadgroup reduction
        {
            id<MTLComputePipelineState> pso = get_pipeline(kernel_name);
            id<MTLCommandBuffer> cb = [q commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
            [enc setComputePipelineState:pso];

            int buf_idx = 0;
            if (a_ptr) {
                [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:buf_idx++];
            }
            if (b_ptr) {
                [enc setBuffer:(__bridge id<MTLBuffer>)(b_ptr) offset:0 atIndex:buf_idx++];
            }
            if (c_ptr) {
                [enc setBuffer:(__bridge id<MTLBuffer>)(c_ptr) offset:0 atIndex:buf_idx++];
            }
            [enc setBuffer:buf offset:0 atIndex:buf_idx++];
            [enc setBytes:&n length:sizeof(int) atIndex:buf_idx];

            [enc dispatchThreadgroups:MTLSizeMake(nblocks, 1, 1)
                threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
            [enc endEncoding];
            [cb commit];
            [cb waitUntilCompleted];
        }

        // Phase 2: reduce partial sums
        if (nblocks > 1) {
            id<MTLComputePipelineState> pso = get_pipeline(@"reduce_kernel");
            id<MTLCommandBuffer> cb = [q commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
            [enc setComputePipelineState:pso];
            [enc setBuffer:buf offset:0 atIndex:0];
            int nb = (int)nblocks;
            [enc setBytes:&nb length:sizeof(int) atIndex:1];
            NSUInteger tg2 = nblocks < METAL_NTHREADS ? nblocks : METAL_NTHREADS;
            [enc dispatchThreadgroups:MTLSizeMake(1, 1, 1)
                threadsPerThreadgroup:MTLSizeMake(tg2, 1, 1)];
            [enc endEncoding];
            [cb commit];
            [cb waitUntilCompleted];
        }

        return ((float *)[buf contents])[0];
    }
}

real metal_vlsc3(void *u_ptr, void *v_ptr, void *w_ptr, int *n,
                   void *strm) {
    if (*n < 1) return 0.0;
    return (real)metal_reduce_sum(u_ptr, v_ptr, w_ptr, *n, strm,
                                    @"glsc3_kernel");
}

real metal_glsc3(void *a_ptr, void *b_ptr, void *c_ptr, int *n,
                   void *strm) {
    if (*n < 1) return 0.0;
    return (real)metal_reduce_sum(a_ptr, b_ptr, c_ptr, *n, strm,
                                    @"glsc3_kernel");
}

real metal_glsc2(void *a_ptr, void *b_ptr, int *n, void *strm) {
    if (*n < 1) return 0.0;
    return (real)metal_reduce_sum(a_ptr, b_ptr, NULL, *n, strm,
                                    @"glsc2_kernel");
}

real metal_glsubnorm2(void *a_ptr, void *b_ptr, int *n, void *strm) {
    if (*n < 1) return 0.0;
    return (real)metal_reduce_sum(a_ptr, b_ptr, NULL, *n, strm,
                                    @"glsubnorm2_kernel");
}

real metal_glsum(void *a_ptr, int *n, void *strm) {
    if (*n < 1) return 0.0;
    return (real)metal_reduce_sum(a_ptr, NULL, NULL, *n, strm,
                                    @"glsum_kernel");
}

static float metal_reduce_extremum(void *a_ptr, int n, void *strm,
                                    NSString *phase1_name,
                                    NSString *phase2_name) {
    @autoreleasepool {
        id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
        NSUInteger nn = (NSUInteger)n;
        NSUInteger tg_size = nn < METAL_NTHREADS ? nn : METAL_NTHREADS;
        NSUInteger nblocks = (nn + tg_size - 1) / tg_size;

        id<MTLBuffer> buf = get_reduction_buf((int)nblocks);

        // Phase 1
        {
            id<MTLComputePipelineState> pso = get_pipeline(phase1_name);
            id<MTLCommandBuffer> cb = [q commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
            [enc setComputePipelineState:pso];
            [enc setBuffer:(__bridge id<MTLBuffer>)(a_ptr) offset:0 atIndex:0];
            [enc setBuffer:buf offset:0 atIndex:1];
            [enc setBytes:&n length:sizeof(int) atIndex:2];
            [enc dispatchThreadgroups:MTLSizeMake(nblocks, 1, 1)
                threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
            [enc endEncoding];
            [cb commit];
            [cb waitUntilCompleted];
        }

        // Phase 2
        if (nblocks > 1) {
            id<MTLComputePipelineState> pso = get_pipeline(phase2_name);
            id<MTLCommandBuffer> cb = [q commandBuffer];
            id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
            [enc setComputePipelineState:pso];
            [enc setBuffer:buf offset:0 atIndex:0];
            int nb = (int)nblocks;
            [enc setBytes:&nb length:sizeof(int) atIndex:1];
            NSUInteger tg2 = nblocks < METAL_NTHREADS ? nblocks : METAL_NTHREADS;
            [enc dispatchThreadgroups:MTLSizeMake(1, 1, 1)
                threadsPerThreadgroup:MTLSizeMake(tg2, 1, 1)];
            [enc endEncoding];
            [cb commit];
            [cb waitUntilCompleted];
        }

        return ((float *)[buf contents])[0];
    }
}

real metal_glmax(void *a_ptr, real *ninf, int *n, void *strm) {
    if (*n < 1) return *ninf;
    return (real)metal_reduce_extremum(a_ptr, *n, strm,
                                         @"glmax_kernel",
                                         @"reduce_max_kernel");
}

real metal_glmin(void *a_ptr, real *pinf, int *n, void *strm) {
    if (*n < 1) return *pinf;
    return (real)metal_reduce_extremum(a_ptr, *n, strm,
                                         @"glmin_kernel",
                                         @"reduce_min_kernel");
}

void metal_glsc3_many(real *h, void *w_ptr, void *v_d_d_ptr,
                      void *mult_ptr, int *j, int *n, void *strm) {
    @autoreleasepool {
        id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
        NSUInteger nn = (NSUInteger)*n;
        NSUInteger tg_size = nn < METAL_NTHREADS ? nn : METAL_NTHREADS;
        NSUInteger nblocks = (nn + tg_size - 1) / tg_size;
        int ncols = *j;

        id<MTLBuffer> buf = get_reduction_buf((int)nblocks);

        /*
         * v_d_d_ptr is a Metal buffer of ncols device pointers.
         * Dispatch one reduction per column using shared memory access.
         */
        id<MTLBuffer> ptrs_buf = (__bridge id<MTLBuffer>)(v_d_d_ptr);
        void **ptrs = (void **)[ptrs_buf contents];

        id<MTLComputePipelineState> pso =
            get_pipeline(@"glsc3_many_kernel");

        for (int col = 0; col < ncols; col++) {
            void *v_col = ptrs[col];
            {
                id<MTLCommandBuffer> cb = [q commandBuffer];
                id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
                [enc setComputePipelineState:pso];
                [enc setBuffer:(__bridge id<MTLBuffer>)(w_ptr)
                        offset:0 atIndex:0];
                [enc setBuffer:(__bridge id<MTLBuffer>)(v_col)
                        offset:0 atIndex:1];
                [enc setBuffer:(__bridge id<MTLBuffer>)(mult_ptr)
                        offset:0 atIndex:2];
                [enc setBuffer:buf offset:0 atIndex:3];
                [enc setBytes:n length:sizeof(int) atIndex:4];
                [enc dispatchThreadgroups:MTLSizeMake(nblocks, 1, 1)
                    threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
                [enc endEncoding];
                [cb commit];
                [cb waitUntilCompleted];
            }

            float *buf_data = (float *)[buf contents];
            float sum = 0.0f;
            for (NSUInteger b = 0; b < nblocks; b++) {
                sum += buf_data[b];
            }
            h[col] = (real)sum;
        }
    }
}
