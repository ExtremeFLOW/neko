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

#ifdef __APPLE__
#include "bc_utils.h"

void metal_zero_dirichlet_apply_scalar(void *msk, void *x,
                                       int *m, void *strm) {
    if (*m < 2) return;
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    int fm = *m;
    bc_dispatch_1d(q, bc_get_pipeline(@"zero_dirichlet_apply_scalar_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(msk) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(x)   offset:0 atIndex:1];
            [enc setBytes:&fm length:sizeof(int) atIndex:2];
        }, (NSUInteger)(*m - 1));
}

void metal_zero_dirichlet_apply_vector(void *msk, void *x, void *y, void *z,
                                       int *m, void *strm) {
    if (*m < 2) return;
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    int fm = *m;
    bc_dispatch_1d(q, bc_get_pipeline(@"zero_dirichlet_apply_vector_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(msk) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(x)   offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(y)   offset:0 atIndex:2];
            [enc setBuffer:(__bridge id<MTLBuffer>)(z)   offset:0 atIndex:3];
            [enc setBytes:&fm length:sizeof(int) atIndex:4];
        }, (NSUInteger)(*m - 1));
}

#endif /* __APPLE__ */
