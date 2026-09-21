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

void metal_symmetry_aligned_apply_vector(void *xmsk, void *ymsk, void *zmsk,
                                 void *x, void *y, void *z,
                                 int *m, int *n, int *l, void *strm) {
    int max_len = *m;
    if (*n > max_len) max_len = *n;
    if (*l > max_len) max_len = *l;
    if (max_len < 1) return;

    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    int fm = *m, fn = *n, fl = *l;
    bc_dispatch_1d(q, bc_get_pipeline(@"symmetry_aligned_apply_vector_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(xmsk) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(ymsk) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(zmsk) offset:0 atIndex:2];
            [enc setBuffer:(__bridge id<MTLBuffer>)(x)    offset:0 atIndex:3];
            [enc setBuffer:(__bridge id<MTLBuffer>)(y)    offset:0 atIndex:4];
            [enc setBuffer:(__bridge id<MTLBuffer>)(z)    offset:0 atIndex:5];
            [enc setBytes:&fm length:sizeof(int) atIndex:6];
            [enc setBytes:&fn length:sizeof(int) atIndex:7];
            [enc setBytes:&fl length:sizeof(int) atIndex:8];
        }, (NSUInteger)max_len);
}

#endif /* __APPLE__ */
