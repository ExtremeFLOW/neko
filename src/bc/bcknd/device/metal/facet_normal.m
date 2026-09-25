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

void metal_facet_normal_apply_surfvec(void *msk, void *facet,
                                      void *x, void *y, void *z,
                                      void *u, void *v, void *w,
                                      void *nx, void *ny, void *nz,
                                      void *area, int *lx, int *m) {
    if (*m < 2) return;
    /* No stream parameter — use the global command queue. */
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(glb_cmd_queue);
    int flx = *lx, fm = *m;
    bc_dispatch_1d(q, bc_get_pipeline(@"facet_normal_apply_surfvec_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(msk)   offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(facet) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(x)     offset:0 atIndex:2];
            [enc setBuffer:(__bridge id<MTLBuffer>)(y)     offset:0 atIndex:3];
            [enc setBuffer:(__bridge id<MTLBuffer>)(z)     offset:0 atIndex:4];
            [enc setBuffer:(__bridge id<MTLBuffer>)(u)     offset:0 atIndex:5];
            [enc setBuffer:(__bridge id<MTLBuffer>)(v)     offset:0 atIndex:6];
            [enc setBuffer:(__bridge id<MTLBuffer>)(w)     offset:0 atIndex:7];
            [enc setBuffer:(__bridge id<MTLBuffer>)(nx)    offset:0 atIndex:8];
            [enc setBuffer:(__bridge id<MTLBuffer>)(ny)    offset:0 atIndex:9];
            [enc setBuffer:(__bridge id<MTLBuffer>)(nz)    offset:0 atIndex:10];
            [enc setBuffer:(__bridge id<MTLBuffer>)(area)  offset:0 atIndex:11];
            [enc setBytes:&flx length:sizeof(int) atIndex:12];
            [enc setBytes:&fm  length:sizeof(int) atIndex:13];
        }, (NSUInteger)(*m - 1));
}

#endif /* __APPLE__ */
