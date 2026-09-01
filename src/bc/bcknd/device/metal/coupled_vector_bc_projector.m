/*
 Copyright (c) 2026, The Neko Authors
 All rights reserved.
*/

#ifdef __APPLE__
#include "bc_utils.h"

void metal_coupled_vector_bc_projector_apply(
    void *mixed_msk, void *x, void *y, void *z, void *constraint_n,
    void *constraint_t1, void *constraint_t2, void *n, void *t1, void *t2,
    int *m, void *strm) {
    if (*m <= 0) return;

    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    int fm = *m;
    bc_dispatch_1d(q,
        bc_get_pipeline(@"coupled_vector_bc_projector_apply_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(mixed_msk) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(x) offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(y) offset:0 atIndex:2];
            [enc setBuffer:(__bridge id<MTLBuffer>)(z) offset:0 atIndex:3];
            [enc setBuffer:(__bridge id<MTLBuffer>)(constraint_n) offset:0 atIndex:4];
            [enc setBuffer:(__bridge id<MTLBuffer>)(constraint_t1) offset:0 atIndex:5];
            [enc setBuffer:(__bridge id<MTLBuffer>)(constraint_t2) offset:0 atIndex:6];
            [enc setBuffer:(__bridge id<MTLBuffer>)(n) offset:0 atIndex:7];
            [enc setBuffer:(__bridge id<MTLBuffer>)(t1) offset:0 atIndex:8];
            [enc setBuffer:(__bridge id<MTLBuffer>)(t2) offset:0 atIndex:9];
            [enc setBytes:&fm length:sizeof(int) atIndex:10];
        }, (NSUInteger)*m);
}

#endif /* __APPLE__ */
