/*
 Copyright (c) 2026, The Neko Authors
 All rights reserved.
*/

#ifdef __APPLE__
#include "bc_utils.h"

void metal_constrain_mixed_bc_zero(
    void *mixed_msk, void *x, void *y, void *z,
    int *constraint_n, int *constraint_t1, int *constraint_t2,
    void *n, void *t1, void *t2, int *m, void *strm) {
    if (*m < 1) return;
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    int cn = *constraint_n, ct1 = *constraint_t1, ct2 = *constraint_t2;
    int fm = *m;
    bc_dispatch_1d(q, bc_get_pipeline(@"constrain_mixed_bc_zero_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(mixed_msk) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(x)  offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(y)  offset:0 atIndex:2];
            [enc setBuffer:(__bridge id<MTLBuffer>)(z)  offset:0 atIndex:3];
            [enc setBytes:&cn  length:sizeof(int) atIndex:4];
            [enc setBytes:&ct1 length:sizeof(int) atIndex:5];
            [enc setBytes:&ct2 length:sizeof(int) atIndex:6];
            [enc setBuffer:(__bridge id<MTLBuffer>)(n)  offset:0 atIndex:7];
            [enc setBuffer:(__bridge id<MTLBuffer>)(t1) offset:0 atIndex:8];
            [enc setBuffer:(__bridge id<MTLBuffer>)(t2) offset:0 atIndex:9];
            [enc setBytes:&fm length:sizeof(int) atIndex:10];
        }, (NSUInteger)*m);
}

void metal_constrain_mixed_bc_set(
    void *mixed_msk, void *x, void *y, void *z,
    int *constraint_n, int *constraint_t1, int *constraint_t2,
    void *n, void *t1, void *t2,
    void *values_n, void *values_t1, void *values_t2, int *m, void *strm) {
    if (*m < 1) return;
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    int cn = *constraint_n, ct1 = *constraint_t1, ct2 = *constraint_t2;
    int fm = *m;
    bc_dispatch_1d(q, bc_get_pipeline(@"constrain_mixed_bc_set_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(mixed_msk) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(x)  offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(y)  offset:0 atIndex:2];
            [enc setBuffer:(__bridge id<MTLBuffer>)(z)  offset:0 atIndex:3];
            [enc setBytes:&cn  length:sizeof(int) atIndex:4];
            [enc setBytes:&ct1 length:sizeof(int) atIndex:5];
            [enc setBytes:&ct2 length:sizeof(int) atIndex:6];
            [enc setBuffer:(__bridge id<MTLBuffer>)(n)  offset:0 atIndex:7];
            [enc setBuffer:(__bridge id<MTLBuffer>)(t1) offset:0 atIndex:8];
            [enc setBuffer:(__bridge id<MTLBuffer>)(t2) offset:0 atIndex:9];
            [enc setBuffer:(__bridge id<MTLBuffer>)(values_n)  offset:0 atIndex:10];
            [enc setBuffer:(__bridge id<MTLBuffer>)(values_t1) offset:0 atIndex:11];
            [enc setBuffer:(__bridge id<MTLBuffer>)(values_t2) offset:0 atIndex:12];
            [enc setBytes:&fm length:sizeof(int) atIndex:13];
        }, (NSUInteger)*m);
}

void metal_constrain_mixed_bc_set_const(
    void *mixed_msk, void *x, void *y, void *z,
    int *constraint_n, int *constraint_t1, int *constraint_t2,
    void *n, void *t1, void *t2,
    float *value_n, float *value_t1, float *value_t2, int *m, void *strm) {
    if (*m < 1) return;
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)(strm);
    int cn = *constraint_n, ct1 = *constraint_t1, ct2 = *constraint_t2;
    float vn = *value_n, vt1 = *value_t1, vt2 = *value_t2;
    int fm = *m;
    bc_dispatch_1d(q, bc_get_pipeline(@"constrain_mixed_bc_set_const_kernel"),
        ^(id<MTLComputeCommandEncoder> enc) {
            [enc setBuffer:(__bridge id<MTLBuffer>)(mixed_msk) offset:0 atIndex:0];
            [enc setBuffer:(__bridge id<MTLBuffer>)(x)  offset:0 atIndex:1];
            [enc setBuffer:(__bridge id<MTLBuffer>)(y)  offset:0 atIndex:2];
            [enc setBuffer:(__bridge id<MTLBuffer>)(z)  offset:0 atIndex:3];
            [enc setBytes:&cn  length:sizeof(int) atIndex:4];
            [enc setBytes:&ct1 length:sizeof(int) atIndex:5];
            [enc setBytes:&ct2 length:sizeof(int) atIndex:6];
            [enc setBuffer:(__bridge id<MTLBuffer>)(n)  offset:0 atIndex:7];
            [enc setBuffer:(__bridge id<MTLBuffer>)(t1) offset:0 atIndex:8];
            [enc setBuffer:(__bridge id<MTLBuffer>)(t2) offset:0 atIndex:9];
            [enc setBytes:&vn  length:sizeof(float) atIndex:10];
            [enc setBytes:&vt1 length:sizeof(float) atIndex:11];
            [enc setBytes:&vt2 length:sizeof(float) atIndex:12];
            [enc setBytes:&fm  length:sizeof(int)   atIndex:13];
        }, (NSUInteger)*m);
}

#endif /* __APPLE__ */
