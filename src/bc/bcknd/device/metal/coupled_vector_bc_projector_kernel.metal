/*
 Copyright (c) 2026, The Neko Authors
 All rights reserved.
*/

#include <metal_stdlib>
using namespace metal;

kernel void coupled_vector_bc_projector_apply_kernel(
    device const int *mixed_msk [[ buffer(0) ]],
    device float *x [[ buffer(1) ]],
    device float *y [[ buffer(2) ]],
    device float *z [[ buffer(3) ]],
    device const int *constraint_n [[ buffer(4) ]],
    device const int *constraint_t1 [[ buffer(5) ]],
    device const int *constraint_t2 [[ buffer(6) ]],
    device const float *n [[ buffer(7) ]],
    device const float *t1 [[ buffer(8) ]],
    device const float *t2 [[ buffer(9) ]],
    constant int &m [[ buffer(10) ]],
    uint tid [[ thread_position_in_grid ]])
{
    const int i = (int)tid;
    if (i >= m) return;

    const int k = mixed_msk[i];
    const int off = 3 * i;
    const float u1 = x[k];
    const float u2 = y[k];
    const float u3 = z[k];

    float uloc_n = u1 * n[off] + u2 * n[off + 1] + u3 * n[off + 2];
    float uloc_t1 = u1 * t1[off] + u2 * t1[off + 1] + u3 * t1[off + 2];
    float uloc_t2 = u1 * t2[off] + u2 * t2[off + 1] + u3 * t2[off + 2];

    if (constraint_n[i] != 0) uloc_n = 0.0f;
    if (constraint_t1[i] != 0) uloc_t1 = 0.0f;
    if (constraint_t2[i] != 0) uloc_t2 = 0.0f;

    x[k] = uloc_n * n[off] + uloc_t1 * t1[off] + uloc_t2 * t2[off];
    y[k] = uloc_n * n[off + 1] + uloc_t1 * t1[off + 1] + uloc_t2 * t2[off + 1];
    z[k] = uloc_n * n[off + 2] + uloc_t1 * t1[off + 2] + uloc_t2 * t2[off + 2];
}
