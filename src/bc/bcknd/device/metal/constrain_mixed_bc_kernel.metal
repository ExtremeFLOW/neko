/*
 Copyright (c) 2026, The Neko Authors
 All rights reserved.
*/

#include <metal_stdlib>
using namespace metal;

/*
 * Constrain a vector field at mixed-boundary nodes in the local (n, t1, t2)
 * basis. Ported from the CUDA/HIP kernels; `mixed_msk` is a resolved 0-based
 * index list of length m, and n/t1/t2 are nodewise basis vectors stored
 * interleaved at offset 3*i.
 */

kernel void constrain_mixed_bc_zero_kernel(
    device const int *mixed_msk [[ buffer(0) ]],
    device float *x [[ buffer(1) ]],
    device float *y [[ buffer(2) ]],
    device float *z [[ buffer(3) ]],
    constant int &constraint_n [[ buffer(4) ]],
    constant int &constraint_t1 [[ buffer(5) ]],
    constant int &constraint_t2 [[ buffer(6) ]],
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

    if (constraint_n != 0) uloc_n = 0.0f;
    if (constraint_t1 != 0) uloc_t1 = 0.0f;
    if (constraint_t2 != 0) uloc_t2 = 0.0f;

    x[k] = uloc_n * n[off] + uloc_t1 * t1[off] + uloc_t2 * t2[off];
    y[k] = uloc_n * n[off + 1] + uloc_t1 * t1[off + 1] + uloc_t2 * t2[off + 1];
    z[k] = uloc_n * n[off + 2] + uloc_t1 * t1[off + 2] + uloc_t2 * t2[off + 2];
}

kernel void constrain_mixed_bc_set_kernel(
    device const int *mixed_msk [[ buffer(0) ]],
    device float *x [[ buffer(1) ]],
    device float *y [[ buffer(2) ]],
    device float *z [[ buffer(3) ]],
    constant int &constraint_n [[ buffer(4) ]],
    constant int &constraint_t1 [[ buffer(5) ]],
    constant int &constraint_t2 [[ buffer(6) ]],
    device const float *n [[ buffer(7) ]],
    device const float *t1 [[ buffer(8) ]],
    device const float *t2 [[ buffer(9) ]],
    device const float *values_n [[ buffer(10) ]],
    device const float *values_t1 [[ buffer(11) ]],
    device const float *values_t2 [[ buffer(12) ]],
    constant int &m [[ buffer(13) ]],
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

    const float glb_1 = values_n[i];
    const float glb_2 = values_t1[i];
    const float glb_3 = values_t2[i];
    const float glb_n = glb_1 * n[off] + glb_2 * n[off + 1]
                      + glb_3 * n[off + 2];
    const float glb_t1 = glb_1 * t1[off] + glb_2 * t1[off + 1]
                       + glb_3 * t1[off + 2];
    const float glb_t2 = glb_1 * t2[off] + glb_2 * t2[off + 1]
                       + glb_3 * t2[off + 2];

    if (constraint_n != 0) uloc_n = glb_n;
    if (constraint_t1 != 0) uloc_t1 = glb_t1;
    if (constraint_t2 != 0) uloc_t2 = glb_t2;

    x[k] = uloc_n * n[off] + uloc_t1 * t1[off] + uloc_t2 * t2[off];
    y[k] = uloc_n * n[off + 1] + uloc_t1 * t1[off + 1] + uloc_t2 * t2[off + 1];
    z[k] = uloc_n * n[off + 2] + uloc_t1 * t1[off + 2] + uloc_t2 * t2[off + 2];
}

kernel void constrain_mixed_bc_set_const_kernel(
    device const int *mixed_msk [[ buffer(0) ]],
    device float *x [[ buffer(1) ]],
    device float *y [[ buffer(2) ]],
    device float *z [[ buffer(3) ]],
    constant int &constraint_n [[ buffer(4) ]],
    constant int &constraint_t1 [[ buffer(5) ]],
    constant int &constraint_t2 [[ buffer(6) ]],
    device const float *n [[ buffer(7) ]],
    device const float *t1 [[ buffer(8) ]],
    device const float *t2 [[ buffer(9) ]],
    constant float &value_n [[ buffer(10) ]],
    constant float &value_t1 [[ buffer(11) ]],
    constant float &value_t2 [[ buffer(12) ]],
    constant int &m [[ buffer(13) ]],
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

    const float glb_n = value_n * n[off] + value_t1 * n[off + 1]
                      + value_t2 * n[off + 2];
    const float glb_t1 = value_n * t1[off] + value_t1 * t1[off + 1]
                       + value_t2 * t1[off + 2];
    const float glb_t2 = value_n * t2[off] + value_t1 * t2[off + 1]
                       + value_t2 * t2[off + 2];

    if (constraint_n != 0) uloc_n = glb_n;
    if (constraint_t1 != 0) uloc_t1 = glb_t1;
    if (constraint_t2 != 0) uloc_t2 = glb_t2;

    x[k] = uloc_n * n[off] + uloc_t1 * t1[off] + uloc_t2 * t2[off];
    y[k] = uloc_n * n[off + 1] + uloc_t1 * t1[off + 1] + uloc_t2 * t2[off + 1];
    z[k] = uloc_n * n[off + 2] + uloc_t1 * t1[off + 2] + uloc_t2 * t2[off + 2];
}
