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

#include <metal_stdlib>
using namespace metal;

// ============================================================================
// Element-wise kernels
// ============================================================================

kernel void cmult_kernel(device float *a[[ buffer(0) ]],
                         constant float &c[[ buffer(1) ]],
                         constant int &n[[ buffer(2) ]],
                         uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = c * a[idx];
}

kernel void cmult2_kernel(device float *a[[ buffer(0) ]],
                          device const float *b[[ buffer(1) ]],
                          constant float &c[[ buffer(2) ]],
                          constant int &n[[ buffer(3) ]],
                          uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = c * b[idx];
}

kernel void cdiv_kernel(device float *a[[ buffer(0) ]],
                        constant float &c[[ buffer(1) ]],
                        constant int &n[[ buffer(2) ]],
                        uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = c / a[idx];
}

kernel void cdiv2_kernel(device float *a[[ buffer(0) ]],
                         device const float *b[[ buffer(1) ]],
                         constant float &c[[ buffer(2) ]],
                         constant int &n[[ buffer(3) ]],
                         uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = c / b[idx];
}

kernel void cadd_kernel(device float *a[[ buffer(0) ]],
                        constant float &c[[ buffer(1) ]],
                        constant int &n[[ buffer(2) ]],
                        uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = a[idx] + c;
}

kernel void cadd2_kernel(device float *a[[ buffer(0) ]],
                         device const float *b[[ buffer(1) ]],
                         constant float &c[[ buffer(2) ]],
                         constant int &n[[ buffer(3) ]],
                         uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = b[idx] + c;
}

kernel void cfill_kernel(device float *a[[ buffer(0) ]],
                         constant float &c[[ buffer(1) ]],
                         constant int &n[[ buffer(2) ]],
                         uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = c;
}

kernel void cfill_mask_kernel(device float *a[[ buffer(0) ]],
                              constant float &c[[ buffer(1) ]],
                              device const int *mask[[ buffer(2) ]],
                              constant int &n_mask[[ buffer(3) ]],
                              uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n_mask) return;
    a[mask[idx]] = c;
}

kernel void masked_copy_kernel(device float *a[[ buffer(0) ]],
                               device const float *b[[ buffer(1) ]],
                               device const int *mask[[ buffer(2) ]],
                               constant int &n_mask[[ buffer(3) ]],
                               uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n_mask) return;
    a[mask[idx + 1] - 1] = b[mask[idx + 1] - 1];
}

kernel void masked_gather_copy_kernel(device float *a[[ buffer(0) ]],
                                      device const float *b[[ buffer(1) ]],
                                      device const int *mask[[ buffer(2) ]],
                                      constant int &n_mask[[ buffer(3) ]],
                                      uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n_mask) return;
    a[idx] = b[mask[idx + 1] - 1];
}

kernel void masked_gather_copy_aligned_kernel(
    device float *a[[ buffer(0) ]],
    device const float *b[[ buffer(1) ]],
    device const int *mask[[ buffer(2) ]],
    constant int &n_mask[[ buffer(3) ]],
    uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n_mask) return;
    a[idx] = b[mask[idx]];
}

kernel void masked_scatter_copy_kernel(device float *a[[ buffer(0) ]],
                                       device const float *b[[ buffer(1) ]],
                                       device const int *mask[[ buffer(2) ]],
                                       constant int &n_mask[[ buffer(3) ]],
                                       uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n_mask) return;
    a[mask[idx + 1] - 1] = b[idx];
}

kernel void masked_scatter_copy_aligned_kernel(
    device float *a[[ buffer(0) ]],
    device const float *b[[ buffer(1) ]],
    device const int *mask[[ buffer(2) ]],
    constant int &n_mask[[ buffer(3) ]],
    uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n_mask) return;
    a[mask[idx]] = b[idx];
}

/* Map a 1-based linear dof index to its (i, j, k, e) location, 1-based. */
static inline void face_gather_nonlinear_index(thread int *index, int idx,
                                               int lx, int ly, int lz) {
    const int idx2 = idx - 1;
    index[3] = idx2 / (lx * ly * lz);
    index[2] = (idx2 - (lx * ly * lz) * index[3]) / (lx * ly);
    index[1] = (idx2 - (lx * ly * lz) * index[3] - (lx * ly) * index[2]) / lx;
    index[0] = (idx2 - (lx * ly * lz) * index[3] - (lx * ly) * index[2]) -
        lx * index[1];
    index[0]++;
    index[1]++;
    index[2]++;
    index[3]++;
}

static inline int face_gather_idx(int i, int j, int k, int l,
                                  int n1, int n2, int nf) {
    return ((i) + (n1) * (((j) - 1) + (n2) * (((k) - 1) + (nf) * (((l) - 1))))) - 1;
}

kernel void face_masked_gather_copy_kernel(
    device float *a[[ buffer(0) ]],
    device const float *b[[ buffer(1) ]],
    device const int *mask[[ buffer(2) ]],
    device const int *facet[[ buffer(3) ]],
    constant int &n1[[ buffer(4) ]],
    constant int &n2[[ buffer(5) ]],
    constant int &lx[[ buffer(6) ]],
    constant int &ly[[ buffer(7) ]],
    constant int &lz[[ buffer(8) ]],
    constant int &n_mask[[ buffer(9) ]],
    uint idx [[ thread_position_in_grid ]]) {
    const int m = (int)idx;
    if (m >= n_mask) return;

    int index[4];
    const int f = facet[m + 1];
    face_gather_nonlinear_index(index, mask[m + 1], lx, ly, lz);

    switch (f) {
    case 1:
    case 2:
        a[m] = b[face_gather_idx(index[1], index[2], f, index[3], n1, n2, 6)];
        break;
    case 3:
    case 4:
        a[m] = b[face_gather_idx(index[0], index[2], f, index[3], n1, n2, 6)];
        break;
    case 5:
    case 6:
        a[m] = b[face_gather_idx(index[0], index[1], f, index[3], n1, n2, 6)];
        break;
    }
}

kernel void cwrap_kernel(device float *a[[ buffer(0) ]],
                         constant float &min_val[[ buffer(1) ]],
                         constant float &max_val[[ buffer(2) ]],
                         constant int &n[[ buffer(3) ]],
                         uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    const float l = max_val - min_val;
    a[idx] = min_val + fmod(fmod(a[idx] - min_val, l) + l, l);
}

kernel void masked_atomic_reduction_kernel(
    device atomic_float *a[[ buffer(0) ]],
    device const float *b[[ buffer(1) ]],
    device const int *mask[[ buffer(2) ]],
    constant int &n_mask[[ buffer(3) ]],
    uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n_mask) return;
    atomic_fetch_add_explicit(&a[mask[idx + 1] - 1], b[idx],
                              memory_order_relaxed);
}

kernel void add2_kernel(device float *a[[ buffer(0) ]],
                        device const float *b[[ buffer(1) ]],
                        constant int &n[[ buffer(2) ]],
                        uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = a[idx] + b[idx];
}

kernel void add3_kernel(device float *a[[ buffer(0) ]],
                        device const float *b[[ buffer(1) ]],
                        device const float *c[[ buffer(2) ]],
                        constant int &n[[ buffer(3) ]],
                        uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = b[idx] + c[idx];
}

kernel void add4_kernel(device float *a[[ buffer(0) ]],
                        device const float *b[[ buffer(1) ]],
                        device const float *c[[ buffer(2) ]],
                        device const float *d[[ buffer(3) ]],
                        constant int &n[[ buffer(4) ]],
                        uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = b[idx] + c[idx] + d[idx];
}

kernel void add2s1_kernel(device float *a[[ buffer(0) ]],
                          device const float *b[[ buffer(1) ]],
                          constant float &c1[[ buffer(2) ]],
                          constant int &n[[ buffer(3) ]],
                          uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = c1 * a[idx] + b[idx];
}

kernel void add2s2_kernel(device float *a[[ buffer(0) ]],
                          device const float *b[[ buffer(1) ]],
                          constant float &c1[[ buffer(2) ]],
                          constant int &n[[ buffer(3) ]],
                          uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = a[idx] + c1 * b[idx];
}

kernel void add2s2_many_kernel(device float *x[[ buffer(0) ]],
                               device const float *p[[ buffer(1) ]],
                               constant float &alpha[[ buffer(2) ]],
                               constant int &n[[ buffer(3) ]],
                               uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    x[idx] += alpha * p[idx];
}

kernel void addsqr2s2_kernel(device float *a[[ buffer(0) ]],
                             device const float *b[[ buffer(1) ]],
                             constant float &c1[[ buffer(2) ]],
                             constant int &n[[ buffer(3) ]],
                             uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = a[idx] + c1 * (b[idx] * b[idx]);
}

kernel void add3s2_kernel(device float *a[[ buffer(0) ]],
                          device const float *b[[ buffer(1) ]],
                          device const float *c[[ buffer(2) ]],
                          constant float &c1[[ buffer(3) ]],
                          constant float &c2[[ buffer(4) ]],
                          constant int &n[[ buffer(5) ]],
                          uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = c1 * b[idx] + c2 * c[idx];
}

kernel void add4s3_kernel(device float *a[[ buffer(0) ]],
                          device const float *b[[ buffer(1) ]],
                          device const float *c[[ buffer(2) ]],
                          device const float *d[[ buffer(3) ]],
                          constant float &c1[[ buffer(4) ]],
                          constant float &c2[[ buffer(5) ]],
                          constant float &c3[[ buffer(6) ]],
                          constant int &n[[ buffer(7) ]],
                          uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = c1 * b[idx] + c2 * c[idx] + c3 * d[idx];
}

kernel void add5s4_kernel(device float *a[[ buffer(0) ]],
                          device const float *b[[ buffer(1) ]],
                          device const float *c[[ buffer(2) ]],
                          device const float *d[[ buffer(3) ]],
                          device const float *e[[ buffer(4) ]],
                          constant float &c1[[ buffer(5) ]],
                          constant float &c2[[ buffer(6) ]],
                          constant float &c3[[ buffer(7) ]],
                          constant float &c4[[ buffer(8) ]],
                          constant int &n[[ buffer(9) ]],
                          uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = a[idx] + c1 * b[idx] + c2 * c[idx] + c3 * d[idx] + c4 * e[idx];
}

kernel void invcol1_kernel(device float *a[[ buffer(0) ]],
                           constant int &n[[ buffer(1) ]],
                           uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = 1.0f / a[idx];
}

kernel void invcol2_kernel(device float *a[[ buffer(0) ]],
                           device const float *b[[ buffer(1) ]],
                           constant int &n[[ buffer(2) ]],
                           uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = a[idx] / b[idx];
}

kernel void invcol3_kernel(device float *a[[ buffer(0) ]],
                           device const float *b[[ buffer(1) ]],
                           device const float *c[[ buffer(2) ]],
                           constant int &n[[ buffer(3) ]],
                           uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = b[idx] / c[idx];
}

kernel void col2_kernel(device float *a[[ buffer(0) ]],
                        device const float *b[[ buffer(1) ]],
                        constant int &n[[ buffer(2) ]],
                        uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = a[idx] * b[idx];
}

kernel void col3_kernel(device float *a[[ buffer(0) ]],
                        device const float *b[[ buffer(1) ]],
                        device const float *c[[ buffer(2) ]],
                        constant int &n[[ buffer(3) ]],
                        uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = b[idx] * c[idx];
}

kernel void subcol3_kernel(device float *a[[ buffer(0) ]],
                           device const float *b[[ buffer(1) ]],
                           device const float *c[[ buffer(2) ]],
                           constant int &n[[ buffer(3) ]],
                           uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = a[idx] - b[idx] * c[idx];
}

kernel void sub2_kernel(device float *a[[ buffer(0) ]],
                        device const float *b[[ buffer(1) ]],
                        constant int &n[[ buffer(2) ]],
                        uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = a[idx] - b[idx];
}

kernel void sub3_kernel(device float *a[[ buffer(0) ]],
                        device const float *b[[ buffer(1) ]],
                        device const float *c[[ buffer(2) ]],
                        constant int &n[[ buffer(3) ]],
                        uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = b[idx] - c[idx];
}

kernel void addcol3_kernel(device float *a[[ buffer(0) ]],
                           device const float *b[[ buffer(1) ]],
                           device const float *c[[ buffer(2) ]],
                           constant int &n[[ buffer(3) ]],
                           uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = a[idx] + b[idx] * c[idx];
}

kernel void addcol4_kernel(device float *a[[ buffer(0) ]],
                           device const float *b[[ buffer(1) ]],
                           device const float *c[[ buffer(2) ]],
                           device const float *d[[ buffer(3) ]],
                           constant int &n[[ buffer(4) ]],
                           uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = a[idx] + b[idx] * c[idx] * d[idx];
}

kernel void addcol3s2_kernel(device float *a[[ buffer(0) ]],
                             device const float *b[[ buffer(1) ]],
                             device const float *c[[ buffer(2) ]],
                             constant float &s[[ buffer(3) ]],
                             constant int &n[[ buffer(4) ]],
                             uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = a[idx] + s * b[idx] * c[idx];
}

kernel void vdot3_kernel(device float *dot[[ buffer(0) ]],
                         device const float *u1[[ buffer(1) ]],
                         device const float *u2[[ buffer(2) ]],
                         device const float *u3[[ buffer(3) ]],
                         device const float *v1[[ buffer(4) ]],
                         device const float *v2[[ buffer(5) ]],
                         device const float *v3[[ buffer(6) ]],
                         constant int &n[[ buffer(7) ]],
                         uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    dot[idx] = u1[idx] * v1[idx] + u2[idx] * v2[idx] + u3[idx] * v3[idx];
}

kernel void vcross_kernel(device float *u1[[ buffer(0) ]],
                          device float *u2[[ buffer(1) ]],
                          device float *u3[[ buffer(2) ]],
                          device const float *v1[[ buffer(3) ]],
                          device const float *v2[[ buffer(4) ]],
                          device const float *v3[[ buffer(5) ]],
                          device const float *w1[[ buffer(6) ]],
                          device const float *w2[[ buffer(7) ]],
                          device const float *w3[[ buffer(8) ]],
                          constant int &n[[ buffer(9) ]],
                          uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    u1[idx] = v2[idx] * w3[idx] - v3[idx] * w2[idx];
    u2[idx] = v3[idx] * w1[idx] - v1[idx] * w3[idx];
    u3[idx] = v1[idx] * w2[idx] - v2[idx] * w1[idx];
}

kernel void absval_kernel(device float *a[[ buffer(0) ]],
                          constant int &n[[ buffer(1) ]],
                          uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = fabs(a[idx]);
}

kernel void iadd_kernel(device int *a[[ buffer(0) ]],
                        constant int &c[[ buffer(1) ]],
                        constant int &n[[ buffer(2) ]],
                        uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = a[idx] + c;
}

// ============================================================================
// Pointwise max/min kernels
// ============================================================================

kernel void pwmax_vec2_kernel(device float *a[[ buffer(0) ]],
                              device const float *b[[ buffer(1) ]],
                              constant int &n[[ buffer(2) ]],
                              uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = max(a[idx], b[idx]);
}

kernel void pwmax_vec3_kernel(device float *a[[ buffer(0) ]],
                              device const float *b[[ buffer(1) ]],
                              device const float *c[[ buffer(2) ]],
                              constant int &n[[ buffer(3) ]],
                              uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = max(b[idx], c[idx]);
}

kernel void pwmax_sca2_kernel(device float *a[[ buffer(0) ]],
                              constant float &c[[ buffer(1) ]],
                              constant int &n[[ buffer(2) ]],
                              uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = max(a[idx], c);
}

kernel void pwmax_sca3_kernel(device float *a[[ buffer(0) ]],
                              device const float *b[[ buffer(1) ]],
                              constant float &c[[ buffer(2) ]],
                              constant int &n[[ buffer(3) ]],
                              uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = max(b[idx], c);
}

kernel void pwmin_vec2_kernel(device float *a[[ buffer(0) ]],
                              device const float *b[[ buffer(1) ]],
                              constant int &n[[ buffer(2) ]],
                              uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = min(a[idx], b[idx]);
}

kernel void pwmin_vec3_kernel(device float *a[[ buffer(0) ]],
                              device const float *b[[ buffer(1) ]],
                              device const float *c[[ buffer(2) ]],
                              constant int &n[[ buffer(3) ]],
                              uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = min(b[idx], c[idx]);
}

kernel void pwmin_sca2_kernel(device float *a[[ buffer(0) ]],
                              constant float &c[[ buffer(1) ]],
                              constant int &n[[ buffer(2) ]],
                              uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = min(a[idx], c);
}

kernel void pwmin_sca3_kernel(device float *a[[ buffer(0) ]],
                              device const float *b[[ buffer(1) ]],
                              constant float &c[[ buffer(2) ]],
                              constant int &n[[ buffer(3) ]],
                              uint idx [[ thread_position_in_grid ]]) {
    if (idx >= (uint)n) return;
    a[idx] = min(b[idx], c);
}

// ============================================================================
// Reduction kernels  (use SIMD group intrinsics)
// ============================================================================

kernel void glsc3_kernel(device const float *a[[ buffer(0) ]],
                         device const float *b[[ buffer(1) ]],
                         device const float *c[[ buffer(2) ]],
                         device float *buf[[ buffer(3) ]],
                         constant int &n[[ buffer(4) ]],
                         uint gid [[ thread_position_in_grid ]],
                         uint tid [[ thread_index_in_threadgroup ]],
                         uint tg_id [[ threadgroup_position_in_grid ]],
                         uint tg_size [[ threads_per_threadgroup ]],
                         uint simd_lane [[ thread_index_in_simdgroup ]],
                         uint simd_id [[ simdgroup_index_in_threadgroup ]],
                         uint num_tg [[ threadgroups_per_grid ]]) {
    float sum = 0.0f;
    for (uint i = gid; i < (uint)n; i += tg_size * num_tg) {
        sum += a[i] * b[i] * c[i];
    }

    sum = simd_sum(sum);

    threadgroup float shared[32];
    if (simd_lane == 0) shared[simd_id] = sum;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    uint num_simd = (tg_size + 31) / 32;
    if (simd_id == 0) {
        sum = (simd_lane < num_simd) ? shared[simd_lane] : 0.0f;
        sum = simd_sum(sum);
        if (simd_lane == 0) buf[tg_id] = sum;
    }
}

kernel void glsc3_many_kernel(device const float *w[[ buffer(0) ]],
                              device const float *v[[ buffer(1) ]],
                              device const float *mult[[ buffer(2) ]],
                              device float *buf[[ buffer(3) ]],
                              constant int &n[[ buffer(4) ]],
                              uint gid [[ thread_position_in_grid ]],
                              uint tid [[ thread_index_in_threadgroup ]],
                              uint tg_id [[ threadgroup_position_in_grid ]],
                              uint tg_size [[ threads_per_threadgroup ]],
                              uint simd_lane [[ thread_index_in_simdgroup ]],
                              uint simd_id [[ simdgroup_index_in_threadgroup ]],
                              uint num_tg [[ threadgroups_per_grid ]]) {
    float sum = 0.0f;
    for (uint i = gid; i < (uint)n; i += tg_size * num_tg) {
        sum += w[i] * v[i] * mult[i];
    }

    sum = simd_sum(sum);

    threadgroup float shared[32];
    if (simd_lane == 0) shared[simd_id] = sum;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    uint num_simd = (tg_size + 31) / 32;
    if (simd_id == 0) {
        sum = (simd_lane < num_simd) ? shared[simd_lane] : 0.0f;
        sum = simd_sum(sum);
        if (simd_lane == 0) buf[tg_id] = sum;
    }
}

kernel void glsc2_kernel(device const float *a[[ buffer(0) ]],
                         device const float *b[[ buffer(1) ]],
                         device float *buf[[ buffer(2) ]],
                         constant int &n[[ buffer(3) ]],
                         uint gid [[ thread_position_in_grid ]],
                         uint tid [[ thread_index_in_threadgroup ]],
                         uint tg_id [[ threadgroup_position_in_grid ]],
                         uint tg_size [[ threads_per_threadgroup ]],
                         uint simd_lane [[ thread_index_in_simdgroup ]],
                         uint simd_id [[ simdgroup_index_in_threadgroup ]],
                         uint num_tg [[ threadgroups_per_grid ]]) {
    float sum = 0.0f;
    for (uint i = gid; i < (uint)n; i += tg_size * num_tg) {
        sum += a[i] * b[i];
    }

    sum = simd_sum(sum);

    threadgroup float shared[32];
    if (simd_lane == 0) shared[simd_id] = sum;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    uint num_simd = (tg_size + 31) / 32;
    if (simd_id == 0) {
        sum = (simd_lane < num_simd) ? shared[simd_lane] : 0.0f;
        sum = simd_sum(sum);
        if (simd_lane == 0) buf[tg_id] = sum;
    }
}

kernel void glsubnorm2_kernel(device const float *a[[ buffer(0) ]],
                              device const float *b[[ buffer(1) ]],
                              device float *buf[[ buffer(2) ]],
                              constant int &n[[ buffer(3) ]],
                              uint gid [[ thread_position_in_grid ]],
                              uint tid [[ thread_index_in_threadgroup ]],
                              uint tg_id [[ threadgroup_position_in_grid ]],
                              uint tg_size [[ threads_per_threadgroup ]],
                              uint simd_lane [[ thread_index_in_simdgroup ]],
                              uint simd_id [[ simdgroup_index_in_threadgroup ]],
                              uint num_tg [[ threadgroups_per_grid ]]) {
    float sum = 0.0f;
    for (uint i = gid; i < (uint)n; i += tg_size * num_tg) {
        float d = a[i] - b[i];
        sum += d * d;
    }

    sum = simd_sum(sum);

    threadgroup float shared[32];
    if (simd_lane == 0) shared[simd_id] = sum;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    uint num_simd = (tg_size + 31) / 32;
    if (simd_id == 0) {
        sum = (simd_lane < num_simd) ? shared[simd_lane] : 0.0f;
        sum = simd_sum(sum);
        if (simd_lane == 0) buf[tg_id] = sum;
    }
}

kernel void glsum_kernel(device const float *a[[ buffer(0) ]],
                         device float *buf[[ buffer(1) ]],
                         constant int &n[[ buffer(2) ]],
                         uint gid [[ thread_position_in_grid ]],
                         uint tid [[ thread_index_in_threadgroup ]],
                         uint tg_id [[ threadgroup_position_in_grid ]],
                         uint tg_size [[ threads_per_threadgroup ]],
                         uint simd_lane [[ thread_index_in_simdgroup ]],
                         uint simd_id [[ simdgroup_index_in_threadgroup ]],
                         uint num_tg [[ threadgroups_per_grid ]]) {
    float sum = 0.0f;
    for (uint i = gid; i < (uint)n; i += tg_size * num_tg) {
        sum += a[i];
    }

    sum = simd_sum(sum);

    threadgroup float shared[32];
    if (simd_lane == 0) shared[simd_id] = sum;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    uint num_simd = (tg_size + 31) / 32;
    if (simd_id == 0) {
        sum = (simd_lane < num_simd) ? shared[simd_lane] : 0.0f;
        sum = simd_sum(sum);
        if (simd_lane == 0) buf[tg_id] = sum;
    }
}

kernel void glmax_kernel(device const float *a[[ buffer(0) ]],
                         device float *buf[[ buffer(1) ]],
                         constant int &n[[ buffer(2) ]],
                         uint gid [[ thread_position_in_grid ]],
                         uint tid [[ thread_index_in_threadgroup ]],
                         uint tg_id [[ threadgroup_position_in_grid ]],
                         uint tg_size [[ threads_per_threadgroup ]],
                         uint simd_lane [[ thread_index_in_simdgroup ]],
                         uint simd_id [[ simdgroup_index_in_threadgroup ]],
                         uint num_tg [[ threadgroups_per_grid ]]) {
    float val = -HUGE_VALF;
    for (uint i = gid; i < (uint)n; i += tg_size * num_tg) {
        val = max(val, a[i]);
    }

    val = simd_max(val);

    threadgroup float shared[32];
    if (simd_lane == 0) shared[simd_id] = val;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    uint num_simd = (tg_size + 31) / 32;
    if (simd_id == 0) {
        val = (simd_lane < num_simd) ? shared[simd_lane] : -HUGE_VALF;
        val = simd_max(val);
        if (simd_lane == 0) buf[tg_id] = val;
    }
}

kernel void glmin_kernel(device const float *a[[ buffer(0) ]],
                         device float *buf[[ buffer(1) ]],
                         constant int &n[[ buffer(2) ]],
                         uint gid [[ thread_position_in_grid ]],
                         uint tid [[ thread_index_in_threadgroup ]],
                         uint tg_id [[ threadgroup_position_in_grid ]],
                         uint tg_size [[ threads_per_threadgroup ]],
                         uint simd_lane [[ thread_index_in_simdgroup ]],
                         uint simd_id [[ simdgroup_index_in_threadgroup ]],
                         uint num_tg [[ threadgroups_per_grid ]]) {
    float val = HUGE_VALF;
    for (uint i = gid; i < (uint)n; i += tg_size * num_tg) {
        val = min(val, a[i]);
    }

    val = simd_min(val);

    threadgroup float shared[32];
    if (simd_lane == 0) shared[simd_id] = val;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    uint num_simd = (tg_size + 31) / 32;
    if (simd_id == 0) {
        val = (simd_lane < num_simd) ? shared[simd_lane] : HUGE_VALF;
        val = simd_min(val);
        if (simd_lane == 0) buf[tg_id] = val;
    }
}

// Second-pass reduction: sum partial results from threadgroups
kernel void reduce_kernel(device float *buf[[ buffer(0) ]],
                          constant int &n[[ buffer(1) ]],
                          uint gid [[ thread_position_in_grid ]],
                          uint tid [[ thread_index_in_threadgroup ]],
                          uint tg_size [[ threads_per_threadgroup ]],
                          uint simd_lane [[ thread_index_in_simdgroup ]],
                          uint simd_id [[ simdgroup_index_in_threadgroup ]]) {
    /* Grid-stride: nblocks may exceed the single threadgroup's thread
       count, and every partial sum must still be folded in. */
    float sum = 0.0f;
    for (uint i = gid; i < (uint)n; i += tg_size) {
        sum += buf[i];
    }
    sum = simd_sum(sum);

    threadgroup float shared[32];
    if (simd_lane == 0) shared[simd_id] = sum;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    uint num_simd = (tg_size + 31) / 32;
    if (simd_id == 0) {
        sum = (simd_lane < num_simd) ? shared[simd_lane] : 0.0f;
        sum = simd_sum(sum);
        if (simd_lane == 0) buf[0] = sum;
    }
}

kernel void reduce_max_kernel(device float *buf[[ buffer(0) ]],
                              constant int &n[[ buffer(1) ]],
                              uint gid [[ thread_position_in_grid ]],
                              uint tid [[ thread_index_in_threadgroup ]],
                              uint tg_size [[ threads_per_threadgroup ]],
                              uint simd_lane [[ thread_index_in_simdgroup ]],
                              uint simd_id [[ simdgroup_index_in_threadgroup ]]) {
    float val = -HUGE_VALF;
    for (uint i = gid; i < (uint)n; i += tg_size) {
        val = max(val, buf[i]);
    }
    val = simd_max(val);

    threadgroup float shared[32];
    if (simd_lane == 0) shared[simd_id] = val;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    uint num_simd = (tg_size + 31) / 32;
    if (simd_id == 0) {
        val = (simd_lane < num_simd) ? shared[simd_lane] : -HUGE_VALF;
        val = simd_max(val);
        if (simd_lane == 0) buf[0] = val;
    }
}

kernel void reduce_min_kernel(device float *buf[[ buffer(0) ]],
                              constant int &n[[ buffer(1) ]],
                              uint gid [[ thread_position_in_grid ]],
                              uint tid [[ thread_index_in_threadgroup ]],
                              uint tg_size [[ threads_per_threadgroup ]],
                              uint simd_lane [[ thread_index_in_simdgroup ]],
                              uint simd_id [[ simdgroup_index_in_threadgroup ]]) {
    float val = HUGE_VALF;
    for (uint i = gid; i < (uint)n; i += tg_size) {
        val = min(val, buf[i]);
    }
    val = simd_min(val);

    threadgroup float shared[32];
    if (simd_lane == 0) shared[simd_id] = val;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    uint num_simd = (tg_size + 31) / 32;
    if (simd_id == 0) {
        val = (simd_lane < num_simd) ? shared[simd_lane] : HUGE_VALF;
        val = simd_min(val);
        if (simd_lane == 0) buf[0] = val;
    }
}
