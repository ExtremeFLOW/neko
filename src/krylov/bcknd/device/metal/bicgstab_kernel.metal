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

/**
 * Metal compute kernels for the BiCGStab method.
 *
 * Each kernel fuses one of the vector updates with the reduction that
 * follows it, mirroring the fused loops of the CPU implementation. Every
 * Metal dispatch is a separate command buffer that the host waits on, so
 * collapsing the updates and the reductions removes round trips as well as
 * memory traffic.
 *
 * The reduction kernels leave one partial sum per threadgroup in @a buf and
 * the host performs the final summation; a two value reduction writes its
 * second set of partials at @a buf[nblocks + tg_id].
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/**
 * Reduce @a v across the threadgroup and leave the result in lane 0 of the
 * first SIMD group.
 */
static inline float bicgstab_tg_sum(float v,
                                    threadgroup float *shared,
                                    uint tg_size,
                                    uint simd_lane,
                                    uint simd_id) {
  v = simd_sum(v);
  if (simd_lane == 0) shared[simd_id] = v;
  threadgroup_barrier(mem_flags::mem_threadgroup);

  const uint num_simd = (tg_size + 31) / 32;
  if (simd_id == 0) {
    v = (simd_lane < num_simd) ? shared[simd_lane] : 0.0f;
    v = simd_sum(v);
  }
  return v;
}

/**
 * Search direction update \f$ p = r + \beta (p - \omega v) \f$.
 */
kernel void bicgstab_update_p_kernel(device float *p[[ buffer(0) ]],
                                     device const float *r[[ buffer(1) ]],
                                     device const float *v[[ buffer(2) ]],
                                     constant float &beta[[ buffer(3) ]],
                                     constant float &omega[[ buffer(4) ]],
                                     constant int &n[[ buffer(5) ]],
                                     uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  p[idx] = r[idx] + beta * (p[idx] - omega * v[idx]);
}

/**
 * Weighted inner product and squared norm in a single pass,
 * \f$ (a^T M b, b^T M b) \f$.
 */
kernel void bicgstab_product_and_norm_kernel(
    device const float *a[[ buffer(0) ]],
    device const float *b[[ buffer(1) ]],
    device const float *mult[[ buffer(2) ]],
    device float *buf[[ buffer(3) ]],
    constant int &n[[ buffer(4) ]],
    constant int &nblocks[[ buffer(5) ]],
    uint gid [[ thread_position_in_grid ]],
    uint tg_id [[ threadgroup_position_in_grid ]],
    uint tg_size [[ threads_per_threadgroup ]],
    uint simd_lane [[ thread_index_in_simdgroup ]],
    uint simd_id [[ simdgroup_index_in_threadgroup ]],
    uint num_tg [[ threadgroups_per_grid ]]) {

  float prod = 0.0f;
  float nrm = 0.0f;
  for (uint i = gid; i < (uint)n; i += tg_size * num_tg) {
    const float bi = b[i];
    const float w = mult[i] * bi;
    prod += a[i] * w;
    nrm += bi * w;
  }

  threadgroup float shared_prod[32];
  threadgroup float shared_nrm[32];
  prod = bicgstab_tg_sum(prod, shared_prod, tg_size, simd_lane, simd_id);
  nrm = bicgstab_tg_sum(nrm, shared_nrm, tg_size, simd_lane, simd_id);

  if (simd_id == 0 && simd_lane == 0) {
    buf[tg_id] = prod;
    buf[(uint)nblocks + tg_id] = nrm;
  }
}

/**
 * BiCGStab part 1: \f$ s = r - \alpha v \f$, returning \f$ s^T M s \f$.
 */
kernel void bicgstab_part1_kernel(
    device float *s[[ buffer(0) ]],
    device const float *r[[ buffer(1) ]],
    device const float *v[[ buffer(2) ]],
    device const float *mult[[ buffer(3) ]],
    device float *buf[[ buffer(4) ]],
    constant float &alpha[[ buffer(5) ]],
    constant int &n[[ buffer(6) ]],
    uint gid [[ thread_position_in_grid ]],
    uint tg_id [[ threadgroup_position_in_grid ]],
    uint tg_size [[ threads_per_threadgroup ]],
    uint simd_lane [[ thread_index_in_simdgroup ]],
    uint simd_id [[ simdgroup_index_in_threadgroup ]],
    uint num_tg [[ threadgroups_per_grid ]]) {

  float sts = 0.0f;
  for (uint i = gid; i < (uint)n; i += tg_size * num_tg) {
    const float si = r[i] - alpha * v[i];
    s[i] = si;
    sts += si * mult[i] * si;
  }

  threadgroup float shared[32];
  sts = bicgstab_tg_sum(sts, shared, tg_size, simd_lane, simd_id);

  if (simd_id == 0 && simd_lane == 0) {
    buf[tg_id] = sts;
  }
}

/**
 * BiCGStab part 2: \f$ x = x + \alpha \hat{p} + \omega \hat{s} \f$ and
 * \f$ r = s - \omega t \f$, returning \f$ (r^T M r, f^T M r) \f$.
 *
 * The second value is the rho inner product of the next iteration, which
 * the recurrence would otherwise reduce separately at the top of the loop.
 */
kernel void bicgstab_part2_kernel(
    device float *x[[ buffer(0) ]],
    device float *r[[ buffer(1) ]],
    device const float *p_hat[[ buffer(2) ]],
    device const float *s_hat[[ buffer(3) ]],
    device const float *s[[ buffer(4) ]],
    device const float *t[[ buffer(5) ]],
    device const float *f[[ buffer(6) ]],
    device const float *mult[[ buffer(7) ]],
    device float *buf[[ buffer(8) ]],
    constant float &alpha[[ buffer(9) ]],
    constant float &omega[[ buffer(10) ]],
    constant int &n[[ buffer(11) ]],
    constant int &nblocks[[ buffer(12) ]],
    uint gid [[ thread_position_in_grid ]],
    uint tg_id [[ threadgroup_position_in_grid ]],
    uint tg_size [[ threads_per_threadgroup ]],
    uint simd_lane [[ thread_index_in_simdgroup ]],
    uint simd_id [[ simdgroup_index_in_threadgroup ]],
    uint num_tg [[ threadgroups_per_grid ]]) {

  float rtr = 0.0f;
  float ftr = 0.0f;
  for (uint i = gid; i < (uint)n; i += tg_size * num_tg) {
    x[i] = x[i] + alpha * p_hat[i] + omega * s_hat[i];
    const float ri = s[i] - omega * t[i];
    r[i] = ri;
    const float w = mult[i] * ri;
    rtr += ri * w;
    ftr += f[i] * w;
  }

  threadgroup float shared_rtr[32];
  threadgroup float shared_ftr[32];
  rtr = bicgstab_tg_sum(rtr, shared_rtr, tg_size, simd_lane, simd_id);
  ftr = bicgstab_tg_sum(ftr, shared_ftr, tg_size, simd_lane, simd_id);

  if (simd_id == 0 && simd_lane == 0) {
    buf[tg_id] = rtr;
    buf[(uint)nblocks + tg_id] = ftr;
  }
}
