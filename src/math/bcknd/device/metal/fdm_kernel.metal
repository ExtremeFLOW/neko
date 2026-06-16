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
 * Metal compute kernels for the fast diagonalization method (FDM).
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/*
 * fdm_do_fast — FDM tensor solve, one element per threadgroup
 */

template <int NL>
void fdm_do_fast_impl(device float *e,
                      device float *r,
                      device const float *s,
                      device const float *d,
                      threadgroup float *shwork,
                      threadgroup float *shwork2,
                      threadgroup float *A,
                      threadgroup float *Bt,
                      threadgroup float *Ct,
                      uint eid, uint tid, uint tpg) {

  const int idx = int(tid);
  const int str = int(tpg);
  const int el = int(eid);

  if (idx < NL * NL) {
    A[idx]  = s[idx + NL * NL + el * NL * NL * 3 * 2];
    Bt[idx] = s[idx + 2 * NL * NL + el * NL * NL * 3 * 2];
    Ct[idx] = s[idx + 2 * 2 * NL * NL + el * NL * NL * 3 * 2];
  }

  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int ii = idx; ii < NL * NL * NL; ii += str) {
    float tmp = 0.0f;
    const int j = ii / NL;
    const int i = ii - j * NL;
    for (int l = 0; l < NL; l++) {
      tmp += A[i + l * NL] * r[l + NL * j + el * NL * NL * NL];
    }
    shwork[ii] = tmp;
  }

  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int ijk = idx; ijk < NL * NL * NL; ijk += str) {
    const int jk = ijk / NL;
    const int i = ijk - jk * NL;
    const int k = jk / NL;
    const int j = jk - k * NL;
    float tmp = 0.0f;
    const int ik2 = i + k * NL * NL;
    for (int l = 0; l < NL; l++) {
      tmp += Bt[l + j * NL] * shwork[l * NL + ik2];
    }
    shwork2[ijk] = tmp;
  }

  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int ijk = idx; ijk < NL * NL * NL; ijk += str) {
    const int jk = ijk / NL;
    const int i = ijk - jk * NL;
    const int k = jk / NL;
    const int j = jk - k * NL;
    float tmp = 0.0f;
    const int ij2 = i + j * NL;
    for (int l = 0; l < NL; l++) {
      tmp += Ct[l + k * NL] * shwork2[ij2 + l * NL * NL];
    }
    r[ijk + el * NL * NL * NL] = tmp * d[ijk + el * NL * NL * NL];
  }

  threadgroup_barrier(mem_flags::mem_threadgroup);

  if (idx < NL * NL) {
    A[idx]  = s[idx + el * NL * NL * 3 * 2];
    Bt[idx] = s[idx + NL * NL + 2 * NL * NL + el * NL * NL * 3 * 2];
    Ct[idx] = s[idx + NL * NL + 2 * 2 * NL * NL + el * NL * NL * 3 * 2];
  }

  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int ii = idx; ii < NL * NL * NL; ii += str) {
    float tmp = 0.0f;
    const int j = ii / NL;
    const int i = ii - j * NL;
    for (int l = 0; l < NL; l++) {
      tmp += A[i + l * NL] * r[l + NL * j + el * NL * NL * NL];
    }
    shwork[ii] = tmp;
  }

  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int ijk = idx; ijk < NL * NL * NL; ijk += str) {
    const int jk = ijk / NL;
    const int i = ijk - jk * NL;
    const int k = jk / NL;
    const int j = jk - k * NL;
    float tmp = 0.0f;
    const int ik2 = i + k * NL * NL;
    for (int l = 0; l < NL; l++) {
      tmp += Bt[l + j * NL] * shwork[l * NL + ik2];
    }
    shwork2[ijk] = tmp;
  }

  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int ijk = idx; ijk < NL * NL * NL; ijk += str) {
    const int jk = ijk / NL;
    const int i = ijk - jk * NL;
    const int k = jk / NL;
    const int j = jk - k * NL;
    float tmp = 0.0f;
    const int ij2 = i + j * NL;
    for (int l = 0; l < NL; l++) {
      tmp += Ct[l + k * NL] * shwork2[ij2 + l * NL * NL];
    }
    e[ijk + el * NL * NL * NL] = tmp;
  }
}

#define INSTANTIATE_FDM_DO_FAST(NL)                                             \
kernel void fdm_do_fast_kernel_nl##NL(                                         \
    device float *e[[ buffer(0) ]],                                            \
    device float *r[[ buffer(1) ]],                                            \
    device const float *s[[ buffer(2) ]],                                      \
    device const float *d[[ buffer(3) ]],                                      \
    uint eid [[ threadgroup_position_in_grid ]],                               \
    uint tid [[ thread_index_in_threadgroup ]],                                \
    uint tpg [[ threads_per_threadgroup ]]) {                                  \
  threadgroup float shwork[NL * NL * NL];                                      \
  threadgroup float shwork2[NL * NL * NL];                                     \
  threadgroup float A[NL * NL];                                                \
  threadgroup float Bt[NL * NL];                                               \
  threadgroup float Ct[NL * NL];                                               \
  fdm_do_fast_impl<NL>(e, r, s, d, shwork, shwork2, A, Bt, Ct,                \
                       eid, tid, tpg);                                         \
}

INSTANTIATE_FDM_DO_FAST(2)
INSTANTIATE_FDM_DO_FAST(3)
INSTANTIATE_FDM_DO_FAST(4)
INSTANTIATE_FDM_DO_FAST(5)
INSTANTIATE_FDM_DO_FAST(6)
INSTANTIATE_FDM_DO_FAST(7)
INSTANTIATE_FDM_DO_FAST(8)
INSTANTIATE_FDM_DO_FAST(9)
INSTANTIATE_FDM_DO_FAST(10)
INSTANTIATE_FDM_DO_FAST(11)
INSTANTIATE_FDM_DO_FAST(12)
INSTANTIATE_FDM_DO_FAST(13)
INSTANTIATE_FDM_DO_FAST(14)
INSTANTIATE_FDM_DO_FAST(15)
