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
 * Metal compute kernels for tensor product application (tnsr3d).
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/*
 * tnsr3d — v = A x Bt x Ct applied to u, one element per group
 */

void tnsr3d_impl(device float *v,
                 const int nv,
                 device const float *u,
                 const int nu,
                 device const float *A,
                 device const float *Bt,
                 device const float *Ct,
                 threadgroup float *shwork,
                 threadgroup float *shwork2,
                 uint eid, uint tid, uint tpg) {

  const int idx = int(tid);
  const int str = int(tpg);
  const int e = int(eid);

  for (int ii = idx; ii < nu*nu*nv; ii += str) {
    float tmp = 0.0f;
    const int j = ii / nv;
    const int i = ii - j*nv;
    for (int l = 0; l < nu; l++) {
      tmp += A[i+l*nv] * u[l+nu*j+e*nu*nu*nu];
    }
    shwork[ii] = tmp;
  }

  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int ijk = idx; ijk < nu*nv*nv; ijk += str) {
    const int jk = ijk / nv;
    const int i = ijk - jk * nv;
    const int k = jk / nv;
    const int j = jk - k * nv;
    float tmp = 0.0f;
    const int ik2 = i + k*nv*nu;
    for (int l = 0; l < nu; l++) {
      tmp += Bt[l+j*nu] * shwork[l*nv+ik2];
    }
    shwork2[ijk] = tmp;
  }

  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int ijk = idx; ijk < nv*nv*nv; ijk += str) {
    const int jk = ijk / nv;
    const int i = ijk - jk * nv;
    const int k = jk / nv;
    const int j = jk - k * nv;
    float tmp = 0.0f;
    const int ij2 = i + j*nv;
    for (int l = 0; l < nu; l++) {
      tmp += Ct[l+k*nu] * shwork2[ij2 + l*nv*nv];
    }
    v[ijk+e*nv*nv*nv] = tmp;
  }
}

#define INSTANTIATE_TNSR3D(N)                                                   \
kernel void tnsr3d_kernel_n##N(                                                 \
    device float *v[[ buffer(0) ]],                                             \
    constant int &nv[[ buffer(1) ]],                                            \
    device const float *u[[ buffer(2) ]],                                       \
    constant int &nu[[ buffer(3) ]],                                            \
    device const float *A[[ buffer(4) ]],                                       \
    device const float *Bt[[ buffer(5) ]],                                      \
    device const float *Ct[[ buffer(6) ]],                                      \
    uint eid [[ threadgroup_position_in_grid ]],                                \
    uint tid [[ thread_index_in_threadgroup ]],                                 \
    uint tpg [[ threads_per_threadgroup ]]) {                                   \
  threadgroup float shwork[N*N*N];                                              \
  threadgroup float shwork2[N*N*N];                                             \
  tnsr3d_impl(v, nv, u, nu, A, Bt, Ct, shwork, shwork2, eid, tid, tpg);        \
}

INSTANTIATE_TNSR3D(2)
INSTANTIATE_TNSR3D(3)
INSTANTIATE_TNSR3D(4)
INSTANTIATE_TNSR3D(5)
INSTANTIATE_TNSR3D(6)
INSTANTIATE_TNSR3D(7)
INSTANTIATE_TNSR3D(8)
INSTANTIATE_TNSR3D(9)
INSTANTIATE_TNSR3D(10)
INSTANTIATE_TNSR3D(11)
INSTANTIATE_TNSR3D(12)
INSTANTIATE_TNSR3D(13)
INSTANTIATE_TNSR3D(14)
INSTANTIATE_TNSR3D(15)
INSTANTIATE_TNSR3D(16)

/*
 * tnsr3d_el — per-point interpolation matrices, one point per group
 */

void tnsr3d_el_impl(device float *v,
                    const int nv,
                    device const float *u,
                    const int nu,
                    device const float *A,
                    device const float *Bt,
                    device const float *Ct,
                    device const int *elements,
                    threadgroup float *shwork,
                    threadgroup float *shwork2,
                    uint pid, uint tid, uint tpg) {

  const int idx = int(tid);
  const int str = int(tpg);
  const int pt = int(pid);
  const int e = elements[pt];

  for (int ii = idx; ii < nu*nu*nv; ii += str) {
    float tmp = 0.0f;
    const int j = ii / nv;
    const int i = ii - j*nv;
    for (int l = 0; l < nu; l++) {
      tmp += A[i+l*nv+pt*nv*nu] * u[l+nu*j+e*nu*nu*nu];
    }
    shwork[ii] = tmp;
  }

  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int ijk = idx; ijk < nu*nv*nv; ijk += str) {
    const int jk = ijk / nv;
    const int i = ijk - jk * nv;
    const int k = jk / nv;
    const int j = jk - k * nv;
    float tmp = 0.0f;
    const int ik2 = i + k*nv*nu;
    for (int l = 0; l < nu; l++) {
      tmp += Bt[l+j*nu+pt*nv*nu] * shwork[l*nv+ik2];
    }
    shwork2[ijk] = tmp;
  }

  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int ijk = idx; ijk < nv*nv*nv; ijk += str) {
    const int jk = ijk / nv;
    const int i = ijk - jk * nv;
    const int k = jk / nv;
    const int j = jk - k * nv;
    float tmp = 0.0f;
    const int ij2 = i + j*nv;
    for (int l = 0; l < nu; l++) {
      tmp += Ct[l+k*nu+pt*nv*nu] * shwork2[ij2 + l*nv*nv];
    }
    v[ijk+pt*nv*nv*nv] = tmp;
  }
}

#define INSTANTIATE_TNSR3D_EL(N)                                                \
kernel void tnsr3d_el_kernel_n##N(                                              \
    device float *v[[ buffer(0) ]],                                             \
    constant int &nv[[ buffer(1) ]],                                            \
    device const float *u[[ buffer(2) ]],                                       \
    constant int &nu[[ buffer(3) ]],                                            \
    device const float *A[[ buffer(4) ]],                                       \
    device const float *Bt[[ buffer(5) ]],                                      \
    device const float *Ct[[ buffer(6) ]],                                      \
    device const int *elements[[ buffer(7) ]],                                  \
    uint pid [[ threadgroup_position_in_grid ]],                                \
    uint tid [[ thread_index_in_threadgroup ]],                                 \
    uint tpg [[ threads_per_threadgroup ]]) {                                   \
  threadgroup float shwork[N*N*N];                                              \
  threadgroup float shwork2[N*N*N];                                             \
  tnsr3d_el_impl(v, nv, u, nu, A, Bt, Ct, elements,                            \
                 shwork, shwork2, pid, tid, tpg);                              \
}

INSTANTIATE_TNSR3D_EL(2)
INSTANTIATE_TNSR3D_EL(3)
INSTANTIATE_TNSR3D_EL(4)
INSTANTIATE_TNSR3D_EL(5)
INSTANTIATE_TNSR3D_EL(6)
INSTANTIATE_TNSR3D_EL(7)
INSTANTIATE_TNSR3D_EL(8)
INSTANTIATE_TNSR3D_EL(9)
INSTANTIATE_TNSR3D_EL(10)
INSTANTIATE_TNSR3D_EL(11)
INSTANTIATE_TNSR3D_EL(12)
INSTANTIATE_TNSR3D_EL(13)
INSTANTIATE_TNSR3D_EL(14)
INSTANTIATE_TNSR3D_EL(15)
INSTANTIATE_TNSR3D_EL(16)
