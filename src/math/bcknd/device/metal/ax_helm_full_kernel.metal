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
 * Metal compute kernels for the stress formulation Ax (full stress tensor).
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/*
 * Template helper – stress formulation kstep
 */

template<int LX>
inline void ax_helm_stress_full_impl(
    device float *au,
    device float *av,
    device float *aw,
    const device float *u,
    const device float *v,
    const device float *w_in,
    threadgroup float *shdx,
    threadgroup float *shdy,
    threadgroup float *shdz,
    threadgroup float *shu,
    threadgroup float *shur,
    threadgroup float *shus,
    threadgroup float *shv,
    threadgroup float *shvr,
    threadgroup float *shvs,
    threadgroup float *shw,
    threadgroup float *shwr,
    threadgroup float *shws,
    const device float *dx,
    const device float *dy,
    const device float *dz,
    const device float *h1,
    const device float *drdx,
    const device float *drdy,
    const device float *drdz,
    const device float *dsdx,
    const device float *dsdy,
    const device float *dsdz,
    const device float *dtdx,
    const device float *dtdy,
    const device float *dtdz,
    const device float *jacinv,
    const device float *weight3,
    uint gid,
    uint i,
    uint j)
{
  const int e   = gid;
  const int ij  = i + j * LX;
  const int ele = e * LX * LX * LX;

  float ru[LX], rv[LX], rw[LX];
  float ruw[LX], rvw[LX], rww[LX];
  float rut, rvt, rwt;

  shdx[ij] = dx[ij];
  shdy[ij] = dy[ij];
  shdz[ij] = dz[ij];

  for (int k = 0; k < LX; ++k) {
    ru[k]  = u[ij + k * LX * LX + ele];
    ruw[k] = 0.0f;
    rv[k]  = v[ij + k * LX * LX + ele];
    rvw[k] = 0.0f;
    rw[k]  = w_in[ij + k * LX * LX + ele];
    rww[k] = 0.0f;
  }

  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int k = 0; k < LX; ++k) {
    const int ijk = ij + k * LX * LX;
    const float drdx_local = drdx[ijk + ele];
    const float drdy_local = drdy[ijk + ele];
    const float drdz_local = drdz[ijk + ele];
    const float dsdx_local = dsdx[ijk + ele];
    const float dsdy_local = dsdy[ijk + ele];
    const float dsdz_local = dsdz[ijk + ele];
    const float dtdx_local = dtdx[ijk + ele];
    const float dtdy_local = dtdy[ijk + ele];
    const float dtdz_local = dtdz[ijk + ele];
    const float dj = h1[ijk + ele] * weight3[ijk] * jacinv[ijk + ele];

    float uttmp = 0.0f;
    float vttmp = 0.0f;
    float wttmp = 0.0f;
    shu[ij] = ru[k];
    shv[ij] = rv[k];
    shw[ij] = rw[k];
    for (int l = 0; l < LX; l++) {
      uttmp += shdz[k + l * LX] * ru[l];
      vttmp += shdz[k + l * LX] * rv[l];
      wttmp += shdz[k + l * LX] * rw[l];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    float urtmp = 0.0f;
    float ustmp = 0.0f;

    float vrtmp = 0.0f;
    float vstmp = 0.0f;

    float wrtmp = 0.0f;
    float wstmp = 0.0f;

    for (int l = 0; l < LX; l++) {
      urtmp += shdx[i + l * LX] * shu[l + j * LX];
      ustmp += shdy[j + l * LX] * shu[i + l * LX];

      vrtmp += shdx[i + l * LX] * shv[l + j * LX];
      vstmp += shdy[j + l * LX] * shv[i + l * LX];

      wrtmp += shdx[i + l * LX] * shw[l + j * LX];
      wstmp += shdy[j + l * LX] * shw[i + l * LX];
    }

    const float u1 = urtmp * drdx_local
      + ustmp * dsdx_local
      + uttmp * dtdx_local;
    const float u2 = urtmp * drdy_local
      + ustmp * dsdy_local
      + uttmp * dtdy_local;
    const float u3 = urtmp * drdz_local
      + ustmp * dsdz_local
      + uttmp * dtdz_local;

    const float v1 = vrtmp * drdx_local
      + vstmp * dsdx_local
      + vttmp * dtdx_local;
    const float v2 = vrtmp * drdy_local
      + vstmp * dsdy_local
      + vttmp * dtdy_local;
    const float v3 = vrtmp * drdz_local
      + vstmp * dsdz_local
      + vttmp * dtdz_local;

    const float w1 = wrtmp * drdx_local
      + wstmp * dsdx_local
      + wttmp * dtdx_local;
    const float w2 = wrtmp * drdy_local
      + wstmp * dsdy_local
      + wttmp * dtdy_local;
    const float w3 = wrtmp * drdz_local
      + wstmp * dsdz_local
      + wttmp * dtdz_local;

    const float s11 = dj * (u1 + u1);
    const float s22 = dj * (v2 + v2);
    const float s33 = dj * (w3 + w3);
    const float s12 = dj * (u2 + v1);
    const float s13 = dj * (u3 + w1);
    const float s23 = dj * (v3 + w2);

    shur[ij] = drdx_local * s11
      + drdy_local * s12
      + drdz_local * s13;
    shus[ij] = dsdx_local * s11
      + dsdy_local * s12
      + dsdz_local * s13;
    rut = dtdx_local * s11
      + dtdy_local * s12
      + dtdz_local * s13;

    shvr[ij] = drdx_local * s12
      + drdy_local * s22
      + drdz_local * s23;
    shvs[ij] = dsdx_local * s12
      + dsdy_local * s22
      + dsdz_local * s23;
    rvt = dtdx_local * s12
      + dtdy_local * s22
      + dtdz_local * s23;

    shwr[ij] = drdx_local * s13
      + drdy_local * s23
      + drdz_local * s33;
    shws[ij] = dsdx_local * s13
      + dsdy_local * s23
      + dsdz_local * s33;
    rwt = dtdx_local * s13
      + dtdy_local * s23
      + dtdz_local * s33;

    threadgroup_barrier(mem_flags::mem_threadgroup);

    float uwijke = 0.0f;
    float vwijke = 0.0f;
    float wwijke = 0.0f;

    for (int l = 0; l < LX; l++) {
      uwijke += shur[l + j * LX] * shdx[l + i * LX];
      ruw[l] += rut * shdz[k + l * LX];
      uwijke += shus[i + l * LX] * shdy[l + j * LX];

      vwijke += shvr[l + j * LX] * shdx[l + i * LX];
      rvw[l] += rvt * shdz[k + l * LX];
      vwijke += shvs[i + l * LX] * shdy[l + j * LX];

      wwijke += shwr[l + j * LX] * shdx[l + i * LX];
      rww[l] += rwt * shdz[k + l * LX];
      wwijke += shws[i + l * LX] * shdy[l + j * LX];
    }
    ruw[k] += uwijke;
    rvw[k] += vwijke;
    rww[k] += wwijke;
  }

  for (int k = 0; k < LX; ++k) {
    au[ij + k * LX * LX + ele] = ruw[k];
    av[ij + k * LX * LX + ele] = rvw[k];
    aw[ij + k * LX * LX + ele] = rww[k];
  }
}

#define INSTANTIATE_AX_HELM_STRESS_FULL(LX)                                   \
kernel void                                                                    \
ax_helm_stress_kernel_full_lx##LX(                                             \
    device float *au[[ buffer(0) ]],                                           \
    device float *av[[ buffer(1) ]],                                           \
    device float *aw[[ buffer(2) ]],                                           \
    const device float *u[[ buffer(3) ]],                                      \
    const device float *v[[ buffer(4) ]],                                      \
    const device float *w[[ buffer(5) ]],                                      \
    const device float *dx[[ buffer(6) ]],                                     \
    const device float *dy[[ buffer(7) ]],                                     \
    const device float *dz[[ buffer(8) ]],                                     \
    const device float *h1[[ buffer(9) ]],                                     \
    const device float *drdx[[ buffer(10) ]],                                  \
    const device float *drdy[[ buffer(11) ]],                                  \
    const device float *drdz[[ buffer(12) ]],                                  \
    const device float *dsdx[[ buffer(13) ]],                                  \
    const device float *dsdy[[ buffer(14) ]],                                  \
    const device float *dsdz[[ buffer(15) ]],                                  \
    const device float *dtdx[[ buffer(16) ]],                                  \
    const device float *dtdy[[ buffer(17) ]],                                  \
    const device float *dtdz[[ buffer(18) ]],                                  \
    const device float *jacinv[[ buffer(19) ]],                                \
    const device float *weight3[[ buffer(20) ]],                               \
    uint2 gid2 [[ threadgroup_position_in_grid ]],                             \
    uint2 tid [[ thread_position_in_threadgroup ]])                            \
{                                                                              \
  const uint gid = gid2.x;                                                    \
  threadgroup float shdx[LX * LX];                                            \
  threadgroup float shdy[LX * LX];                                            \
  threadgroup float shdz[LX * LX];                                            \
  threadgroup float shu [LX * LX];                                            \
  threadgroup float shur[LX * LX];                                            \
  threadgroup float shus[LX * LX];                                            \
  threadgroup float shv [LX * LX];                                            \
  threadgroup float shvr[LX * LX];                                            \
  threadgroup float shvs[LX * LX];                                            \
  threadgroup float shw [LX * LX];                                            \
  threadgroup float shwr[LX * LX];                                            \
  threadgroup float shws[LX * LX];                                            \
  ax_helm_stress_full_impl<LX>(                                               \
      au, av, aw, u, v, w,                                                    \
      shdx, shdy, shdz,                                                       \
      shu, shur, shus, shv, shvr, shvs, shw, shwr, shws,                     \
      dx, dy, dz, h1,                                                        \
      drdx, drdy, drdz, dsdx, dsdy, dsdz, dtdx, dtdy, dtdz,                  \
      jacinv, weight3,                                                       \
      gid, tid.x, tid.y);                                                     \
}

/* ================================================================== */
/*  Stress part2 – pointwise, no shared memory, no template needed     */
/* ================================================================== */

kernel void ax_helm_stress_kernel_vector_part2(
    device float *au[[ buffer(0) ]],
    device float *av[[ buffer(1) ]],
    device float *aw[[ buffer(2) ]],
    const device float *u[[ buffer(3) ]],
    const device float *v[[ buffer(4) ]],
    const device float *w[[ buffer(5) ]],
    const device float *h2[[ buffer(6) ]],
    const device float *B[[ buffer(7) ]],
    constant int &n[[ buffer(8) ]],
    uint idx [[ thread_position_in_grid ]],
    uint str [[ threads_per_grid ]])
{
  for (uint i = idx; i < (uint)n; i += str) {
    au[i] = au[i] + h2[i] * B[i] * u[i];
    av[i] = av[i] + h2[i] * B[i] * v[i];
    aw[i] = aw[i] + h2[i] * B[i] * w[i];
  }
}

/* ================================================================== */
/*  Instantiations for all supported polynomial orders                 */
/* ================================================================== */

INSTANTIATE_AX_HELM_STRESS_FULL(2)
INSTANTIATE_AX_HELM_STRESS_FULL(3)
INSTANTIATE_AX_HELM_STRESS_FULL(4)
INSTANTIATE_AX_HELM_STRESS_FULL(5)
INSTANTIATE_AX_HELM_STRESS_FULL(6)
INSTANTIATE_AX_HELM_STRESS_FULL(7)
INSTANTIATE_AX_HELM_STRESS_FULL(8)
INSTANTIATE_AX_HELM_STRESS_FULL(9)
INSTANTIATE_AX_HELM_STRESS_FULL(10)
INSTANTIATE_AX_HELM_STRESS_FULL(11)
INSTANTIATE_AX_HELM_STRESS_FULL(12)
INSTANTIATE_AX_HELM_STRESS_FULL(13)
INSTANTIATE_AX_HELM_STRESS_FULL(14)
INSTANTIATE_AX_HELM_STRESS_FULL(15)
INSTANTIATE_AX_HELM_STRESS_FULL(16)
