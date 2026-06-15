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
 * Metal compute kernels for Ax helm (Helmholtz operator).
 *
 * Uses C++ templates for the inner computation; thin kernel entry points
 * are generated via macros because MSL `kernel` functions themselves
 * cannot be templates.
 *
 * Apple GPUs only support FP32, so T = float throughout.
 */

#include <metal_stdlib>
using namespace metal;

/*
 * Template helpers – kstep (non-padded)
 */

template<int LX>
inline void ax_helm_kstep_impl(
    device float *w,
    const device float *u,
    threadgroup float *shdx,
    threadgroup float *shdy,
    threadgroup float *shdz,
    threadgroup float *shu,
    threadgroup float *shur,
    threadgroup float *shus,
    const device float *dx,
    const device float *dy,
    const device float *dz,
    const device float *h1,
    const device float *g11,
    const device float *g22,
    const device float *g33,
    const device float *g12,
    const device float *g13,
    const device float *g23,
    uint gid,
    uint i,
    uint j)
{
  const int e   = gid;
  const int ij  = i + j * LX;
  const int ele = e * LX * LX * LX;

  float ru[LX];
  float rw[LX];
  float rut;

  shdx[ij] = dx[ij];
  shdy[ij] = dy[ij];
  shdz[ij] = dz[ij];

  for (int k = 0; k < LX; ++k) {
    ru[k] = u[ij + k * LX * LX + ele];
    rw[k] = 0.0f;
  }

  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int k = 0; k < LX; ++k) {
    const int ijk = ij + k * LX * LX;
    const float G00 = g11[ijk + ele];
    const float G11 = g22[ijk + ele];
    const float G22 = g33[ijk + ele];
    const float G01 = g12[ijk + ele];
    const float G02 = g13[ijk + ele];
    const float G12 = g23[ijk + ele];
    const float H1  = h1[ijk + ele];

    float ttmp = 0.0f;
    shu[ij] = ru[k];
    for (int l = 0; l < LX; l++) {
      ttmp += shdz[k + l * LX] * ru[l];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    float rtmp = 0.0f;
    float stmp = 0.0f;
    for (int l = 0; l < LX; l++) {
      rtmp += shdx[i + l * LX] * shu[l + j * LX];
      stmp += shdy[j + l * LX] * shu[i + l * LX];
    }

    shur[ij] = H1 * (G00 * rtmp + G01 * stmp + G02 * ttmp);
    shus[ij] = H1 * (G01 * rtmp + G11 * stmp + G12 * ttmp);
    rut      = H1 * (G02 * rtmp + G12 * stmp + G22 * ttmp);

    threadgroup_barrier(mem_flags::mem_threadgroup);

    float wijke = 0.0f;
    for (int l = 0; l < LX; l++) {
      wijke += shur[l + j * LX] * shdx[l + i * LX];
      rw[l] += rut * shdz[k + l * LX];
      wijke += shus[i + l * LX] * shdy[l + j * LX];
    }
    rw[k] += wijke;
  }

  for (int k = 0; k < LX; ++k) {
    w[ij + k * LX * LX + ele] = rw[k];
  }
}

/*
 * Template helpers – kstep padded
 */

template<int LX>
inline void ax_helm_kstep_padded_impl(
    device float *w,
    const device float *u,
    threadgroup float *shdx,
    threadgroup float *shdy,
    threadgroup float *shdz,
    threadgroup float *shu,
    threadgroup float *shur,
    threadgroup float *shus,
    const device float *dx,
    const device float *dy,
    const device float *dz,
    const device float *h1,
    const device float *g11,
    const device float *g22,
    const device float *g33,
    const device float *g12,
    const device float *g13,
    const device float *g23,
    uint gid,
    uint i,
    uint j)
{
  const int e    = gid;
  const int ij   = i + j * LX;
  const int ij_p = i + j * (LX + 1);
  const int ele  = e * LX * LX * LX;

  float ru[LX];
  float rw[LX];
  float rut;

  shdx[ij_p] = dx[ij];
  shdy[ij_p] = dy[ij];
  shdz[ij_p] = dz[ij];

  for (int k = 0; k < LX; ++k) {
    ru[k] = u[ij + k * LX * LX + ele];
    rw[k] = 0.0f;
  }

  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int k = 0; k < LX; ++k) {
    const int ijk = ij + k * LX * LX;
    const float G00 = g11[ijk + ele];
    const float G11 = g22[ijk + ele];
    const float G22 = g33[ijk + ele];
    const float G01 = g12[ijk + ele];
    const float G02 = g13[ijk + ele];
    const float G12 = g23[ijk + ele];
    const float H1  = h1[ijk + ele];

    float ttmp = 0.0f;
    shu[ij_p] = ru[k];
    for (int l = 0; l < LX; l++) {
      ttmp += shdz[k + l * (LX + 1)] * ru[l];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    float rtmp = 0.0f;
    float stmp = 0.0f;
    for (int l = 0; l < LX; l++) {
      rtmp += shdx[i + l * (LX + 1)] * shu[l + j * (LX + 1)];
      stmp += shdy[j + l * (LX + 1)] * shu[i + l * (LX + 1)];
    }

    shur[ij]   = H1 * (G00 * rtmp + G01 * stmp + G02 * ttmp);
    shus[ij_p] = H1 * (G01 * rtmp + G11 * stmp + G12 * ttmp);
    rut        = H1 * (G02 * rtmp + G12 * stmp + G22 * ttmp);

    threadgroup_barrier(mem_flags::mem_threadgroup);

    float wijke = 0.0f;
    for (int l = 0; l < LX; l++) {
      wijke += shur[l + j * LX]       * shdx[l + i * (LX + 1)];
      rw[l] += rut * shdz[k + l * (LX + 1)];
      wijke += shus[i + l * (LX + 1)] * shdy[l + j * (LX + 1)];
    }
    rw[k] += wijke;
  }

  for (int k = 0; k < LX; ++k) {
    w[ij + k * LX * LX + ele] = rw[k];
  }
}

/*
 * Template helpers – vector kstep (non-padded)
 */

template<int LX>
inline void ax_helm_vector_kstep_impl(
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
    const device float *g11,
    const device float *g22,
    const device float *g33,
    const device float *g12,
    const device float *g13,
    const device float *g23,
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
    const float G00 = g11[ijk + ele];
    const float G11 = g22[ijk + ele];
    const float G22 = g33[ijk + ele];
    const float G01 = g12[ijk + ele];
    const float G02 = g13[ijk + ele];
    const float G12 = g23[ijk + ele];
    const float H1  = h1[ijk + ele];

    float uttmp = 0.0f, vttmp = 0.0f, wttmp = 0.0f;
    shu[ij] = ru[k];
    shv[ij] = rv[k];
    shw[ij] = rw[k];
    for (int l = 0; l < LX; l++) {
      uttmp += shdz[k + l * LX] * ru[l];
      vttmp += shdz[k + l * LX] * rv[l];
      wttmp += shdz[k + l * LX] * rw[l];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    float urtmp = 0.0f, ustmp = 0.0f;
    float vrtmp = 0.0f, vstmp = 0.0f;
    float wrtmp = 0.0f, wstmp = 0.0f;
    for (int l = 0; l < LX; l++) {
      urtmp += shdx[i + l * LX] * shu[l + j * LX];
      ustmp += shdy[j + l * LX] * shu[i + l * LX];
      vrtmp += shdx[i + l * LX] * shv[l + j * LX];
      vstmp += shdy[j + l * LX] * shv[i + l * LX];
      wrtmp += shdx[i + l * LX] * shw[l + j * LX];
      wstmp += shdy[j + l * LX] * shw[i + l * LX];
    }

    shur[ij] = H1 * (G00 * urtmp + G01 * ustmp + G02 * uttmp);
    shus[ij] = H1 * (G01 * urtmp + G11 * ustmp + G12 * uttmp);
    rut      = H1 * (G02 * urtmp + G12 * ustmp + G22 * uttmp);

    shvr[ij] = H1 * (G00 * vrtmp + G01 * vstmp + G02 * vttmp);
    shvs[ij] = H1 * (G01 * vrtmp + G11 * vstmp + G12 * vttmp);
    rvt      = H1 * (G02 * vrtmp + G12 * vstmp + G22 * vttmp);

    shwr[ij] = H1 * (G00 * wrtmp + G01 * wstmp + G02 * wttmp);
    shws[ij] = H1 * (G01 * wrtmp + G11 * wstmp + G12 * wttmp);
    rwt      = H1 * (G02 * wrtmp + G12 * wstmp + G22 * wttmp);

    threadgroup_barrier(mem_flags::mem_threadgroup);

    float uwijke = 0.0f, vwijke = 0.0f, wwijke = 0.0f;
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

/*
 * Template helpers – vector kstep padded
 */

template<int LX>
inline void ax_helm_vector_kstep_padded_impl(
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
    const device float *g11,
    const device float *g22,
    const device float *g33,
    const device float *g12,
    const device float *g13,
    const device float *g23,
    uint gid,
    uint i,
    uint j)
{
  const int e    = gid;
  const int ij   = i + j * LX;
  const int ij_p = i + j * (LX + 1);
  const int ele  = e * LX * LX * LX;

  float ru[LX], rv[LX], rw[LX];
  float ruw[LX], rvw[LX], rww[LX];
  float rut, rvt, rwt;

  shdx[ij_p] = dx[ij];
  shdy[ij_p] = dy[ij];
  shdz[ij_p] = dz[ij];

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
    const float G00 = g11[ijk + ele];
    const float G11 = g22[ijk + ele];
    const float G22 = g33[ijk + ele];
    const float G01 = g12[ijk + ele];
    const float G02 = g13[ijk + ele];
    const float G12 = g23[ijk + ele];
    const float H1  = h1[ijk + ele];

    float uttmp = 0.0f, vttmp = 0.0f, wttmp = 0.0f;
    shu[ij_p] = ru[k];
    shv[ij_p] = rv[k];
    shw[ij_p] = rw[k];
    for (int l = 0; l < LX; l++) {
      uttmp += shdz[k + l * (LX + 1)] * ru[l];
      vttmp += shdz[k + l * (LX + 1)] * rv[l];
      wttmp += shdz[k + l * (LX + 1)] * rw[l];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    float urtmp = 0.0f, ustmp = 0.0f;
    float vrtmp = 0.0f, vstmp = 0.0f;
    float wrtmp = 0.0f, wstmp = 0.0f;
    for (int l = 0; l < LX; l++) {
      urtmp += shdx[i + l * (LX + 1)] * shu[l + j * (LX + 1)];
      ustmp += shdy[j + l * (LX + 1)] * shu[i + l * (LX + 1)];
      vrtmp += shdx[i + l * (LX + 1)] * shv[l + j * (LX + 1)];
      vstmp += shdy[j + l * (LX + 1)] * shv[i + l * (LX + 1)];
      wrtmp += shdx[i + l * (LX + 1)] * shw[l + j * (LX + 1)];
      wstmp += shdy[j + l * (LX + 1)] * shw[i + l * (LX + 1)];
    }

    shur[ij]   = H1 * (G00 * urtmp + G01 * ustmp + G02 * uttmp);
    shus[ij_p] = H1 * (G01 * urtmp + G11 * ustmp + G12 * uttmp);
    rut        = H1 * (G02 * urtmp + G12 * ustmp + G22 * uttmp);

    shvr[ij]   = H1 * (G00 * vrtmp + G01 * vstmp + G02 * vttmp);
    shvs[ij_p] = H1 * (G01 * vrtmp + G11 * vstmp + G12 * vttmp);
    rvt        = H1 * (G02 * vrtmp + G12 * vstmp + G22 * vttmp);

    shwr[ij]   = H1 * (G00 * wrtmp + G01 * wstmp + G02 * wttmp);
    shws[ij_p] = H1 * (G01 * wrtmp + G11 * wstmp + G12 * wttmp);
    rwt        = H1 * (G02 * wrtmp + G12 * wstmp + G22 * wttmp);

    threadgroup_barrier(mem_flags::mem_threadgroup);

    float uwijke = 0.0f, vwijke = 0.0f, wwijke = 0.0f;
    for (int l = 0; l < LX; l++) {
      uwijke += shur[l + j * LX]       * shdx[l + i * (LX + 1)];
      ruw[l] += rut * shdz[k + l * (LX + 1)];
      uwijke += shus[i + l * (LX + 1)] * shdy[l + j * (LX + 1)];

      vwijke += shvr[l + j * LX]       * shdx[l + i * (LX + 1)];
      rvw[l] += rvt * shdz[k + l * (LX + 1)];
      vwijke += shvs[i + l * (LX + 1)] * shdy[l + j * (LX + 1)];

      wwijke += shwr[l + j * LX]       * shdx[l + i * (LX + 1)];
      rww[l] += rwt * shdz[k + l * (LX + 1)];
      wwijke += shws[i + l * (LX + 1)] * shdy[l + j * (LX + 1)];
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

/* ================================================================== */
/*  Kernel entry-point macros                                          */
/*  MSL `kernel` functions cannot be templates, so we stamp out thin   */
/*  wrappers that allocate threadgroup memory and forward to the       */
/*  template implementations above.                                    */
/* ================================================================== */

/* ---- Scalar kstep (non-padded) ----------------------------------- */

#define INSTANTIATE_AX_HELM_KSTEP(LX)                                        \
kernel void                                                                   \
ax_helm_kernel_kstep_lx##LX(                                                  \
    device float *w[[ buffer(0) ]],                                           \
    const device float *u[[ buffer(1) ]],                                     \
    const device float *dx[[ buffer(2) ]],                                    \
    const device float *dy[[ buffer(3) ]],                                    \
    const device float *dz[[ buffer(4) ]],                                    \
    const device float *h1[[ buffer(5) ]],                                    \
    const device float *g11[[ buffer(6) ]],                                   \
    const device float *g22[[ buffer(7) ]],                                   \
    const device float *g33[[ buffer(8) ]],                                   \
    const device float *g12[[ buffer(9) ]],                                   \
    const device float *g13[[ buffer(10) ]],                                  \
    const device float *g23[[ buffer(11) ]],                                  \
    uint2 gid2 [[ threadgroup_position_in_grid ]],                            \
    uint2 tid [[ thread_position_in_threadgroup ]])                           \
{                                                                             \
  const uint gid = gid2.x;                                                   \
  threadgroup float shdx[LX * LX];                                           \
  threadgroup float shdy[LX * LX];                                           \
  threadgroup float shdz[LX * LX];                                           \
  threadgroup float shu [LX * LX];                                           \
  threadgroup float shur[LX * LX];                                           \
  threadgroup float shus[LX * LX];                                           \
  ax_helm_kstep_impl<LX>(                                                    \
      w, u, shdx, shdy, shdz, shu, shur, shus,                              \
      dx, dy, dz, h1, g11, g22, g33, g12, g13, g23,                         \
      gid, tid.x, tid.y);                                                    \
}

/* ---- Scalar kstep padded ----------------------------------------- */

#define INSTANTIATE_AX_HELM_KSTEP_PADDED(LX)                                 \
kernel void                                                                   \
ax_helm_kernel_kstep_padded_lx##LX(                                           \
    device float *w[[ buffer(0) ]],                                           \
    const device float *u[[ buffer(1) ]],                                     \
    const device float *dx[[ buffer(2) ]],                                    \
    const device float *dy[[ buffer(3) ]],                                    \
    const device float *dz[[ buffer(4) ]],                                    \
    const device float *h1[[ buffer(5) ]],                                    \
    const device float *g11[[ buffer(6) ]],                                   \
    const device float *g22[[ buffer(7) ]],                                   \
    const device float *g33[[ buffer(8) ]],                                   \
    const device float *g12[[ buffer(9) ]],                                   \
    const device float *g13[[ buffer(10) ]],                                  \
    const device float *g23[[ buffer(11) ]],                                  \
    uint2 gid2 [[ threadgroup_position_in_grid ]],                            \
    uint2 tid [[ thread_position_in_threadgroup ]])                           \
{                                                                             \
  const uint gid = gid2.x;                                                   \
  threadgroup float shdx[LX * (LX + 1)];                                     \
  threadgroup float shdy[LX * (LX + 1)];                                     \
  threadgroup float shdz[LX * (LX + 1)];                                     \
  threadgroup float shu [LX * (LX + 1)];                                     \
  threadgroup float shur[LX * LX];                                           \
  threadgroup float shus[LX * (LX + 1)];                                     \
  ax_helm_kstep_padded_impl<LX>(                                             \
      w, u, shdx, shdy, shdz, shu, shur, shus,                              \
      dx, dy, dz, h1, g11, g22, g33, g12, g13, g23,                         \
      gid, tid.x, tid.y);                                                    \
}

/* ---- Vector kstep (non-padded) ----------------------------------- */

#define INSTANTIATE_AX_HELM_VECTOR_KSTEP(LX)                                 \
kernel void                                                                   \
ax_helm_kernel_vector_kstep_lx##LX(                                           \
    device float *au[[ buffer(0) ]],                                          \
    device float *av[[ buffer(1) ]],                                          \
    device float *aw[[ buffer(2) ]],                                          \
    const device float *u[[ buffer(3) ]],                                     \
    const device float *v[[ buffer(4) ]],                                     \
    const device float *w[[ buffer(5) ]],                                     \
    const device float *dx[[ buffer(6) ]],                                    \
    const device float *dy[[ buffer(7) ]],                                    \
    const device float *dz[[ buffer(8) ]],                                    \
    const device float *h1[[ buffer(9) ]],                                    \
    const device float *g11[[ buffer(10) ]],                                  \
    const device float *g22[[ buffer(11) ]],                                  \
    const device float *g33[[ buffer(12) ]],                                  \
    const device float *g12[[ buffer(13) ]],                                  \
    const device float *g13[[ buffer(14) ]],                                  \
    const device float *g23[[ buffer(15) ]],                                  \
    uint2 gid2 [[ threadgroup_position_in_grid ]],                            \
    uint2 tid [[ thread_position_in_threadgroup ]])                           \
{                                                                             \
  const uint gid = gid2.x;                                                   \
  threadgroup float shdx[LX * LX];                                           \
  threadgroup float shdy[LX * LX];                                           \
  threadgroup float shdz[LX * LX];                                           \
  threadgroup float shu [LX * LX];                                           \
  threadgroup float shur[LX * LX];                                           \
  threadgroup float shus[LX * LX];                                           \
  threadgroup float shv [LX * LX];                                           \
  threadgroup float shvr[LX * LX];                                           \
  threadgroup float shvs[LX * LX];                                           \
  threadgroup float shw [LX * LX];                                           \
  threadgroup float shwr[LX * LX];                                           \
  threadgroup float shws[LX * LX];                                           \
  ax_helm_vector_kstep_impl<LX>(                                             \
      au, av, aw, u, v, w,                                                   \
      shdx, shdy, shdz,                                                      \
      shu, shur, shus, shv, shvr, shvs, shw, shwr, shws,                    \
      dx, dy, dz, h1, g11, g22, g33, g12, g13, g23,                         \
      gid, tid.x, tid.y);                                                    \
}

/* ---- Vector kstep padded ----------------------------------------- */

#define INSTANTIATE_AX_HELM_VECTOR_KSTEP_PADDED(LX)                          \
kernel void                                                                   \
ax_helm_kernel_vector_kstep_padded_lx##LX(                                    \
    device float *au[[ buffer(0) ]],                                          \
    device float *av[[ buffer(1) ]],                                          \
    device float *aw[[ buffer(2) ]],                                          \
    const device float *u[[ buffer(3) ]],                                     \
    const device float *v[[ buffer(4) ]],                                     \
    const device float *w[[ buffer(5) ]],                                     \
    const device float *dx[[ buffer(6) ]],                                    \
    const device float *dy[[ buffer(7) ]],                                    \
    const device float *dz[[ buffer(8) ]],                                    \
    const device float *h1[[ buffer(9) ]],                                    \
    const device float *g11[[ buffer(10) ]],                                  \
    const device float *g22[[ buffer(11) ]],                                  \
    const device float *g33[[ buffer(12) ]],                                  \
    const device float *g12[[ buffer(13) ]],                                  \
    const device float *g13[[ buffer(14) ]],                                  \
    const device float *g23[[ buffer(15) ]],                                  \
    uint2 gid2 [[ threadgroup_position_in_grid ]],                            \
    uint2 tid [[ thread_position_in_threadgroup ]])                           \
{                                                                             \
  const uint gid = gid2.x;                                                   \
  threadgroup float shdx[LX * (LX + 1)];                                     \
  threadgroup float shdy[LX * (LX + 1)];                                     \
  threadgroup float shdz[LX * (LX + 1)];                                     \
  threadgroup float shu [LX * (LX + 1)];                                     \
  threadgroup float shur[LX * LX];                                           \
  threadgroup float shus[LX * (LX + 1)];                                     \
  threadgroup float shv [LX * (LX + 1)];                                     \
  threadgroup float shvr[LX * LX];                                           \
  threadgroup float shvs[LX * (LX + 1)];                                     \
  threadgroup float shw [LX * (LX + 1)];                                     \
  threadgroup float shwr[LX * LX];                                           \
  threadgroup float shws[LX * (LX + 1)];                                     \
  ax_helm_vector_kstep_padded_impl<LX>(                                      \
      au, av, aw, u, v, w,                                                   \
      shdx, shdy, shdz,                                                      \
      shu, shur, shus, shv, shvr, shvs, shw, shwr, shws,                    \
      dx, dy, dz, h1, g11, g22, g33, g12, g13, g23,                         \
      gid, tid.x, tid.y);                                                    \
}

/* ================================================================== */
/*  Vector part2 – pointwise, no shared memory, no template needed     */
/* ================================================================== */

kernel void ax_helm_kernel_vector_part2(
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

/* Non-padded: LX = 2..16 */
INSTANTIATE_AX_HELM_KSTEP(2)
INSTANTIATE_AX_HELM_KSTEP(3)
INSTANTIATE_AX_HELM_KSTEP(4)
INSTANTIATE_AX_HELM_KSTEP(5)
INSTANTIATE_AX_HELM_KSTEP(6)
INSTANTIATE_AX_HELM_KSTEP(7)
INSTANTIATE_AX_HELM_KSTEP(8)
INSTANTIATE_AX_HELM_KSTEP(9)
INSTANTIATE_AX_HELM_KSTEP(10)
INSTANTIATE_AX_HELM_KSTEP(11)
INSTANTIATE_AX_HELM_KSTEP(12)
INSTANTIATE_AX_HELM_KSTEP(13)
INSTANTIATE_AX_HELM_KSTEP(14)
INSTANTIATE_AX_HELM_KSTEP(15)
INSTANTIATE_AX_HELM_KSTEP(16)

/* Padded: LX = 4, 8, 16 (powers of 2 that benefit from padding) */
INSTANTIATE_AX_HELM_KSTEP_PADDED(4)
INSTANTIATE_AX_HELM_KSTEP_PADDED(8)
INSTANTIATE_AX_HELM_KSTEP_PADDED(16)

/* Vector non-padded: LX = 2..16 */
INSTANTIATE_AX_HELM_VECTOR_KSTEP(2)
INSTANTIATE_AX_HELM_VECTOR_KSTEP(3)
INSTANTIATE_AX_HELM_VECTOR_KSTEP(4)
INSTANTIATE_AX_HELM_VECTOR_KSTEP(5)
INSTANTIATE_AX_HELM_VECTOR_KSTEP(6)
INSTANTIATE_AX_HELM_VECTOR_KSTEP(7)
INSTANTIATE_AX_HELM_VECTOR_KSTEP(8)
INSTANTIATE_AX_HELM_VECTOR_KSTEP(9)
INSTANTIATE_AX_HELM_VECTOR_KSTEP(10)
INSTANTIATE_AX_HELM_VECTOR_KSTEP(11)
INSTANTIATE_AX_HELM_VECTOR_KSTEP(12)
INSTANTIATE_AX_HELM_VECTOR_KSTEP(13)
INSTANTIATE_AX_HELM_VECTOR_KSTEP(14)
INSTANTIATE_AX_HELM_VECTOR_KSTEP(15)
INSTANTIATE_AX_HELM_VECTOR_KSTEP(16)

/* Vector padded: LX = 4, 8, 16 */
INSTANTIATE_AX_HELM_VECTOR_KSTEP_PADDED(4)
INSTANTIATE_AX_HELM_VECTOR_KSTEP_PADDED(8)
INSTANTIATE_AX_HELM_VECTOR_KSTEP_PADDED(16)
