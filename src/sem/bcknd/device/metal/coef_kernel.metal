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
 * Metal compute kernels for SEM coefficient generation.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/*
 * Geometry kernel — compute G11..G33 and multiply by w3
 */

template <int LX>
void coef_generate_geo_impl(device float *G11,
                            device float *G12,
                            device float *G13,
                            device float *G22,
                            device float *G23,
                            device float *G33,
                            device const float *drdx,
                            device const float *drdy,
                            device const float *drdz,
                            device const float *dsdx,
                            device const float *dsdy,
                            device const float *dsdz,
                            device const float *dtdx,
                            device const float *dtdy,
                            device const float *dtdz,
                            device const float *jacinv,
                            device const float *w3,
                            constant int &gdim,
                            threadgroup float *shw3,
                            uint eid,
                            uint tid,
                            uint tpg) {

  constexpr int N3 = LX * LX * LX;
  const int e = eid;

  /* First pass: compute Gij and load w3 into shared memory */
  for (int j = tid; j < N3; j += tpg) {
    const int i = j + e * N3;
    G11[i] = (drdx[i]*drdx[i] + drdy[i]*drdy[i] + drdz[i]*drdz[i]) * jacinv[i];
    G22[i] = (dsdx[i]*dsdx[i] + dsdy[i]*dsdy[i] + dsdz[i]*dsdz[i]) * jacinv[i];
    G33[i] = (dtdx[i]*dtdx[i] + dtdy[i]*dtdy[i] + dtdz[i]*dtdz[i]) * jacinv[i];

    G12[i] = (drdx[i]*dsdx[i] + drdy[i]*dsdy[i] + drdz[i]*dsdz[i]) * jacinv[i];
    G13[i] = (drdx[i]*dtdx[i] + drdy[i]*dtdy[i] + drdz[i]*dtdz[i]) * jacinv[i];
    G23[i] = (dsdx[i]*dtdx[i] + dsdy[i]*dtdy[i] + dsdz[i]*dtdz[i]) * jacinv[i];

    shw3[j] = w3[j];
  }

  threadgroup_barrier(mem_flags::mem_threadgroup);

  /* Second pass: multiply by quadrature weights */
  for (int j = tid; j < N3; j += tpg) {
    const int idx = j + e * N3;
    G11[idx] *= shw3[j];
    G12[idx] *= shw3[j];
    G13[idx] *= shw3[j];
    G22[idx] *= shw3[j];
    G23[idx] *= shw3[j];
    G33[idx] *= shw3[j];
  }
}

/*
 * dxyz kernel — compute dxdr, dydr, dzdr, etc. via matrix multiply
 */

template <int LX>
void coef_generate_dxyz_impl(device float *dxdr,
                             device float *dydr,
                             device float *dzdr,
                             device float *dxds,
                             device float *dyds,
                             device float *dzds,
                             device float *dxdt,
                             device float *dydt,
                             device float *dzdt,
                             device const float *dx,
                             device const float *dy,
                             device const float *dz,
                             device const float *x,
                             device const float *y,
                             device const float *z,
                             threadgroup float *shdx,
                             threadgroup float *shdy,
                             threadgroup float *shdz,
                             threadgroup float *shu,
                             uint eid,
                             uint tid,
                             uint tpg) {

  constexpr int N3 = LX * LX * LX;
  const int e = eid;

  if (tid < (uint)(LX * LX)) {
    shdx[tid] = dx[tid];
    shdy[tid] = dy[tid];
    shdz[tid] = dz[tid];
  }

  /* --- x coordinate --- */
  for (int j = tid; j < N3; j += tpg) {
    shu[j] = x[j + e * N3];
  }
  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int j = tid; j < N3; j += tpg) {
    const int jk = j / LX;
    const int ii = j - jk * LX;
    const int kk = jk / LX;
    const int jj = jk - kk * LX;
    if (ii < LX && jj < LX && kk < LX) {
      float rtmp = 0.0f, stmp = 0.0f, ttmp = 0.0f;
      for (int l = 0; l < LX; l++) {
        rtmp += shdx[ii + l * LX] * shu[l + jj * LX + kk * LX * LX];
        stmp += shdy[jj + l * LX] * shu[ii + l * LX + kk * LX * LX];
        ttmp += shdz[kk + l * LX] * shu[ii + jj * LX + l * LX * LX];
      }
      dxdr[j + e * N3] = rtmp;
      dxds[j + e * N3] = stmp;
      dxdt[j + e * N3] = ttmp;
    }
  }
  threadgroup_barrier(mem_flags::mem_threadgroup);

  /* --- y coordinate --- */
  for (int j = tid; j < N3; j += tpg) {
    shu[j] = y[j + e * N3];
  }
  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int j = tid; j < N3; j += tpg) {
    const int jk = j / LX;
    const int ii = j - jk * LX;
    const int kk = jk / LX;
    const int jj = jk - kk * LX;
    if (ii < LX && jj < LX && kk < LX) {
      float rtmp = 0.0f, stmp = 0.0f, ttmp = 0.0f;
      for (int l = 0; l < LX; l++) {
        rtmp += shdx[ii + l * LX] * shu[l + jj * LX + kk * LX * LX];
        stmp += shdy[jj + l * LX] * shu[ii + l * LX + kk * LX * LX];
        ttmp += shdz[kk + l * LX] * shu[ii + jj * LX + l * LX * LX];
      }
      dydr[j + e * N3] = rtmp;
      dyds[j + e * N3] = stmp;
      dydt[j + e * N3] = ttmp;
    }
  }
  threadgroup_barrier(mem_flags::mem_threadgroup);

  /* --- z coordinate --- */
  for (int j = tid; j < N3; j += tpg) {
    shu[j] = z[j + e * N3];
  }
  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int j = tid; j < N3; j += tpg) {
    const int jk = j / LX;
    const int ii = j - jk * LX;
    const int kk = jk / LX;
    const int jj = jk - kk * LX;
    if (ii < LX && jj < LX && kk < LX) {
      float rtmp = 0.0f, stmp = 0.0f, ttmp = 0.0f;
      for (int l = 0; l < LX; l++) {
        rtmp += shdx[ii + l * LX] * shu[l + jj * LX + kk * LX * LX];
        stmp += shdy[jj + l * LX] * shu[ii + l * LX + kk * LX * LX];
        ttmp += shdz[kk + l * LX] * shu[ii + jj * LX + l * LX * LX];
      }
      dzdr[j + e * N3] = rtmp;
      dzds[j + e * N3] = stmp;
      dzdt[j + e * N3] = ttmp;
    }
  }
}

/*
 * drst kernel — compute Jacobian and inverse mapping
 */

kernel void coef_generate_drst_kernel(device float *jac[[ buffer(0) ]],
                                      device float *jacinv[[ buffer(1) ]],
                                      device float *drdx[[ buffer(2) ]],
                                      device float *drdy[[ buffer(3) ]],
                                      device float *drdz[[ buffer(4) ]],
                                      device float *dsdx[[ buffer(5) ]],
                                      device float *dsdy[[ buffer(6) ]],
                                      device float *dsdz[[ buffer(7) ]],
                                      device float *dtdx[[ buffer(8) ]],
                                      device float *dtdy[[ buffer(9) ]],
                                      device float *dtdz[[ buffer(10) ]],
                                      device const float *dxdr[[ buffer(11) ]],
                                      device const float *dydr[[ buffer(12) ]],
                                      device const float *dzdr[[ buffer(13) ]],
                                      device const float *dxds[[ buffer(14) ]],
                                      device const float *dyds[[ buffer(15) ]],
                                      device const float *dzds[[ buffer(16) ]],
                                      device const float *dxdt[[ buffer(17) ]],
                                      device const float *dydt[[ buffer(18) ]],
                                      device const float *dzdt[[ buffer(19) ]],
                                      constant int &n[[ buffer(20) ]],
                                      uint tid [[ thread_position_in_grid ]],
                                      uint stride [[ threads_per_grid ]]) {

  for (int i = tid; i < n; i += int(stride)) {
    jac[i] = (dxdr[i] * dyds[i] * dzdt[i])
           + (dxdt[i] * dydr[i] * dzds[i])
           + (dxds[i] * dydt[i] * dzdr[i])
           - (dxdr[i] * dydt[i] * dzds[i])
           - (dxds[i] * dydr[i] * dzdt[i])
           - (dxdt[i] * dyds[i] * dzdr[i]);
    jacinv[i] = 1.0f / jac[i];

    drdx[i] = dyds[i]*dzdt[i] - dydt[i]*dzds[i];
    drdy[i] = dxdt[i]*dzds[i] - dxds[i]*dzdt[i];
    drdz[i] = dxds[i]*dydt[i] - dxdt[i]*dyds[i];
    dsdx[i] = dydt[i]*dzdr[i] - dydr[i]*dzdt[i];
    dsdy[i] = dxdr[i]*dzdt[i] - dxdt[i]*dzdr[i];
    dsdz[i] = dxdt[i]*dydr[i] - dxdr[i]*dydt[i];
    dtdx[i] = dydr[i]*dzds[i] - dyds[i]*dzdr[i];
    dtdy[i] = dxds[i]*dzdr[i] - dxdr[i]*dzds[i];
    dtdz[i] = dxdr[i]*dyds[i] - dxds[i]*dydr[i];
  }
}

/*
 * Instantiation macros for geo and dxyz kernels (LX = 2..16)
 */

#define INSTANTIATE_GEO(LX)                                                     \
kernel void coef_generate_geo_kernel_lx##LX(                                    \
    device float *G11[[ buffer(0) ]],                                           \
    device float *G12[[ buffer(1) ]],                                           \
    device float *G13[[ buffer(2) ]],                                           \
    device float *G22[[ buffer(3) ]],                                           \
    device float *G23[[ buffer(4) ]],                                           \
    device float *G33[[ buffer(5) ]],                                           \
    device const float *drdx[[ buffer(6) ]],                                    \
    device const float *drdy[[ buffer(7) ]],                                    \
    device const float *drdz[[ buffer(8) ]],                                    \
    device const float *dsdx[[ buffer(9) ]],                                    \
    device const float *dsdy[[ buffer(10) ]],                                   \
    device const float *dsdz[[ buffer(11) ]],                                   \
    device const float *dtdx[[ buffer(12) ]],                                   \
    device const float *dtdy[[ buffer(13) ]],                                   \
    device const float *dtdz[[ buffer(14) ]],                                   \
    device const float *jacinv[[ buffer(15) ]],                                 \
    device const float *w3[[ buffer(16) ]],                                     \
    constant int &gdim[[ buffer(17) ]],                                         \
    uint eid [[ threadgroup_position_in_grid ]],                                \
    uint tid [[ thread_index_in_threadgroup ]],                                 \
    uint tpg [[ threads_per_threadgroup ]]) {                                   \
  threadgroup float shw3[LX * LX * LX];                                         \
  coef_generate_geo_impl<LX>(G11, G12, G13, G22, G23, G33,                     \
                             drdx, drdy, drdz, dsdx, dsdy, dsdz,               \
                             dtdx, dtdy, dtdz, jacinv, w3, gdim,               \
                             shw3, eid, tid, tpg);                             \
}

#define INSTANTIATE_DXYZ(LX)                                                    \
kernel void coef_generate_dxyz_kernel_lx##LX(                                   \
    device float *dxdr[[ buffer(0) ]],                                          \
    device float *dydr[[ buffer(1) ]],                                          \
    device float *dzdr[[ buffer(2) ]],                                          \
    device float *dxds[[ buffer(3) ]],                                          \
    device float *dyds[[ buffer(4) ]],                                          \
    device float *dzds[[ buffer(5) ]],                                          \
    device float *dxdt[[ buffer(6) ]],                                          \
    device float *dydt[[ buffer(7) ]],                                          \
    device float *dzdt[[ buffer(8) ]],                                          \
    device const float *dx[[ buffer(9) ]],                                      \
    device const float *dy[[ buffer(10) ]],                                     \
    device const float *dz[[ buffer(11) ]],                                     \
    device const float *x[[ buffer(12) ]],                                      \
    device const float *y[[ buffer(13) ]],                                      \
    device const float *z[[ buffer(14) ]],                                      \
    uint eid [[ threadgroup_position_in_grid ]],                                \
    uint tid [[ thread_index_in_threadgroup ]],                                 \
    uint tpg [[ threads_per_threadgroup ]]) {                                   \
  threadgroup float shdx[LX * LX];                                              \
  threadgroup float shdy[LX * LX];                                              \
  threadgroup float shdz[LX * LX];                                              \
  threadgroup float shu[LX * LX * LX];                                          \
  coef_generate_dxyz_impl<LX>(dxdr, dydr, dzdr, dxds, dyds, dzds,              \
                              dxdt, dydt, dzdt, dx, dy, dz, x, y, z,           \
                              shdx, shdy, shdz, shu,                           \
                              eid, tid, tpg);                                   \
}

INSTANTIATE_GEO(2)
INSTANTIATE_GEO(3)
INSTANTIATE_GEO(4)
INSTANTIATE_GEO(5)
INSTANTIATE_GEO(6)
INSTANTIATE_GEO(7)
INSTANTIATE_GEO(8)
INSTANTIATE_GEO(9)
INSTANTIATE_GEO(10)
INSTANTIATE_GEO(11)
INSTANTIATE_GEO(12)
INSTANTIATE_GEO(13)
INSTANTIATE_GEO(14)
INSTANTIATE_GEO(15)
INSTANTIATE_GEO(16)

INSTANTIATE_DXYZ(2)
INSTANTIATE_DXYZ(3)
INSTANTIATE_DXYZ(4)
INSTANTIATE_DXYZ(5)
INSTANTIATE_DXYZ(6)
INSTANTIATE_DXYZ(7)
INSTANTIATE_DXYZ(8)
INSTANTIATE_DXYZ(9)
INSTANTIATE_DXYZ(10)
INSTANTIATE_DXYZ(11)
INSTANTIATE_DXYZ(12)
INSTANTIATE_DXYZ(13)
INSTANTIATE_DXYZ(14)
INSTANTIATE_DXYZ(15)
INSTANTIATE_DXYZ(16)

/*
 * Mass matrix kernel — B = Binv = jac * w3
 */

kernel void coef_generate_mass_kernel(device float *B[[ buffer(0) ]],
                                      device float *Binv[[ buffer(1) ]],
                                      device const float *jac[[ buffer(2) ]],
                                      device const float *w3[[ buffer(3) ]],
                                      constant int &lxyz[[ buffer(4) ]],
                                      constant int &nel[[ buffer(5) ]],
                                      uint tid [[ thread_position_in_grid ]],
                                      uint stride [[ threads_per_grid ]]) {

  const int n = lxyz * nel;
  for (int idx = int(tid); idx < n; idx += int(stride)) {
    const int local_idx = idx - (idx / lxyz) * lxyz;
    const float mass_val = jac[idx] * w3[local_idx];

    B[idx] = mass_val;
    Binv[idx] = mass_val;
  }
}

/*
 * Facet area and unit normal kernel
 */

template <int LX>
void coef_generate_area_and_normal_impl(device float *area,
                                        device float *nx,
                                        device float *ny,
                                        device float *nz,
                                        device const float *dxdr,
                                        device const float *dydr,
                                        device const float *dzdr,
                                        device const float *dxds,
                                        device const float *dyds,
                                        device const float *dzds,
                                        device const float *dxdt,
                                        device const float *dydt,
                                        device const float *dzdt,
                                        device const float *wx,
                                        device const float *wy,
                                        device const float *wz,
                                        const float eps,
                                        uint eid, uint tid, uint tpg) {

  const int e = int(eid);
  const int lxyz = LX * LX * LX;
  const int lxy = LX * LX;

  for (int ijk = int(tid); ijk < lxyz; ijk += int(tpg)) {

    const int jk = ijk / LX;
    const int i = ijk - jk * LX;
    const int k = jk / LX;
    const int j = jk - k * LX;

    const int offset = ijk + (e * lxyz);
    const int face_offset = e * lxy * 6;

    int f, out_idx;
    float tx, ty, tz, dot, weight, length, sgn;

    /* ds x dt */
    if (i == 0 || i == LX - 1) {
      tx = dyds[offset] * dzdt[offset] - dzds[offset] * dydt[offset];
      ty = dzds[offset] * dxdt[offset] - dxds[offset] * dzdt[offset];
      tz = dxds[offset] * dydt[offset] - dyds[offset] * dxdt[offset];

      dot = tx*tx + ty*ty + tz*tz;
      length = sqrt(dot);
      weight = wy[j] * wz[k];

      if (i == 0) {
        f = 0;
        sgn = -1.0f;
      } else {
        f = 1;
        sgn = 1.0f;
      }

      out_idx = j + (k * LX) + (f * lxy) + face_offset;

      area[out_idx] = length * weight;

      if (length > eps) {
        nx[out_idx] = (tx / length) * sgn;
        ny[out_idx] = (ty / length) * sgn;
        nz[out_idx] = (tz / length) * sgn;
      } else {
        nx[out_idx] = tx * sgn;
        ny[out_idx] = ty * sgn;
        nz[out_idx] = tz * sgn;
      }
    }

    /* dr x dt */
    if (j == 0 || j == LX - 1) {
      tx = dydr[offset] * dzdt[offset] - dzdr[offset] * dydt[offset];
      ty = dzdr[offset] * dxdt[offset] - dxdr[offset] * dzdt[offset];
      tz = dxdr[offset] * dydt[offset] - dydr[offset] * dxdt[offset];

      dot = tx*tx + ty*ty + tz*tz;
      length = sqrt(dot);
      weight = wx[i] * wz[k];

      if (j == 0) {
        f = 2;
        sgn = 1.0f;
      } else {
        f = 3;
        sgn = -1.0f;
      }

      out_idx = i + (k * LX) + (f * lxy) + face_offset;

      area[out_idx] = length * weight;

      if (length > eps) {
        nx[out_idx] = (tx / length) * sgn;
        ny[out_idx] = (ty / length) * sgn;
        nz[out_idx] = (tz / length) * sgn;
      } else {
        nx[out_idx] = tx * sgn;
        ny[out_idx] = ty * sgn;
        nz[out_idx] = tz * sgn;
      }
    }

    /* dr x ds */
    if (k == 0 || k == LX - 1) {
      tx = dydr[offset] * dzds[offset] - dzdr[offset] * dyds[offset];
      ty = dzdr[offset] * dxds[offset] - dxdr[offset] * dzds[offset];
      tz = dxdr[offset] * dyds[offset] - dydr[offset] * dxds[offset];

      dot = tx*tx + ty*ty + tz*tz;
      length = sqrt(dot);
      weight = wx[i] * wy[j];

      if (k == 0) {
        f = 4;
        sgn = -1.0f;
      } else {
        f = 5;
        sgn = 1.0f;
      }

      out_idx = i + (j * LX) + (f * lxy) + face_offset;

      area[out_idx] = length * weight;

      if (length > eps) {
        nx[out_idx] = (tx / length) * sgn;
        ny[out_idx] = (ty / length) * sgn;
        nz[out_idx] = (tz / length) * sgn;
      } else {
        nx[out_idx] = tx * sgn;
        ny[out_idx] = ty * sgn;
        nz[out_idx] = tz * sgn;
      }
    }
  }
}

#define INSTANTIATE_AREA_NORMAL(LX)                                             \
kernel void coef_generate_area_and_normal_kernel_lx##LX(                        \
    device float *area[[ buffer(0) ]],                                          \
    device float *nx[[ buffer(1) ]],                                            \
    device float *ny[[ buffer(2) ]],                                            \
    device float *nz[[ buffer(3) ]],                                            \
    device const float *dxdr[[ buffer(4) ]],                                    \
    device const float *dydr[[ buffer(5) ]],                                    \
    device const float *dzdr[[ buffer(6) ]],                                    \
    device const float *dxds[[ buffer(7) ]],                                    \
    device const float *dyds[[ buffer(8) ]],                                    \
    device const float *dzds[[ buffer(9) ]],                                    \
    device const float *dxdt[[ buffer(10) ]],                                   \
    device const float *dydt[[ buffer(11) ]],                                   \
    device const float *dzdt[[ buffer(12) ]],                                   \
    device const float *wx[[ buffer(13) ]],                                     \
    device const float *wy[[ buffer(14) ]],                                     \
    device const float *wz[[ buffer(15) ]],                                     \
    constant float &eps[[ buffer(16) ]],                                        \
    uint eid [[ threadgroup_position_in_grid ]],                                \
    uint tid [[ thread_index_in_threadgroup ]],                                 \
    uint tpg [[ threads_per_threadgroup ]]) {                                   \
  coef_generate_area_and_normal_impl<LX>(area, nx, ny, nz,                     \
                                         dxdr, dydr, dzdr,                     \
                                         dxds, dyds, dzds,                     \
                                         dxdt, dydt, dzdt,                     \
                                         wx, wy, wz, eps, eid, tid, tpg);      \
}

INSTANTIATE_AREA_NORMAL(2)
INSTANTIATE_AREA_NORMAL(3)
INSTANTIATE_AREA_NORMAL(4)
INSTANTIATE_AREA_NORMAL(5)
INSTANTIATE_AREA_NORMAL(6)
INSTANTIATE_AREA_NORMAL(7)
INSTANTIATE_AREA_NORMAL(8)
INSTANTIATE_AREA_NORMAL(9)
INSTANTIATE_AREA_NORMAL(10)
INSTANTIATE_AREA_NORMAL(11)
INSTANTIATE_AREA_NORMAL(12)
INSTANTIATE_AREA_NORMAL(13)
INSTANTIATE_AREA_NORMAL(14)
INSTANTIATE_AREA_NORMAL(15)
INSTANTIATE_AREA_NORMAL(16)
