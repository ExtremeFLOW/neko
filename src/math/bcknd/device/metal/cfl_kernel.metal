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
 * Metal compute kernels for CFL computation.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 * @note Threadgroup memory (4 x LX^3 floats) limits LX to 12.
 */

#include <metal_stdlib>
using namespace metal;

#define CFL_CHUNKS 256

/*
 * CFL kernel — per-element max, one element per threadgroup
 */

template <int LX>
void cfl_impl(const float dt,
              device const float *u,
              device const float *v,
              device const float *w,
              device const float *drdx,
              device const float *dsdx,
              device const float *dtdx,
              device const float *drdy,
              device const float *dsdy,
              device const float *dtdy,
              device const float *drdz,
              device const float *dsdz,
              device const float *dtdz,
              device const float *dr_inv,
              device const float *ds_inv,
              device const float *dt_inv,
              device const float *jacinv,
              device float *cfl_h,
              threadgroup float *shu,
              threadgroup float *shv,
              threadgroup float *shw,
              threadgroup float *shjacinv,
              threadgroup float *shdr_inv,
              threadgroup float *shds_inv,
              threadgroup float *shdt_inv,
              threadgroup float *shcfl,
              uint eid, uint tid, uint tpg) {

  int i, j, k;

  const int e = int(eid);
  const int iii = int(tid);
  const int nchunks = (LX * LX * LX - 1) / CFL_CHUNKS + 1;

  if (iii < LX) {
    shdr_inv[iii] = dr_inv[iii];
    shds_inv[iii] = ds_inv[iii];
    shdt_inv[iii] = dt_inv[iii];
  }

  j = iii;
  while (j < (LX * LX * LX)) {
    shu[j] = u[j + e * LX * LX * LX];
    shv[j] = v[j + e * LX * LX * LX];
    shw[j] = w[j + e * LX * LX * LX];

    shjacinv[j] = jacinv[j + e * LX * LX * LX];

    j = j + CFL_CHUNKS;
  }

  threadgroup_barrier(mem_flags::mem_threadgroup);

  float cfl_tmp = 0.0f;
  for (int n = 0; n < nchunks; n++) {
    const int ijk = iii + n * CFL_CHUNKS;
    const int jk = ijk / LX;
    i = ijk - jk * LX;
    k = jk / LX;
    j = jk - k * LX;
    if (i < LX && j < LX && k < LX) {
      const float cflr = fabs(dt * ((shu[ijk] * drdx[ijk + e * LX * LX * LX]
                                     + shv[ijk] * drdy[ijk + e * LX * LX * LX]
                                     + shw[ijk] * drdz[ijk + e * LX * LX * LX]
                                     ) * shjacinv[ijk]) * shdr_inv[i]);
      const float cfls = fabs(dt * ((shu[ijk] * dsdx[ijk + e * LX * LX * LX]
                                     + shv[ijk] * dsdy[ijk + e * LX * LX * LX]
                                     + shw[ijk] * dsdz[ijk + e * LX * LX * LX]
                                     ) * shjacinv[ijk]) * shds_inv[j]);
      const float cflt = fabs(dt * ((shu[ijk] * dtdx[ijk + e * LX * LX * LX]
                                     + shv[ijk] * dtdy[ijk + e * LX * LX * LX]
                                     + shw[ijk] * dtdz[ijk + e * LX * LX * LX]
                                     ) * shjacinv[ijk]) * shdt_inv[k]);

      cfl_tmp = fmax(cflr + cfls + cflt, cfl_tmp);
    }
  }
  shcfl[iii] = cfl_tmp;

  threadgroup_barrier(mem_flags::mem_threadgroup);

  i = int(tpg) >> 1;
  while (i != 0) {
    if (iii < i) {
      shcfl[iii] = fmax(shcfl[iii], shcfl[iii + i]);
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    i = i >> 1;
  }

  if (iii == 0) {
    cfl_h[e] = shcfl[0];
  }
}

#define INSTANTIATE_CFL(LX)                                                     \
kernel void cfl_kernel_lx##LX(                                                  \
    constant float &dt[[ buffer(0) ]],                                          \
    device const float *u[[ buffer(1) ]],                                       \
    device const float *v[[ buffer(2) ]],                                       \
    device const float *w[[ buffer(3) ]],                                       \
    device const float *drdx[[ buffer(4) ]],                                    \
    device const float *dsdx[[ buffer(5) ]],                                    \
    device const float *dtdx[[ buffer(6) ]],                                    \
    device const float *drdy[[ buffer(7) ]],                                    \
    device const float *dsdy[[ buffer(8) ]],                                    \
    device const float *dtdy[[ buffer(9) ]],                                    \
    device const float *drdz[[ buffer(10) ]],                                   \
    device const float *dsdz[[ buffer(11) ]],                                   \
    device const float *dtdz[[ buffer(12) ]],                                   \
    device const float *dr_inv[[ buffer(13) ]],                                 \
    device const float *ds_inv[[ buffer(14) ]],                                 \
    device const float *dt_inv[[ buffer(15) ]],                                 \
    device const float *jacinv[[ buffer(16) ]],                                 \
    device float *cfl_h[[ buffer(17) ]],                                        \
    uint eid [[ threadgroup_position_in_grid ]],                                \
    uint tid [[ thread_index_in_threadgroup ]],                                 \
    uint tpg [[ threads_per_threadgroup ]]) {                                   \
  threadgroup float shu[LX * LX * LX];                                          \
  threadgroup float shv[LX * LX * LX];                                          \
  threadgroup float shw[LX * LX * LX];                                          \
  threadgroup float shjacinv[LX * LX * LX];                                     \
  threadgroup float shdr_inv[LX];                                               \
  threadgroup float shds_inv[LX];                                               \
  threadgroup float shdt_inv[LX];                                               \
  threadgroup float shcfl[CFL_CHUNKS];                                          \
  cfl_impl<LX>(dt, u, v, w, drdx, dsdx, dtdx, drdy, dsdy, dtdy,                \
               drdz, dsdz, dtdz, dr_inv, ds_inv, dt_inv, jacinv,               \
               cfl_h, shu, shv, shw, shjacinv,                                 \
               shdr_inv, shds_inv, shdt_inv, shcfl, eid, tid, tpg);            \
}

INSTANTIATE_CFL(2)
INSTANTIATE_CFL(3)
INSTANTIATE_CFL(4)
INSTANTIATE_CFL(5)
INSTANTIATE_CFL(6)
INSTANTIATE_CFL(7)
INSTANTIATE_CFL(8)
INSTANTIATE_CFL(9)
INSTANTIATE_CFL(10)
INSTANTIATE_CFL(11)
INSTANTIATE_CFL(12)
