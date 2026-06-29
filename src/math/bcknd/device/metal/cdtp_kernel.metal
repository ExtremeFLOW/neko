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
 * Metal compute kernels for the transposed derivative operator (D^T x).
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/*
 * cdtp — kstep formulation, one element per (LX x LX) threadgroup
 */

template <int LX>
void cdtp_kstep_impl(device float *dtx,
                     device const float *x,
                     device const float *dr,
                     device const float *ds,
                     device const float *dt,
                     device const float *dxt,
                     device const float *dyt,
                     device const float *dzt,
                     device const float *w3,
                     threadgroup float *shtar,
                     threadgroup float *shtas,
                     threadgroup float *shdxt,
                     threadgroup float *shdyt,
                     threadgroup float *shdzt,
                     uint e, uint i, uint j) {

  const int ij = i + j * LX;
  const int ele = e * LX * LX * LX;

  shdxt[ij] = dxt[ij];
  shdyt[ij] = dyt[ij];
  shdzt[ij] = dzt[ij];

  float rtar[LX];
  float rtas[LX];
  float rtat[LX];

  for (int k = 0; k < LX; ++k) {
    const float wx = x[ij + k * LX * LX + ele] * w3[ij + k * LX * LX];
    rtar[k] = wx * dr[ij + k * LX * LX + ele];
    rtas[k] = wx * ds[ij + k * LX * LX + ele];
    rtat[k] = wx * dt[ij + k * LX * LX + ele];
  }

  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int k = 0; k < LX; ++k) {
    const int ijk = ij + k * LX * LX;
    float ttmp = 0.0f;
    shtar[ij] = rtar[k];
    shtas[ij] = rtas[k];
    for (int l = 0; l < LX; l++) {
      ttmp += shdzt[k + l * LX] * rtat[l];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    float rtmp = 0.0f;
    float stmp = 0.0f;
    for (int l = 0; l < LX; l++) {
      rtmp += shdxt[i + l * LX] * shtar[l + j * LX];
      stmp += shdyt[j + l * LX] * shtas[i + l * LX];
    }

    dtx[ijk + ele] = (rtmp + stmp + ttmp);

    threadgroup_barrier(mem_flags::mem_threadgroup);
  }
}

#define INSTANTIATE_CDTP(LX)                                                    \
kernel void cdtp_kernel_lx##LX(                                                \
    device float *dtx[[ buffer(0) ]],                                          \
    device const float *x[[ buffer(1) ]],                                      \
    device const float *dr[[ buffer(2) ]],                                     \
    device const float *ds[[ buffer(3) ]],                                     \
    device const float *dt[[ buffer(4) ]],                                     \
    device const float *dxt[[ buffer(5) ]],                                    \
    device const float *dyt[[ buffer(6) ]],                                    \
    device const float *dzt[[ buffer(7) ]],                                    \
    device const float *w3[[ buffer(8) ]],                                     \
    uint2 gid2 [[ threadgroup_position_in_grid ]],                             \
    uint2 tid [[ thread_position_in_threadgroup ]]) {                          \
  threadgroup float shtar[LX * LX];                                            \
  threadgroup float shtas[LX * LX];                                            \
  threadgroup float shdxt[LX * LX];                                            \
  threadgroup float shdyt[LX * LX];                                            \
  threadgroup float shdzt[LX * LX];                                            \
  cdtp_kstep_impl<LX>(dtx, x, dr, ds, dt, dxt, dyt, dzt, w3,                  \
                      shtar, shtas, shdxt, shdyt, shdzt,                       \
                      gid2.x, tid.x, tid.y);                                   \
}

INSTANTIATE_CDTP(2)
INSTANTIATE_CDTP(3)
INSTANTIATE_CDTP(4)
INSTANTIATE_CDTP(5)
INSTANTIATE_CDTP(6)
INSTANTIATE_CDTP(7)
INSTANTIATE_CDTP(8)
INSTANTIATE_CDTP(9)
INSTANTIATE_CDTP(10)
INSTANTIATE_CDTP(11)
INSTANTIATE_CDTP(12)
INSTANTIATE_CDTP(13)
INSTANTIATE_CDTP(14)
INSTANTIATE_CDTP(15)
INSTANTIATE_CDTP(16)
