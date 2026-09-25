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
 * Metal compute kernels for the convective term (conv1).
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/*
 * conv1 — kstep formulation, one element per (LX x LX) threadgroup
 */

template <int LX>
void conv1_kstep_impl(device float *du,
                      device const float *u,
                      device const float *vx,
                      device const float *vy,
                      device const float *vz,
                      device const float *dx,
                      device const float *dy,
                      device const float *dz,
                      device const float *drdx,
                      device const float *dsdx,
                      device const float *dtdx,
                      device const float *drdy,
                      device const float *dsdy,
                      device const float *dtdy,
                      device const float *drdz,
                      device const float *dsdz,
                      device const float *dtdz,
                      device const float *jacinv,
                      threadgroup float *shu,
                      threadgroup float *shdx,
                      threadgroup float *shdy,
                      threadgroup float *shdz,
                      uint e, uint i, uint j) {

  const int ij = i + j * LX;
  const int ele = e * LX * LX * LX;

  shdx[ij] = dx[ij];
  shdy[ij] = dy[ij];
  shdz[ij] = dz[ij];

  float ru[LX];
  float rvx[LX];
  float rvy[LX];
  float rvz[LX];
  float rjacinv[LX];

  for (int k = 0; k < LX; ++k) {
    ru[k] = u[ij + k * LX * LX + ele];
    rvx[k] = vx[ij + k * LX * LX + ele];
    rvy[k] = vy[ij + k * LX * LX + ele];
    rvz[k] = vz[ij + k * LX * LX + ele];
    rjacinv[k] = jacinv[ij + k * LX * LX + ele];
  }

  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int k = 0; k < LX; ++k) {
    const int ijk = ij + k * LX * LX;
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

    du[ijk + ele] = rjacinv[k] *
      (rvx[k] * (drdx[ijk + ele] * rtmp
                 + dsdx[ijk + ele] * stmp
                 + dtdx[ijk + ele] * ttmp)
       + rvy[k] * (drdy[ijk + ele] * rtmp
                   + dsdy[ijk + ele] * stmp
                   + dtdy[ijk + ele] * ttmp)
       + rvz[k] * (drdz[ijk + ele] * rtmp
                   + dsdz[ijk + ele] * stmp
                   + dtdz[ijk + ele] * ttmp));
    threadgroup_barrier(mem_flags::mem_threadgroup);
  }
}

#define INSTANTIATE_CONV1(LX)                                                   \
kernel void conv1_kernel_lx##LX(                                              \
    device float *du[[ buffer(0) ]],                                           \
    device const float *u[[ buffer(1) ]],                                      \
    device const float *vx[[ buffer(2) ]],                                     \
    device const float *vy[[ buffer(3) ]],                                     \
    device const float *vz[[ buffer(4) ]],                                     \
    device const float *dx[[ buffer(5) ]],                                     \
    device const float *dy[[ buffer(6) ]],                                     \
    device const float *dz[[ buffer(7) ]],                                     \
    device const float *drdx[[ buffer(8) ]],                                   \
    device const float *dsdx[[ buffer(9) ]],                                   \
    device const float *dtdx[[ buffer(10) ]],                                  \
    device const float *drdy[[ buffer(11) ]],                                  \
    device const float *dsdy[[ buffer(12) ]],                                  \
    device const float *dtdy[[ buffer(13) ]],                                  \
    device const float *drdz[[ buffer(14) ]],                                  \
    device const float *dsdz[[ buffer(15) ]],                                  \
    device const float *dtdz[[ buffer(16) ]],                                  \
    device const float *jacinv[[ buffer(17) ]],                                \
    uint2 gid2 [[ threadgroup_position_in_grid ]],                             \
    uint2 tid [[ thread_position_in_threadgroup ]]) {                          \
  threadgroup float shu[LX * LX];                                              \
  threadgroup float shdx[LX * LX];                                             \
  threadgroup float shdy[LX * LX];                                             \
  threadgroup float shdz[LX * LX];                                             \
  conv1_kstep_impl<LX>(du, u, vx, vy, vz, dx, dy, dz,                         \
                       drdx, dsdx, dtdx, drdy, dsdy, dtdy,                     \
                       drdz, dsdz, dtdz, jacinv,                               \
                       shu, shdx, shdy, shdz, gid2.x, tid.x, tid.y);          \
}

INSTANTIATE_CONV1(2)
INSTANTIATE_CONV1(3)
INSTANTIATE_CONV1(4)
INSTANTIATE_CONV1(5)
INSTANTIATE_CONV1(6)
INSTANTIATE_CONV1(7)
INSTANTIATE_CONV1(8)
INSTANTIATE_CONV1(9)
INSTANTIATE_CONV1(10)
INSTANTIATE_CONV1(11)
INSTANTIATE_CONV1(12)
INSTANTIATE_CONV1(13)
INSTANTIATE_CONV1(14)
INSTANTIATE_CONV1(15)
INSTANTIATE_CONV1(16)
