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
 * Metal compute kernels for the scalar convective term (convect_scalar).
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/*
 * convect_scalar — kstep, one element per (LX x LX) threadgroup
 */

template <int LX>
void convect_scalar_kstep_impl(device float *du,
                               device const float *u,
                               device const float *cr,
                               device const float *cs,
                               device const float *ct,
                               device const float *dx,
                               device const float *dy,
                               device const float *dz,
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
  float rcr[LX];
  float rcs[LX];
  float rct[LX];

  for (int k = 0; k < LX; ++k) {
    ru[k] = u[ij + k * LX * LX + ele];
    rcr[k] = cr[ij + k * LX * LX + ele];
    rcs[k] = cs[ij + k * LX * LX + ele];
    rct[k] = ct[ij + k * LX * LX + ele];
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

    du[ijk + ele] = rcr[k] * rtmp + rcs[k] * stmp + rct[k] * ttmp;
    threadgroup_barrier(mem_flags::mem_threadgroup);
  }
}

#define INSTANTIATE_CONVECT_SCALAR(LX)                                          \
kernel void convect_scalar_kernel_lx##LX(                                      \
    device float *du[[ buffer(0) ]],                                           \
    device const float *u[[ buffer(1) ]],                                      \
    device const float *cr[[ buffer(2) ]],                                     \
    device const float *cs[[ buffer(3) ]],                                     \
    device const float *ct[[ buffer(4) ]],                                     \
    device const float *dx[[ buffer(5) ]],                                     \
    device const float *dy[[ buffer(6) ]],                                     \
    device const float *dz[[ buffer(7) ]],                                     \
    uint2 gid2 [[ threadgroup_position_in_grid ]],                             \
    uint2 tid [[ thread_position_in_threadgroup ]]) {                          \
  threadgroup float shu[LX * LX];                                              \
  threadgroup float shdx[LX * LX];                                             \
  threadgroup float shdy[LX * LX];                                             \
  threadgroup float shdz[LX * LX];                                             \
  convect_scalar_kstep_impl<LX>(du, u, cr, cs, ct, dx, dy, dz,                \
                                shu, shdx, shdy, shdz,                         \
                                gid2.x, tid.x, tid.y);                         \
}

INSTANTIATE_CONVECT_SCALAR(2)
INSTANTIATE_CONVECT_SCALAR(3)
INSTANTIATE_CONVECT_SCALAR(4)
INSTANTIATE_CONVECT_SCALAR(5)
INSTANTIATE_CONVECT_SCALAR(6)
INSTANTIATE_CONVECT_SCALAR(7)
INSTANTIATE_CONVECT_SCALAR(8)
INSTANTIATE_CONVECT_SCALAR(9)
INSTANTIATE_CONVECT_SCALAR(10)
INSTANTIATE_CONVECT_SCALAR(11)
INSTANTIATE_CONVECT_SCALAR(12)
INSTANTIATE_CONVECT_SCALAR(13)
INSTANTIATE_CONVECT_SCALAR(14)
INSTANTIATE_CONVECT_SCALAR(15)
INSTANTIATE_CONVECT_SCALAR(16)
