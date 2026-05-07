#ifndef __MATH_CONV1_KERNEL_CL__
#define __MATH_CONV1_KERNEL_CL__
/*
 Copyright (c) 2021-2025, The Neko Authors
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
 * Device kernel for convective terms
 */

#define DEFINE_CONV1_KERNEL(LX, CHUNKS)                                        \
__kernel void conv1_kernel_lx##LX(__global real * __restrict__ du,             \
                                  __global const real * __restrict__ u,        \
                                  __global const real * __restrict__ vx,       \
                                  __global const real * __restrict__ vy,       \
                                  __global const real * __restrict__ vz,       \
                                  __global const real * __restrict__ dx,       \
                                  __global const real * __restrict__ dy,       \
                                  __global const real * __restrict__ dz,       \
                                  __global const real * __restrict__ drdx,     \
                                  __global const real * __restrict__ dsdx,     \
                                  __global const real * __restrict__ dtdx,     \
                                  __global const real * __restrict__ drdy,     \
                                  __global const real * __restrict__ dsdy,     \
                                  __global const real * __restrict__ dtdy,     \
                                  __global const real * __restrict__ drdz,     \
                                  __global const real * __restrict__ dsdz,     \
                                  __global const real * __restrict__ dtdz,     \
                                  __global const real * __restrict__ jacinv) { \
                                                                               \
  __local real shu[LX * LX * LX];                                              \
                                                                               \
  __local real shvx[LX * LX * LX];                                             \
  __local real shvy[LX * LX * LX];                                             \
  __local real shvz[LX * LX * LX];                                             \
                                                                               \
  __local real shdx[LX * LX];                                                  \
  __local real shdy[LX * LX];                                                  \
  __local real shdz[LX * LX];                                                  \
                                                                               \
  __local real shjacinv[LX * LX * LX];                                         \
                                                                               \
                                                                               \
                                                                               \
  const int e = get_group_id(0);                                               \
  const int iii = get_local_id(0);                                             \
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;                         \
  const int ele = e*LX*LX*LX;                                                  \
                                                                               \
  if (iii < (LX * LX)) {                                                       \
    shdx[iii] = dx[iii];                                                       \
    shdy[iii] = dy[iii];                                                       \
    shdz[iii] = dz[iii];                                                       \
  }                                                                            \
                                                                               \
  int l = iii;                                                                 \
  while(l < (LX * LX * LX)) {                                                  \
    shu[l] = u[l + ele];                                                       \
                                                                               \
    shvx[l] = vx[l + ele];                                                     \
    shvy[l] = vy[l + ele];                                                     \
    shvz[l] = vz[l + ele];                                                     \
                                                                               \
    shjacinv[l] = jacinv[l + ele];                                             \
                                                                               \
    l = l + CHUNKS;                                                            \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for (int n = 0; n < nchunks; n++) {                                          \
    const int ijk = iii + n * CHUNKS;                                          \
    const int jk = ijk / LX;                                                   \
    const int i = ijk - jk * LX;                                               \
    const int k = jk / LX;                                                     \
    const int j = jk - k * LX;                                                 \
    if ( i < LX && j < LX && k < LX) {                                         \
      real rtmp = 0.0;                                                         \
      real stmp = 0.0;                                                         \
      real ttmp = 0.0;                                                         \
      for (int l = 0; l < LX; l++) {                                           \
        rtmp += shdx[i + l * LX] * shu[l + j * LX + k * LX * LX];              \
        stmp += shdy[j + l * LX] * shu[i + l * LX + k * LX * LX];              \
        ttmp += shdz[k + l * LX] * shu[i + j * LX + l * LX * LX];              \
      }                                                                        \
                                                                               \
      du[ijk + e * LX * LX * LX] = shjacinv[ijk] *                             \
        (shvx[ijk] * (drdx[ijk + ele] * rtmp                                   \
                      + dsdx[ijk + ele] * stmp                                 \
                      + dtdx[ijk + ele] * ttmp)                                \
         + shvy[ijk] * (drdy[ijk + ele] * rtmp                                 \
                        + dsdy[ijk + ele] * stmp                               \
                        + dtdy[ijk + ele] * ttmp)                              \
         + shvz[ijk] * (drdz[ijk + ele] * rtmp                                 \
                        + dsdz[ijk + ele] * stmp                               \
                        + dtdz[ijk + ele] * ttmp));                            \
    }                                                                          \
  }                                                                            \
}

DEFINE_CONV1_KERNEL(2, 256)
DEFINE_CONV1_KERNEL(3, 256)
DEFINE_CONV1_KERNEL(4, 256)
DEFINE_CONV1_KERNEL(5, 256)
DEFINE_CONV1_KERNEL(6, 256)
DEFINE_CONV1_KERNEL(7, 256)
DEFINE_CONV1_KERNEL(8, 256)
DEFINE_CONV1_KERNEL(9, 256)
DEFINE_CONV1_KERNEL(10, 256)
DEFINE_CONV1_KERNEL(11, 256)

#define DEFINE_CONV1_KERNEL_KSTEP(LX)                                          \
__kernel                                                                       \
void conv1_kernel_kstep_lx##LX(__global real * __restrict__ du,                \
                               __global const real * __restrict__ u,           \
                               __global const real * __restrict__ vx,          \
                               __global const real * __restrict__ vy,          \
                               __global const real * __restrict__ vz,          \
                               __global const real * __restrict__ dx,          \
                               __global const real * __restrict__ dy,          \
                               __global const real * __restrict__ dz,          \
                               __global const real * __restrict__ drdx,        \
                               __global const real * __restrict__ dsdx,        \
                               __global const real * __restrict__ dtdx,        \
                               __global const real * __restrict__ drdy,        \
                               __global const real * __restrict__ dsdy,        \
                               __global const real * __restrict__ dtdy,        \
                               __global const real * __restrict__ drdz,        \
                               __global const real * __restrict__ dsdz,        \
                               __global const real * __restrict__ dtdz,        \
                               __global const real * __restrict__ jacinv) {    \
                                                                               \
  __local real shu[LX * LX];                                                   \
                                                                               \
  __local real shdx[LX*LX];                                                    \
  __local real shdy[LX*LX];                                                    \
  __local real shdz[LX*LX];                                                    \
                                                                               \
  const int e = get_group_id(0);                                               \
  const int j = get_local_id(1);                                               \
  const int i = get_local_id(0);                                               \
  const int ij = i + j * LX;                                                   \
  const int ele = e*LX*LX*LX;                                                  \
                                                                               \
  shdx[ij] = dx[ij];                                                           \
  shdy[ij] = dy[ij];                                                           \
  shdz[ij] = dz[ij];                                                           \
                                                                               \
  real ru[LX];                                                                 \
  real rvx[LX];                                                                \
  real rvy[LX];                                                                \
  real rvz[LX];                                                                \
  real rjacinv[LX];                                                            \
                                                                               \
  for (int k = 0; k < LX; ++k) {                                               \
    ru[k] = u[ij + k*LX*LX + ele];                                             \
    rvx[k] = vx[ij + k*LX*LX + ele];                                           \
    rvy[k] = vy[ij + k*LX*LX + ele];                                           \
    rvz[k] = vz[ij + k*LX*LX + ele];                                           \
    rjacinv[k] = jacinv[ij + k*LX*LX + ele];                                   \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for (int k = 0; k < LX; ++k) {                                               \
    const int ijk = ij + k*LX*LX;                                              \
    real ttmp = 0.0;                                                           \
    shu[ij] = ru[k];                                                           \
    for (int l = 0; l < LX; l++) {                                             \
      ttmp += shdz[k+l*LX] * ru[l];                                            \
    }                                                                          \
    barrier(CLK_LOCAL_MEM_FENCE);                                              \
                                                                               \
    real rtmp = 0.0;                                                           \
    real stmp = 0.0;                                                           \
                                                                               \
    for (int l = 0; l < LX; l++) {                                             \
      rtmp += shdx[i+l*LX] * shu[l+j*LX];                                      \
      stmp += shdy[j+l*LX] * shu[i+l*LX];                                      \
    }                                                                          \
                                                                               \
    du[ijk + ele] = rjacinv[k] *                                               \
	(rvx[k] * (drdx[ijk + ele] * rtmp                                            \
                   + dsdx[ijk + ele] * stmp                                    \
                   + dtdx[ijk + ele] * ttmp)                                   \
	 + rvy[k] * (drdy[ijk + ele] * rtmp                                          \
                     + dsdy[ijk + ele] * stmp                                  \
                     + dtdy[ijk + ele] * ttmp)                                 \
	 + rvz[k] * (drdz[ijk + ele] * rtmp                                          \
                     + dsdz[ijk + ele] * stmp                                  \
                     + dtdz[ijk + ele] * ttmp));                               \
    barrier(CLK_LOCAL_MEM_FENCE);                                              \
  }                                                                            \
}

DEFINE_CONV1_KERNEL_KSTEP(2)
DEFINE_CONV1_KERNEL_KSTEP(3)
DEFINE_CONV1_KERNEL_KSTEP(4)
DEFINE_CONV1_KERNEL_KSTEP(5)
DEFINE_CONV1_KERNEL_KSTEP(6)
DEFINE_CONV1_KERNEL_KSTEP(7)
DEFINE_CONV1_KERNEL_KSTEP(8)
DEFINE_CONV1_KERNEL_KSTEP(9)
DEFINE_CONV1_KERNEL_KSTEP(10)
DEFINE_CONV1_KERNEL_KSTEP(11)
DEFINE_CONV1_KERNEL_KSTEP(12)
DEFINE_CONV1_KERNEL_KSTEP(13)
DEFINE_CONV1_KERNEL_KSTEP(14)
DEFINE_CONV1_KERNEL_KSTEP(15)
DEFINE_CONV1_KERNEL_KSTEP(16)

#endif // __MATH_CONV1_KERNEL_CL__
