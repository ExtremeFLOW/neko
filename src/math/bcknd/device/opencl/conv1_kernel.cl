#ifndef __MATH_CONV1_KERNEL_CL__
#define __MATH_CONV1_KERNEL_CL__
/*
 Copyright (c) 2021-2022, The Neko Authors
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
  int i,j,k;                                                                   \
                                                                               \
  const int e = get_group_id(0);                                               \
  const int iii = get_local_id(0);                                             \
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;                         \
                                                                               \
  if (iii < (LX * LX)) {                                                       \
    shdx[iii] = dx[iii];                                                       \
    shdy[iii] = dy[iii];                                                       \
    shdz[iii] = dz[iii];                                                       \
  }                                                                            \
                                                                               \
  j = iii;                                                                     \
  while(j < (LX * LX * LX)) {                                                  \
    shu[j] = u[j + e * LX * LX * LX];                                          \
                                                                               \
    shvx[j] = vx[j + e * LX * LX * LX];                                        \
    shvy[j] = vy[j + e * LX * LX * LX];                                        \
    shvz[j] = vz[j + e * LX * LX * LX];                                        \
                                                                               \
    shjacinv[j] = jacinv[j + e * LX * LX * LX];                                \
                                                                               \
    j = j + CHUNKS;                                                            \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for (int n = 0; n < nchunks; n++) {                                          \
    const int ijk = iii + n * CHUNKS;                                          \
    const int jk = ijk / LX;                                                   \
    i = ijk - jk * LX;                                                         \
    k = jk / LX;                                                               \
    j = jk - k * LX;                                                           \
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
        (shvx[ijk] * (drdx[ijk + e * LX * LX * LX] * rtmp                      \
                      + dsdx[ijk + e * LX * LX * LX] * stmp                    \
                      + dtdx[ijk + e * LX * LX * LX] * ttmp)                   \
         + shvy[ijk] * (drdy[ijk + e * LX * LX * LX] * rtmp                    \
                        + dsdy[ijk + e * LX * LX * LX] * stmp                  \
                        + dtdy[ijk + e * LX * LX * LX] * ttmp)                 \
         + shvz[ijk] * (drdz[ijk + e * LX * LX * LX] * rtmp                    \
                        + dsdz[ijk + e * LX * LX * LX] * stmp                  \
                        + dtdz[ijk + e * LX * LX * LX] * ttmp));               \
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


#endif // __MATH_CONV1_KERNEL_CL__
