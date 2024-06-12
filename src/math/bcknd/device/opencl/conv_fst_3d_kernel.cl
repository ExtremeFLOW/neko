#ifndef __MATH_CONV_FST_3D_KERNEL_CL__
#define __MATH_CONV_FST_3D_KERNEL_CL__
/*
 Copyright (c) 2021-2024, The Neko Authors
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
 * Device kernel for conv_fst_3d
 */

#define DEFINE_CONV_FST_3D_KERNEL(LX, CHUNKS)                                      \
__kernel void conv_fst_3d_kernel_lx##LX(__global real * __restrict__ ud,           \
                                        __global const real * __restrict__ u,      \
                                        __global const real * __restrict__ cr,     \
                                        __global const real * __restrict__ cs,     \
                                        __global const real * __restrict__ ct,     \
                                        __global const real * __restrict__ dx,     \
                                        __global const real * __restrict__ dy,     \
                                        __global const real * __restrict__ dz) {   \
                                                                                   \
  __local real shu[LX * LX * LX];                                                  \
                                                                                   \
  __local real shcr[LX * LX * LX];                                                 \
  __local real shcs[LX * LX * LX];                                                 \
  __local real shct[LX * LX * LX];                                                 \
                                                                                   \
  __local real shdx[LX * LX];                                                      \
  __local real shdy[LX * LX];                                                      \
  __local real shdz[LX * LX];                                                      \
                                                                                   \
  int i,j,k;                                                                       \
                                                                                   \
  const int e = get_group_id(0);                                                   \
  const int iii = get_local_id(0);                                                 \
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;                             \
                                                                                   \
  if (iii < (LX * LX)) {                                                           \
    shdx[iii] = dx[iii];                                                           \
    shdy[iii] = dy[iii];                                                           \
    shdz[iii] = dz[iii];                                                           \
  }                                                                                \
                                                                                   \
  j = iii;                                                                         \
  while(j < (LX * LX * LX)) {                                                      \
    shu[j] = u[j + e * LX * LX * LX];                                              \
                                                                                   \
    shcr[j] = cr[j + e * LX * LX * LX];                                            \
    shcs[j] = cs[j + e * LX * LX * LX];                                            \
    shct[j] = ct[j + e * LX * LX * LX];                                            \
                                                                                   \
    j = j + CHUNKS;                                                                \
  }                                                                                \
                                                                                   \
  barrier(CLK_LOCAL_MEM_FENCE);                                                    \
                                                                                   \
  for (int n = 0; n < nchunks; n++) {                                              \
    const int ijk = iii + n * CHUNKS;                                              \
    const int jk = ijk / LX;                                                       \
    i = ijk - jk * LX;                                                             \
    k = jk / LX;                                                                   \
    j = jk - k * LX;                                                               \
    if ( i < LX && j < LX && k < LX) {                                             \
      real rtmp = 0.0;                                                             \
      real stmp = 0.0;                                                             \
      real ttmp = 0.0;                                                             \
      for (int l = 0; l < LX; l++) {                                               \
        rtmp += shdx[i + l * LX] * shu[l + j * LX + k * LX * LX];                  \
        stmp += shdy[j + l * LX] * shu[i + l * LX + k * LX * LX];                  \
        ttmp += shdz[k + l * LX] * shu[i + j * LX + l * LX * LX];                  \
      }                                                                            \
                                                                                   \
      ud[ijk + e * LX * LX * LX] =                                                 \
          (shcr[ijk] * rtmp + shcs[ijk] * stmp + shct[ijk] * ttmp);                \
                                                                                   \
    }                                                                              \
  }                                                                                \
}        

DEFINE_CONV_FST_3D_KERNEL(18, 256)
DEFINE_CONV_FST_3D_KERNEL(17, 256)
DEFINE_CONV_FST_3D_KERNEL(16, 256)
DEFINE_CONV_FST_3D_KERNEL(15, 256)
DEFINE_CONV_FST_3D_KERNEL(14, 256)
DEFINE_CONV_FST_3D_KERNEL(13, 256)
DEFINE_CONV_FST_3D_KERNEL(12, 256)
DEFINE_CONV_FST_3D_KERNEL(11, 256)
DEFINE_CONV_FST_3D_KERNEL(10, 256)
DEFINE_CONV_FST_3D_KERNEL(9, 256)
DEFINE_CONV_FST_3D_KERNEL(8, 256)
DEFINE_CONV_FST_3D_KERNEL(7, 256)
DEFINE_CONV_FST_3D_KERNEL(6, 256)
DEFINE_CONV_FST_3D_KERNEL(5, 256)
DEFINE_CONV_FST_3D_KERNEL(4, 256)
DEFINE_CONV_FST_3D_KERNEL(3, 256)
DEFINE_CONV_FST_3D_KERNEL(2, 256)


#endif // __MATH_CONV_FST_3D_KERNEL_CL__