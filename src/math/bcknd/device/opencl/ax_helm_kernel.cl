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
 * Device kernel for Ax helm
 */

#define DEFINE_AX_HELM_KERNEL(LX, CHUNKS)                                      \
__kernel void ax_helm_kernel_lx##LX(__global real * __restrict__ w,            \
                                    __global const real * __restrict__ u,      \
                                    __global const real * __restrict__ dx,     \
                                    __global const real * __restrict__ dy,     \
                                    __global const real * __restrict__ dz,     \
                                    __global const real * __restrict__ dxt,    \
                                    __global const real * __restrict__ dyt,    \
                                    __global const real * __restrict__ dzt,    \
                                    __global const real * __restrict__ h1,     \
                                    __global const real * __restrict__ g11,    \
                                    __global const real * __restrict__ g22,    \
                                    __global const real * __restrict__ g33,    \
                                    __global const real * __restrict__ g12,    \
                                    __global const real * __restrict__ g13,    \
                                    __global const real * __restrict__ g23) {  \
                                                                               \
  __local real shdx[LX*LX];                                                    \
  __local real shdy[LX*LX];                                                    \
  __local real shdzt[LX*LX];                                                   \
                                                                               \
  __local real shdxt[LX*LX];                                                   \
  __local real shdyt[LX*LX];                                                   \
  __local real shdz[LX*LX];                                                    \
                                                                               \
  __local real shu[LX*LX*LX];                                                  \
  __local real shur[LX*LX*LX];                                                 \
  __local real shus[LX*LX*LX];                                                 \
  __local real shut[LX*LX*LX];                                                 \
                                                                               \
  int l,i,j,k,n;                                                               \
                                                                               \
  const int e = get_group_id(0);                                               \
  const int iii = get_local_id(0);                                             \
  const int nchunks = (LX * LX * LX - 1)/CHUNKS + 1;                           \
                                                                               \
  if (iii<LX*LX) {                                                             \
    shdx[iii] = dx[iii];                                                       \
    shdy[iii] = dy[iii];                                                       \
    shdz[iii] = dz[iii];                                                       \
  }                                                                            \
  i = iii;                                                                     \
  while (i < LX * LX * LX){                                                    \
    shu[i] = u[i+e*LX*LX*LX];                                                  \
    i = i + CHUNKS;                                                            \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  if (iii<LX*LX){                                                              \
    shdxt[iii] = dxt[iii];                                                     \
    shdyt[iii] = dyt[iii];                                                     \
    shdzt[iii] = dzt[iii];                                                     \
  }                                                                            \
                                                                               \
  for (n=0; n<nchunks; n++){                                                   \
    const int ijk = iii+n*CHUNKS;                                              \
    const int jk = ijk/LX;                                                     \
    i = ijk-jk*LX;                                                             \
    k = jk/LX;                                                                 \
    j = jk-k*LX;                                                               \
    if (i<LX && j<LX && k<LX && ijk < LX*LX*LX) {                              \
      real rtmp = 0.0;                                                         \
      real stmp = 0.0;                                                         \
      real ttmp = 0.0;                                                         \
      for (l = 0; l<LX; l++){                                                  \
        rtmp = rtmp + shdx[i+l*LX] * shu[l+j*LX+k*LX*LX];                      \
        stmp = stmp + shdy[j+l*LX] * shu[i+l*LX+k*LX*LX];                      \
        ttmp = ttmp + shdz[k+l*LX] * shu[i+j*LX+l*LX*LX];                      \
      }                                                                        \
      shur[ijk] = h1[ijk+e*LX*LX*LX]                                           \
                * (g11[ijk+e*LX*LX*LX] * rtmp                                  \
                   + g12[ijk+e*LX*LX*LX] * stmp                                \
                   + g13[ijk+e*LX*LX*LX] * ttmp);                              \
      shus[ijk] = h1[ijk+e*LX*LX*LX]                                           \
                * (g12[ijk+e*LX*LX*LX] * rtmp                                  \
                   + g22[ijk+e*LX*LX*LX] * stmp                                \
                   + g23[ijk+e*LX*LX*LX] * ttmp);                              \
      shut[ijk] = h1[ijk+e*LX*LX*LX]                                           \
                * (g13[ijk+e*LX*LX*LX] * rtmp                                  \
                   + g23[ijk+e*LX*LX*LX] * stmp                                \
                   + g33[ijk+e*LX*LX*LX] * ttmp);                              \
    }                                                                          \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for (n=0; n<nchunks; n++){                                                   \
    const int ijk = iii+n*CHUNKS;                                              \
    const int jk = ijk/LX;                                                     \
    i = ijk-jk*LX;                                                             \
    k = jk/LX;                                                                 \
    j = jk-k*LX;                                                               \
    if (i<LX && j<LX && k<LX && ijk <LX*LX*LX) {                               \
      real wijke = 0.0;                                                        \
      for (l = 0; l<LX; l++){                                                  \
        wijke = wijke                                                          \
              + shdxt[i+l*LX] * shur[l+j*LX+k*LX*LX]                           \
              + shdyt[j+l*LX] * shus[i+l*LX+k*LX*LX]                           \
              + shdzt[k+l*LX] * shut[i+j*LX+l*LX*LX];                          \
      }                                                                        \
      w[ijk+e*LX*LX*LX] = wijke;                                               \
    }                                                                          \
  }                                                                            \
}

DEFINE_AX_HELM_KERNEL(1, 256)
DEFINE_AX_HELM_KERNEL(2, 256)
DEFINE_AX_HELM_KERNEL(3, 256)
DEFINE_AX_HELM_KERNEL(4, 256)
DEFINE_AX_HELM_KERNEL(5, 256)
DEFINE_AX_HELM_KERNEL(6, 256)
DEFINE_AX_HELM_KERNEL(7, 256)
DEFINE_AX_HELM_KERNEL(8, 256)
DEFINE_AX_HELM_KERNEL(9, 256)

