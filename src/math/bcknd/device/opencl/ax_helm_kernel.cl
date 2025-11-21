#ifndef __MATH_AX_HELM_KERNEL_CL__
#define __MATH_AX_HELM_KERNEL_CL__
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

DEFINE_AX_HELM_KERNEL(2, 256)
DEFINE_AX_HELM_KERNEL(3, 256)
DEFINE_AX_HELM_KERNEL(4, 256)
DEFINE_AX_HELM_KERNEL(5, 256)
DEFINE_AX_HELM_KERNEL(6, 256)
DEFINE_AX_HELM_KERNEL(7, 256)
DEFINE_AX_HELM_KERNEL(8, 256)
DEFINE_AX_HELM_KERNEL(9, 256)
DEFINE_AX_HELM_KERNEL(10, 256)
DEFINE_AX_HELM_KERNEL(11, 256)
DEFINE_AX_HELM_KERNEL(12, 256)


#define DEFINE_AX_HELM_KERNEL_KSTEP(LX)                                        \
__kernel                                                                       \
void ax_helm_kernel_kstep_lx##LX(__global real * __restrict__ w,               \
                                 __global const real * __restrict__ u,         \
                                 __global const real * __restrict__ dx,        \
                                 __global const real * __restrict__ dy,        \
                                 __global const real * __restrict__ dz,        \
                                 __global const real * __restrict__ h1,        \
                                 __global const real * __restrict__ g11,       \
                                 __global const real * __restrict__ g22,       \
                                 __global const real * __restrict__ g33,       \
                                 __global const real * __restrict__ g12,       \
                                 __global const real * __restrict__ g13,       \
                                 __global const real * __restrict__ g23) {     \
                                                                               \
  __local real shdx[LX * LX];                                                  \
  __local real shdy[LX * LX];                                                  \
  __local real shdz[LX * LX];                                                  \
                                                                               \
  __local real shu[LX * LX];                                                   \
  __local real shur[LX * LX];                                                  \
  __local real shus[LX * LX];                                                  \
                                                                               \
  real ru[LX];                                                                 \
  real rw[LX];                                                                 \
  real rut;                                                                    \
                                                                               \
  const int e = get_group_id(0);                                               \
  const int j = get_local_id(1);                                               \
  const int i = get_local_id(0);                                               \
  const int ij = i + j*LX;                                                     \
  const int ele = e*LX*LX*LX;                                                  \
                                                                               \
  shdx[ij] = dx[ij];                                                           \
  shdy[ij] = dy[ij];                                                           \
  shdz[ij] = dz[ij];                                                           \
                                                                               \
  for(int k = 0; k < LX; ++k){                                                 \
    ru[k] = u[ij + k*LX*LX + ele];                                             \
    rw[k] = 0.0;                                                               \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for (int k = 0; k < LX; ++k){                                                \
    const int ijk = ij + k*LX*LX;                                              \
    const real G00 = g11[ijk+ele];                                             \
    const real G11 = g22[ijk+ele];                                             \
    const real G22 = g33[ijk+ele];                                             \
    const real G01 = g12[ijk+ele];                                             \
    const real G02 = g13[ijk+ele];                                             \
    const real G12 = g23[ijk+ele];                                             \
    const real H1  = h1[ijk+ele];                                              \
    real ttmp = 0.0;                                                           \
    shu[ij] = ru[k];                                                           \
    for (int l = 0; l < LX; l++){                                              \
      ttmp += shdz[k+l*LX] * ru[l];                                            \
    }                                                                          \
    barrier(CLK_LOCAL_MEM_FENCE);                                              \
                                                                               \
    real rtmp = 0.0;                                                           \
    real stmp = 0.0;                                                           \
                                                                               \
    for (int l = 0; l < LX; l++){                                              \
      rtmp += shdx[i+l*LX] * shu[l+j*LX];                                      \
      stmp += shdy[j+l*LX] * shu[i+l*LX];                                      \
    }                                                                          \
    shur[ij] = H1                                                              \
             * (G00 * rtmp                                                     \
                + G01 * stmp                                                   \
                + G02 * ttmp);                                                 \
    shus[ij] = H1                                                              \
             * (G01 * rtmp                                                     \
                + G11 * stmp                                                   \
                + G12 * ttmp);                                                 \
    rut      = H1                                                              \
             * (G02 * rtmp                                                     \
                + G12 * stmp                                                   \
                + G22 * ttmp);                                                 \
                                                                               \
    barrier(CLK_LOCAL_MEM_FENCE);                                              \
                                                                               \
    real wijke = 0.0;                                                          \
                                                                               \
    for (int l = 0; l < LX; l++){                                              \
      wijke += shur[l+j*LX] * shdx[l+i*LX];                                    \
      rw[l] += rut * shdz[k+l*LX];                                             \
      wijke += shus[i+l*LX] * shdy[l + j*LX];                                  \
    }                                                                          \
    rw[k] += wijke;                                                            \
  }                                                                            \
                                                                               \
  for (int k = 0; k < LX; ++k){                                                \
    w[ij + k*LX*LX + ele] = rw[k];                                             \
  }                                                                            \
}

DEFINE_AX_HELM_KERNEL_KSTEP(2)
DEFINE_AX_HELM_KERNEL_KSTEP(3)
DEFINE_AX_HELM_KERNEL_KSTEP(4)
DEFINE_AX_HELM_KERNEL_KSTEP(5)
DEFINE_AX_HELM_KERNEL_KSTEP(6)
DEFINE_AX_HELM_KERNEL_KSTEP(7)
DEFINE_AX_HELM_KERNEL_KSTEP(8)
DEFINE_AX_HELM_KERNEL_KSTEP(9)
DEFINE_AX_HELM_KERNEL_KSTEP(10)
DEFINE_AX_HELM_KERNEL_KSTEP(11)
DEFINE_AX_HELM_KERNEL_KSTEP(12)
DEFINE_AX_HELM_KERNEL_KSTEP(13)
DEFINE_AX_HELM_KERNEL_KSTEP(14)
DEFINE_AX_HELM_KERNEL_KSTEP(15)
DEFINE_AX_HELM_KERNEL_KSTEP(16)


/*
 * Vector versions
 */


#define DEFINE_AX_HELM_KERNEL_VECTOR_KSTEP(LX)                                 \
__kernel void                                                                  \
ax_helm_kernel_vector_kstep_lx##LX(__global real * __restrict__ au,            \
                                   __global real * __restrict__ av,            \
                                   __global real * __restrict__ aw,            \
                                   __global const real * __restrict__ u,       \
                                   __global const real * __restrict__ v,       \
                                   __global const real * __restrict__ w,       \
                                   __global const real * __restrict__ dx,      \
                                   __global const real * __restrict__ dy,      \
                                   __global const real * __restrict__ dz,      \
                                   __global const real * __restrict__ h1,      \
                                   __global const real * __restrict__ g11,     \
                                   __global const real * __restrict__ g22,     \
                                   __global const real * __restrict__ g33,     \
                                   __global const real * __restrict__ g12,     \
                                   __global const real * __restrict__ g13,     \
                                   __global const real * __restrict__ g23) {   \
                                                                               \
  __local real shdx[LX * LX];                                                  \
  __local real shdy[LX * LX];                                                  \
  __local real shdz[LX * LX];                                                  \
                                                                               \
  __local real shu[LX * LX];                                                   \
  __local real shur[LX * LX];                                                  \
  __local real shus[LX * LX];                                                  \
                                                                               \
  __local real shv[LX * LX];                                                   \
  __local real shvr[LX * LX];                                                  \
  __local real shvs[LX * LX];                                                  \
                                                                               \
  __local real shw[LX * LX];                                                   \
  __local real shwr[LX * LX];                                                  \
  __local real shws[LX * LX];                                                  \
                                                                               \
  real ru[LX];                                                                 \
  real rv[LX];                                                                 \
  real rw[LX];                                                                 \
                                                                               \
  real ruw[LX];                                                                \
  real rvw[LX];                                                                \
  real rww[LX];                                                                \
                                                                               \
  real rut;                                                                    \
  real rvt;                                                                    \
  real rwt;                                                                    \
                                                                               \
  const int e = get_group_id(0);                                               \
  const int j = get_local_id(1);                                               \
  const int i = get_local_id(0);                                               \
  const int ij = i + j*LX;                                                     \
  const int ele = e*LX*LX*LX;                                                  \
                                                                               \
  shdx[ij] = dx[ij];                                                           \
  shdy[ij] = dy[ij];                                                           \
  shdz[ij] = dz[ij];                                                           \
                                                                               \
  for(int k = 0; k < LX; ++k){                                                 \
    ru[k] = u[ij + k*LX*LX + ele];                                             \
    ruw[k] = 0.0;                                                              \
                                                                               \
    rv[k] = v[ij + k*LX*LX + ele];                                             \
    rvw[k] = 0.0;                                                              \
                                                                               \
    rw[k] = w[ij + k*LX*LX + ele];                                             \
    rww[k] = 0.0;                                                              \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for (int k = 0; k < LX; ++k){                                                \
    const int ijk = ij + k*LX*LX;                                              \
    const real G00 = g11[ijk+ele];                                             \
    const real G11 = g22[ijk+ele];                                             \
    const real G22 = g33[ijk+ele];                                             \
    const real G01 = g12[ijk+ele];                                             \
    const real G02 = g13[ijk+ele];                                             \
    const real G12 = g23[ijk+ele];                                             \
    const real H1  = h1[ijk+ele];                                              \
    real uttmp = 0.0;                                                          \
    real vttmp = 0.0;                                                          \
    real wttmp = 0.0;                                                          \
    shu[ij] = ru[k];                                                           \
    shv[ij] = rv[k];                                                           \
    shw[ij] = rw[k];                                                           \
    for (int l = 0; l < LX; l++){                                              \
      uttmp += shdz[k+l*LX] * ru[l];                                           \
      vttmp += shdz[k+l*LX] * rv[l];                                           \
      wttmp += shdz[k+l*LX] * rw[l];                                           \
    }                                                                          \
    barrier(CLK_LOCAL_MEM_FENCE);                                              \
                                                                               \
    real urtmp = 0.0;                                                          \
    real ustmp = 0.0;                                                          \
                                                                               \
    real vrtmp = 0.0;                                                          \
    real vstmp = 0.0;                                                          \
                                                                               \
    real wrtmp = 0.0;                                                          \
    real wstmp = 0.0;                                                          \
                                                                               \
    for (int l = 0; l < LX; l++){                                              \
      urtmp += shdx[i+l*LX] * shu[l+j*LX];                                     \
      ustmp += shdy[j+l*LX] * shu[i+l*LX];                                     \
                                                                               \
      vrtmp += shdx[i+l*LX] * shv[l+j*LX];                                     \
      vstmp += shdy[j+l*LX] * shv[i+l*LX];                                     \
                                                                               \
      wrtmp += shdx[i+l*LX] * shw[l+j*LX];                                     \
      wstmp += shdy[j+l*LX] * shw[i+l*LX];                                     \
    }                                                                          \
                                                                               \
    shur[ij] = H1                                                              \
             * (G00 * urtmp                                                    \
                + G01 * ustmp                                                  \
                + G02 * uttmp);                                                \
    shus[ij] = H1                                                              \
             * (G01 * urtmp                                                    \
                + G11 * ustmp                                                  \
                + G12 * uttmp);                                                \
    rut      = H1                                                              \
             * (G02 * urtmp                                                    \
                + G12 * ustmp                                                  \
                + G22 * uttmp);                                                \
                                                                               \
    shvr[ij] = H1                                                              \
             * (G00 * vrtmp                                                    \
                + G01 * vstmp                                                  \
                + G02 * vttmp);                                                \
    shvs[ij] = H1                                                              \
             * (G01 * vrtmp                                                    \
                + G11 * vstmp                                                  \
                + G12 * vttmp);                                                \
    rvt      = H1                                                              \
             * (G02 * vrtmp                                                    \
                + G12 * vstmp                                                  \
                + G22 * vttmp);                                                \
                                                                               \
    shwr[ij] = H1                                                              \
             * (G00 * wrtmp                                                    \
                + G01 * wstmp                                                  \
                + G02 * wttmp);                                                \
    shws[ij] = H1                                                              \
             * (G01 * wrtmp                                                    \
                + G11 * wstmp                                                  \
                + G12 * wttmp);                                                \
    rwt      = H1                                                              \
             * (G02 * wrtmp                                                    \
                + G12 * wstmp                                                  \
                + G22 * wttmp);                                                \
                                                                               \
    barrier(CLK_LOCAL_MEM_FENCE);                                              \
                                                                               \
    real uwijke = 0.0;                                                         \
    real vwijke = 0.0;                                                         \
    real wwijke = 0.0;                                                         \
                                                                               \
    for (int l = 0; l < LX; l++){                                              \
      uwijke += shur[l+j*LX] * shdx[l+i*LX];                                   \
      ruw[l] += rut * shdz[k+l*LX];                                            \
      uwijke += shus[i+l*LX] * shdy[l + j*LX];                                 \
                                                                               \
      vwijke += shvr[l+j*LX] * shdx[l+i*LX];                                   \
      rvw[l] += rvt * shdz[k+l*LX];                                            \
      vwijke += shvs[i+l*LX] * shdy[l + j*LX];                                 \
                                                                               \
      wwijke += shwr[l+j*LX] * shdx[l+i*LX];                                   \
      rww[l] += rwt * shdz[k+l*LX];                                            \
      wwijke += shws[i+l*LX] * shdy[l + j*LX];                                 \
    }                                                                          \
    ruw[k] += uwijke;                                                          \
    rvw[k] += vwijke;                                                          \
    rww[k] += wwijke;                                                          \
  }                                                                            \
                                                                               \
  for (int k = 0; k < LX; ++k){                                                \
   au[ij + k*LX*LX + ele] = ruw[k];                                            \
   av[ij + k*LX*LX + ele] = rvw[k];                                            \
   aw[ij + k*LX*LX + ele] = rww[k];                                            \
  }                                                                            \
}

DEFINE_AX_HELM_KERNEL_VECTOR_KSTEP(2)
DEFINE_AX_HELM_KERNEL_VECTOR_KSTEP(3)
DEFINE_AX_HELM_KERNEL_VECTOR_KSTEP(4)
DEFINE_AX_HELM_KERNEL_VECTOR_KSTEP(5)
DEFINE_AX_HELM_KERNEL_VECTOR_KSTEP(6)
DEFINE_AX_HELM_KERNEL_VECTOR_KSTEP(7)
DEFINE_AX_HELM_KERNEL_VECTOR_KSTEP(8)
DEFINE_AX_HELM_KERNEL_VECTOR_KSTEP(9)
DEFINE_AX_HELM_KERNEL_VECTOR_KSTEP(10)
DEFINE_AX_HELM_KERNEL_VECTOR_KSTEP(11)
DEFINE_AX_HELM_KERNEL_VECTOR_KSTEP(12)
DEFINE_AX_HELM_KERNEL_VECTOR_KSTEP(13)
DEFINE_AX_HELM_KERNEL_VECTOR_KSTEP(14)
DEFINE_AX_HELM_KERNEL_VECTOR_KSTEP(15)
DEFINE_AX_HELM_KERNEL_VECTOR_KSTEP(16)

#endif // __MATH_AX_HELM_KERNEL_CL__
