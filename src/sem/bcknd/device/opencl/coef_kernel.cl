/*
 Copyright (c) 2022-2023, The Neko Authors
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
 * Device kernel for coef geometry
 */
#define DEFINE_GENERATE_GEO_KERNEL(LX, CHUNKS)                                 \
__kernel                                                                       \
void coef_generate_geo_kernel_lx##LX(__global real * __restrict__ G11,         \
                                     __global real * __restrict__ G12,         \
                                     __global real * __restrict__ G13,         \
                                     __global real * __restrict__ G22,         \
                                     __global real * __restrict__ G23,         \
                                     __global real * __restrict__ G33,         \
                                     __global const real * __restrict__ drdx,  \
                                     __global const real * __restrict__ drdy,  \
                                     __global const real * __restrict__ drdz,  \
                                     __global const real * __restrict__ dsdx,  \
                                     __global const real * __restrict__ dsdy,  \
                                     __global const real * __restrict__ dsdz,  \
                                     __global const real * __restrict__ dtdx,  \
                                     __global const real * __restrict__ dtdy,  \
                                     __global const real * __restrict__ dtdz,  \
                                     __global const real * __restrict__ jacinv,\
                                     __global const real * __restrict__ w3,    \
                                     const int gdim) {                         \
                                                                               \
  int i,j,k;                                                                   \
                                                                               \
  const int e = get_group_id(0);                                               \
  const int iii = get_local_id(0);                                             \
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;                         \
                                                                               \
  __local real shw3[LX * LX * LX];                                             \
                                                                               \
  if (iii < (LX * LX * LX)) {                                                  \
    shw3[iii] = w3[iii];                                                       \
  }                                                                            \
                                                                               \
  j = iii;                                                                     \
  while( j < (LX * LX * LX)) {                                                 \
    const int i = j + e * LX * LX * LX;                                        \
    G11[i] = (drdx[i]*drdx[i] + drdy[i]*drdy[i] + drdz[i]*drdz[i]) * jacinv[i];\
    G22[i] = (dsdx[i]*dsdx[i] + dsdy[i]*dsdy[i] + dsdz[i]*dsdz[i]) * jacinv[i];\
    G33[i] = (dtdx[i]*dtdx[i] + dtdy[i]*dtdy[i] + dtdz[i]*dtdz[i]) * jacinv[i];\
                                                                               \
    G12[i] = (drdx[i]*dsdx[i] + drdy[i]*dsdy[i] + drdz[i]*dsdz[i]) * jacinv[i];\
    G13[i] = (drdx[i]*dtdx[i] + drdy[i]*dtdy[i] + drdz[i]*dtdz[i]) * jacinv[i];\
    G23[i] = (dsdx[i]*dtdx[i] + dsdy[i]*dtdy[i] + dsdz[i]*dtdz[i]) * jacinv[i];\
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
      G11[ijk + e * LX * LX * LX] *= shw3[ijk];                                \
      G12[ijk + e * LX * LX * LX] *= shw3[ijk];                                \
      G13[ijk + e * LX * LX * LX] *= shw3[ijk];                                \
      G22[ijk + e * LX * LX * LX] *= shw3[ijk];                                \
      G23[ijk + e * LX * LX * LX] *= shw3[ijk];                                \
      G33[ijk + e * LX * LX * LX] *= shw3[ijk];                                \
    }                                                                          \
  }                                                                            \
}

DEFINE_GENERATE_GEO_KERNEL(2, 256)
DEFINE_GENERATE_GEO_KERNEL(3, 256)
DEFINE_GENERATE_GEO_KERNEL(4, 256)
DEFINE_GENERATE_GEO_KERNEL(5, 256)
DEFINE_GENERATE_GEO_KERNEL(6, 256)
DEFINE_GENERATE_GEO_KERNEL(7, 256)
DEFINE_GENERATE_GEO_KERNEL(8, 256)
DEFINE_GENERATE_GEO_KERNEL(9, 256)
DEFINE_GENERATE_GEO_KERNEL(10, 256)
DEFINE_GENERATE_GEO_KERNEL(11, 256)
DEFINE_GENERATE_GEO_KERNEL(12, 256)
DEFINE_GENERATE_GEO_KERNEL(13, 256)
DEFINE_GENERATE_GEO_KERNEL(14, 256)
DEFINE_GENERATE_GEO_KERNEL(15, 256)
DEFINE_GENERATE_GEO_KERNEL(16, 256)

/**
 * Device kernel for coef dxyz
 */
#define DEFINE_GENERATE_DXYZ_KERNEL(LX, CHUNKS)                                \
__kernel                                                                       \
void coef_generate_dxyz_kernel_lx##LX(__global real * __restrict__ dxdr,       \
                                      __global real * __restrict__ dydr,       \
                                      __global real * __restrict__ dzdr,       \
                                      __global real * __restrict__ dxds,       \
                                      __global real * __restrict__ dyds,       \
                                      __global real * __restrict__ dzds,       \
                                      __global real * __restrict__ dxdt,       \
                                      __global real * __restrict__ dydt,       \
                                      __global real * __restrict__ dzdt,       \
                                      __global const real * __restrict__ dx,   \
                                      __global const real * __restrict__ dy,   \
                                      __global const real * __restrict__ dz,   \
                                      __global const real * __restrict__ x,    \
                                      __global const real * __restrict__ y,    \
                                      __global const real * __restrict__ z) {  \
                                                                               \
  int i,j,k;                                                                   \
                                                                               \
  const int e = get_group_id(0);                                               \
  const int iii = get_local_id(0);                                             \
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;                         \
                                                                               \
  __local real shdx[LX * LX];                                                  \
  __local real shdy[LX * LX];                                                  \
  __local real shdz[LX * LX];                                                  \
                                                                               \
  __local real  shu[LX * LX * LX];                                             \
                                                                               \
  if (iii < (LX * LX)) {                                                       \
    shdx[iii] = dx[iii];                                                       \
    shdy[iii] = dy[iii];                                                       \
    shdz[iii] = dz[iii];                                                       \
  }                                                                            \
                                                                               \
  j = iii;                                                                     \
  while(j < (LX * LX * LX)) {                                                  \
    shu[iii] = x[j + e * LX * LX * LX];                                        \
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
      dxdr[ijk + e * LX * LX * LX] = rtmp;                                     \
      dxds[ijk + e * LX * LX * LX] = stmp;                                     \
      dxdt[ijk + e * LX * LX * LX] = ttmp;                                     \
    }                                                                          \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  j = iii;                                                                     \
  while(j < (LX * LX * LX)) {                                                  \
    shu[iii] = y[j + e * LX * LX * LX];                                        \
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
      dydr[ijk + e * LX * LX * LX] = rtmp;                                     \
      dyds[ijk + e * LX * LX * LX] = stmp;                                     \
      dydt[ijk + e * LX * LX * LX] = ttmp;                                     \
    }                                                                          \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  j = iii;                                                                     \
  while(j < (LX * LX * LX)) {                                                  \
    shu[iii] = z[j + e * LX * LX * LX];                                        \
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
      dzdr[ijk + e * LX * LX * LX] = rtmp;                                     \
      dzds[ijk + e * LX * LX * LX] = stmp;                                     \
      dzdt[ijk + e * LX * LX * LX] = ttmp;                                     \
    }                                                                          \
  }                                                                            \
}

DEFINE_GENERATE_DXYZ_KERNEL(2, 256)
DEFINE_GENERATE_DXYZ_KERNEL(3, 256)
DEFINE_GENERATE_DXYZ_KERNEL(4, 256)
DEFINE_GENERATE_DXYZ_KERNEL(5, 256)
DEFINE_GENERATE_DXYZ_KERNEL(6, 256)
DEFINE_GENERATE_DXYZ_KERNEL(7, 256)
DEFINE_GENERATE_DXYZ_KERNEL(8, 256)
DEFINE_GENERATE_DXYZ_KERNEL(9, 256)
DEFINE_GENERATE_DXYZ_KERNEL(10, 256)
DEFINE_GENERATE_DXYZ_KERNEL(11, 256)
DEFINE_GENERATE_DXYZ_KERNEL(12, 256)
DEFINE_GENERATE_DXYZ_KERNEL(13, 256)
DEFINE_GENERATE_DXYZ_KERNEL(14, 256)
DEFINE_GENERATE_DXYZ_KERNEL(15, 256)
DEFINE_GENERATE_DXYZ_KERNEL(16, 256)

/**
 * Device kernel for coef drst
 */
__kernel void coef_generate_drst_kernel(__global real * __restrict__ jac,
                                        __global real * __restrict__ jacinv,
                                        __global real * __restrict__ drdx,
                                        __global real * __restrict__ drdy,
                                        __global real * __restrict__ drdz,
                                        __global real * __restrict__ dsdx,
                                        __global real * __restrict__ dsdy,
                                        __global real * __restrict__ dsdz,
                                        __global real * __restrict__ dtdx,
                                        __global real * __restrict__ dtdy,
                                        __global real * __restrict__ dtdz,
                                        __global const real * __restrict__ dxdr,
                                        __global const real * __restrict__ dydr,
                                        __global const real * __restrict__ dzdr,
                                        __global const real * __restrict__ dxds,
                                        __global const real * __restrict__ dyds,
                                        __global const real * __restrict__ dzds,
                                        __global const real * __restrict__ dxdt,
                                        __global const real * __restrict__ dydt,
                                        __global const real * __restrict__ dzdt,
                                        const int n) {
                                
  const int idx = get_global_id(0);
  const int str = get_local_size(0) * get_num_groups(0);
  const real one = 1.0;

  for (int i = idx; i < n; i += str) {
    jac[i] = (dxdr[i] * dyds[i] * dzdt[i])
           + (dxdt[i] * dydr[i] * dzds[i])
           + (dxds[i] * dydt[i] * dzdr[i])
           - (dxdr[i] * dydt[i] * dzds[i])
           - (dxds[i] * dydr[i] * dzdt[i])
           - (dxdt[i] * dyds[i] * dzdr[i]);
    jacinv[i] = one / jac[i];    

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
