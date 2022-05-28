/*
 Copyright (c) 2022, The Neko Authors
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
 * Device kernel for CFL
 */

#define DEFINE_CFL_KERNEL(LX, CHUNKS)                                          \
__kernel void cfl_kernel_lx##LX(const real dt,                                 \
				__global const real * __restrict__ u,	       \
				__global const real * __restrict__ v,	       \
				__global const real * __restrict__ w,	       \
				__global const real * __restrict__ drdx,       \
				__global const real * __restrict__ dsdx,       \
				__global const real * __restrict__ dtdx,       \
				__global const real * __restrict__ drdy,       \
				__global const real * __restrict__ dsdy,       \
				__global const real * __restrict__ dtdy,       \
				__global const real * __restrict__ drdz,       \
				__global const real * __restrict__ dsdz,       \
				__global const real * __restrict__ dtdz,       \
				__global const real * __restrict__ dr_inv,     \
				__global const real * __restrict__ ds_inv,     \
				__global const real * __restrict__ dt_inv,     \
				__global const real * __restrict__ jacinv,     \
				__global real * __restrict__ cfl_h) {	       \
                                                                               \
  int i,j,k;								       \
                                                                               \
  const int e = get_group_id(0);                                               \
  const int iii = get_local_id(0);                                             \
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;                         \
                                                                               \
  __local real shu[LX * LX * LX];                                              \
  __local real shv[LX * LX * LX];                                              \
  __local real shw[LX * LX * LX];                                              \
                                                                               \
  __local real shdr_inv[LX];                                                   \
  __local real shds_inv[LX];                                                   \
  __local real shdt_inv[LX];                                                   \
                                                                               \
  __local real shjacinv[LX * LX * LX];                                         \
                                                                               \
  __local real shcfl[256];                                                     \
                                                                               \
  if (iii < LX) {                                                              \
    shdr_inv[iii] = dr_inv[iii];                                               \
    shds_inv[iii] = ds_inv[iii];                                               \
    shdt_inv[iii] = dt_inv[iii];                                               \
  }                                                                            \
                                                                               \
  j = iii;                                                                     \
  while(j < (LX * LX * LX)) {                                                  \
    shu[j] = u[j + e * LX * LX * LX];                                          \
    shv[j] = v[j + e * LX * LX * LX];                                          \
    shw[j] = w[j + e * LX * LX * LX];                                          \
                                                                               \
    shjacinv[j] = jacinv[j + e * LX * LX * LX];                                \
                                                                               \
    j = j + CHUNKS;                                                            \
  }                                                                            \
                                                                               \
  shcfl[iii] = 0.0;                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  real cfl_tmp = 0.0;                                                          \
  for (int n = 0; n < nchunks; n++) {                                          \
    const int ijk = iii + n * CHUNKS;                                          \
    const int jk = ijk / LX;                                                   \
    i = ijk - jk * LX;                                                         \
    k = jk / LX;                                                               \
    j = jk - k * LX;                                                           \
    if ( i < LX && j < LX && k < LX) {                                         \
      const real cflr = fabs(dt * (( shu[ijk] * drdx[ijk + e * LX * LX * LX]   \
                                     + shv[ijk] * drdy[ijk + e * LX * LX * LX] \
                                     * shw[ijk] * drdz[ijk + e * LX * LX * LX] \
                                     ) * shjacinv[ijk]) * shdr_inv[i]);        \
      const real cfls = fabs(dt * (( shu[ijk] * dsdx[ijk + e * LX * LX * LX]   \
                                     + shv[ijk] * dsdy[ijk + e * LX * LX * LX] \
                                     + shw[ijk] * dsdz[ijk + e * LX * LX * LX] \
                                     ) * shjacinv[ijk]) * shds_inv[j]);        \
      const real cflt = fabs( dt * ( ( shu[ijk] * dtdx[ijk + e * LX * LX * LX] \
                                      + shv[ijk] * dtdy[ijk + e * LX * LX * LX]\
				      + shw[ijk] * dtdz[ijk + e * LX * LX * LX]\
				      ) * shjacinv[ijk]) * shdt_inv[k]);       \
                                                                               \
      cfl_tmp = fmax(cflr + cfls + cflt, cfl_tmp);                             \
                                                                               \
    }                                                                          \
  }                                                                            \
  shcfl[iii] = cfl_tmp;                                                        \
                                                                               \
  i = (get_local_size(0)) >> 1;                                                \
  while (i != 0) {                                                             \
    if (iii < i) {                                                             \
      shcfl[iii] = fmax(shcfl[iii], shcfl[iii + i]);                           \
    }                                                                          \
    barrier(CLK_LOCAL_MEM_FENCE);                                              \
    i = i>>1;                                                                  \
  }                                                                            \
                                                                               \
  if (get_local_id(0) == 0) {                                                  \
    cfl_h[get_group_id(0)] = shcfl[0];                                         \
  }                                                                            \
}                                                                             

DEFINE_CFL_KERNEL(2, 256)
DEFINE_CFL_KERNEL(3, 256)
DEFINE_CFL_KERNEL(4, 256)
DEFINE_CFL_KERNEL(5, 256)
DEFINE_CFL_KERNEL(6, 256)
DEFINE_CFL_KERNEL(7, 256)
DEFINE_CFL_KERNEL(8, 256)
DEFINE_CFL_KERNEL(9, 256)

