#ifndef __MATH_SET_CONVECT_RST_KERNEL_CL__
#define __MATH_SET_CONVECT_RST_KERNEL_CL__
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
 * Device kernel for velocity gradients
 */

#define DEFINE_SET_CONVECT_RST_KERNEL(LX, CHUNKS)                              \
__kernel void set_convect_rst_kernel_lx##LX(__global real * __restrict__ cr,   \
                                   __global real * __restrict__ cs,            \
                                   __global real * __restrict__ ct,            \
                                   __global const real * __restrict__ cx,      \
                                   __global const real * __restrict__ cy,      \
                                   __global const real * __restrict__ cz,      \
                                   __global const real * __restrict__ drdx,    \
                                   __global const real * __restrict__ dsdx,    \
                                   __global const real * __restrict__ dtdx,    \
                                   __global const real * __restrict__ drdy,    \
                                   __global const real * __restrict__ dsdy,    \
                                   __global const real * __restrict__ dtdy,    \
                                   __global const real * __restrict__ drdz,    \
                                   __global const real * __restrict__ dsdz,    \
                                   __global const real * __restrict__ dtdz,    \
                                   __global const real * __restrict__ w3) {    \
                                                                               \
                                                                               \
  int i,j,k;                                                                   \
                                                                               \
  const int e = get_group_id(0);                                               \
  const int iii = get_local_id(0);                                             \
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;                         \
                                                                               \
                                                                               \
  j = iii;                                                                     \
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
                                                                               \
      cr[ijk + e * LX * LX * LX] = w3[ijk]                                     \
        * (drdx[ijk + e * LX * LX * LX] * cx[ijk + e * LX * LX * LX]           \
           + drdy[ijk + e * LX * LX * LX] * cy[ijk + e * LX * LX * LX]         \
           + drdz[ijk + e * LX * LX * LX] * cz[ijk + e * LX * LX * LX]);       \
                                                                               \
      cs[ijk + e * LX * LX * LX] = w3[ijk]                                     \
        * (dsdx[ijk + e * LX * LX * LX] * cx[ijk + e * LX * LX * LX]           \
           + dsdy[ijk + e * LX * LX * LX] * cy[ijk + e * LX * LX * LX]         \
           + dsdz[ijk + e * LX * LX * LX] * cz[ijk + e * LX * LX * LX]);       \
                                                                               \
      ct[ijk + e * LX * LX * LX] = w3[ijk]                                     \
        * (dtdx[ijk + e * LX * LX * LX] * cx[ijk + e * LX * LX * LX]           \
           + dtdy[ijk + e * LX * LX * LX] * cy[ijk + e * LX * LX * LX]         \
           + dtdz[ijk + e * LX * LX * LX] * cz[ijk + e * LX * LX * LX]);       \
                                                                               \
    }                                                                          \
  }                                                                            \
}

DEFINE_SET_CONVECT_RST_KERNEL(12, 256)
DEFINE_SET_CONVECT_RST_KERNEL(11, 256)
DEFINE_SET_CONVECT_RST_KERNEL(10, 256)
DEFINE_SET_CONVECT_RST_KERNEL(9, 256)
DEFINE_SET_CONVECT_RST_KERNEL(8, 256)
DEFINE_SET_CONVECT_RST_KERNEL(7, 256)
DEFINE_SET_CONVECT_RST_KERNEL(6, 256)
DEFINE_SET_CONVECT_RST_KERNEL(5, 256)
DEFINE_SET_CONVECT_RST_KERNEL(4, 256)
DEFINE_SET_CONVECT_RST_KERNEL(3, 256)
DEFINE_SET_CONVECT_RST_KERNEL(2, 256)

#endif // __MATH_SET_CONVECT_RST_KERNEL_CL__
