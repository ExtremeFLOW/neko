#ifndef __MATH_CONVECT_SCALAR_KERNEL_H__
#define __MATH_CONVECT_SCALAR_KERNEL_H__
/*
 Copyright (c) 2021-2023, The Neko Authors
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

template< typename T, const int LX, const int CHUNKS >
__global__ void convect_scalar_kernel_1d(T * __restrict__ du,
                                const T * __restrict__ u,
                                const T * __restrict__ cr,
                                const T * __restrict__ cs,
                                const T * __restrict__ ct,
                                const T * __restrict__ dx,
                                const T * __restrict__ dy,
                                const T * __restrict__ dz) {

  __shared__ T shu[LX * LX * LX];

  __shared__ T shcr[LX * LX * LX];
  __shared__ T shcs[LX * LX * LX];
  __shared__ T shct[LX * LX * LX];

  __shared__ T shdx[LX * LX];
  __shared__ T shdy[LX * LX];
  __shared__ T shdz[LX * LX];

  const int e = blockIdx.x;
  const int iii = threadIdx.x;
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;
  const int ele = e*LX*LX*LX;

  if (iii < (LX * LX)) {
    shdx[iii] = dx[iii];
    shdy[iii] = dy[iii];
    shdz[iii] = dz[iii];
  }

  int l = iii;
  while(l < (LX * LX * LX)) {
    shu[l] = u[l + ele];

    shcr[l] = cr[l + ele];
    shcs[l] = cs[l + ele];
    shct[l] = ct[l + ele];

    l = l + CHUNKS;
  }

  __syncthreads();

  for (int n = 0; n < nchunks; n++) {
    const int ijk = iii + n * CHUNKS;
    const int jk = ijk / LX;
    const int i = ijk - jk * LX;
    const int k = jk / LX;
    const int j = jk - k * LX;
    if ( i < LX && j < LX && k < LX) {
      T rtmp = 0.0;
      T stmp = 0.0;
      T ttmp = 0.0;
      for (int l = 0; l < LX; l++) {
	rtmp += shdx[i + l * LX] * shu[l + j * LX + k * LX * LX];
	stmp += shdy[j + l * LX] * shu[i + l * LX + k * LX * LX];
	ttmp += shdz[k + l * LX] * shu[i + j * LX + l * LX * LX];
      }

      du[ijk + e * LX * LX * LX] = shcr[ijk] * rtmp
	                               + shcs[ijk] * stmp
	                               + shct[ijk] * ttmp;
    }
  }
}

template< typename T, const int LX >
__global__ void __launch_bounds__(LX*LX,3)
  convect_scalar_kernel_kstep(T * __restrict__ du,
                     const T * __restrict__ u,
                     const T * __restrict__ cr,
                     const T * __restrict__ cs,
                     const T * __restrict__ ct,
                     const T * __restrict__ dx,
                     const T * __restrict__ dy,
                     const T * __restrict__ dz) {

  __shared__ T shu[LX * LX];

  __shared__ T shdx[LX*LX];
  __shared__ T shdy[LX*LX];
  __shared__ T shdz[LX*LX];

  const int e = blockIdx.x;
  const int j = threadIdx.y;
  const int i = threadIdx.x;
  const int ij = i + j * LX;
  const int ele = e*LX*LX*LX;

  shdx[ij] = dx[ij];
  shdy[ij] = dy[ij];
  shdz[ij] = dz[ij];

  T ru[LX];
  T rcr[LX];
  T rcs[LX];
  T rct[LX];

#pragma unroll LX
  for (int k = 0; k < LX; ++k) {
    ru[k] = u[ij + k*LX*LX + ele];
    rcr[k] = cr[ij + k*LX*LX + ele];
    rcs[k] = cs[ij + k*LX*LX + ele];
    rct[k] = ct[ij + k*LX*LX + ele];
  }

  __syncthreads();

#pragma unroll
  for (int k = 0; k < LX; ++k) {
    const int ijk = ij + k*LX*LX;
    T ttmp = 0.0;
    shu[ij] = ru[k];
    for (int l = 0; l < LX; l++) {
      ttmp += shdz[k+l*LX] * ru[l];
    }
    __syncthreads();

    T rtmp = 0.0;
    T stmp = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++) {
      rtmp += shdx[i+l*LX] * shu[l+j*LX];
      stmp += shdy[j+l*LX] * shu[i+l*LX];
    }

    du[ijk + ele] = rcr[k] * rtmp
	                + rcs[k] * stmp
	                + rct[k] * ttmp;
    __syncthreads();
  }
}


#endif // __MATH_CONVECT_SCALAR_KERNEL_H__
