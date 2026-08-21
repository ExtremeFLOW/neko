#ifndef __MATH_CDTP_KERNEL_H__
#define __MATH_CDTP_KERNEL_H__

#include "elem_block.h"
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
 * Device kernel for \f$ D^T x \f$
 */

template< typename T, const int LX, const int CHUNKS >
__global__ void cdtp_kernel_1d(T * __restrict__ dtx,
                               const T * __restrict__ x,
                               const T * __restrict__ dr,
                               const T * __restrict__ ds,
                               const T * __restrict__ dt,
                               const T * __restrict__ dxt,
                               const T * __restrict__ dyt,
                               const T * __restrict__ dzt,
                               const T * __restrict__ w3) { 
  
  __shared__ T shdxt[LX * LX];
  __shared__ T shdyt[LX * LX];
  __shared__ T shdzt[LX * LX];
  
  __shared__ T shtar[LX * LX * LX];
  __shared__ T shtas[LX * LX * LX];
  __shared__ T shtat[LX * LX * LX];
  
  const int e = blockIdx.x;
  const int iii = threadIdx.x;
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;

  if (iii < (LX * LX)) {
    shdxt[iii] = dxt[iii];
    shdyt[iii] = dyt[iii];
    shdzt[iii] = dzt[iii];
  }

  int l = iii;
  while(l < (LX * LX * LX)) {
    T wx = x[l + e * LX * LX * LX] * w3[l]; 

    shtar[l] = wx*dr[l + e * LX * LX * LX];
    shtas[l] = wx*ds[l + e * LX * LX * LX];
    shtat[l] = wx*dt[l + e * LX * LX * LX];

    l = l + CHUNKS;
  }

  __syncthreads();
  for (int n = 0; n < nchunks; n++) {
    const int ijk = iii + n * CHUNKS;
    const int jk = ijk / LX;
    const int i = ijk - jk * LX;
    const int k = jk / LX;
    const int j = jk - k * LX;
    if ( i < LX && j < LX && k < LX && ijk < LX*LX*LX) {
      T rtmp = 0.0;
      T stmp = 0.0;
      T ttmp = 0.0;
      for (int l = 0; l < LX; l++) {
	    rtmp += shdxt[i + l * LX] * shtar[l+j*LX+k*LX*LX];	
	    stmp += shdyt[j + l * LX] * shtas[i+l*LX + k*LX*LX];
	    ttmp += shdzt[k + l * LX] * shtat[i + j*LX + l*LX*LX];
      }
      dtx[ijk + e * LX * LX * LX] = ( rtmp + stmp + ttmp );
      
    }
  }
}

template< typename T, const int LX, const int EB >
__global__ void NEKO_EB_BOUNDS(LX*LX*EB)
  cdtp_kernel_kstep(T * __restrict__ dtx,
                    const T * __restrict__ x,
                    const T * __restrict__ dr,
                    const T * __restrict__ ds,
                    const T * __restrict__ dt,
                    const T * __restrict__ dxt,
                    const T * __restrict__ dyt,
                    const T * __restrict__ dzt,
                    const T * __restrict__ w3,
                    const int nelv) { 
  
  __shared__ T shdxt[LX * LX];
  __shared__ T shdyt[LX * LX];
  __shared__ T shdzt[LX * LX];
  
  __shared__ T shtar[EB * LX * LX];
  __shared__ T shtas[EB * LX * LX];

  static_assert(sizeof(shdxt) +
                sizeof(shdyt) +
                sizeof(shdzt) +
                sizeof(shtar) +
                sizeof(shtas)
                <= NEKO_EB_MAX_SMEM,
                "kstep block exceeds the shared memory budget");

  T rtar[LX];
  T rtas[LX];
  T rtat[LX];

  const int eb = (EB == 1) ? 0 : threadIdx.z;
  const int e_blk = blockIdx.x * EB + eb;
  /* Threads past the last element still have to reach the barriers in
     the k loop, so clamp their reads and drop their stores rather than
     returning early. At EB == 1 this all constant folds away */
  const bool active = (EB == 1) ? true : (e_blk < nelv);
  const int e = active ? e_blk : (nelv - 1);
  const int sh = eb * LX * LX;
  const int j = threadIdx.y;
  const int i = threadIdx.x;
  const int ij = i + j * LX;
  const int ele = e*LX*LX*LX;
  
  if (eb == 0) {
    shdxt[ij] = dxt[ij];
    shdyt[ij] = dyt[ij];
    shdzt[ij] = dzt[ij];
  }


#pragma unroll LX
  for (int k = 0; k < LX; ++k) {
    T wx = x[ij + k*LX*LX + ele] * w3[ij + k*LX*LX];

    rtar[k] = wx *dr[ij + k*LX*LX + ele];
    rtas[k] = wx *ds[ij + k*LX*LX + ele];            
    rtat[k] = wx *dt[ij + k*LX*LX + ele];
  }
  
  __syncthreads();

#pragma unroll
  for (int k = 0; k < LX; ++k) {
    const int ijk = ij + k*LX*LX;
    T ttmp = 0.0;
    shtar[sh + ij] = rtar[k];
    shtas[sh + ij] = rtas[k];
#pragma unroll
    for (int l = 0; l < LX; l++) {
      ttmp += shdzt[k+l*LX] * rtat[l];
    }
    __syncthreads();

    T rtmp = 0.0;
    T stmp = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++) {
      rtmp += shdxt[i+l*LX] * shtar[sh + l+j*LX];
      stmp += shdyt[j+l*LX] * shtas[sh + i+l*LX];
    }

    if (active) {
      dtx[ijk + ele] = ( rtmp + stmp + ttmp );
    }

    __syncthreads();
  }
}


#endif // __MATH_CDTP_KERNEL_H__
