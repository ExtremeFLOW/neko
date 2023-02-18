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
                               const T * __restrict__ B,
                               const T * __restrict__ jac) { 
  
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
     // We can probably avoid this division by use of jacinv instead.
    T wx = (x[l + e * LX * LX * LX] * B[l + e * LX * LX * LX]) /
      jac[l + e * LX * LX * LX];

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

template< typename T, const int LX >
__global__ void __launch_bounds__(LX*LX,3)
  cdtp_kernel_kstep(T * __restrict__ dtx,
                    const T * __restrict__ x,
                    const T * __restrict__ dr,
                    const T * __restrict__ ds,
                    const T * __restrict__ dt,
                    const T * __restrict__ dxt,
                    const T * __restrict__ dyt,
                    const T * __restrict__ dzt,
                    const T * __restrict__ B,
                    const T * __restrict__ jac) { 
  
  __shared__ T shdxt[LX * LX];
  __shared__ T shdyt[LX * LX];
  __shared__ T shdzt[LX * LX];
  
  __shared__ T shtar[LX * LX];
  __shared__ T shtas[LX * LX];

  T rtar[LX];
  T rtas[LX];
  T rtat[LX];
  T rjac[LX];

  const int e = blockIdx.x;
  const int j = threadIdx.y;
  const int i = threadIdx.x;
  const int ij = i + j * LX;
  const int ele = e*LX*LX*LX;
  
  shdxt[ij] = dxt[ij];
  shdyt[ij] = dyt[ij];
  shdzt[ij] = dzt[ij];


#pragma unroll LX
  for (int k = 0; k < LX; ++k) {
    T wx = (x[ij + k*LX*LX + ele] * B[ij + k*LX*LX + ele]) /
      jac[ij + k*LX*LX + ele];

    rtar[k] = wx *dr[ij + k*LX*LX + ele];
    rtas[k] = wx *ds[ij + k*LX*LX + ele];            
    rtat[k] = wx *dt[ij + k*LX*LX + ele];
  }
  
  __syncthreads();

#pragma unroll
  for (int k = 0; k < LX; ++k) {
    const int ijk = ij + k*LX*LX;
    T ttmp = 0.0;
    shtar[ij] = rtar[k];
    shtas[ij] = rtas[k];
    for (int l = 0; l < LX; l++) {
      ttmp += shdzt[k+l*LX] * rtat[l];
    }
    __syncthreads();

    T rtmp = 0.0;
    T stmp = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++) {
      rtmp += shdxt[i+l*LX] * shtar[l+j*LX];
      stmp += shdyt[j+l*LX] * shtas[i+l*LX];
    }

    dtx[ijk + ele] = ( rtmp + stmp + ttmp );

    __syncthreads();
  }
}

