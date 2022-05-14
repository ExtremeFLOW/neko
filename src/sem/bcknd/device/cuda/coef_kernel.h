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
 * Device kernel for coef geometry
 */
template< typename T, const int LX, const int CHUNKS >
__global__ void coef_generate_geo_kernel(T * __restrict__ G11,
					 T * __restrict__ G12, 
					 T * __restrict__ G13,
					 T * __restrict__ G22,
					 T * __restrict__ G23, 
					 T * __restrict__ G33,
					 const T * __restrict__ drdx, 
					 const T * __restrict__ drdy, 
					 const T * __restrict__ drdz,
					 const T * __restrict__ dsdx, 
					 const T * __restrict__ dsdy, 
					 const T * __restrict__ dsdz, 
					 const T * __restrict__ dtdx,
					 const T * __restrict__ dtdy, 
					 const T * __restrict__ dtdz,
					 const T * __restrict__ jacinv, 
					 const T * __restrict__ w3,
					 const int gdim) {

  int i,j,k;
  
  const int e = blockIdx.x;
  const int iii = threadIdx.x;
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  __shared__ T shw3[LX * LX * LX];

  if (iii < (LX * LX * LX)) {
    shw3[iii] = w3[iii];
  }

  j = iii;
  while( j < (LX * LX * LX)) {
    const int i = j + e * LX * LX * LX;
    G11[i] = (drdx[i]*drdx[i] + drdy[i]*drdy[i] + drdz[i]*drdz[i]) * jacinv[i];
    G22[i] = (dsdx[i]*dsdx[i] + dsdy[i]*dsdy[i] + dsdz[i]*dsdz[i]) * jacinv[i];
    G33[i] = (dtdx[i]*dtdx[i] + dtdy[i]*dtdy[i] + dtdz[i]*dtdz[i]) * jacinv[i];

    G12[i] = (drdx[i]*dsdx[i] + drdy[i]*dsdy[i] + drdz[i]*dsdz[i]) * jacinv[i];
    G13[i] = (drdx[i]*dtdx[i] + drdy[i]*dtdy[i] + drdz[i]*dtdz[i]) * jacinv[i];
    G23[i] = (dsdx[i]*dtdx[i] + dsdy[i]*dtdy[i] + dsdz[i]*dtdz[i]) * jacinv[i];
    j = j + CHUNKS;
  }

  __syncthreads();

  for (int n = 0; n < nchunks; n++) {
    const int ijk = iii + n * CHUNKS;
    const int jk = ijk / LX;
    i = ijk - jk * LX;
    k = jk / LX;
    j = jk - k * LX;
    if ( i < LX && j < LX && k < LX) {
      G11[ijk + e * LX * LX * LX] *= shw3[ijk];
      G12[ijk + e * LX * LX * LX] *= shw3[ijk];
      G13[ijk + e * LX * LX * LX] *= shw3[ijk];
      G22[ijk + e * LX * LX * LX] *= shw3[ijk];
      G23[ijk + e * LX * LX * LX] *= shw3[ijk];
      G33[ijk + e * LX * LX * LX] *= shw3[ijk];
    }
  }
}

