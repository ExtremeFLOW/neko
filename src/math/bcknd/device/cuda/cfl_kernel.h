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

template< typename T, const int LX, const int CHUNKS >
__global__ void cfl_kernel(const T dt,
			   const T * __restrict__ u,
			   const T * __restrict__ v,
			   const T * __restrict__ w,
			   const T * __restrict__ drdx,
			   const T * __restrict__ dsdx,
			   const T * __restrict__ dtdx,
			   const T * __restrict__ drdy,
			   const T * __restrict__ dsdy,
			   const T * __restrict__ dtdy,
			   const T * __restrict__ drdz,
			   const T * __restrict__ dsdz,
			   const T * __restrict__ dtdz,
			   const T * __restrict__ dr_inv,
			   const T * __restrict__ ds_inv, 
			   const T * __restrict__ dt_inv,
			   const T * __restrict__ jacinv,
			   T * __restrict__ cfl_h) { 

  int i,j,k;
  
  const int e = blockIdx.x;
  const int iii = threadIdx.x;
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;

  __shared__ T shu[LX * LX * LX];
  __shared__ T shv[LX * LX * LX];
  __shared__ T shw[LX * LX * LX];

  __shared__ T shdr_inv[LX];
  __shared__ T shds_inv[LX];
  __shared__ T shdt_inv[LX];

  __shared__ T shjacinv[LX * LX * LX];

  __shared__ T shcfl[1024];

  if (iii < LX) {
    shdr_inv[iii] = dr_inv[iii];
    shds_inv[iii] = ds_inv[iii];
    shdt_inv[iii] = dt_inv[iii];
  }

  j = iii;
  while(j < (LX * LX * LX)) {
    shu[j] = u[j + e * LX * LX * LX];
    shv[j] = v[j + e * LX * LX * LX];
    shw[j] = w[j + e * LX * LX * LX];
    
    shjacinv[j] = jacinv[j + e * LX * LX * LX];

    j = j + CHUNKS;
  }
  
  shcfl[iii] = 0.0;

  __syncthreads();

  T cfl_tmp = 0.0;
  for (int n = 0; n < nchunks; n++) {
    const int ijk = iii + n * CHUNKS;
    const int jk = ijk / LX;
    i = ijk - jk * LX;
    k = jk / LX;
    j = jk - k * LX;
    if ( i < LX && j < LX && k < LX) {
      const T cflr = fabs( dt * ( ( shu[ijk] * drdx[ijk + e * LX * LX * LX]
                                    + shv[ijk] * drdy[ijk + e * LX * LX * LX]
                                    * shw[ijk] * drdz[ijk + e * LX * LX * LX] 
                                    ) * shjacinv[ijk]) * shdr_inv[i]);
      const T cfls = fabs( dt * ( ( shu[ijk] * dsdx[ijk + e * LX * LX * LX]
                                    + shv[ijk] * dsdy[ijk + e * LX * LX * LX]
                                    + shw[ijk] * dsdz[ijk + e * LX * LX * LX] 
                                    ) * shjacinv[ijk]) * shds_inv[j]);
      const T cflt = fabs( dt * ( ( shu[ijk] * dtdx[ijk + e * LX * LX * LX]
                                    + shv[ijk] * dtdy[ijk + e * LX * LX * LX]
                                    + shw[ijk] * dtdz[ijk + e * LX * LX * LX] 
                                    ) * shjacinv[ijk]) * shdt_inv[k]);

      cfl_tmp = fmax(cflr + cfls + cflt, cfl_tmp);

    }
  }
  shcfl[iii] = cfl_tmp;
  
  i = blockDim.x >> 1;
  while (i != 0) {
    if (iii < i) {
      shcfl[iii] = fmax(shcfl[iii], shcfl[iii + i]);
    }
    __syncthreads();
    i = i>>1;
  }

  if (threadIdx.x == 0) {
    cfl_h[blockIdx.x] = shcfl[0];
  }
}

