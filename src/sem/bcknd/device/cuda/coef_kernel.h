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

  __shared__ T shw3[LX * LX * LX];

  j = iii;
  while( j < (LX * LX * LX)) {
    const int i = j + e * LX * LX * LX;
    G11[i] = (drdx[i]*drdx[i] + drdy[i]*drdy[i] + drdz[i]*drdz[i]) * jacinv[i];
    G22[i] = (dsdx[i]*dsdx[i] + dsdy[i]*dsdy[i] + dsdz[i]*dsdz[i]) * jacinv[i];
    G33[i] = (dtdx[i]*dtdx[i] + dtdy[i]*dtdy[i] + dtdz[i]*dtdz[i]) * jacinv[i];

    G12[i] = (drdx[i]*dsdx[i] + drdy[i]*dsdy[i] + drdz[i]*dsdz[i]) * jacinv[i];
    G13[i] = (drdx[i]*dtdx[i] + drdy[i]*dtdy[i] + drdz[i]*dtdz[i]) * jacinv[i];
    G23[i] = (dsdx[i]*dtdx[i] + dsdy[i]*dtdy[i] + dsdz[i]*dtdz[i]) * jacinv[i];

    shw3[j] = w3[j];

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

/**
 * Device kernel for coef dxyz
 */
template< typename T, const int LX, const int CHUNKS >
__global__ void coef_generate_dxyz_kernel(T * __restrict__ dxdr, 
					  T * __restrict__ dydr, 
					  T * __restrict__ dzdr, 
					  T * __restrict__ dxds, 
					  T * __restrict__ dyds, 
					  T * __restrict__ dzds, 
					  T * __restrict__ dxdt, 
					  T * __restrict__ dydt, 
					  T * __restrict__ dzdt,
					  const T * __restrict__ dx, 
					  const T * __restrict__ dy, 
					  const T * __restrict__ dz, 
					  const T * __restrict__ x, 
					  const T * __restrict__ y, 
					  const T * __restrict__ z) {

  int i,j,k;
  
  const int e = blockIdx.x;
  const int iii = threadIdx.x;
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;

  __shared__ T shdx[LX * LX];
  __shared__ T shdy[LX * LX];
  __shared__ T shdz[LX * LX];

  __shared__ T shu[LX * LX * LX];

  if (iii < (LX * LX)) {
    shdx[iii] = dx[iii];
    shdy[iii] = dy[iii];
    shdz[iii] = dz[iii];
  }

  j = iii;
  while(j < (LX * LX * LX)) {
    shu[j] = x[j + e * LX * LX * LX];
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
      T rtmp = 0.0;
      T stmp = 0.0;
      T ttmp = 0.0;
      for (int l = 0; l < LX; l++) {		
	rtmp += shdx[i + l * LX] * shu[l + j * LX + k * LX * LX];	
	stmp += shdy[j + l * LX] * shu[i + l * LX + k * LX * LX];
	ttmp += shdz[k + l * LX] * shu[i + j * LX + l * LX * LX];	
      } 
      dxdr[ijk + e * LX * LX * LX] = rtmp;
      dxds[ijk + e * LX * LX * LX] = stmp;
      dxdt[ijk + e * LX * LX * LX] = ttmp;
    }
  }

  __syncthreads();

  j = iii;
  while(j < (LX * LX * LX)) {
    shu[j] = y[j + e * LX * LX * LX];
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
      T rtmp = 0.0;
      T stmp = 0.0;
      T ttmp = 0.0;
      for (int l = 0; l < LX; l++) {		
	rtmp += shdx[i + l * LX] * shu[l + j * LX + k * LX * LX];	
	stmp += shdy[j + l * LX] * shu[i + l * LX + k * LX * LX];
	ttmp += shdz[k + l * LX] * shu[i + j * LX + l * LX * LX];	
      } 
      dydr[ijk + e * LX * LX * LX] = rtmp;
      dyds[ijk + e * LX * LX * LX] = stmp;
      dydt[ijk + e * LX * LX * LX] = ttmp;
    }
  }

  __syncthreads();

  j = iii;
  while(j < (LX * LX * LX)) {
    shu[j] = z[j + e * LX * LX * LX];
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
      T rtmp = 0.0;
      T stmp = 0.0;
      T ttmp = 0.0;
      for (int l = 0; l < LX; l++) {		
	rtmp += shdx[i + l * LX] * shu[l + j * LX + k * LX * LX];	
	stmp += shdy[j + l * LX] * shu[i + l * LX + k * LX * LX];
	ttmp += shdz[k + l * LX] * shu[i + j * LX + l * LX * LX];	
      } 
      dzdr[ijk + e * LX * LX * LX] = rtmp;
      dzds[ijk + e * LX * LX * LX] = stmp;
      dzdt[ijk + e * LX * LX * LX] = ttmp;
    }
  }
}

/**
 * Device kernel for coef drst
 */
template< typename T >
__global__ void coef_generate_drst_kernel(T * __restrict__ jac,
					  T * __restrict__ jacinv,
					  T * __restrict__ drdx,
					  T * __restrict__ drdy,
					  T * __restrict__ drdz,
					  T * __restrict__ dsdx,
					  T * __restrict__ dsdy,
					  T * __restrict__ dsdz,
					  T * __restrict__ dtdx,
					  T * __restrict__ dtdy,
					  T * __restrict__ dtdz,
					  const T * __restrict__ dxdr, 
					  const T * __restrict__ dydr, 
					  const T * __restrict__ dzdr, 
					  const T * __restrict__ dxds, 
					  const T * __restrict__ dyds, 
					  const T * __restrict__ dzds, 
					  const T * __restrict__ dxdt, 
					  const T * __restrict__ dydt, 
					  const T * __restrict__ dzdt,
					  const int n) {
				
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  const T one = 1.0;

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
