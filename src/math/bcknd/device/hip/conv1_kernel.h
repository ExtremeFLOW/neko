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
__global__ void conv1_kernel_1d(T * __restrict__ du,
                                const T * __restrict__ u,
                                const T * __restrict__ vx,
                                const T * __restrict__ vy,
                                const T * __restrict__ vz,
                                const T * __restrict__ dx,
                                const T * __restrict__ dy,
                                const T * __restrict__ dz,
                                const T * __restrict__ drdx,
                                const T * __restrict__ dsdx,
                                const T * __restrict__ dtdx,
                                const T * __restrict__ drdy,
                                const T * __restrict__ dsdy,
                                const T * __restrict__ dtdy,
                                const T * __restrict__ drdz,
                                const T * __restrict__ dsdz,
                                const T * __restrict__ dtdz,
                                const T * __restrict__ jacinv) { 

  __shared__ T shu[LX * LX * LX];

  __shared__ T shvx[LX * LX * LX];
  __shared__ T shvy[LX * LX * LX];
  __shared__ T shvz[LX * LX * LX];
  
  __shared__ T shdx[LX * LX];
  __shared__ T shdy[LX * LX];
  __shared__ T shdz[LX * LX];
  
  __shared__ T shjacinv[LX * LX * LX];

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

    shvx[l] = vx[l + ele];
    shvy[l] = vy[l + ele];
    shvz[l] = vz[l + ele];
    
    shjacinv[l] = jacinv[l + ele];

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
      
      du[ijk + e * LX * LX * LX] = shjacinv[ijk] *
	(shvx[ijk] * (drdx[ijk + ele] * rtmp
		      + dsdx[ijk + ele] * stmp
		      + dtdx[ijk + ele] * ttmp)
	 + shvy[ijk] * (drdy[ijk + ele] * rtmp
			+ dsdy[ijk + ele] * stmp
			+ dtdy[ijk + ele] * ttmp)
	 + shvz[ijk] * (drdz[ijk + ele] * rtmp
			+ dsdz[ijk + ele] * stmp
			+ dtdz[ijk + ele] * ttmp));
    }
  }  
}

template< typename T, const int LX >
__global__ void __launch_bounds__(LX*LX,3)
  conv1_kernel_kstep(T * __restrict__ du,
                     const T * __restrict__ u,
                     const T * __restrict__ vx,
                     const T * __restrict__ vy,
                     const T * __restrict__ vz,
                     const T * __restrict__ dx,
                     const T * __restrict__ dy,
                     const T * __restrict__ dz,
                     const T * __restrict__ drdx,
                     const T * __restrict__ dsdx,
                     const T * __restrict__ dtdx,
                     const T * __restrict__ drdy,
                     const T * __restrict__ dsdy,
                     const T * __restrict__ dtdy,
                     const T * __restrict__ drdz,
                     const T * __restrict__ dsdz,
                     const T * __restrict__ dtdz,
                     const T * __restrict__ jacinv) { 

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
  T rvx[LX];
  T rvy[LX];
  T rvz[LX];
  T rjacinv[LX];
  
#pragma unroll LX
  for (int k = 0; k < LX; ++k) {
    ru[k] = u[ij + k*LX*LX + ele];
    rvx[k] = vx[ij + k*LX*LX + ele];
    rvy[k] = vy[ij + k*LX*LX + ele];
    rvz[k] = vz[ij + k*LX*LX + ele];
    rjacinv[k] = jacinv[ij + k*LX*LX + ele];
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

    du[ijk + ele] = rjacinv[k] *
	(rvx[k] * (drdx[ijk + ele] * rtmp
                   + dsdx[ijk + ele] * stmp
                   + dtdx[ijk + ele] * ttmp)
	 + rvy[k] * (drdy[ijk + ele] * rtmp
                     + dsdy[ijk + ele] * stmp
                     + dtdy[ijk + ele] * ttmp)
	 + rvz[k] * (drdz[ijk + ele] * rtmp
                     + dsdz[ijk + ele] * stmp
                     + dtdz[ijk + ele] * ttmp));
    __syncthreads();
  }
}

