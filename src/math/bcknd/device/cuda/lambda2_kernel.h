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

#include <math/bcknd/device/cuda/opgrad_kernel.h>

template< typename T>
__inline__ __device__ T eigen_val_calc(T grad11,T grad12,T grad13,T grad21,T grad22,T grad23, T grad31, T grad32, T grad33){
    T s11 = grad11;
    T s22 = grad22;
    T s33 = grad33;
    T s12 = 0.5*(grad12+grad21);
    T s13 = 0.5*(grad13+grad31);
    T s23 = 0.5*(grad23+grad32);

    T o12 = 0.5*(grad12-grad21);
    T o13 = 0.5*(grad13-grad31);
    T o23 = 0.5*(grad23-grad32);

    T a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13;
    T a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23;
    T a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23;
        
    T a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23;
    T a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13;
    T a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23;
               
               
    T B = -(a11 + a22 + a33);
    T C = -(a12*a12 + a13*a13 + a23*a23 - a11 * a22 - a11 * a33 - a22 * a33);
    T D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 - a22 * a13*a13 - a33 * a12*a12  +  a11 * a22 * a33);
                     
                     
    T q = (3.0 * C - B*B) / 9.0;
    T r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0;

    T theta = acos( r / sqrt(-q*q*q) ); 
    T pi = 4.0*atan(1.0);
    T eigen1 = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0;
    T eigen2 = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0;
    T eigen3 = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0;
                                 
    if (eigen1 <= eigen2 && eigen2 <= eigen3)
        return eigen2;
    else if (eigen3 <= eigen2 && eigen2 <= eigen1)
        return eigen2;
    else if (eigen1 <= eigen3 && eigen3 <= eigen2)
        return eigen3;
    else if (eigen2 <= eigen3 && eigen3 <= eigen1)
        return eigen3;
    else if (eigen2 <= eigen1 && eigen1 <= eigen3)
        return eigen1;
    else if (eigen3 <= eigen1 && eigen1 <= eigen2)
        return eigen1;
    return 0.0;
}


template< typename T, const int LX, const int CHUNKS >
__global__ void lambda2_kernel_1d(T * __restrict__ lambda2,
                                 const T * __restrict__ u,
                                 const T * __restrict__ v,
                                 const T * __restrict__ w,
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
  __shared__ T shv[LX * LX * LX];
  __shared__ T shw[LX * LX * LX];

  __shared__ T shdx[LX * LX];
  __shared__ T shdy[LX * LX];
  __shared__ T shdz[LX * LX];


  
  
  int i,j,k;
  
  const int e = blockIdx.x;
  const int ele = blockIdx.x*LX*LX*LX;
  const int iii = threadIdx.x;
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;


  if (iii < (LX * LX)) {
    shdx[iii] = dx[iii];
    shdy[iii] = dy[iii];
    shdz[iii] = dz[iii];
  }

  j = iii;
  while(j < (LX * LX * LX)) {
    shu[j] = u[j + ele];
    shv[j] = v[j + ele];
    shw[j] = w[j + ele];
    j = j + CHUNKS;
  }
  
  __syncthreads();
  
  for (int n = 0; n < nchunks; n++) {
    const int ijk = iii + n * CHUNKS;
    const int jk = ijk / LX;
    i = ijk - jk * LX;
    k = jk / LX;
    j = jk - k * LX;
    if ( i < LX && j < LX && k < LX ) {
      T rtmpu = 0.0;
      T stmpu = 0.0;
      T ttmpu = 0.0;
      
      T rtmpv = 0.0;
      T stmpv = 0.0;
      T ttmpv = 0.0;

      T rtmpw = 0.0;
      T stmpw = 0.0;
      T ttmpw = 0.0;
      for (int l = 0; l < LX; l++) {		
	      rtmpu += shdx[i + l * LX] * shu[l + j * LX + k * LX * LX];	
	      stmpu += shdy[j + l * LX] * shu[i + l * LX + k * LX * LX];
	      ttmpu += shdz[k + l * LX] * shu[i + j * LX + l * LX * LX];
	      
          rtmpv += shdx[i + l * LX] * shv[l + j * LX + k * LX * LX];	
	      stmpv += shdy[j + l * LX] * shv[i + l * LX + k * LX * LX];
	      ttmpv += shdz[k + l * LX] * shv[i + j * LX + l * LX * LX];
	      
          rtmpw += shdx[i + l * LX] * shw[l + j * LX + k * LX * LX];	
	      stmpw += shdy[j + l * LX] * shw[i + l * LX + k * LX * LX];
	      ttmpw += shdz[k + l * LX] * shw[i + j * LX + l * LX * LX];
      }

      T jinv = jacinv[ijk + ele];

      T grad11 = jinv 
	  * (drdx[ijk + ele] * rtmpu
	   + dsdx[ijk + ele] * stmpu
	   + dtdx[ijk + ele] * ttmpu);

      T grad12 = jinv 
	  * (drdy[ijk + ele] * rtmpu
	   + dsdy[ijk + ele] * stmpu
	   + dtdy[ijk + ele] * ttmpu);
      
      T grad13 = jinv 
	  * (drdz[ijk + ele] * rtmpu
	   + dsdz[ijk + ele] * stmpu
	   + dtdz[ijk + ele] * ttmpu);

      T grad21 = jinv 
	  * (drdx[ijk + ele] * rtmpv
	   + dsdx[ijk + ele] * stmpv
	   + dtdx[ijk + ele] * ttmpv);

      T grad22 = jinv 
	  * (drdy[ijk + ele] * rtmpv
	   + dsdy[ijk + ele] * stmpv
	   + dtdy[ijk + ele] * ttmpv);
      
      T grad23 = jinv 
	  * (drdz[ijk + ele] * rtmpv
	   + dsdz[ijk + ele] * stmpv
	   + dtdz[ijk + ele] * ttmpv);


      T grad31 = jinv
	  * (drdx[ijk + ele] * rtmpw
	   + dsdx[ijk + ele] * stmpw
	   + dtdx[ijk + ele] * ttmpw);

      T grad32 = jinv
	  * (drdy[ijk + ele] * rtmpw
	   + dsdy[ijk + ele] * stmpw
	   + dtdy[ijk + ele] * ttmpw);
      
      T grad33 = jinv
	  * (drdz[ijk + ele] * rtmpw
	   + dsdz[ijk + ele] * stmpw
	   + dtdz[ijk + ele] * ttmpw);
      lambda2[ijk + e*LX*LX*LX] = eigen_val_calc<T>( grad11, grad12, grad13, grad21, grad22, grad23, grad31, grad32, grad33);
    }
  } 
  
}

template< typename T, const int LX >
__global__ void __launch_bounds__(LX*LX,3)
lambda2_kernel_kstep(T * __restrict__ lambda2,
                    const T * __restrict__ u,
                    const T * __restrict__ v,
                    const T * __restrict__ w,
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
  __shared__ T shv[LX * LX];
  __shared__ T shw[LX * LX];

  __shared__ T shdx[LX * LX];
  __shared__ T shdy[LX * LX];
  __shared__ T shdz[LX * LX];
    
  const int e = blockIdx.x;
  const int j = threadIdx.y;
  const int i = threadIdx.x;
  const int ij = i + j * LX;
  const int ele = e*LX*LX*LX;
  
  shdx[ij] = dx[ij];
  shdy[ij] = dy[ij];
  shdz[ij] = dz[ij];

  T ru[LX];
  T rv[LX];
  T rw[LX];
  
#pragma unroll LX
  for (int k = 0; k < LX; ++k) {
    ru[k] = u[ij + k*LX*LX + ele];
    rv[k] = v[ij + k*LX*LX + ele];
    rw[k] = w[ij + k*LX*LX + ele];
  }
    
  __syncthreads();

  #pragma unroll
  for (int k = 0; k < LX; ++k) {
    const int ijk = ij + k*LX*LX;
    const T jinv = jacinv[ijk+ele];
    T ttmpu = 0.0;
    T ttmpv = 0.0;
    T ttmpw = 0.0;
    shu[ij] = ru[k];
    shv[ij] = rv[k];
    shw[ij] = rw[k];
    for (int l = 0; l < LX; l++) {
      ttmpu += shdz[k+l*LX] * ru[l];
      ttmpv += shdz[k+l*LX] * rv[l];
      ttmpw += shdz[k+l*LX] * rw[l];
    }
    __syncthreads();

    T rtmpu = 0.0;
    T stmpu = 0.0;
    T rtmpv = 0.0;
    T stmpv = 0.0;
    T rtmpw = 0.0;
    T stmpw = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++) {
      rtmpu += shdx[i+l*LX] * shu[l+j*LX];
      stmpu += shdy[j+l*LX] * shu[i+l*LX];
      rtmpv += shdx[i+l*LX] * shv[l+j*LX];
      stmpv += shdy[j+l*LX] * shv[i+l*LX];
      rtmpw += shdx[i+l*LX] * shw[l+j*LX];
      stmpw += shdy[j+l*LX] * shw[i+l*LX];
    }

    T grad11 = jinv * (drdx[ijk + ele] * rtmpu
                          + dsdx[ijk + ele] * stmpu
                          + dtdx[ijk + ele] * ttmpu);

    T grad12 = jinv * (drdy[ijk + ele] * rtmpu
                          + dsdy[ijk + ele] * stmpu
                          + dtdy[ijk + ele] * ttmpu);

    T grad13 = jinv * (drdz[ijk + ele] * rtmpu
                          + dsdz[ijk + ele] * stmpu
                          + dtdz[ijk + ele] * ttmpu);
    T grad21 = jinv * (drdx[ijk + ele] * rtmpv
                          + dsdx[ijk + ele] * stmpv
                          + dtdx[ijk + ele] * ttmpv);

    T grad22 = jinv * (drdy[ijk + ele] * rtmpv
                          + dsdy[ijk + ele] * stmpv
                          + dtdy[ijk + ele] * ttmpv);

    T grad23 = jinv * (drdz[ijk + ele] * rtmpv
                          + dsdz[ijk + ele] * stmpv
                          + dtdz[ijk + ele] * ttmpv);
    T grad31 = jinv * (drdx[ijk + ele] * rtmpw
                          + dsdx[ijk + ele] * stmpw
                          + dtdx[ijk + ele] * ttmpw);

    T grad32 = jinv * (drdy[ijk + ele] * rtmpw
                          + dsdy[ijk + ele] * stmpw
                          + dtdy[ijk + ele] * ttmpw);

    T grad33 = jinv * (drdz[ijk + ele] * rtmpw
                          + dsdz[ijk + ele] * stmpw
                          + dtdz[ijk + ele] * ttmpw);
    lambda2[ijk + ele] = eigen_val_calc<T>( grad11, grad12, grad13, grad21, grad22, grad23, grad31, grad32, grad33);
    __syncthreads();
  }
}

