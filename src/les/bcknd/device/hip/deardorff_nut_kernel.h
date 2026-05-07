#ifndef __COMMON_DEARDORFF_NUT_KERNEL_H__
#define __COMMON_DEARDORFF_NUT_KERNEL_H__
/*
 Copyright (c) 2025, The Neko Authors
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
 * Device kernel for deardorff_nut_compute
 */
#include <cmath>
#include <algorithm>
template< typename T>
__global__ void deardorff_nut_compute(T *__restrict__ TKE,
                                    const T *__restrict__ dTdx,
                                    const T *__restrict__ dTdy,
                                    const T *__restrict__ dTdz,
                                    const T * __restrict__ a11,
                                    const T * __restrict__ a12,
                                    const T * __restrict__ a13,
                                    const T * __restrict__ a21,
                                    const T * __restrict__ a22,
                                    const T * __restrict__ a23,
                                    const T * __restrict__ a31,
                                    const T * __restrict__ a32,
                                    const T * __restrict__ a33,
                                    const T * __restrict__ delta,
                                    T * __restrict__ nut,
                                    T * __restrict__ temperature_alphat,
                                    T * __restrict__ TKE_alphat,
                                    T * __restrict__ TKE_source,
                                    const T c_k,
                                    const T T0,
                                    const T g1,
                                    const T g2,
                                    const T g3,
                                    const T eps,
                                    const int n){

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  T s11, s22, s33;
  T s12, s13, s23;
  T N2, l;
  T shear, buoyancy, dissipation;

  for (int i = idx; i < n; i += str) {
    const T dTdx_r = dTdx[i];
    const T dTdy_r = dTdy[i];
    const T dTdz_r = dTdz[i];
    const T a11_r = a11[i];
    const T a12_r = a12[i];
    const T a13_r = a13[i];
    const T a21_r = a21[i];
    const T a22_r = a22[i];
    const T a23_r = a23[i];
    const T a31_r = a31[i];
    const T a32_r = a32[i];
    const T a33_r = a33[i];
    const T delta_r = delta[i];

    if (TKE[i] < eps) {
      TKE[i] = eps;
    }
    const T TKE_r = TKE[i];

    N2 = (dTdx_r * g1 + dTdy_r * g2 + dTdz_r * g3) / T0;
    if (N2 > 0.0) {
      l = 0.76 * sqrt(TKE_r / N2);
      l = min(l, delta_r);
    }
    else {
      l = delta_r;
    }

    nut[i] = c_k * l * sqrt(TKE_r);
    temperature_alphat[i] = (1.0 + 2.0 * l/delta_r) * nut[i];
    TKE_alphat[i] = 2.0 * nut[i];

    s11 = a11_r + a11_r;
    s22 = a22_r + a22_r;
    s33 = a33_r + a33_r;
    s12 = a12_r + a21_r;
    s13 = a13_r + a31_r;
    s23 = a23_r + a32_r;

    shear = nut[i] * (s11*a11_r + s12*a12_r + s13*a13_r
                    + s12*a21_r + s22*a22_r + s23*a23_r
                    + s13*a31_r + s23*a32_r + s33*a33_r);
    buoyancy = - (dTdx_r * g1 + 
                  dTdy_r * g2 + 
                  dTdz_r * g3) * temperature_alphat[i] / T0;
    dissipation = - (0.19 + 0.74 * l/delta_r) * sqrt(TKE_r*TKE_r*TKE_r) / l;

    TKE_source[i] = shear + buoyancy + dissipation;
  }
}
#endif // __COMMON_DEARDORFF_NUT_KERNEL_H__