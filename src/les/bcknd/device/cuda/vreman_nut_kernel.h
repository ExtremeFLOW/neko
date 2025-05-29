#ifndef __COMMON_VREMAN_NUT_KERNEL_H__
#define __COMMON_VREMAN_NUT_KERNEL_H__
/*
 Copyright (c) 2024, The Neko Authors
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
 * Device kernel for vreman_nut_compute
 */
#include <cmath>
#include <algorithm>
template< typename T>
__global__ void vreman_nut_compute(const T * __restrict__ a11,
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
                                      const T * __restrict__ mult,
                                      const T c,
                                      const T eps,
                                      const int n){

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  T beta11, beta22, beta33;
  T beta12, beta13, beta23;
  T b_beta, aijaij;

  for (int i = idx; i < n; i += str) {
    
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

    beta11 = a11_r*a11_r + a21_r*a21_r + a31_r*a31_r;
    beta22 = a12_r*a12_r + a22_r*a22_r + a32_r*a32_r;
    beta33 = a13_r*a13_r + a23_r*a23_r + a33_r*a33_r;
    beta12 = a11_r*a12_r + a21_r*a22_r + a31_r*a32_r;
    beta13 = a11_r*a13_r + a21_r*a23_r + a31_r*a33_r;
    beta23 = a12_r*a13_r + a22_r*a23_r + a32_r*a33_r;

    b_beta = beta11*beta22 - beta12*beta12
           + beta11*beta33 - beta13*beta13
           + beta22*beta33 - beta23*beta23;
    
    b_beta = max(0.0, b_beta);

    aijaij = beta11 + beta22 + beta33;

    nut[i] = c * delta_r * delta_r
             * sqrt(b_beta / (aijaij + eps)) * mult[i];
  }
}
#endif // __COMMON_VREMAN_NUT_KERNEL_H__