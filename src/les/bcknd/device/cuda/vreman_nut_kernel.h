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
__global__ void hip_vreman_nut_compute(T * __restrict__ a11,
                                      T * __restrict__ a12,
                                      T * __restrict__ a13,
                                      T * __restrict__ a21,
                                      T * __restrict__ a22,
                                      T * __restrict__ a23,
                                      T * __restrict__ a31,
                                      T * __restrict__ a32,
                                      T * __restrict__ a33,
                                      T * __restrict__ delta,
                                      T * __restrict__ nut,
                                      T * __restrict__ mult,
                                      const T c,
                                      const T eps,
                                      const int n){

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    T beta11, beta22, beta33;
    T beta12, beta13, beta23;
    T b_beta, aijaij;

    a11[i] = a11[i] * mult[i];
    a12[i] = a12[i] * mult[i];
    a13[i] = a13[i] * mult[i];
    a21[i] = a21[i] * mult[i];
    a22[i] = a22[i] * mult[i];
    a23[i] = a23[i] * mult[i];
    a31[i] = a31[i] * mult[i];
    a32[i] = a32[i] * mult[i];
    a33[i] = a33[i] * mult[i];

    beta11 = a11[i]*a11[i] + a21[i]*a21[i] + a31[i]*a31[i];
    beta22 = a12[i]*a12[i] + a22[i]*a22[i] + a32[i]*a32[i];
    beta33 = a13[i]*a13[i] + a23[i]*a23[i] + a33[i]*a33[i];
    beta12 = a11[i]*a12[i] + a21[i]*a22[i] + a31[i]*a32[i];
    beta13 = a11[i]*a13[i] + a21[i]*a23[i] + a31[i]*a33[i];
    beta23 = a12[i]*a13[i] + a22[i]*a23[i] + a32[i]*a33[i];

    b_beta = beta11*beta22 - beta12*beta12
           + beta11*beta33 - beta13*beta13
           + beta22*beta33 - beta23*beta23;
    
    b_beta = max(0.0, b_beta);

    aijaij = beta11 + beta22 + beta33;

    nut[i] = c * delta[i] * delta[i]
             * sqrt(b_beta / (aijaij + eps));
  }
}
#endif // __COMMON_VREMAN_NUT_KERNEL_H__