#ifndef __COMMON_SMAGORINSKY_NUT_KERNEL_H__
#define __COMMON_SMAGORINSKY_NUT_KERNEL_H__
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
 * Device kernel for smagorinsky_nut_compute
 */
#include <cmath>
#include <algorithm>
template< typename T>
__global__ void smagorinsky_nut_compute(T * __restrict__ s11,
                                      T * __restrict__ s22,
                                      T * __restrict__ s33,
                                      T * __restrict__ s12,
                                      T * __restrict__ s13,
                                      T * __restrict__ s23,
                                      T * __restrict__ delta,
                                      T * __restrict__ nut,
                                      T * __restrict__ mult,
                                      const T c_s,
                                      const int n){

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    T s_abs;

    s11[i] = s11[i] * mult[i];
    s22[i] = s22[i] * mult[i];
    s33[i] = s33[i] * mult[i];
    s12[i] = s12[i] * mult[i];
    s13[i] = s13[i] * mult[i];
    s23[i] = s23[i] * mult[i];

    s_abs = sqrt(2.0 * (s11[i] * s11[i] +
                        s22[i] * s22[i] +
                        s33[i] * s33[i]) + 
                 4.0 * (s12[i] * s12[i] +
                        s13[i] * s13[i] +
                        s23[i] * s23[i]));

    nut[i] = c_s * c_s * delta[i] * delta[i] * s_abs;
  }
}
#endif // __COMMON_SMAGORINSKY_NUT_KERNEL_H__