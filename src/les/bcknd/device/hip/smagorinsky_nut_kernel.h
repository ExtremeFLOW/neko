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
__global__ void smagorinsky_nut_compute(const T * __restrict__ s11,
                                      const T * __restrict__ s22,
                                      const T * __restrict__ s33,
                                      const T * __restrict__ s12,
                                      const T * __restrict__ s13,
                                      const T * __restrict__ s23,
                                      const T * __restrict__ delta,
                                      T * __restrict__ nut,
                                      const T * __restrict__ mult,
                                      const T c_s,
                                      const int n){

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  T s_abs;

  for (int i = idx; i < n; i += str) {

    const T s11_r = s11[i];
    const T s22_r = s22[i];
    const T s33_r = s33[i];
    const T s12_r = s12[i];
    const T s13_r = s13[i];
    const T s23_r = s23[i];
    const T delta_r = delta[i];
    
    s_abs = sqrt(2.0 * (s11_r * s11_r +
                        s22_r * s22_r +
                        s33_r * s33_r) + 
                 4.0 * (s12_r * s12_r +
                        s13_r * s13_r +
                        s23_r * s23_r));

    nut[i] = c_s * c_s * delta_r * delta_r * s_abs * mult[i];
  }
}
#endif // __COMMON_SMAGORINSKY_NUT_KERNEL_H__