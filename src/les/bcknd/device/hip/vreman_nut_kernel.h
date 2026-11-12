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

/**
 * Richardson-number correction of Lilly (1962), Ri = N^2 / |S|^2, for the
 * Vreman model, applied to an already computed element-local nut field
 * (before the final gather-scatter average).
 * The velocity gradients a_ij are gather-scatter sums, so they are scaled
 * by mult to obtain the averaged gradient before forming |S|^2.
 * gdT is the interface-averaged g . grad(T), with g pointing upwards
 * (reversed from the case-file input at init), so Ri is single-valued.
 */
template<typename T>
__global__ void vreman_stability(const T * __restrict__ a11,
                                 const T * __restrict__ a12,
                                 const T * __restrict__ a13,
                                 const T * __restrict__ a21,
                                 const T * __restrict__ a22,
                                 const T * __restrict__ a23,
                                 const T * __restrict__ a31,
                                 const T * __restrict__ a32,
                                 const T * __restrict__ a33,
                                 const T * __restrict__ gdT,
                                 T * nut,
                                 const T * __restrict__ mult,
                                 const T ri_c, const T ref_temp,
                                 const T eps, const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    const T m = mult[i];

    // N^2 = (g . grad T) / T0 with g pointing upwards, > 0 when stable
    const T N2 = gdT[i] / ref_temp;

    // |S|^2 = 2 S_ij S_ij of the averaged gradient (the a_ij are
    // gather-scatter sums, hence m*m)
    const T s12 = a12[i] + a21[i];
    const T s13 = a13[i] + a31[i];
    const T s23 = a23[i] + a32[i];
    const T s_sq = ((T) 2.0 * (a11[i]*a11[i] + a22[i]*a22[i] + a33[i]*a33[i])
                    + s12*s12 + s13*s13 + s23*s23) * m * m;

    // Richardson number of Lilly (1962), Ri = N^2 / |S|^2
    const T ri = N2 / (s_sq + eps);
    if (ri <= ri_c) {
      T v = (T) 1.0 - ri / ri_c;
      if (v < (T) 0.0) v = (T) 0.0;
      nut[i] = nut[i] * sqrt(v);
    } else {
      nut[i] = eps;
    }
  }
}
#endif // __COMMON_VREMAN_NUT_KERNEL_H__
