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

template<typename T>
__global__ void vreman_nut_compute_buoy(const T * __restrict__ a11,
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
                                        const int n,
                                        const T * __restrict__ dTdx,
                                        const T * __restrict__ dTdy,
                                        const T * __restrict__ dTdz,
                                        const T n1, const T n2, const T n3,
                                        const T g1, const T g2, const T g3,
                                        const T ri_c, const T ref_temp){
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  T beta11, beta22, beta33;
  T beta12, beta13, beta23;
  T b_beta, aijaij;
  T nut0;

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

    nut0 = c * delta_r * delta_r
             * sqrt(b_beta / (aijaij + eps)) * mult[i];

        // Scalar gradient for buoyancy
        T buoyancy = (g1*dTdx[i] + g2*dTdy[i] + g3*dTdz[i]) / ref_temp;

        // Directional derivative along n
        T du_n1 = a11_r*n1 + a12_r*n2 + a13_r*n3;
        T du_n2 = a21_r*n1 + a22_r*n2 + a23_r*n3;
        T du_n3 = a31_r*n1 + a32_r*n2 + a33_r*n3;

        T du_parallel = du_n1*n1 + du_n2*n2 + du_n3*n3;

        T sh1 = du_n1 - du_parallel*n1;
        T sh2 = du_n2 - du_parallel*n2;
        T sh3 = du_n3 - du_parallel*n3;

        T shear_sq = sh1*sh1 + sh2*sh2 + sh3*sh3;

        T ri = buoyancy / (shear_sq + eps);

        T out;
        if (ri <= ri_c) {
            T v = (T)1 - ri/ri_c;
            if (v < (T)0) v = (T)0;
            out = nut0 * sqrt(v);
        } else {
            out = eps;
        }

        nut[i] = out;
    }
}
#endif // __COMMON_VREMAN_NUT_KERNEL_H__
