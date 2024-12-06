#ifndef __COMMON_DYNAMIC_SMAGORINSKY_NUT_KERNEL_H__
#define __COMMON_DYNAMIC_SMAGORINSKY_NUT_KERNEL_H__
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
 * Device kernel for dynamic_smagorinsky
 */
#include <cmath>
#include <algorithm>
template< typename T>
__global__ void s_abs_compute(T * __restrict__ s_abs,
                              T * __restrict__ s11,
                              T * __restrict__ s22,
                              T * __restrict__ s33,
                              T * __restrict__ s12,
                              T * __restrict__ s13,
                              T * __restrict__ s23,
                              T * __restrict__ mult,
                              const int n){

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    s11[i] = s11[i] * mult[i];
    s22[i] = s22[i] * mult[i];
    s33[i] = s33[i] * mult[i];
    s12[i] = s12[i] * mult[i];
    s13[i] = s13[i] * mult[i];
    s23[i] = s23[i] * mult[i];

    s_abs[i] = sqrt(2.0 * (s11[i] * s11[i] +
                        s22[i] * s22[i] +
                        s33[i] * s33[i]) + 
                 4.0 * (s12[i] * s12[i] +
                        s13[i] * s13[i] +
                        s23[i] * s23[i]));

  }
}

template< typename T>
__global__ void lij_compute_part1(T * __restrict__ l11,
                                  T * __restrict__ l22,
                                  T * __restrict__ l33,
                                  T * __restrict__ l12,
                                  T * __restrict__ l13,
                                  T * __restrict__ l23,
                                  T * __restrict__ u,
                                  T * __restrict__ v,
                                  T * __restrict__ w,
                                  T * __restrict__ fu,
                                  T * __restrict__ fv,
                                  T * __restrict__ fw,
                                  T * __restrict__ fuu,
                                  T * __restrict__ fvv,
                                  T * __restrict__ fww,
                                  T * __restrict__ fuv,
                                  T * __restrict__ fuw,
                                  T * __restrict__ fvw,
                                  const int n){

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    l11[i] = fu[i] * fu[i];
    l22[i] = fv[i] * fv[i];
    l33[i] = fw[i] * fw[i];
    l12[i] = fu[i] * fv[i];
    l13[i] = fu[i] * fw[i];
    l23[i] = fv[i] * fw[i];

    fuu[i] = u[i] * u[i];
    fvv[i] = v[i] * v[i];
    fww[i] = w[i] * w[i];
    fuv[i] = u[i] * v[i];
    fuw[i] = u[i] * w[i];
    fvw[i] = v[i] * w[i];
  }
}

template< typename T>
__global__ void lij_compute_part2(T * __restrict__ l11,
                                  T * __restrict__ l22,
                                  T * __restrict__ l33,
                                  T * __restrict__ l12,
                                  T * __restrict__ l13,
                                  T * __restrict__ l23,
                                  T * __restrict__ fuu,
                                  T * __restrict__ fvv,
                                  T * __restrict__ fww,
                                  T * __restrict__ fuv,
                                  T * __restrict__ fuw,
                                  T * __restrict__ fvw,
                                  const int n){

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    l11[i] -= fuu[i];
    l22[i] -= fvv[i];
    l33[i] -= fww[i];
    l12[i] -= fuv[i];
    l13[i] -= fuw[i];
    l23[i] -= fvw[i];
  }
}

template< typename T>
__global__ void mij_compute_part1(T * __restrict__ m11,
                                  T * __restrict__ m22,
                                  T * __restrict__ m33,
                                  T * __restrict__ m12,
                                  T * __restrict__ m13,
                                  T * __restrict__ m23,
                                  T * __restrict__ s_abs,
                                  T * __restrict__ s11,
                                  T * __restrict__ s22,
                                  T * __restrict__ s33,
                                  T * __restrict__ s12,
                                  T * __restrict__ s13,
                                  T * __restrict__ s23,
                                  T * __restrict__ fs_abs,
                                  T * __restrict__ fs11,
                                  T * __restrict__ fs22,
                                  T * __restrict__ fs33,
                                  T * __restrict__ fs12,
                                  T * __restrict__ fs13,
                                  T * __restrict__ fs23,
                                  T * __restrict__ fsabss11,
                                  T * __restrict__ fsabss22,
                                  T * __restrict__ fsabss33,
                                  T * __restrict__ fsabss12,
                                  T * __restrict__ fsabss13,
                                  T * __restrict__ fsabss23,
                                  const T delta_ratio2,
                                  const int n){

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    m11[i] = delta_ratio2 * (fs_abs[i] * fs11[i]);
    m22[i] = delta_ratio2 * (fs_abs[i] * fs22[i]);
    m33[i] = delta_ratio2 * (fs_abs[i] * fs33[i]);
    m12[i] = delta_ratio2 * (fs_abs[i] * fs12[i]);
    m13[i] = delta_ratio2 * (fs_abs[i] * fs13[i]);
    m23[i] = delta_ratio2 * (fs_abs[i] * fs23[i]);

    fsabss11[i] = s_abs[i] * s11[i];
    fsabss22[i] = s_abs[i] * s22[i];
    fsabss33[i] = s_abs[i] * s33[i];
    fsabss12[i] = s_abs[i] * s12[i];
    fsabss13[i] = s_abs[i] * s13[i];
    fsabss23[i] = s_abs[i] * s23[i];
  }
}

template< typename T>
__global__ void mij_nut_compute_part2(T * __restrict__ m11,
                                      T * __restrict__ m22,
                                      T * __restrict__ m33,
                                      T * __restrict__ m12,
                                      T * __restrict__ m13,
                                      T * __restrict__ m23,
                                      T * __restrict__ l11,
                                      T * __restrict__ l22,
                                      T * __restrict__ l33,
                                      T * __restrict__ l12,
                                      T * __restrict__ l13,
                                      T * __restrict__ l23,
                                      T * __restrict__ fsabss11,
                                      T * __restrict__ fsabss22,
                                      T * __restrict__ fsabss33,
                                      T * __restrict__ fsabss12,
                                      T * __restrict__ fsabss13,
                                      T * __restrict__ fsabss23,
                                      T * __restrict__ num,
                                      T * __restrict__ den,
                                      T * __restrict__ c_dyn,
                                      T * __restrict__ delta,
                                      T * __restrict__ s_abs,
                                      T * __restrict__ nut,
                                      const T alpha,
                                      const int n){

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    T num_curr, den_curr;
    T delta2 = delta[i] * delta[i];

    m11[i] = delta2 * (m11[i] - fsabss11[i]);
    m22[i] = delta2 * (m22[i] - fsabss22[i]);
    m33[i] = delta2 * (m33[i] - fsabss33[i]);
    m12[i] = delta2 * (m12[i] - fsabss12[i]);
    m13[i] = delta2 * (m13[i] - fsabss13[i]);
    m23[i] = delta2 * (m23[i] - fsabss23[i]);

    num_curr = m11[i] * l11[i]
             + m22[i] * l22[i]
             + m33[i] * l33[i]
      + 2.0 * (m12[i] * l12[i] 
             + m13[i] * l13[i]
             + m23[i] * l23[i]);
    
    den_curr = m11[i] * m11[i]
             + m22[i] * m22[i]
             + m33[i] * m33[i]
      + 2.0 * (m12[i] * m12[i] 
             + m13[i] * m13[i]
             + m23[i] * m23[i]);

    num[i] = alpha * num[i] + (1.0 - alpha) * num_curr;
    den[i] = alpha * den[i] + (1.0 - alpha) * den_curr;

    if (den[i] > 0.0) {
        c_dyn[i] = 0.5 * num[i] / den[i];
    } else {
        c_dyn[i] = 0.0;
    }

    c_dyn[i] = max(c_dyn[i], 0.0);
    nut[i] = c_dyn[i] * delta2 * s_abs[i];

  }
}
#endif // __COMMON_DYNAMIC_SMAGORINSKY_NUT_KERNEL_H__