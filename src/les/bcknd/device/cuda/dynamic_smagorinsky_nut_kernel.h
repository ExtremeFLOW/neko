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
                              const T * __restrict__ s11,
                              const T * __restrict__ s22,
                              const T * __restrict__ s33,
                              const T * __restrict__ s12,
                              const T * __restrict__ s13,
                              const T * __restrict__ s23,
                              const int n){

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    const T s11_r = s11[i];
    const T s22_r = s22[i];
    const T s33_r = s33[i];
    const T s12_r = s12[i];
    const T s13_r = s13[i];
    const T s23_r = s23[i];

    s_abs[i] = sqrt(2.0 * (s11_r * s11_r +
                        s22_r * s22_r +
                        s33_r * s33_r) + 
                 4.0 * (s12_r * s12_r +
                        s13_r * s13_r +
                        s23_r * s23_r));

  }
}

template< typename T>
__global__ void lij_compute_part1(T * __restrict__ l11,
                                  T * __restrict__ l22,
                                  T * __restrict__ l33,
                                  T * __restrict__ l12,
                                  T * __restrict__ l13,
                                  T * __restrict__ l23,
                                  const T * __restrict__ u,
                                  const T * __restrict__ v,
                                  const T * __restrict__ w,
                                  const T * __restrict__ fu,
                                  const T * __restrict__ fv,
                                  const T * __restrict__ fw,
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
    const T u_r = u[i];
    const T v_r = v[i];
    const T w_r = w[i];
    const T fu_r = fu[i];
    const T fv_r = fv[i];
    const T fw_r = fw[i];

    l11[i] = fu_r * fu_r;
    l22[i] = fv_r * fv_r;
    l33[i] = fw_r * fw_r;
    l12[i] = fu_r * fv_r;
    l13[i] = fu_r * fw_r;
    l23[i] = fv_r * fw_r;

    fuu[i] = u_r * u_r;
    fvv[i] = v_r * v_r;
    fww[i] = w_r * w_r;
    fuv[i] = u_r * v_r;
    fuw[i] = u_r * w_r;
    fvw[i] = v_r * w_r;
  }
}

template< typename T>
__global__ void lij_compute_part2(T * __restrict__ l11,
                                  T * __restrict__ l22,
                                  T * __restrict__ l33,
                                  T * __restrict__ l12,
                                  T * __restrict__ l13,
                                  T * __restrict__ l23,
                                  const T * __restrict__ fuu,
                                  const T * __restrict__ fvv,
                                  const T * __restrict__ fww,
                                  const T * __restrict__ fuv,
                                  const T * __restrict__ fuw,
                                  const T * __restrict__ fvw,
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
                                  const T * __restrict__ s_abs,
                                  const T * __restrict__ s11,
                                  const T * __restrict__ s22,
                                  const T * __restrict__ s33,
                                  const T * __restrict__ s12,
                                  const T * __restrict__ s13,
                                  const T * __restrict__ s23,
                                  const T * __restrict__ fs_abs,
                                  const T * __restrict__ fs11,
                                  const T * __restrict__ fs22,
                                  const T * __restrict__ fs33,
                                  const T * __restrict__ fs12,
                                  const T * __restrict__ fs13,
                                  const T * __restrict__ fs23,
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
    const T s_abs_r = s_abs[i];
    const T fs_abs_r = fs_abs[i];

    m11[i] = delta_ratio2 * (fs_abs_r * fs11[i]);
    m22[i] = delta_ratio2 * (fs_abs_r * fs22[i]);
    m33[i] = delta_ratio2 * (fs_abs_r * fs33[i]);
    m12[i] = delta_ratio2 * (fs_abs_r * fs12[i]);
    m13[i] = delta_ratio2 * (fs_abs_r * fs13[i]);
    m23[i] = delta_ratio2 * (fs_abs_r * fs23[i]);

    fsabss11[i] = s_abs_r * s11[i];
    fsabss22[i] = s_abs_r * s22[i];
    fsabss33[i] = s_abs_r * s33[i];
    fsabss12[i] = s_abs_r * s12[i];
    fsabss13[i] = s_abs_r * s13[i];
    fsabss23[i] = s_abs_r * s23[i];
  }
}

template< typename T>
__global__ void mij_nut_compute_part2(const T * __restrict__ m11,
                                      const T * __restrict__ m22,
                                      const T * __restrict__ m33,
                                      const T * __restrict__ m12,
                                      const T * __restrict__ m13,
                                      const T * __restrict__ m23,
                                      T * __restrict__ l11,
                                      T * __restrict__ l22,
                                      T * __restrict__ l33,
                                      T * __restrict__ l12,
                                      T * __restrict__ l13,
                                      T * __restrict__ l23,
                                      const T * __restrict__ fsabss11,
                                      const T * __restrict__ fsabss22,
                                      const T * __restrict__ fsabss33,
                                      const T * __restrict__ fsabss12,
                                      const T * __restrict__ fsabss13,
                                      const T * __restrict__ fsabss23,
                                      T * __restrict__ num,
                                      T * __restrict__ den,
                                      T * __restrict__ c_dyn,
                                      const T * __restrict__ delta,
                                      const T * __restrict__ s_abs,
                                      T * __restrict__ nut,
                                      const T alpha,
                                      const int n){

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  T num_curr, den_curr;

  for (int i = idx; i < n; i += str) {
    const T delta_r = delta[i];
    const T delta2 =  delta_r * delta_r;
    T c_dyn_r;

    const T m11_r = delta2 * (m11[i] - fsabss11[i]);
    const T m22_r = delta2 * (m22[i] - fsabss22[i]);
    const T m33_r = delta2 * (m33[i] - fsabss33[i]);
    const T m12_r = delta2 * (m12[i] - fsabss12[i]);
    const T m13_r = delta2 * (m13[i] - fsabss13[i]);
    const T m23_r = delta2 * (m23[i] - fsabss23[i]);

    num_curr = m11_r * l11[i]
             + m22_r * l22[i]
             + m33_r * l33[i]
      + 2.0 * (m12_r * l12[i]
             + m13_r * l13[i]
             + m23_r * l23[i]);
    
    den_curr = m11_r * m11_r
             + m22_r * m22_r
             + m33_r * m33_r
      + 2.0 * (m12_r * m12_r 
             + m13_r * m13_r
             + m23_r * m23_r);

    const T num_r = alpha * num[i] + (1.0 - alpha) * num_curr;
    const T den_r = alpha * den[i] + (1.0 - alpha) * den_curr;

    num[i] = num_r;
    den[i] = den_r;

    if (den_r > 0.0) {
        c_dyn_r = 0.5 * num_r / den_r;
    } else {
        c_dyn_r = 0.0;
    }

    c_dyn_r = max(c_dyn_r, 0.0);
    c_dyn[i] = c_dyn_r;
    nut[i] = c_dyn_r * delta2 * s_abs[i];

  }
}
#endif // __COMMON_DYNAMIC_SMAGORINSKY_NUT_KERNEL_H__