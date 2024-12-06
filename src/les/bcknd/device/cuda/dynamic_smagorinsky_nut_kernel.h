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
#endif // __COMMON_DYNAMIC_SMAGORINSKY_NUT_KERNEL_H__