#ifndef __COMMON_WALE_NUT_KERNEL_H__
#define __COMMON_WALE_NUT_KERNEL_H__
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
 * Device kernel for wale_nut_compute
 */
#include <cmath>
#include <algorithm>
template< typename T>
__global__ void wale_nut_compute(const T * __restrict__ g11,
                                 const T * __restrict__ g12,
                                 const T * __restrict__ g13,
                                 const T * __restrict__ g21,
                                 const T * __restrict__ g22,
                                 const T * __restrict__ g23,
                                 const T * __restrict__ g31,
                                 const T * __restrict__ g32,
                                 const T * __restrict__ g33,
                                 const T * __restrict__ delta,
                                 T * __restrict__ nut,
                                 const T * __restrict__ mult,
                                 const T c,
                                 const T eps,
                                 const int n){

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  T s11, s22, s33;
  T s12, s13, s23;
  T gsqr_11, gsqr_12, gsqr_13;
  T gsqr_21, gsqr_22, gsqr_23;
  T gsqr_31, gsqr_32, gsqr_33;
  T sd11, sd22, sd33;
  T sd12, sd13, sd23;
  T Sdij_Sdij, Sij_Sij, OP_wale;

  for (int i = idx; i < n; i += str) {
    
    const T g11_r = g11[i];
    const T g12_r = g12[i];
    const T g13_r = g13[i];
    const T g21_r = g21[i];
    const T g22_r = g22[i];
    const T g23_r = g23[i];
    const T g31_r = g31[i];
    const T g32_r = g32[i];
    const T g33_r = g33[i];
    const T delta_r = delta[i];

    s11 = g11_r;
    s22 = g22_r;
    s33 = g33_r;
    s12 = 0.5 * (g12_r + g21_r);
    s13 = 0.5 * (g13_r + g31_r);
    s23 = 0.5 * (g23_r + g32_r);

    gsqr_11 = g11_r * g11_r 
            + g12_r * g21_r 
            + g13_r * g31_r;
    gsqr_12 = g11_r * g12_r 
            + g12_r * g22_r 
            + g13_r * g32_r;
    gsqr_13 = g11_r * g13_r 
            + g12_r * g23_r 
            + g13_r * g33_r;
    gsqr_21 = g21_r * g11_r 
            + g22_r * g21_r 
            + g23_r * g31_r;
    gsqr_22 = g21_r * g12_r 
            + g22_r * g22_r 
            + g23_r * g32_r;
    gsqr_23 = g21_r * g13_r 
            + g22_r * g23_r 
            + g23_r * g33_r;
    gsqr_31 = g31_r * g11_r 
            + g32_r * g21_r 
            + g33_r * g31_r;
    gsqr_32 = g31_r * g12_r 
            + g32_r * g22_r 
            + g33_r * g32_r;
    gsqr_33 = g31_r * g13_r 
            + g32_r * g23_r 
            + g33_r * g33_r;
    
    sd11 = gsqr_11 - ( (gsqr_11 + gsqr_22 + gsqr_33) / 3.0);
    sd22 = gsqr_22 - ( (gsqr_11 + gsqr_22 + gsqr_33) / 3.0);
    sd33 = gsqr_33 - ( (gsqr_11 + gsqr_22 + gsqr_33) / 3.0);
    sd12 = 0.5 * (gsqr_12 + gsqr_21);
    sd13 = 0.5 * (gsqr_13 + gsqr_31);
    sd23 = 0.5 * (gsqr_23 + gsqr_32);
    
    Sdij_Sdij = sd11*sd11 + sd22*sd22 + sd33*sd33 +
         2.0 * (sd12*sd12 + sd13*sd13 + sd23*sd23);
    Sij_Sij = s11*s11 + s22*s22 + s33*s33 +
         2.0 * (s12*s12 + s13*s13 + s23*s23);
    OP_wale = sqrt(Sdij_Sdij*Sdij_Sdij*Sdij_Sdij) / 
         max((sqrt(Sij_Sij*Sij_Sij*Sij_Sij*Sij_Sij*Sij_Sij) 
            + sqrt(sqrt(Sdij_Sdij*Sdij_Sdij*Sdij_Sdij*Sdij_Sdij*Sdij_Sdij))), 
            eps);


    nut[i] = c * c * delta_r * delta_r
             * OP_wale * mult[i];
  }
}
#endif // __COMMON_WALE_NUT_KERNEL_H__