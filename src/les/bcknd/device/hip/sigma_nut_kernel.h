#ifndef __COMMON_SIGMA_NUT_KERNEL_H__
#define __COMMON_SIGMA_NUT_KERNEL_H__
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
 * Device kernel for sigma_nut_compute
 */
#include <cmath>
#include <algorithm>
template< typename T>
__global__ void sigma_nut_compute(const T * __restrict__ g11,
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

  T sigG11, sigG12, sigG13, sigG22, sigG23, sigG33;
  T sigma1, sigma2, sigma3;
  T Invariant1, Invariant2, Invariant3;
  T alpha1, alpha2, alpha3;
  T Dsigma;
  const T pi_3 = 4.0/3.0*atan(1.0);
  T tmp1;

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

    sigG11 = g11_r*g11_r + g21_r*g21_r + g31_r*g31_r;
    sigG22 = g12_r*g12_r + g22_r*g22_r + g32_r*g32_r;
    sigG33 = g13_r*g13_r + g23_r*g23_r + g33_r*g33_r;
    sigG12 = g11_r*g12_r + \
             g21_r*g22_r + \
             g31_r*g32_r;
    sigG13 = g11_r*g13_r + \
             g21_r*g23_r + \
             g31_r*g33_r;
    sigG23 = g12_r*g13_r + \
             g22_r*g23_r + \
             g32_r*g33_r;

    if (abs(sigG11) < eps) {
        sigG11 = 0.0;
    }
    if (abs(sigG12) < eps) {
        sigG12 = 0.0;
    }
    if (abs(sigG13) < eps) {
        sigG13 = 0.0;
    }
    if (abs(sigG22) < eps) {
        sigG22 = 0.0;
    }
    if (abs(sigG23) < eps) {
        sigG23 = 0.0;
    }
    if (abs(sigG33) < eps) {
        sigG33 = 0.0;
    }

    if (abs(sigG12*sigG12 + \
            sigG13*sigG13 + sigG23*sigG23) < eps) {
        sigma1 = sqrt(max(max(max(sigG11, sigG22), sigG33), 0.0));
        sigma3 = sqrt(max(min(min(sigG11, sigG22), sigG33), 0.0));
        Invariant1 = sigG11 + sigG22 + sigG33;
        sigma2 = sqrt(abs(Invariant1 - sigma1*sigma1 - sigma3*sigma3));
    } else { 
        Invariant1 = sigG11 + sigG22 + sigG33;
        Invariant2 = sigG11*sigG22 + sigG11*sigG33 + sigG22*sigG33 - \
                (sigG12*sigG12 + sigG13*sigG13 + sigG23*sigG23);
        Invariant3 = sigG11*sigG22*sigG33 + \
                    2.0*sigG12*sigG13*sigG23 - \
                (sigG11*sigG23*sigG23 + sigG22*sigG13*sigG13 + \
                sigG33*sigG12*sigG12);

        Invariant1 = max(Invariant1, 0.0);
        Invariant2 = max(Invariant2, 0.0);
        Invariant3 = max(Invariant3, 0.0);

        alpha1 = Invariant1*Invariant1/9.0 - Invariant2/3.0;

        alpha1 = max(alpha1, 0.0);

        alpha2 = Invariant1*Invariant1*Invariant1/27.0 - \
            Invariant1*Invariant2/6.0 + Invariant3/2.0;

        tmp1 = alpha2/sqrt(alpha1 * alpha1 * alpha1);

        if (tmp1 <= -1.0) {
            sigma1 = sqrt(max(Invariant1/3.0 + sqrt(alpha1), 0.0));
            sigma2 = sigma1;
            sigma3 = sqrt(Invariant1/3.0 - 2.0*sqrt(alpha1));

        } else if (tmp1 >= 1.0) {
            sigma1 = sqrt(max(Invariant1/3.0 + 2.0*sqrt(alpha1), 0.0));
            sigma2 = sqrt(Invariant1/3.0 - sqrt(alpha1));
            sigma3 = sigma2;

        } else {
            alpha3 = acos(tmp1)/3.0;

            if (abs(Invariant3) < eps) {
                sigma1 = sqrt(max(Invariant1/3.0 + \
                            2.0*sqrt(alpha1)*cos(alpha3), 0.0));
                sigma2 = sqrt(abs(Invariant1 - sigma1*sigma1));
                sigma3 = 0.0;
            } else {
                sigma1 = sqrt(max(Invariant1/3.0 + \
                            2.0*sqrt(alpha1)*cos(alpha3), 0.0));
                sigma2 = sqrt(Invariant1/3.0 - \
                            2.0*sqrt(alpha1)*cos(pi_3 + alpha3));
                sigma3 = sqrt(abs(Invariant1 - \
                                sigma1*sigma1-sigma2*sigma2));
            }
        }
    }

    if (sigma1 > 0.0) {
        Dsigma = sigma3*(sigma1 - sigma2)*(sigma2 - sigma3)/(sigma1*sigma1);
    } else {
        Dsigma = 0.0;
    }

    Dsigma = max(Dsigma, 0.0);

    nut[i] = (c*delta_r)*(c*delta_r) * Dsigma * mult[i];
  }
}
#endif // __COMMON_SIGMA_NUT_KERNEL_H__