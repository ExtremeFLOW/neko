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
__global__ void hip_sigma_nut_compute(T * __restrict__ g11,
                                      T * __restrict__ g12,
                                      T * __restrict__ g13,
                                      T * __restrict__ g21,
                                      T * __restrict__ g22,
                                      T * __restrict__ g23,
                                      T * __restrict__ g31,
                                      T * __restrict__ g32,
                                      T * __restrict__ g33,
                                      T * __restrict__ delta,
                                      T * __restrict__ nut,
                                      T * __restrict__ mult,
                                      const T c,
                                      const T eps,
                                      const int n){

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    T sigG11, sigG12, sigG13, sigG22, sigG23, sigG33;
    T sigma1, sigma2, sigma3;
    T Invariant1, Invariant2, Invariant3;
    T alpha1, alpha2, alpha3;
    T Dsigma;
    T pi_3;
    T tmp1;
    
    pi_3 = 4.0/3.0*atan(1.0);

    g11[i] = g11[i] * mult[i];
    g12[i] = g12[i] * mult[i];
    g13[i] = g13[i] * mult[i];
    g21[i] = g21[i] * mult[i];
    g22[i] = g22[i] * mult[i];
    g23[i] = g23[i] * mult[i];
    g31[i] = g31[i] * mult[i];
    g32[i] = g32[i] * mult[i];
    g33[i] = g33[i] * mult[i];

    sigG11 = g11[i]*g11[i] + g21[i]*g21[i] + g31[i]*g31[i];
    sigG22 = g12[i]*g12[i] + g22[i]*g22[i] + g32[i]*g32[i];
    sigG33 = g13[i]*g13[i] + g23[i]*g23[i] + g33[i]*g33[i];
    sigG12 = g11[i]*g12[i] + \
             g21[i]*g22[i] + \
             g31[i]*g32[i];
    sigG13 = g11[i]*g13[i] + \
             g21[i]*g23[i] + \
             g31[i]*g33[i];
    sigG23 = g12[i]*g13[i] + \
             g22[i]*g23[i] + \
             g32[i]*g33[i];

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

        tmp1 = alpha2/pow(alpha1,(3.0/2.0));

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

    nut[i] = (c*delta[i])*(c*delta[i]) * Dsigma;
  }
}
#endif // __COMMON_SIGMA_NUT_KERNEL_H__