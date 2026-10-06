#ifndef __LES_SIGMA_NUT_KERNEL_CL__
#define __LES_SIGMA_NUT_KERNEL_CL__
/*
 Copyright (c) 2026, The Neko Authors
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
__kernel void sigma_nut_compute_kernel(__global const real * __restrict__ g11,
                                       __global const real * __restrict__ g12,
                                       __global const real * __restrict__ g13,
                                       __global const real * __restrict__ g21,
                                       __global const real * __restrict__ g22,
                                       __global const real * __restrict__ g23,
                                       __global const real * __restrict__ g31,
                                       __global const real * __restrict__ g32,
                                       __global const real * __restrict__ g33,
                                       __global const real * __restrict__ delta,
                                       __global real * __restrict__ nut,
                                       __global const real * __restrict__ mult,
                                       const real c,
                                       const real eps,
                                       const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  const real zero = (real) 0.0;
  const real one = (real) 1.0;
  const real two = (real) 2.0;
  const real pi_3 = (real) 4.0 / (real) 3.0 * atan(one);

  real sigG11, sigG12, sigG13, sigG22, sigG23, sigG33;
  real sigma1, sigma2, sigma3;
  real Invariant1, Invariant2, Invariant3;
  real alpha1, alpha2, alpha3;
  real Dsigma;
  real tmp1;

  for (int i = idx; i < n; i += str) {

    const real g11_r = g11[i];
    const real g12_r = g12[i];
    const real g13_r = g13[i];
    const real g21_r = g21[i];
    const real g22_r = g22[i];
    const real g23_r = g23[i];
    const real g31_r = g31[i];
    const real g32_r = g32[i];
    const real g33_r = g33[i];
    const real delta_r = delta[i];

    sigG11 = g11_r * g11_r + g21_r * g21_r + g31_r * g31_r;
    sigG22 = g12_r * g12_r + g22_r * g22_r + g32_r * g32_r;
    sigG33 = g13_r * g13_r + g23_r * g23_r + g33_r * g33_r;
    sigG12 = g11_r * g12_r + g21_r * g22_r + g31_r * g32_r;
    sigG13 = g11_r * g13_r + g21_r * g23_r + g31_r * g33_r;
    sigG23 = g12_r * g13_r + g22_r * g23_r + g32_r * g33_r;

    if (fabs(sigG11) < eps) {
      sigG11 = zero;
    }
    if (fabs(sigG12) < eps) {
      sigG12 = zero;
    }
    if (fabs(sigG13) < eps) {
      sigG13 = zero;
    }
    if (fabs(sigG22) < eps) {
      sigG22 = zero;
    }
    if (fabs(sigG23) < eps) {
      sigG23 = zero;
    }
    if (fabs(sigG33) < eps) {
      sigG33 = zero;
    }

    if (fabs(sigG12 * sigG12 +
             sigG13 * sigG13 + sigG23 * sigG23) < eps) {
      sigma1 = sqrt(max(max(max(sigG11, sigG22), sigG33), zero));
      sigma3 = sqrt(max(min(min(sigG11, sigG22), sigG33), zero));
      Invariant1 = sigG11 + sigG22 + sigG33;
      sigma2 = sqrt(fabs(Invariant1 - sigma1 * sigma1 - sigma3 * sigma3));
    } else {
      Invariant1 = sigG11 + sigG22 + sigG33;
      Invariant2 = sigG11 * sigG22 + sigG11 * sigG33 + sigG22 * sigG33 -
        (sigG12 * sigG12 + sigG13 * sigG13 + sigG23 * sigG23);
      Invariant3 = sigG11 * sigG22 * sigG33 +
        two * sigG12 * sigG13 * sigG23 -
        (sigG11 * sigG23 * sigG23 + sigG22 * sigG13 * sigG13 +
         sigG33 * sigG12 * sigG12);

      Invariant1 = max(Invariant1, zero);
      Invariant2 = max(Invariant2, zero);
      Invariant3 = max(Invariant3, zero);

      alpha1 = Invariant1 * Invariant1 / (real) 9.0 - Invariant2 / (real) 3.0;

      alpha1 = max(alpha1, zero);

      alpha2 = Invariant1 * Invariant1 * Invariant1 / (real) 27.0 -
        Invariant1 * Invariant2 / (real) 6.0 + Invariant3 / two;

      tmp1 = alpha2 / sqrt(alpha1 * alpha1 * alpha1);

      if (tmp1 <= -one) {
        sigma1 = sqrt(max(Invariant1 / (real) 3.0 + sqrt(alpha1), zero));
        sigma2 = sigma1;
        sigma3 = sqrt(Invariant1 / (real) 3.0 - two * sqrt(alpha1));

      } else if (tmp1 >= one) {
        sigma1 = sqrt(max(Invariant1 / (real) 3.0 + two * sqrt(alpha1), zero));
        sigma2 = sqrt(Invariant1 / (real) 3.0 - sqrt(alpha1));
        sigma3 = sigma2;

      } else {
        alpha3 = acos(tmp1) / (real) 3.0;

        if (fabs(Invariant3) < eps) {
          sigma1 = sqrt(max(Invariant1 / (real) 3.0 +
                            two * sqrt(alpha1) * cos(alpha3), zero));
          sigma2 = sqrt(fabs(Invariant1 - sigma1 * sigma1));
          sigma3 = zero;
        } else {
          sigma1 = sqrt(max(Invariant1 / (real) 3.0 +
                            two * sqrt(alpha1) * cos(alpha3), zero));
          sigma2 = sqrt(Invariant1 / (real) 3.0 -
                        two * sqrt(alpha1) * cos(pi_3 + alpha3));
          sigma3 = sqrt(fabs(Invariant1 -
                             sigma1 * sigma1 - sigma2 * sigma2));
        }
      }
    }

    if (sigma1 > zero) {
      Dsigma = sigma3 * (sigma1 - sigma2) * (sigma2 - sigma3) /
        (sigma1 * sigma1);
    } else {
      Dsigma = zero;
    }

    Dsigma = max(Dsigma, zero);

    nut[i] = (c * delta_r) * (c * delta_r) * Dsigma * mult[i];
  }
}

#endif // __LES_SIGMA_NUT_KERNEL_CL__
