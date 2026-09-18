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
 * Metal compute kernel for the Sigma eddy viscosity.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

kernel void sigma_nut_compute_kernel(device const float *g11 [[ buffer(0) ]],
                                     device const float *g12 [[ buffer(1) ]],
                                     device const float *g13 [[ buffer(2) ]],
                                     device const float *g21 [[ buffer(3) ]],
                                     device const float *g22 [[ buffer(4) ]],
                                     device const float *g23 [[ buffer(5) ]],
                                     device const float *g31 [[ buffer(6) ]],
                                     device const float *g32 [[ buffer(7) ]],
                                     device const float *g33 [[ buffer(8) ]],
                                     device const float *delta [[ buffer(9) ]],
                                     device float *nut [[ buffer(10) ]],
                                     device const float *mult [[ buffer(11) ]],
                                     constant float &c [[ buffer(12) ]],
                                     constant float &eps [[ buffer(13) ]],
                                     constant int &n [[ buffer(14) ]],
                                     uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint) n) return;

  const int i = (int) idx;

  const float pi_3 = 4.0f / 3.0f * atan(1.0f);

  float sigG11, sigG12, sigG13, sigG22, sigG23, sigG33;
  float sigma1, sigma2, sigma3;
  float Invariant1, Invariant2, Invariant3;
  float alpha1, alpha2, alpha3;
  float Dsigma;
  float tmp1;

  const float g11_r = g11[i];
  const float g12_r = g12[i];
  const float g13_r = g13[i];
  const float g21_r = g21[i];
  const float g22_r = g22[i];
  const float g23_r = g23[i];
  const float g31_r = g31[i];
  const float g32_r = g32[i];
  const float g33_r = g33[i];
  const float delta_r = delta[i];

  sigG11 = g11_r * g11_r + g21_r * g21_r + g31_r * g31_r;
  sigG22 = g12_r * g12_r + g22_r * g22_r + g32_r * g32_r;
  sigG33 = g13_r * g13_r + g23_r * g23_r + g33_r * g33_r;
  sigG12 = g11_r * g12_r + g21_r * g22_r + g31_r * g32_r;
  sigG13 = g11_r * g13_r + g21_r * g23_r + g31_r * g33_r;
  sigG23 = g12_r * g13_r + g22_r * g23_r + g32_r * g33_r;

  if (fabs(sigG11) < eps) {
    sigG11 = 0.0f;
  }
  if (fabs(sigG12) < eps) {
    sigG12 = 0.0f;
  }
  if (fabs(sigG13) < eps) {
    sigG13 = 0.0f;
  }
  if (fabs(sigG22) < eps) {
    sigG22 = 0.0f;
  }
  if (fabs(sigG23) < eps) {
    sigG23 = 0.0f;
  }
  if (fabs(sigG33) < eps) {
    sigG33 = 0.0f;
  }

  if (fabs(sigG12 * sigG12 +
           sigG13 * sigG13 + sigG23 * sigG23) < eps) {
    sigma1 = sqrt(max(max(max(sigG11, sigG22), sigG33), 0.0f));
    sigma3 = sqrt(max(min(min(sigG11, sigG22), sigG33), 0.0f));
    Invariant1 = sigG11 + sigG22 + sigG33;
    sigma2 = sqrt(fabs(Invariant1 - sigma1 * sigma1 - sigma3 * sigma3));
  } else {
    Invariant1 = sigG11 + sigG22 + sigG33;
    Invariant2 = sigG11 * sigG22 + sigG11 * sigG33 + sigG22 * sigG33 -
      (sigG12 * sigG12 + sigG13 * sigG13 + sigG23 * sigG23);
    Invariant3 = sigG11 * sigG22 * sigG33 +
      2.0f * sigG12 * sigG13 * sigG23 -
      (sigG11 * sigG23 * sigG23 + sigG22 * sigG13 * sigG13 +
       sigG33 * sigG12 * sigG12);

    Invariant1 = max(Invariant1, 0.0f);
    Invariant2 = max(Invariant2, 0.0f);
    Invariant3 = max(Invariant3, 0.0f);

    alpha1 = Invariant1 * Invariant1 / 9.0f - Invariant2 / 3.0f;

    alpha1 = max(alpha1, 0.0f);

    alpha2 = Invariant1 * Invariant1 * Invariant1 / 27.0f -
      Invariant1 * Invariant2 / 6.0f + Invariant3 / 2.0f;

    tmp1 = alpha2 / sqrt(alpha1 * alpha1 * alpha1);

    if (tmp1 <= -1.0f) {
      sigma1 = sqrt(max(Invariant1 / 3.0f + sqrt(alpha1), 0.0f));
      sigma2 = sigma1;
      sigma3 = sqrt(Invariant1 / 3.0f - 2.0f * sqrt(alpha1));

    } else if (tmp1 >= 1.0f) {
      sigma1 = sqrt(max(Invariant1 / 3.0f + 2.0f * sqrt(alpha1), 0.0f));
      sigma2 = sqrt(Invariant1 / 3.0f - sqrt(alpha1));
      sigma3 = sigma2;

    } else {
      alpha3 = acos(tmp1) / 3.0f;

      if (fabs(Invariant3) < eps) {
        sigma1 = sqrt(max(Invariant1 / 3.0f +
                          2.0f * sqrt(alpha1) * cos(alpha3), 0.0f));
        sigma2 = sqrt(fabs(Invariant1 - sigma1 * sigma1));
        sigma3 = 0.0f;
      } else {
        sigma1 = sqrt(max(Invariant1 / 3.0f +
                          2.0f * sqrt(alpha1) * cos(alpha3), 0.0f));
        sigma2 = sqrt(Invariant1 / 3.0f -
                      2.0f * sqrt(alpha1) * cos(pi_3 + alpha3));
        sigma3 = sqrt(fabs(Invariant1 -
                           sigma1 * sigma1 - sigma2 * sigma2));
      }
    }
  }

  if (sigma1 > 0.0f) {
    Dsigma = sigma3 * (sigma1 - sigma2) * (sigma2 - sigma3) /
      (sigma1 * sigma1);
  } else {
    Dsigma = 0.0f;
  }

  Dsigma = max(Dsigma, 0.0f);

  nut[i] = (c * delta_r) * (c * delta_r) * Dsigma * mult[i];
}
