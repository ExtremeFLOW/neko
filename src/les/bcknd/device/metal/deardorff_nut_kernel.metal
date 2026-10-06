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
 * Metal compute kernel for the Deardorff eddy viscosity.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

kernel void deardorff_nut_compute_kernel(device float *TKE [[ buffer(0) ]],
                                         device const float *dTdx [[ buffer(1) ]],
                                         device const float *dTdy [[ buffer(2) ]],
                                         device const float *dTdz [[ buffer(3) ]],
                                         device const float *a11 [[ buffer(4) ]],
                                         device const float *a12 [[ buffer(5) ]],
                                         device const float *a13 [[ buffer(6) ]],
                                         device const float *a21 [[ buffer(7) ]],
                                         device const float *a22 [[ buffer(8) ]],
                                         device const float *a23 [[ buffer(9) ]],
                                         device const float *a31 [[ buffer(10) ]],
                                         device const float *a32 [[ buffer(11) ]],
                                         device const float *a33 [[ buffer(12) ]],
                                         device const float *delta [[ buffer(13) ]],
                                         device float *nut [[ buffer(14) ]],
                                         device float *temperature_alphat [[ buffer(15) ]],
                                         device float *TKE_alphat [[ buffer(16) ]],
                                         device float *TKE_source [[ buffer(17) ]],
                                         constant float &c_k [[ buffer(18) ]],
                                         constant float &T0 [[ buffer(19) ]],
                                         constant float &g1 [[ buffer(20) ]],
                                         constant float &g2 [[ buffer(21) ]],
                                         constant float &g3 [[ buffer(22) ]],
                                         constant float &eps [[ buffer(23) ]],
                                         constant int &n [[ buffer(24) ]],
                                         uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint) n) return;

  const int i = (int) idx;

  const float dTdx_r = dTdx[i];
  const float dTdy_r = dTdy[i];
  const float dTdz_r = dTdz[i];
  const float a11_r = a11[i];
  const float a12_r = a12[i];
  const float a13_r = a13[i];
  const float a21_r = a21[i];
  const float a22_r = a22[i];
  const float a23_r = a23[i];
  const float a31_r = a31[i];
  const float a32_r = a32[i];
  const float a33_r = a33[i];
  const float delta_r = delta[i];

  if (TKE[i] < eps) {
    TKE[i] = eps;
  }
  const float TKE_r = TKE[i];

  const float N2 = (dTdx_r * g1 + dTdy_r * g2 + dTdz_r * g3) / T0;
  float l;
  if (N2 > 0.0f) {
    l = 0.76f * sqrt(TKE_r / N2);
    l = min(l, delta_r);
  } else {
    l = delta_r;
  }

  const float nut_r = c_k * l * sqrt(TKE_r);
  const float temperature_alphat_r = (1.0f + 2.0f * l / delta_r) * nut_r;

  nut[i] = nut_r;
  temperature_alphat[i] = temperature_alphat_r;
  TKE_alphat[i] = 2.0f * nut_r;

  const float s11 = a11_r + a11_r;
  const float s22 = a22_r + a22_r;
  const float s33 = a33_r + a33_r;
  const float s12 = a12_r + a21_r;
  const float s13 = a13_r + a31_r;
  const float s23 = a23_r + a32_r;

  const float shear = nut_r * (s11 * a11_r + s12 * a12_r + s13 * a13_r
                               + s12 * a21_r + s22 * a22_r + s23 * a23_r
                               + s13 * a31_r + s23 * a32_r + s33 * a33_r);
  const float buoyancy = - (dTdx_r * g1 +
                            dTdy_r * g2 +
                            dTdz_r * g3) * temperature_alphat_r / T0;
  const float dissipation = - (0.19f + 0.74f * l / delta_r)
    * sqrt(TKE_r * TKE_r * TKE_r) / l;

  TKE_source[i] = shear + buoyancy + dissipation;
}
