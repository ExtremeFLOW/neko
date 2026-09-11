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
 * Metal compute kernels for the Vreman eddy viscosity.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

kernel void vreman_nut_compute_kernel(device const float *a11 [[ buffer(0) ]],
                                      device const float *a12 [[ buffer(1) ]],
                                      device const float *a13 [[ buffer(2) ]],
                                      device const float *a21 [[ buffer(3) ]],
                                      device const float *a22 [[ buffer(4) ]],
                                      device const float *a23 [[ buffer(5) ]],
                                      device const float *a31 [[ buffer(6) ]],
                                      device const float *a32 [[ buffer(7) ]],
                                      device const float *a33 [[ buffer(8) ]],
                                      device const float *delta [[ buffer(9) ]],
                                      device float *nut [[ buffer(10) ]],
                                      device const float *mult [[ buffer(11) ]],
                                      constant float &c [[ buffer(12) ]],
                                      constant float &eps [[ buffer(13) ]],
                                      constant int &n [[ buffer(14) ]],
                                      uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint) n) return;

  const int i = (int) idx;

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

  const float beta11 = a11_r * a11_r + a21_r * a21_r + a31_r * a31_r;
  const float beta22 = a12_r * a12_r + a22_r * a22_r + a32_r * a32_r;
  const float beta33 = a13_r * a13_r + a23_r * a23_r + a33_r * a33_r;
  const float beta12 = a11_r * a12_r + a21_r * a22_r + a31_r * a32_r;
  const float beta13 = a11_r * a13_r + a21_r * a23_r + a31_r * a33_r;
  const float beta23 = a12_r * a13_r + a22_r * a23_r + a32_r * a33_r;

  float b_beta = beta11 * beta22 - beta12 * beta12
    + beta11 * beta33 - beta13 * beta13
    + beta22 * beta33 - beta23 * beta23;

  b_beta = max(0.0f, b_beta);

  const float aijaij = beta11 + beta22 + beta33;

  nut[i] = c * delta_r * delta_r
    * sqrt(b_beta / (aijaij + eps)) * mult[i];
}

kernel void vreman_nut_compute_buoy_kernel(device const float *a11 [[ buffer(0) ]],
                                           device const float *a12 [[ buffer(1) ]],
                                           device const float *a13 [[ buffer(2) ]],
                                           device const float *a21 [[ buffer(3) ]],
                                           device const float *a22 [[ buffer(4) ]],
                                           device const float *a23 [[ buffer(5) ]],
                                           device const float *a31 [[ buffer(6) ]],
                                           device const float *a32 [[ buffer(7) ]],
                                           device const float *a33 [[ buffer(8) ]],
                                           device const float *delta [[ buffer(9) ]],
                                           device float *nut [[ buffer(10) ]],
                                           device const float *mult [[ buffer(11) ]],
                                           constant float &c [[ buffer(12) ]],
                                           constant float &eps [[ buffer(13) ]],
                                           constant int &n [[ buffer(14) ]],
                                           device const float *dTdx [[ buffer(15) ]],
                                           device const float *dTdy [[ buffer(16) ]],
                                           device const float *dTdz [[ buffer(17) ]],
                                           constant float &n1 [[ buffer(18) ]],
                                           constant float &n2 [[ buffer(19) ]],
                                           constant float &n3 [[ buffer(20) ]],
                                           constant float &g1 [[ buffer(21) ]],
                                           constant float &g2 [[ buffer(22) ]],
                                           constant float &g3 [[ buffer(23) ]],
                                           constant float &ri_c [[ buffer(24) ]],
                                           constant float &ref_temp [[ buffer(25) ]],
                                           uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint) n) return;

  const int i = (int) idx;

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

  const float beta11 = a11_r * a11_r + a21_r * a21_r + a31_r * a31_r;
  const float beta22 = a12_r * a12_r + a22_r * a22_r + a32_r * a32_r;
  const float beta33 = a13_r * a13_r + a23_r * a23_r + a33_r * a33_r;
  const float beta12 = a11_r * a12_r + a21_r * a22_r + a31_r * a32_r;
  const float beta13 = a11_r * a13_r + a21_r * a23_r + a31_r * a33_r;
  const float beta23 = a12_r * a13_r + a22_r * a23_r + a32_r * a33_r;

  float b_beta = beta11 * beta22 - beta12 * beta12
    + beta11 * beta33 - beta13 * beta13
    + beta22 * beta33 - beta23 * beta23;

  b_beta = max(0.0f, b_beta);

  const float aijaij = beta11 + beta22 + beta33;

  const float nut0 = c * delta_r * delta_r
    * sqrt(b_beta / (aijaij + eps)) * mult[i];

  /* Scalar gradient for buoyancy */
  const float buoyancy = (g1 * dTdx[i] + g2 * dTdy[i] + g3 * dTdz[i])
    / ref_temp;

  /* Directional derivative along n */
  const float du_n1 = a11_r * n1 + a12_r * n2 + a13_r * n3;
  const float du_n2 = a21_r * n1 + a22_r * n2 + a23_r * n3;
  const float du_n3 = a31_r * n1 + a32_r * n2 + a33_r * n3;

  const float du_parallel = du_n1 * n1 + du_n2 * n2 + du_n3 * n3;

  const float sh1 = du_n1 - du_parallel * n1;
  const float sh2 = du_n2 - du_parallel * n2;
  const float sh3 = du_n3 - du_parallel * n3;

  const float shear_sq = sh1 * sh1 + sh2 * sh2 + sh3 * sh3;

  const float ri = buoyancy / (shear_sq + eps);

  float out;
  if (ri <= ri_c) {
    float v = 1.0f - ri / ri_c;
    if (v < 0.0f) v = 0.0f;
    out = nut0 * sqrt(v);
  } else {
    out = eps;
  }

  nut[i] = out;
}
