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

/*
 * Richardson-number correction of Lilly (1962), Ri = N^2 / |S|^2, for the
 * Vreman model. Gradients are gather-scatter sums, scaled by mult.
 * gdT is the interface-averaged g . grad(T), g pointing upwards
 * (reversed from the case-file input at init).
 */
kernel void vreman_stability_kernel(
                                      device const float *a11 [[ buffer(0) ]],
                                      device const float *a12 [[ buffer(1) ]],
                                      device const float *a13 [[ buffer(2) ]],
                                      device const float *a21 [[ buffer(3) ]],
                                      device const float *a22 [[ buffer(4) ]],
                                      device const float *a23 [[ buffer(5) ]],
                                      device const float *a31 [[ buffer(6) ]],
                                      device const float *a32 [[ buffer(7) ]],
                                      device const float *a33 [[ buffer(8) ]],
                                      device const float *gdT [[ buffer(9) ]],
                                      device float *nut [[ buffer(10) ]],
                                      device const float *mult [[ buffer(11) ]],
                                      constant float &ri_c [[ buffer(12) ]],
                                      constant float &ref_temp [[ buffer(13) ]],
                                      constant float &eps [[ buffer(14) ]],
                                      constant int &n [[ buffer(15) ]],
                                      uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint) n) return;

  const int i = (int) idx;
  const float m = mult[i];
  const float N2 = gdT[i] / ref_temp;

  // |S|^2 = 2 S_ij S_ij of the averaged gradient (sums, hence m*m)
  const float s12 = a12[i] + a21[i];
  const float s13 = a13[i] + a31[i];
  const float s23 = a23[i] + a32[i];
  const float s_sq = (2.0f * (a11[i] * a11[i] + a22[i] * a22[i]
                              + a33[i] * a33[i])
                      + s12 * s12 + s13 * s13 + s23 * s23) * m * m;

  const float ri = N2 / (s_sq + eps);
  if (ri <= ri_c) {
    float v = 1.0f - ri / ri_c;
    if (v < 0.0f) v = 0.0f;
    nut[i] = nut[i] * sqrt(v);
  } else {
    nut[i] = eps;
  }
}
