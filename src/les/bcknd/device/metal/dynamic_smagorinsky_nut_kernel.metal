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
 * Metal compute kernels for the dynamic Smagorinsky model.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

kernel void s_abs_compute_kernel(device float *s_abs [[ buffer(0) ]],
                                 device const float *s11 [[ buffer(1) ]],
                                 device const float *s22 [[ buffer(2) ]],
                                 device const float *s33 [[ buffer(3) ]],
                                 device const float *s12 [[ buffer(4) ]],
                                 device const float *s13 [[ buffer(5) ]],
                                 device const float *s23 [[ buffer(6) ]],
                                 constant int &n [[ buffer(7) ]],
                                 uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint) n) return;

  const int i = (int) idx;

  const float s11_r = s11[i];
  const float s22_r = s22[i];
  const float s33_r = s33[i];
  const float s12_r = s12[i];
  const float s13_r = s13[i];
  const float s23_r = s23[i];

  s_abs[i] = sqrt(2.0f * (s11_r * s11_r +
                          s22_r * s22_r +
                          s33_r * s33_r) +
                  4.0f * (s12_r * s12_r +
                          s13_r * s13_r +
                          s23_r * s23_r));
}

kernel void lij_compute_part1_kernel(device float *l11 [[ buffer(0) ]],
                                     device float *l22 [[ buffer(1) ]],
                                     device float *l33 [[ buffer(2) ]],
                                     device float *l12 [[ buffer(3) ]],
                                     device float *l13 [[ buffer(4) ]],
                                     device float *l23 [[ buffer(5) ]],
                                     device const float *u [[ buffer(6) ]],
                                     device const float *v [[ buffer(7) ]],
                                     device const float *w [[ buffer(8) ]],
                                     device const float *fu [[ buffer(9) ]],
                                     device const float *fv [[ buffer(10) ]],
                                     device const float *fw [[ buffer(11) ]],
                                     device float *fuu [[ buffer(12) ]],
                                     device float *fvv [[ buffer(13) ]],
                                     device float *fww [[ buffer(14) ]],
                                     device float *fuv [[ buffer(15) ]],
                                     device float *fuw [[ buffer(16) ]],
                                     device float *fvw [[ buffer(17) ]],
                                     constant int &n [[ buffer(18) ]],
                                     uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint) n) return;

  const int i = (int) idx;

  const float u_r = u[i];
  const float v_r = v[i];
  const float w_r = w[i];
  const float fu_r = fu[i];
  const float fv_r = fv[i];
  const float fw_r = fw[i];

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

kernel void lij_compute_part2_kernel(device float *l11 [[ buffer(0) ]],
                                     device float *l22 [[ buffer(1) ]],
                                     device float *l33 [[ buffer(2) ]],
                                     device float *l12 [[ buffer(3) ]],
                                     device float *l13 [[ buffer(4) ]],
                                     device float *l23 [[ buffer(5) ]],
                                     device const float *fuu [[ buffer(6) ]],
                                     device const float *fvv [[ buffer(7) ]],
                                     device const float *fww [[ buffer(8) ]],
                                     device const float *fuv [[ buffer(9) ]],
                                     device const float *fuw [[ buffer(10) ]],
                                     device const float *fvw [[ buffer(11) ]],
                                     constant int &n [[ buffer(12) ]],
                                     uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint) n) return;

  const int i = (int) idx;

  l11[i] -= fuu[i];
  l22[i] -= fvv[i];
  l33[i] -= fww[i];
  l12[i] -= fuv[i];
  l13[i] -= fuw[i];
  l23[i] -= fvw[i];
}

kernel void mij_compute_part1_kernel(device float *m11 [[ buffer(0) ]],
                                     device float *m22 [[ buffer(1) ]],
                                     device float *m33 [[ buffer(2) ]],
                                     device float *m12 [[ buffer(3) ]],
                                     device float *m13 [[ buffer(4) ]],
                                     device float *m23 [[ buffer(5) ]],
                                     device const float *s_abs [[ buffer(6) ]],
                                     device const float *s11 [[ buffer(7) ]],
                                     device const float *s22 [[ buffer(8) ]],
                                     device const float *s33 [[ buffer(9) ]],
                                     device const float *s12 [[ buffer(10) ]],
                                     device const float *s13 [[ buffer(11) ]],
                                     device const float *s23 [[ buffer(12) ]],
                                     device const float *fs_abs [[ buffer(13) ]],
                                     device const float *fs11 [[ buffer(14) ]],
                                     device const float *fs22 [[ buffer(15) ]],
                                     device const float *fs33 [[ buffer(16) ]],
                                     device const float *fs12 [[ buffer(17) ]],
                                     device const float *fs13 [[ buffer(18) ]],
                                     device const float *fs23 [[ buffer(19) ]],
                                     device float *fsabss11 [[ buffer(20) ]],
                                     device float *fsabss22 [[ buffer(21) ]],
                                     device float *fsabss33 [[ buffer(22) ]],
                                     device float *fsabss12 [[ buffer(23) ]],
                                     device float *fsabss13 [[ buffer(24) ]],
                                     device float *fsabss23 [[ buffer(25) ]],
                                     constant float &delta_ratio2 [[ buffer(26) ]],
                                     constant int &n [[ buffer(27) ]],
                                     uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint) n) return;

  const int i = (int) idx;

  const float s_abs_r = s_abs[i];
  const float fs_abs_r = fs_abs[i];

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

kernel void mij_nut_compute_part2_kernel(device const float *m11 [[ buffer(0) ]],
                                         device const float *m22 [[ buffer(1) ]],
                                         device const float *m33 [[ buffer(2) ]],
                                         device const float *m12 [[ buffer(3) ]],
                                         device const float *m13 [[ buffer(4) ]],
                                         device const float *m23 [[ buffer(5) ]],
                                         device const float *l11 [[ buffer(6) ]],
                                         device const float *l22 [[ buffer(7) ]],
                                         device const float *l33 [[ buffer(8) ]],
                                         device const float *l12 [[ buffer(9) ]],
                                         device const float *l13 [[ buffer(10) ]],
                                         device const float *l23 [[ buffer(11) ]],
                                         device const float *fsabss11 [[ buffer(12) ]],
                                         device const float *fsabss22 [[ buffer(13) ]],
                                         device const float *fsabss33 [[ buffer(14) ]],
                                         device const float *fsabss12 [[ buffer(15) ]],
                                         device const float *fsabss13 [[ buffer(16) ]],
                                         device const float *fsabss23 [[ buffer(17) ]],
                                         device float *num [[ buffer(18) ]],
                                         device float *den [[ buffer(19) ]],
                                         device float *c_dyn [[ buffer(20) ]],
                                         device const float *delta [[ buffer(21) ]],
                                         device const float *s_abs [[ buffer(22) ]],
                                         device float *nut [[ buffer(23) ]],
                                         constant float &alpha [[ buffer(24) ]],
                                         constant int &n [[ buffer(25) ]],
                                         uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint) n) return;

  const int i = (int) idx;

  const float delta_r = delta[i];
  const float delta2 = delta_r * delta_r;
  float c_dyn_r;

  const float m11_r = delta2 * (m11[i] - fsabss11[i]);
  const float m22_r = delta2 * (m22[i] - fsabss22[i]);
  const float m33_r = delta2 * (m33[i] - fsabss33[i]);
  const float m12_r = delta2 * (m12[i] - fsabss12[i]);
  const float m13_r = delta2 * (m13[i] - fsabss13[i]);
  const float m23_r = delta2 * (m23[i] - fsabss23[i]);

  const float num_curr = m11_r * l11[i]
    + m22_r * l22[i]
    + m33_r * l33[i]
    + 2.0f * (m12_r * l12[i]
              + m13_r * l13[i]
              + m23_r * l23[i]);

  const float den_curr = m11_r * m11_r
    + m22_r * m22_r
    + m33_r * m33_r
    + 2.0f * (m12_r * m12_r
              + m13_r * m13_r
              + m23_r * m23_r);

  const float num_r = alpha * num[i] + (1.0f - alpha) * num_curr;
  const float den_r = alpha * den[i] + (1.0f - alpha) * den_curr;

  num[i] = num_r;
  den[i] = den_r;

  if (den_r > 0.0f) {
    c_dyn_r = 0.5f * num_r / den_r;
  } else {
    c_dyn_r = 0.0f;
  }

  c_dyn_r = max(c_dyn_r, 0.0f);
  c_dyn[i] = c_dyn_r;
  nut[i] = c_dyn_r * delta2 * s_abs[i];
}
