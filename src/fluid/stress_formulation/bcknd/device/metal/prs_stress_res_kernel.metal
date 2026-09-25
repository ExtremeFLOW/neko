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
 * Metal compute kernels for the stress formulation pressure residual.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/**
 * Device kernel for the stress formulation pressure residual, part 1
 */
kernel void prs_stress_res_part1_kernel(device float *ta1[[ buffer(0) ]],
                                        device float *ta2[[ buffer(1) ]],
                                        device float *ta3[[ buffer(2) ]],
                                        device float *wa1[[ buffer(3) ]],
                                        device float *wa2[[ buffer(4) ]],
                                        device float *wa3[[ buffer(5) ]],
                                        device const float *s11[[ buffer(6) ]],
                                        device const float *s22[[ buffer(7) ]],
                                        device const float *s33[[ buffer(8) ]],
                                        device const float *s12[[ buffer(9) ]],
                                        device const float *s13[[ buffer(10) ]],
                                        device const float *s23[[ buffer(11) ]],
                                        device const float *f_u[[ buffer(12) ]],
                                        device const float *f_v[[ buffer(13) ]],
                                        device const float *f_w[[ buffer(14) ]],
                                        device const float *B[[ buffer(15) ]],
                                        device const float *rho[[ buffer(16) ]],
                                        constant int &n[[ buffer(17) ]],
                                        uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  wa1[idx] -= 2.0f * (ta1[idx] * s11[idx]
                      + ta2[idx] * s12[idx]
                      + ta3[idx] * s13[idx]);
  wa2[idx] -= 2.0f * (ta1[idx] * s12[idx]
                      + ta2[idx] * s22[idx]
                      + ta3[idx] * s23[idx]);
  wa3[idx] -= 2.0f * (ta1[idx] * s13[idx]
                      + ta2[idx] * s23[idx]
                      + ta3[idx] * s33[idx]);

  ta1[idx] = (f_u[idx] / rho[idx]) - ((wa1[idx] / rho[idx]) * B[idx]);
  ta2[idx] = (f_v[idx] / rho[idx]) - ((wa2[idx] / rho[idx]) * B[idx]);
  ta3[idx] = (f_w[idx] / rho[idx]) - ((wa3[idx] / rho[idx]) * B[idx]);
}

/**
 * Device kernel for the stress formulation pressure residual, part 3
 */
kernel void prs_stress_res_part3_kernel(device float *p_res[[ buffer(0) ]],
                                        device const float *ta1[[ buffer(1) ]],
                                        device const float *ta2[[ buffer(2) ]],
                                        device const float *ta3[[ buffer(3) ]],
                                        device const float *wa1[[ buffer(4) ]],
                                        device const float *wa2[[ buffer(5) ]],
                                        device const float *wa3[[ buffer(6) ]],
                                        constant float &dtbd[[ buffer(7) ]],
                                        constant int &n[[ buffer(8) ]],
                                        uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  p_res[idx] = p_res[idx] - (dtbd * (ta1[idx] + ta2[idx] + ta3[idx]))
    - (wa1[idx] + wa2[idx] + wa3[idx]);
}
