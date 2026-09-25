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
 * Metal compute kernel for the WALE eddy viscosity.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

kernel void wale_nut_compute_kernel(device const float *g11 [[ buffer(0) ]],
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

  const float s11 = g11_r;
  const float s22 = g22_r;
  const float s33 = g33_r;
  const float s12 = 0.5f * (g12_r + g21_r);
  const float s13 = 0.5f * (g13_r + g31_r);
  const float s23 = 0.5f * (g23_r + g32_r);

  const float gsqr_11 = g11_r * g11_r + g12_r * g21_r + g13_r * g31_r;
  const float gsqr_12 = g11_r * g12_r + g12_r * g22_r + g13_r * g32_r;
  const float gsqr_13 = g11_r * g13_r + g12_r * g23_r + g13_r * g33_r;
  const float gsqr_21 = g21_r * g11_r + g22_r * g21_r + g23_r * g31_r;
  const float gsqr_22 = g21_r * g12_r + g22_r * g22_r + g23_r * g32_r;
  const float gsqr_23 = g21_r * g13_r + g22_r * g23_r + g23_r * g33_r;
  const float gsqr_31 = g31_r * g11_r + g32_r * g21_r + g33_r * g31_r;
  const float gsqr_32 = g31_r * g12_r + g32_r * g22_r + g33_r * g32_r;
  const float gsqr_33 = g31_r * g13_r + g32_r * g23_r + g33_r * g33_r;

  const float trace = (gsqr_11 + gsqr_22 + gsqr_33) / 3.0f;

  const float sd11 = gsqr_11 - trace;
  const float sd22 = gsqr_22 - trace;
  const float sd33 = gsqr_33 - trace;
  const float sd12 = 0.5f * (gsqr_12 + gsqr_21);
  const float sd13 = 0.5f * (gsqr_13 + gsqr_31);
  const float sd23 = 0.5f * (gsqr_23 + gsqr_32);

  const float Sdij_Sdij = sd11 * sd11 + sd22 * sd22 + sd33 * sd33 +
    2.0f * (sd12 * sd12 + sd13 * sd13 + sd23 * sd23);
  const float Sij_Sij = s11 * s11 + s22 * s22 + s33 * s33 +
    2.0f * (s12 * s12 + s13 * s13 + s23 * s23);

  const float OP_wale = sqrt(Sdij_Sdij * Sdij_Sdij * Sdij_Sdij) /
    max((sqrt(Sij_Sij * Sij_Sij * Sij_Sij * Sij_Sij * Sij_Sij)
         + sqrt(sqrt(Sdij_Sdij * Sdij_Sdij * Sdij_Sdij
                     * Sdij_Sdij * Sdij_Sdij))), eps);

  nut[i] = c * c * delta_r * delta_r * OP_wale * mult[i];
}
