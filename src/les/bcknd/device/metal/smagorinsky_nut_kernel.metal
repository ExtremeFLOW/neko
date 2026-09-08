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
 * Metal compute kernel for the Smagorinsky eddy viscosity.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

kernel void smagorinsky_nut_compute_kernel(device const float *s11 [[ buffer(0) ]],
                                           device const float *s22 [[ buffer(1) ]],
                                           device const float *s33 [[ buffer(2) ]],
                                           device const float *s12 [[ buffer(3) ]],
                                           device const float *s13 [[ buffer(4) ]],
                                           device const float *s23 [[ buffer(5) ]],
                                           device const float *delta [[ buffer(6) ]],
                                           device float *nut [[ buffer(7) ]],
                                           device const float *mult [[ buffer(8) ]],
                                           constant float &c_s [[ buffer(9) ]],
                                           constant int &n [[ buffer(10) ]],
                                           uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint) n) return;

  const int i = (int) idx;

  const float s11_r = s11[i];
  const float s22_r = s22[i];
  const float s33_r = s33[i];
  const float s12_r = s12[i];
  const float s13_r = s13[i];
  const float s23_r = s23[i];
  const float delta_r = delta[i];

  const float s_abs = sqrt(2.0f * (s11_r * s11_r +
                                   s22_r * s22_r +
                                   s33_r * s33_r) +
                           4.0f * (s12_r * s12_r +
                                   s13_r * s13_r +
                                   s23_r * s23_r));

  nut[i] = c_s * c_s * delta_r * delta_r * s_abs * mult[i];
}
