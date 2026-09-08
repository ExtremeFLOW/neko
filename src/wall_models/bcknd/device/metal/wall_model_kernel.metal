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
 * Metal compute kernel for the wall shear stress magnitude field.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

kernel void wall_model_compute_mag_field_kernel(device const float *tau_x_d [[ buffer(0) ]],
                                                device const float *tau_y_d [[ buffer(1) ]],
                                                device const float *tau_z_d [[ buffer(2) ]],
                                                device float *tau_field_d [[ buffer(3) ]],
                                                device const int *msk_d [[ buffer(4) ]],
                                                constant int &m [[ buffer(5) ]],
                                                uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint) m) return;

  const int i = (int) idx;

  /* Magnitude of the shear stress vector */
  const float magtau = sqrt(tau_x_d[i] * tau_x_d[i] +
                            tau_y_d[i] * tau_y_d[i] +
                            tau_z_d[i] * tau_z_d[i]);

  /* Store the result in the tau_field array at the masked index */
  tau_field_d[msk_d[i + 1] - 1] = magtau;
}
