#ifndef __WALL_MODELS_WALL_MODEL_KERNEL_CL__
#define __WALL_MODELS_WALL_MODEL_KERNEL_CL__
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
 * Device kernel for wall_model_compute_mag_field
 */
__kernel void wall_model_compute_mag_field_kernel(
    __global const real * __restrict__ tau_x_d,
    __global const real * __restrict__ tau_y_d,
    __global const real * __restrict__ tau_z_d,
    __global real * __restrict__ tau_field_d,
    __global const int * __restrict__ msk_d,
    const int m) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < m; i += str) {
    /* Magnitude of the shear stress vector */
    const real magtau = sqrt(tau_x_d[i] * tau_x_d[i] +
                             tau_y_d[i] * tau_y_d[i] +
                             tau_z_d[i] * tau_z_d[i]);

    /* Store the result in the tau_field array at the masked index */
    tau_field_d[msk_d[i + 1] - 1] = magtau;
  }
}

#endif // __WALL_MODELS_WALL_MODEL_KERNEL_CL__
