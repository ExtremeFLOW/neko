#ifndef __COMMON_WALL_MODEL_KERNEL_H__
#define __COMMON_WALL_MODEL_KERNEL_H__
/*
 Copyright (c) 2025, The Neko Authors
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
 * Device kernel for wall_model_tau_field_compute
 */
#include <cmath>
#include <algorithm>
template<typename T>
__global__ void wall_model_compute_mag_field(const T * __restrict__ tau_x_d,
                                             const T * __restrict__ tau_y_d,
                                             const T * __restrict__ tau_z_d,
                                             T * __restrict__ tau_field_d,
                                             const int * __restrict__ msk_d,
                                             const int m) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int str = blockDim.x * gridDim.x;

    for (int i = idx; i < m; i += str) {
        // Compute the magnitude of the shear stress vector
        T magtau = sqrt(tau_x_d[i] * tau_x_d[i] +
        tau_y_d[i] * tau_y_d[i] +
        tau_z_d[i] * tau_z_d[i]);

        // Store the result in the tau_field array at the masked index
        tau_field_d[msk_d[i]-1] = magtau;
    }
} 
#endif // __COMMON_WALL_model_KERNEL_H__