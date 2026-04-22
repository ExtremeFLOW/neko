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

#ifndef ROUGH_LOG_LAW_KERNEL_H
#define ROUGH_LOG_LAW_KERNEL_H

#include <cmath>
#include <algorithm>

/**
 * CUDA kernel for the rough log-law wall model.
 */
template<typename T>
__global__ void rough_log_law_compute(const T* __restrict__ u_d,
                                      const T* __restrict__ v_d,
                                      const T* __restrict__ w_d,
                                      const int* __restrict__ ind_r_d,
                                      const int* __restrict__ ind_s_d,
                                      const int* __restrict__ ind_t_d,
                                      const int* __restrict__ ind_e_d,
                                      const T* __restrict__ n_x_d,
                                      const T* __restrict__ n_y_d,
                                      const T* __restrict__ n_z_d,
                                      const T* __restrict__ h_d,
                                      T* __restrict__ tau_x_d,
                                      T* __restrict__ tau_y_d,
                                      T* __restrict__ tau_z_d,
                                      const int n_nodes,
                                      const int lx,
                                      const T kappa,
                                      const T rho,
                                      const T B,
                                      const T z0) {

    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int str = blockDim.x * gridDim.x;
    for (int i = idx; i < n_nodes; i += str) {
        // Sample the velocity
        const int index = (ind_e_d[i] - 1) * lx * lx * lx +
                          (ind_t_d[i] - 1) * lx * lx +
                          (ind_s_d[i] - 1) * lx +
                          (ind_r_d[i] - 1);

        T ui = u_d[index];
        T vi = v_d[index];
        T wi = w_d[index];

        // Load normal vectors and wall shear stress values once
        T nx = n_x_d[i];
        T ny = n_y_d[i];
        T nz = n_z_d[i];
        T h = h_d[i];

        // Project on tangential direction
        T normu = ui * nx + vi * ny + wi * nz;

        ui -= normu * nx;
        vi -= normu * ny;
        wi -= normu * nz;

        T magu = sqrt(ui * ui + vi * vi + wi * wi);

        // Compute the wall shear stress using the rough log-law
        T utau = 0.0;
        if (h > z0) {
            utau = (magu - B) * kappa / log(h / z0);
        }

        // Distribute according to the velocity vector
        tau_x_d[i] = -rho *utau * utau * ui / magu;
        tau_y_d[i] = -rho *utau * utau * vi / magu;
        tau_z_d[i] = -rho *utau * utau * wi / magu;
    }
}

#endif // ROUGH_LOG_LAW_KERNEL_H