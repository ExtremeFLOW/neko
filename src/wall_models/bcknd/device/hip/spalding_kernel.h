#ifndef __COMMON_SPALDING_KERNEL_H__
#define __COMMON_SPALDING_KERNEL_H__
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
 * Device kernel for spalding_compute
 */
#include <cmath>
#include <algorithm>
template<typename T>
__device__ T solve(const T u, const T y, const T guess, const T nu,
                   const T kappa, const T B);

/**
 * CUDA kernel for Spalding's wall model.
 */
template<typename T>
__global__ void spalding_compute(const T * __restrict__ u_d,
                                 const T * __restrict__ v_d,
                                 const T * __restrict__ w_d,
                                 const int * __restrict__ ind_r_d,
                                 const int * __restrict__ ind_s_d,
                                 const int * __restrict__ ind_t_d,
                                 const int * __restrict__ ind_e_d,
                                 const T * __restrict__ n_x_d,
                                 const T * __restrict__ n_y_d,
                                 const T * __restrict__ n_z_d,
                                 const T * __restrict__ nu_d,
                                 const T * __restrict__ h_d,
                                 T * __restrict__ tau_x_d,
                                 T * __restrict__ tau_y_d,
                                 T * __restrict__ tau_z_d,
                                 const int n_nodes,
                                 const int lx,
                                 const T kappa,
                                 const T B,
                                 const int tstep) {
                                    
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

        // Get initial guess for Newton solver
        T guess;
        if (tstep == 1) {
            guess = sqrt(magu * nu_d[i] / h);
        } else {
            guess = tau_x_d[i] * tau_x_d[i] +
                    tau_y_d[i] * tau_y_d[i] +
                    tau_z_d[i] * tau_z_d[i];
            guess = sqrt(sqrt(guess));
        }

        // Solve for utau using Newton's method
        T utau = solve(magu, h, guess, nu_d[i], kappa, B);

        // Distribute according to the velocity vector
        tau_x_d[i] = -utau * utau * ui / magu;
        tau_y_d[i] = -utau * utau * vi / magu;
        tau_z_d[i] = -utau * utau * wi / magu;
    }
}

/**
 * Newton solver for the algebraic equation defined by the law on GPU.
 */
template<typename T>
__device__ T solve(const T u, const T y, const T guess, const T nu,
                   const T kappa, const T B) {
    T utau = guess;
    T yp, up, f, df, old, error;
    const int maxiter = 100;

    for (int k = 0; k < maxiter; ++k) {
        up = u / utau;
        yp = y * utau / nu;
        old = utau;

        // Evaluate function and its derivative
        f = (up + exp(-kappa * B) *
                  (exp(kappa * up) - 1.0 - kappa * up -
                   0.5 * (kappa * up) * (kappa * up) -
                   (1.0 / 6.0) * (kappa * up) * (kappa * up) * (kappa * up)) -
             yp);

        df = (-y / nu - u / (utau * utau) -
              kappa * up / utau * exp(-kappa * B) *
                  (exp(kappa * up) - 1.0 - kappa * up -
                   0.5 * (kappa * up) * (kappa * up)));

        // Update solution
        utau -= f / df;

        error = fabs((old - utau) / old);

        if (error < 1e-3) {
            break;
        }
    }

    return utau;
}

#endif // __COMMON_SPALDING_KERNEL_H__