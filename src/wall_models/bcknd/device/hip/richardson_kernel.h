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

#ifndef RICHARDSON_KERNEL_H
#define RICHARDSON_KERNEL_H

#include <cmath>
#include <algorithm>

/*
* Similarity laws and corrections for the STABLE regime:
* Based on Mauritsen et al. 2007
*/

template<typename T>
__device__ T f_tau_stable(T Ri_b)
{
    return 0.17 * (0.25 + 0.75 / (1.0 + 4.0 * Ri_b));
}

template<typename T>
__device__ T f_theta_stable(T Ri_b)
{
    return -0.145 / (1.0 + 4.0 * Ri_b);
}

template<typename T>
__device__ T tau_stable(T magu, T Ri_b, T h, T z0, T kappa)
{
    T log_hz0 = log(h / z0);
    return (magu * magu) / (log_hz0 * log_hz0) * (f_tau_stable<T>(Ri_b) / f_tau_stable<T>(0.0)) * (kappa * kappa);
}

template<typename T>
__device__ T heat_flux_stable(T ti, T ts, T Ri_b, T h, T z0h, T utau, T kappa, T Pr)
{
    return (ti - ts) / log(h / z0h) * (f_theta_stable<T>(Ri_b) / f_theta_stable<T>(0.0)) * kappa * (utau / Pr);
}

/*
* Similarity laws and corrections for the UNSTABLE (convective) regime:
* Based on Louis 1979
*/

template<typename T>
__device__ T f_tau_convective(T Ri_b, T c)
{
    return 1.0 - (2.0 * Ri_b) / (1.0 + c * sqrt(fabs(Ri_b)));
}

template<typename T>
__device__ T f_theta_convective(T Ri_b, T c)
{
    // Functionally identical to f_tau_convective in Louis 1979
    return 1.0 - (2.0 * Ri_b) / (1.0 + c * sqrt(fabs(Ri_b)));
}

template<typename T>
__device__ T tau_convective(T magu, T Ri_b, T h, T z0, T kappa)
{
    T a = kappa / log(h / z0);
    T b = 2.0;
    T c = 7.4 * (a * a) * b * sqrt(h / z0);
    
    return (a * a) * (magu * magu) * f_tau_convective<T>(Ri_b, c);
}

template<typename T>
__device__ T heat_flux_convective(T ti, T ts, T Ri_b, T h, T magu, T z0h, T kappa)
{
    T a = kappa / log(h / z0h);
    T b = 2.0;
    T c = 5.3 * (a * a) * b * sqrt(h / z0h);
    
    return -(a * a) / 0.74 * magu * (ti - ts) * f_theta_convective<T>(Ri_b, c);
}

/*
* Similarity laws and corrections for the NEUTRAL regime:
*/

template<typename T>
__device__ T tau_neutral(T magu, T h, T z0, T kappa)
{
    T val = (kappa * magu) / log(h / z0);
    return val * val;
}

template<typename T>
__device__ T heat_flux_neutral(T ti, T ts, T h, T z0h, T utau, T kappa)
{
    return kappa * utau * (ti - ts) / log(h / z0h);
}

/*
 * CUDA kernel for the Richardson wall model.   
 */
template<typename T, int BC_TYPE>  
__global__ void richardson_compute(
    const T* __restrict__ u_d,     
    const T* __restrict__ v_d,
    const T* __restrict__ w_d,
    const T* __restrict__ temp_d,
    const T* __restrict__ h_d,     
    const T* __restrict__ n_x_d,  
    const T* __restrict__ n_y_d,
    const T* __restrict__ n_z_d,
    const int* __restrict__ ind_r_d,
    const int* __restrict__ ind_s_d,
    const int* __restrict__ ind_t_d,
    const int* __restrict__ ind_e_d,
    T* __restrict__ tau_x_d,       
    T* __restrict__ tau_y_d,
    T* __restrict__ tau_z_d,
    int n_nodes,
    int lx,                     
    T kappa,
    T mu,
    T rho, 
    T g1,
    T g2,
    T g3,
    T Pr,
    T z0,
    T z0h_in,
    T bc_value
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int str = blockDim.x * gridDim.x;
    if(idx >= n_nodes) return;

    const T Ri_threshold = 1e-4;

    // Use 64-bit integer to prevent overflow for large polynomial order / element counts        
    for (int i = idx; i < n_nodes; i += str) {
        const int index = (ind_e_d[i] - 1) * lx * lx * lx +
                          (ind_t_d[i] - 1) * lx * lx +
                          (ind_s_d[i] - 1) * lx +
                          (ind_r_d[i] - 1);  // 'long long' might be needed to avoid
                                             // overflowing 32-bits of 'int'

        T ui = u_d[index];
        T vi = v_d[index];
        T wi = w_d[index];
        T ti = temp_d[index];
        T hi = h_d[i];

        // Extract the local normal vector
        T nx = n_x_d[i];
        T ny = n_y_d[i];
        T nz = n_z_d[i];
        
        // Get the tangential component
        T normu = ui * nx + vi * ny + wi * nz;
        ui -= normu * nx;
        vi -= normu * ny;
        wi -= normu * nz;

        T magu = sqrt(ui*ui + vi*vi + wi*wi);
        magu = fmax(magu, (T)1e-6);

        T utau = kappa * magu / log(hi/z0);

        // Zilitinkevich 1995 correlation for thermal roughness
        T z0h;
        if (z0h_in < 0.0) {
            // Note that this uses previous timestep's utau, hence 
            // lags behind by one dt. usually very negligible
            z0h = z0 * exp(z0h_in * sqrt((utau*z0)/(mu/rho)));  
        } else {
            z0h = z0h_in;
        }

        T ts = 0;
        T q  = 0;

        // Initialize variables based on Boundary Condition
        if constexpr (BC_TYPE == 0) { // Neumann
            q = bc_value;
            ts = ti - (q * Pr * log(hi/z0h)) / (fmax(utau, (T)1e-6) * kappa);
        } else {                      // Dirichlet
            ts = bc_value;
            q  = kappa * utau * (ts - ti) / log(hi/z0h);
        }

        T Ri_b;
        T g_dot_n = fabs(g1*nx + g2*ny + g3*nz);

        // Compute Bulk Richardson Number
        if constexpr (BC_TYPE == 0) {
            Ri_b = -g_dot_n * hi / ti * q / (magu * magu * magu * kappa * kappa);
        } else {
            Ri_b = g_dot_n * hi / ti * (ti - ts) / (magu * magu);
        }

        T tau_mag = 0;

        // Stability regime branching
        if (Ri_b > Ri_threshold) { // Stable
            tau_mag = tau_stable<T>(magu, Ri_b, hi, z0, kappa);
            utau = sqrt(tau_mag);
            if constexpr (BC_TYPE == 1) {
                q = heat_flux_stable<T>(ti, ts, Ri_b, hi, z0h, utau, kappa, Pr);
            }
        } 
        else if (Ri_b < -Ri_threshold) { // Convective
            tau_mag = tau_convective<T>(magu, Ri_b, hi, z0, kappa);
            utau = sqrt(tau_mag);
            if constexpr (BC_TYPE == 1) {
                q = heat_flux_convective<T>(ti, ts, Ri_b, hi, magu, z0h, kappa);
            }
        } 
        else { // Neutral
            tau_mag = tau_neutral<T>(magu, hi, z0, kappa);
            utau = sqrt(tau_mag);
            if constexpr (BC_TYPE == 1) {
                q = heat_flux_neutral<T>(ti, ts, hi, z0h, utau, kappa);
            }
        }

        // Apply spatial distribution
        tau_x_d[i] = -rho*tau_mag * ui / magu;
        tau_y_d[i] = -rho*tau_mag * vi / magu;
        tau_z_d[i] = -rho*tau_mag * wi / magu;
        
        // Note: L_ob calculation is omitted from the kernel as it is a pure
        // diagnostic variable and writing it to global memory would require 
        // passing an extra L_ob_d array pointer if GPU diagnostics are needed.
    }
}

#endif // RICHARDSON_KERNEL_H