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

#ifndef __FLUID_COMPRESSIBLE_OPS_KERNEL_H__
#define __FLUID_COMPRESSIBLE_OPS_KERNEL_H__

#include <device/device_config.h>

/**
 * Device kernel for compute_max_wave_speed
 */
template< typename T>
__global__ void compute_max_wave_speed_kernel(T * __restrict__ max_wave_speed,
                                              const T * __restrict__ u,
                                              const T * __restrict__ v,
                                              const T * __restrict__ w,
                                              const T gamma,
                                              const T * __restrict__ p,
                                              const T * __restrict__ rho,
                                              const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    T vel_mag = sqrt(u[i]*u[i] + v[i]*v[i] + w[i]*w[i]);
    T sound_speed = sqrt(gamma * p[i] / rho[i]);
    max_wave_speed[i] = vel_mag + sound_speed;
  }

}

/**
 * Device kernel for compute_entropy
 */
template< typename T>
__global__ void compute_entropy_kernel(T * __restrict__ S,
                                       const T * __restrict__ p,
                                       const T * __restrict__ rho,
                                       const T gamma,
                                       const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    // S = 1/(gamma-1) * rho * (log(p) - gamma * log(rho))
    T log_p = log(p[i]);
    T log_rho = log(rho[i]);
    S[i] = (1.0 / (gamma - 1.0)) * rho[i] * (log_p - gamma * log_rho);
  }

}

/**
 * Device kernel for update u,v,w
 */
template< typename T>
__global__ void update_uvw_kernel(T* __restrict__ u,
                                  T* __restrict__ v,
                                  T* __restrict__ w,
                                  const T* __restrict__ m_x,
                                  const T* __restrict__ m_y,
                                  const T* __restrict__ m_z,
                                  const T* __restrict__ rho,
                                  const int n) {


  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    u[i] = m_x[i] / rho[i];
    v[i] = m_y[i] / rho[i];
    w[i] = m_z[i] / rho[i];
  }

}


/**
 * Device kernel for update m_x, m_y, m_z, ruvw
 */
template< typename T>
__global__ void update_mxyz_p_ruvw_kernel(T* __restrict__ m_x,
                                          T* __restrict__ m_y,
                                          T* __restrict__ m_z,
                                          T* __restrict__ p,
                                          T* __restrict__ ruvw,
                                          const T* __restrict__ u,
                                          const T* __restrict__ v,
                                          const T* __restrict__ w,
                                          const T* __restrict__ E,
                                          const T* __restrict__ rho,
                                          const T gamma,
                                          const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    m_x[i] = u[i] * rho[i];
    m_y[i] = v[i] * rho[i];
    m_z[i] = w[i] * rho[i];

    /* Update p = (gamma - 1) * (E - 0.5 * rho * (u^2 + v^2 + w^2)) */
    const real tmp = 0.5 * rho[i] * (u[i]*u[i] + v[i]*v[i] + w[i]*w[i]);
    p[i] = (gamma - 1.0) * (E[i] - tmp);
    ruvw[i] = tmp;
  }
}

#define MAX(a,b) (((a)>(b))?(a):(b))

/**
 * Device kernel for update E
 */
template< typename T>
__global__ void update_e_kernel(T* __restrict__ E,
                                T* __restrict__ p,
                                const T* __restrict__ ruvw,
                                const T gamma,
                                const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    /* Ensure pressure is positive */
    p[i] = MAX(p[i], 1e-12);
    /* E = p / (gamma - 1) + 0.5 * rho * (u^2 + v^2 + w^2) */
    E[i] = p[i] * (1.0 / (gamma - 1.0)) + ruvw[i];
  }
}

/**
 * Device kernel for update T = p / (rho * (gamma - 1))
 */
template< typename T>
__global__ void update_temperature_kernel(T* __restrict__ T_field,
                                          const T* __restrict__ p,
                                          const T* __restrict__ rho,
                                          const T gamma,
                                          const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    T_field[i] = p[i] / (rho[i] * (gamma - 1.0));
  }
}

/**
 * Device kernel for preparing physical Navier-Stokes flux terms.
 */
template< typename T>
__global__ void ns_flux_prepare_kernel(T* __restrict__ div_flux,
                                       T* __restrict__ dissipation,
                                       T* __restrict__ h1,
                                       const T* __restrict__ dudx,
                                       const T* __restrict__ dudy,
                                       const T* __restrict__ dudz,
                                       const T* __restrict__ dvdx,
                                       const T* __restrict__ dvdy,
                                       const T* __restrict__ dvdz,
                                       const T* __restrict__ dwdx,
                                       const T* __restrict__ dwdy,
                                       const T* __restrict__ dwdz,
                                       const T* __restrict__ mu,
                                       const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    const T div_u = dudx[i] + dvdy[i] + dwdz[i];
    const T two_thirds = 2.0 / 3.0;
    const T tau_xx = mu[i] * (2.0 * dudx[i] - two_thirds * div_u);
    const T tau_yy = mu[i] * (2.0 * dvdy[i] - two_thirds * div_u);
    const T tau_zz = mu[i] * (2.0 * dwdz[i] - two_thirds * div_u);
    const T tau_xy = mu[i] * (dudy[i] + dvdx[i]);
    const T tau_xz = mu[i] * (dudz[i] + dwdx[i]);
    const T tau_yz = mu[i] * (dvdz[i] + dwdy[i]);

    div_flux[i] = mu[i] * div_u;
    h1[i] = mu[i];
    dissipation[i] = tau_xx * dudx[i] +
      tau_xy * (dudy[i] + dvdx[i]) +
      tau_xz * (dudz[i] + dwdx[i]) +
      tau_yy * dvdy[i] +
      tau_yz * (dvdz[i] + dwdy[i]) +
      tau_zz * dwdz[i];
  }
}

/**
 * Device kernel for finalizing physical Navier-Stokes flux terms.
 */
template< typename T>
__global__ void ns_flux_finalize_kernel(T* __restrict__ visc_m_x,
                                        T* __restrict__ visc_m_y,
                                        T* __restrict__ visc_m_z,
                                        T* __restrict__ visc_E,
                                        T* __restrict__ f_x,
                                        T* __restrict__ f_y,
                                        T* __restrict__ f_z,
                                        const T* __restrict__ opgrad_x,
                                        const T* __restrict__ opgrad_y,
                                        const T* __restrict__ opgrad_z,
                                        const T* __restrict__ u,
                                        const T* __restrict__ v,
                                        const T* __restrict__ w,
                                        const T* __restrict__ B,
                                        const T* __restrict__ dissipation,
                                        const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    f_x[i] -= (2.0 / 3.0) * opgrad_x[i];
    f_y[i] -= (2.0 / 3.0) * opgrad_y[i];
    f_z[i] -= (2.0 / 3.0) * opgrad_z[i];
    visc_m_x[i] += f_x[i];
    visc_m_y[i] += f_y[i];
    visc_m_z[i] += f_z[i];
    visc_E[i] += u[i] * f_x[i] + v[i] * f_y[i] + w[i] * f_z[i] -
      B[i] * dissipation[i];
  }
}

/**
 * Device kernel for preparing conductive energy flux terms.
 */
template< typename T>
__global__ void ns_flux_temperature_kernel(T* __restrict__ div_flux,
                                           T* __restrict__ h1,
                                           const T* __restrict__ p,
                                           const T* __restrict__ rho,
                                           const T* __restrict__ kappa,
                                           const T gamma,
                                           const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    div_flux[i] = p[i] / (rho[i] * (gamma - 1.0));
    h1[i] = kappa[i];
  }
}

#endif // __FLUID_COMPRESSIBLE_OPS_KERNEL_H__
