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

#endif // __FLUID_COMPRESSIBLE_OPS_KERNEL_H__
