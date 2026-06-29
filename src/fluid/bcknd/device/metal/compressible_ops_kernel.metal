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
 * Metal compute kernels for compressible flow operations.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/**
 * Device kernel for updating the velocity components u = m / rho
 */
kernel void update_uvw_kernel(device float *u[[ buffer(0) ]],
                              device float *v[[ buffer(1) ]],
                              device float *w[[ buffer(2) ]],
                              device const float *m_x[[ buffer(3) ]],
                              device const float *m_y[[ buffer(4) ]],
                              device const float *m_z[[ buffer(5) ]],
                              device const float *rho[[ buffer(6) ]],
                              constant int &n[[ buffer(7) ]],
                              uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  u[idx] = m_x[idx] / rho[idx];
  v[idx] = m_y[idx] / rho[idx];
  w[idx] = m_z[idx] / rho[idx];
}

/**
 * Device kernel for updating momentum, pressure and kinetic energy
 */
kernel void update_mxyz_p_ruvw_kernel(device float *m_x[[ buffer(0) ]],
                                      device float *m_y[[ buffer(1) ]],
                                      device float *m_z[[ buffer(2) ]],
                                      device float *p[[ buffer(3) ]],
                                      device float *ruvw[[ buffer(4) ]],
                                      device const float *u[[ buffer(5) ]],
                                      device const float *v[[ buffer(6) ]],
                                      device const float *w[[ buffer(7) ]],
                                      device const float *E[[ buffer(8) ]],
                                      device const float *rho[[ buffer(9) ]],
                                      constant float &gamma[[ buffer(10) ]],
                                      constant int &n[[ buffer(11) ]],
                                      uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  m_x[idx] = u[idx] * rho[idx];
  m_y[idx] = v[idx] * rho[idx];
  m_z[idx] = w[idx] * rho[idx];

  /* Update p = (gamma - 1) * (E - 0.5 * rho * (u^2 + v^2 + w^2)) */
  const float tmp = 0.5f * rho[idx] * (u[idx] * u[idx]
                                       + v[idx] * v[idx]
                                       + w[idx] * w[idx]);
  p[idx] = (gamma - 1.0f) * (E[idx] - tmp);
  ruvw[idx] = tmp;
}

/**
 * Device kernel for updating the total energy
 */
kernel void update_e_kernel(device float *E[[ buffer(0) ]],
                            device float *p[[ buffer(1) ]],
                            device const float *ruvw[[ buffer(2) ]],
                            constant float &gamma[[ buffer(3) ]],
                            constant int &n[[ buffer(4) ]],
                            uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  /* Ensure pressure is positive */
  p[idx] = fmax(p[idx], 1e-12f);
  /* E = p / (gamma - 1) + 0.5 * rho * (u^2 + v^2 + w^2) */
  E[idx] = p[idx] * (1.0f / (gamma - 1.0f)) + ruvw[idx];
}

/**
 * Device kernel for computing the entropy
 * S = 1/(gamma-1) * rho * (log(p) - gamma * log(rho))
 */
kernel void compute_entropy_kernel(device float *S[[ buffer(0) ]],
                                   device const float *p[[ buffer(1) ]],
                                   device const float *rho[[ buffer(2) ]],
                                   constant float &gamma[[ buffer(3) ]],
                                   constant int &n[[ buffer(4) ]],
                                   uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  const float log_p = log(p[idx]);
  const float log_rho = log(rho[idx]);
  S[idx] = (1.0f / (gamma - 1.0f)) * rho[idx] * (log_p - gamma * log_rho);
}

/**
 * Device kernel for computing the maximum wave speed |u| + c
 */
kernel void compute_max_wave_speed_kernel(device float *max_wave_speed[[ buffer(0) ]],
                                          device const float *u[[ buffer(1) ]],
                                          device const float *v[[ buffer(2) ]],
                                          device const float *w[[ buffer(3) ]],
                                          constant float &gamma[[ buffer(4) ]],
                                          device const float *p[[ buffer(5) ]],
                                          device const float *rho[[ buffer(6) ]],
                                          constant int &n[[ buffer(7) ]],
                                          uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  const float vel_mag = sqrt(u[idx] * u[idx]
                             + v[idx] * v[idx]
                             + w[idx] * w[idx]);
  const float sound_speed = sqrt(gamma * p[idx] / rho[idx]);
  max_wave_speed[idx] = vel_mag + sound_speed;
}
