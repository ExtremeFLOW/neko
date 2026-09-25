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
 * Metal compute kernels for compressible flow update operations.
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
 * Device kernel for updating the temperature T = p / (rho * (gamma - 1))
 */
kernel void update_temperature_kernel(device float *T_field[[ buffer(0) ]],
                                      device const float *p[[ buffer(1) ]],
                                      device const float *rho[[ buffer(2) ]],
                                      constant float &gamma[[ buffer(3) ]],
                                      constant int &n[[ buffer(4) ]],
                                      uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  T_field[idx] = p[idx] / (rho[idx] * (gamma - 1.0f));
}

/**
 * Device kernel for preparing physical Navier-Stokes flux terms
 */
kernel void ns_flux_prepare_kernel(device float *div_flux[[ buffer(0) ]],
                                   device float *dissipation[[ buffer(1) ]],
                                   device float *h1[[ buffer(2) ]],
                                   device const float *dudx[[ buffer(3) ]],
                                   device const float *dudy[[ buffer(4) ]],
                                   device const float *dudz[[ buffer(5) ]],
                                   device const float *dvdx[[ buffer(6) ]],
                                   device const float *dvdy[[ buffer(7) ]],
                                   device const float *dvdz[[ buffer(8) ]],
                                   device const float *dwdx[[ buffer(9) ]],
                                   device const float *dwdy[[ buffer(10) ]],
                                   device const float *dwdz[[ buffer(11) ]],
                                   device const float *mu[[ buffer(12) ]],
                                   constant int &n[[ buffer(13) ]],
                                   uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  const float div_u = dudx[idx] + dvdy[idx] + dwdz[idx];
  const float two_thirds = 2.0f / 3.0f;
  const float tau_xx = mu[idx] * (2.0f * dudx[idx] - two_thirds * div_u);
  const float tau_yy = mu[idx] * (2.0f * dvdy[idx] - two_thirds * div_u);
  const float tau_zz = mu[idx] * (2.0f * dwdz[idx] - two_thirds * div_u);
  const float tau_xy = mu[idx] * (dudy[idx] + dvdx[idx]);
  const float tau_xz = mu[idx] * (dudz[idx] + dwdx[idx]);
  const float tau_yz = mu[idx] * (dvdz[idx] + dwdy[idx]);

  div_flux[idx] = mu[idx] * div_u;
  h1[idx] = mu[idx];
  dissipation[idx] = tau_xx * dudx[idx]
    + tau_xy * (dudy[idx] + dvdx[idx])
    + tau_xz * (dudz[idx] + dwdx[idx])
    + tau_yy * dvdy[idx]
    + tau_yz * (dvdz[idx] + dwdy[idx])
    + tau_zz * dwdz[idx];
}

/**
 * Device kernel for finalizing physical Navier-Stokes flux terms
 */
kernel void ns_flux_finalize_kernel(device float *visc_m_x[[ buffer(0) ]],
                                    device float *visc_m_y[[ buffer(1) ]],
                                    device float *visc_m_z[[ buffer(2) ]],
                                    device float *visc_E[[ buffer(3) ]],
                                    device float *f_x[[ buffer(4) ]],
                                    device float *f_y[[ buffer(5) ]],
                                    device float *f_z[[ buffer(6) ]],
                                    device const float *opgrad_x[[ buffer(7) ]],
                                    device const float *opgrad_y[[ buffer(8) ]],
                                    device const float *opgrad_z[[ buffer(9) ]],
                                    device const float *u[[ buffer(10) ]],
                                    device const float *v[[ buffer(11) ]],
                                    device const float *w[[ buffer(12) ]],
                                    device const float *B[[ buffer(13) ]],
                                    device const float *dissipation[[ buffer(14) ]],
                                    constant int &n[[ buffer(15) ]],
                                    uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  f_x[idx] -= (2.0f / 3.0f) * opgrad_x[idx];
  f_y[idx] -= (2.0f / 3.0f) * opgrad_y[idx];
  f_z[idx] -= (2.0f / 3.0f) * opgrad_z[idx];
  visc_m_x[idx] += f_x[idx];
  visc_m_y[idx] += f_y[idx];
  visc_m_z[idx] += f_z[idx];
  visc_E[idx] += u[idx] * f_x[idx] + v[idx] * f_y[idx] + w[idx] * f_z[idx]
    - B[idx] * dissipation[idx];
}

/**
 * Device kernel for preparing conductive energy flux terms
 */
kernel void ns_flux_temperature_kernel(device float *div_flux[[ buffer(0) ]],
                                       device float *h1[[ buffer(1) ]],
                                       device const float *p[[ buffer(2) ]],
                                       device const float *rho[[ buffer(3) ]],
                                       device const float *kappa[[ buffer(4) ]],
                                       constant float &gamma[[ buffer(5) ]],
                                       constant int &n[[ buffer(6) ]],
                                       uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  div_flux[idx] = p[idx] / (rho[idx] * (gamma - 1.0f));
  h1[idx] = kappa[idx];
}
