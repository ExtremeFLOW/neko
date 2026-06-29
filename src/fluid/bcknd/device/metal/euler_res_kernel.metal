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
 * Metal compute kernels for the compressible Euler residual.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/**
 * Device kernel for the viscous part of the Euler residual
 */
kernel void euler_res_part_visc_kernel(device float *rhs[[ buffer(0) ]],
                                       device const float *Binv[[ buffer(1) ]],
                                       device const float *lap_sol[[ buffer(2) ]],
                                       device const float *effective_visc[[ buffer(3) ]],
                                       constant int &n[[ buffer(4) ]],
                                       uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  rhs[idx] = -rhs[idx] - effective_visc[idx] * Binv[idx] * lap_sol[idx];
}

/**
 * Device kernel for the x-momentum flux
 */
kernel void euler_res_part_mx_flux_kernel(device float *f_x[[ buffer(0) ]],
                                          device float *f_y[[ buffer(1) ]],
                                          device float *f_z[[ buffer(2) ]],
                                          device const float *m_x[[ buffer(3) ]],
                                          device const float *m_y[[ buffer(4) ]],
                                          device const float *m_z[[ buffer(5) ]],
                                          device const float *rho_field[[ buffer(6) ]],
                                          device const float *p[[ buffer(7) ]],
                                          constant int &n[[ buffer(8) ]],
                                          uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  f_x[idx] = m_x[idx] * m_x[idx] / rho_field[idx] + p[idx];
  f_y[idx] = m_x[idx] * m_y[idx] / rho_field[idx];
  f_z[idx] = m_x[idx] * m_z[idx] / rho_field[idx];
}

/**
 * Device kernel for the y-momentum flux
 */
kernel void euler_res_part_my_flux_kernel(device float *f_x[[ buffer(0) ]],
                                          device float *f_y[[ buffer(1) ]],
                                          device float *f_z[[ buffer(2) ]],
                                          device const float *m_x[[ buffer(3) ]],
                                          device const float *m_y[[ buffer(4) ]],
                                          device const float *m_z[[ buffer(5) ]],
                                          device const float *rho_field[[ buffer(6) ]],
                                          device const float *p[[ buffer(7) ]],
                                          constant int &n[[ buffer(8) ]],
                                          uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  f_x[idx] = m_y[idx] * m_x[idx] / rho_field[idx];
  f_y[idx] = m_y[idx] * m_y[idx] / rho_field[idx] + p[idx];
  f_z[idx] = m_y[idx] * m_z[idx] / rho_field[idx];
}

/**
 * Device kernel for the z-momentum flux
 */
kernel void euler_res_part_mz_flux_kernel(device float *f_x[[ buffer(0) ]],
                                          device float *f_y[[ buffer(1) ]],
                                          device float *f_z[[ buffer(2) ]],
                                          device const float *m_x[[ buffer(3) ]],
                                          device const float *m_y[[ buffer(4) ]],
                                          device const float *m_z[[ buffer(5) ]],
                                          device const float *rho_field[[ buffer(6) ]],
                                          device const float *p[[ buffer(7) ]],
                                          constant int &n[[ buffer(8) ]],
                                          uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  f_x[idx] = m_z[idx] * m_x[idx] / rho_field[idx];
  f_y[idx] = m_z[idx] * m_y[idx] / rho_field[idx];
  f_z[idx] = m_z[idx] * m_z[idx] / rho_field[idx] + p[idx];
}

/**
 * Device kernel for the energy flux
 */
kernel void euler_res_part_E_flux_kernel(device float *f_x[[ buffer(0) ]],
                                         device float *f_y[[ buffer(1) ]],
                                         device float *f_z[[ buffer(2) ]],
                                         device const float *m_x[[ buffer(3) ]],
                                         device const float *m_y[[ buffer(4) ]],
                                         device const float *m_z[[ buffer(5) ]],
                                         device const float *rho_field[[ buffer(6) ]],
                                         device const float *p[[ buffer(7) ]],
                                         device const float *E[[ buffer(8) ]],
                                         constant int &n[[ buffer(9) ]],
                                         uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  f_x[idx] = (E[idx] + p[idx]) * m_x[idx] / rho_field[idx];
  f_y[idx] = (E[idx] + p[idx]) * m_y[idx] / rho_field[idx];
  f_z[idx] = (E[idx] + p[idx]) * m_z[idx] / rho_field[idx];
}

/**
 * Device kernel for multiplying the residual with the inverse multiplicity
 */
kernel void euler_res_part_coef_mult_kernel(device float *rhs_rho[[ buffer(0) ]],
                                            device float *rhs_m_x[[ buffer(1) ]],
                                            device float *rhs_m_y[[ buffer(2) ]],
                                            device float *rhs_m_z[[ buffer(3) ]],
                                            device float *rhs_E[[ buffer(4) ]],
                                            device const float *mult[[ buffer(5) ]],
                                            constant int &n[[ buffer(6) ]],
                                            uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  rhs_rho[idx] = rhs_rho[idx] * mult[idx];
  rhs_m_x[idx] = rhs_m_x[idx] * mult[idx];
  rhs_m_y[idx] = rhs_m_y[idx] * mult[idx];
  rhs_m_z[idx] = rhs_m_z[idx] * mult[idx];
  rhs_E[idx] = rhs_E[idx] * mult[idx];
}

/**
 * Device kernel for the Runge-Kutta stage summation
 */
kernel void euler_res_part_rk_sum_kernel(device float *rho[[ buffer(0) ]],
                                         device float *m_x[[ buffer(1) ]],
                                         device float *m_y[[ buffer(2) ]],
                                         device float *m_z[[ buffer(3) ]],
                                         device float *E[[ buffer(4) ]],
                                         device const float *k_rho_i[[ buffer(5) ]],
                                         device const float *k_m_x_i[[ buffer(6) ]],
                                         device const float *k_m_y_i[[ buffer(7) ]],
                                         device const float *k_m_z_i[[ buffer(8) ]],
                                         device const float *k_E_i[[ buffer(9) ]],
                                         constant float &dt[[ buffer(10) ]],
                                         constant float &c[[ buffer(11) ]],
                                         constant int &n[[ buffer(12) ]],
                                         uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  rho[idx] = rho[idx] + dt * c * k_rho_i[idx];
  m_x[idx] = m_x[idx] + dt * c * k_m_x_i[idx];
  m_y[idx] = m_y[idx] + dt * c * k_m_y_i[idx];
  m_z[idx] = m_z[idx] + dt * c * k_m_z_i[idx];
  E[idx] = E[idx] + dt * c * k_E_i[idx];
}
