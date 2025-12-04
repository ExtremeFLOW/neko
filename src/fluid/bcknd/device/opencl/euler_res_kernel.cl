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

#ifndef __FLUID_EULER_RES_KERNEL__
#define __FLUID_EULER_RES_KERNEL__

__kernel 
void euler_res_part_visc_kernel(__global real * __restrict__ rhs,
                                __global const real * __restrict__ Binv,
                                __global const real * __restrict__ lap_sol,
                                __global const real * __restrict__ h,
                                const real c_avisc,
                                const int n) {
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);
  
  for (int i = idx; i < n; i += str) {
    rhs[i] =  -rhs[i] - c_avisc * h[i] * Binv[i] * lap_sol[i];
  }
}

__kernel
void euler_res_part_mx_flux_kernel(__global real * __restrict__ f_x,
                                   __global real * __restrict__ f_y,
                                   __global real * __restrict__ f_z,
                                   __global const real * __restrict__ m_x,
                                   __global const real * __restrict__ m_y,
                                   __global const real * __restrict__ m_z,
                                   __global const real * __restrict__ rho_field,
                                   __global const real * __restrict__ p,
                                   const int n) {
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);
  
  for (int i = idx; i < n; i += str) {
    f_x[i] = m_x[i] * m_x[i] / rho_field[i] + p[i];
    f_y[i] = m_x[i] * m_y[i] / rho_field[i];
    f_z[i] = m_x[i] * m_z[i] / rho_field[i];
  }
}

__kernel
void euler_res_part_my_flux_kernel(__global real * __restrict__ f_x,
                                   __global real * __restrict__ f_y,
                                   __global real * __restrict__ f_z,
                                   __global const real * __restrict__ m_x,
                                   __global const real * __restrict__ m_y,
                                   __global const real * __restrict__ m_z,
                                   __global const real * __restrict__ rho_field,
                                   __global const real * __restrict__ p,
                                   const int n) {
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);
  
  for (int i = idx; i < n; i += str) {
    f_x[i] = m_y[i] * m_x[i] / rho_field[i];
    f_y[i] = m_y[i] * m_y[i] / rho_field[i] + p[i];
    f_z[i] = m_y[i] * m_z[i] / rho_field[i];
  }
}

__kernel
void euler_res_part_mz_flux_kernel(__global real * __restrict__ f_x,
                                   __global real * __restrict__ f_y,
                                   __global real * __restrict__ f_z,
                                   __global const real * __restrict__ m_x,
                                   __global const real * __restrict__ m_y,
                                   __global const real * __restrict__ m_z,
                                   __global const real * __restrict__ rho_field,
                                   __global const real * __restrict__ p,
                                   const int n) {
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);
  
  for (int i = idx; i < n; i += str) {
    f_x[i] = m_z[i] * m_x[i] / rho_field[i];
    f_y[i] = m_z[i] * m_y[i] / rho_field[i];
    f_z[i] = m_z[i] * m_z[i] / rho_field[i] + p[i];
  }
}

__kernel
void euler_res_part_E_flux_kernel(__global real * __restrict__ f_x,
                                  __global real * __restrict__ f_y,
                                  __global real * __restrict__ f_z,
                                  __global const real * __restrict__ m_x,
                                  __global const real * __restrict__ m_y,
                                  __global const real * __restrict__ m_z,
                                  __global const real * __restrict__ rho_field,
                                  __global const real * __restrict__ p,
                                  __global const real * __restrict__ E,
                                  const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);  
  
  for (int i = idx; i < n; i += str) {
    f_x[i] = (E[i] + p[i]) * m_x[i] / rho_field[i];
    f_y[i] = (E[i] + p[i]) * m_y[i] / rho_field[i];
    f_z[i] = (E[i] + p[i]) * m_z[i] / rho_field[i];
  }
}

__kernel
void euler_res_part_coef_mult_kernel(__global real * __restrict__ rhs_rho,
                                     __global real * __restrict__ rhs_m_x,
                                     __global real * __restrict__ rhs_m_y,
                                     __global real * __restrict__ rhs_m_z,
                                     __global real * __restrict__ rhs_E,
                                     __global const real * __restrict__ mult,
                                     const int n) {
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);  
  
  for (int i = idx; i < n; i += str) {
    rhs_rho[i] = rhs_rho[i] * mult[i];
    rhs_m_x[i] = rhs_m_x[i] * mult[i];
    rhs_m_y[i] = rhs_m_y[i] * mult[i];
    rhs_m_z[i] = rhs_m_z[i] * mult[i];
    rhs_E[i] = rhs_E[i] * mult[i];
  }
}

__kernel
void euler_res_part_rk_sum_kernel(__global real * __restrict__ rho,
                                  __global real * __restrict__ m_x,
                                  __global real * __restrict__ m_y,
                                  __global real * __restrict__ m_z,
                                  __global real * __restrict__ E,
                                  __global const real * __restrict__ k_rho_i,
                                  __global const real * __restrict__ k_m_x_i,
                                  __global const real * __restrict__ k_m_y_i,
                                  __global const real * __restrict__ k_m_z_i,
                                  __global const real * __restrict__ k_E_i,
                                  const real dt,
                                  const real c,
                                  const int n) {
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);  
  
  for (int i = idx; i < n; i += str) {
    rho[i] = rho[i] + dt * c * k_rho_i[i];
    m_x[i] = m_x[i] + dt * c * k_m_x_i[i];
    m_y[i] = m_y[i] + dt * c * k_m_y_i[i];
    m_z[i] = m_z[i] + dt * c * k_m_z_i[i];
    E[i] = E[i] + dt * c * k_E_i[i];
  }
}

#endif // __FLUID_EULER_RES_KERNEL__
