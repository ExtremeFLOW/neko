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

#ifndef __OPENCL_COMPRESSIBLE_OPS_UPDATE_KERNEL__
#define __OPENCL_COMPRESSIBLE_OPS_UPDATE_KERNEL__

/**
 * Device kernel for update u,v,w
 */
__kernel void update_uvw_kernel(__global real* __restrict__ u,
                                __global real* __restrict__ v,
                                __global real* __restrict__ w,
                                __global const real* __restrict__ m_x,
                                __global const real* __restrict__ m_y,
                                __global const real* __restrict__ m_z,
                                __global const real* __restrict__ rho,
                                const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    u[i] = m_x[i] / rho[i];
    v[i] = m_y[i] / rho[i];
    w[i] = m_z[i] / rho[i];
  }

}


/**
 * Device kernel for update m_x, m_y, m_z, ruvw
 */
__kernel void update_mxyz_p_ruvw_kernel(__global real* __restrict__ m_x,
                                        __global real* __restrict__ m_y,
                                        __global real* __restrict__ m_z,
                                        __global real* __restrict__ p,
                                        __global real* __restrict__ ruvw,
                                        __global const real* __restrict__ u,
                                        __global const real* __restrict__ v,
                                        __global const real* __restrict__ w,
                                        __global const real* __restrict__ E,
                                        __global const real* __restrict__ rho,
                                        const real gamma,
                                        const int n) {
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);
  
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
__kernel void update_e_kernel(__global real* __restrict__ E,
                              __global real* __restrict__ p,
                              __global const real* __restrict__ ruvw,
                              const real gamma,
                              const int n) {
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);
  
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
__kernel void update_temperature_kernel(__global real* __restrict__ T,
                                        __global const real* __restrict__ p,
                                        __global const real* __restrict__ rho,
                                        const real gamma,
                                        const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    T[i] = p[i] / (rho[i] * (gamma - 1.0));
  }
}

/**
 * Device kernel for preparing physical Navier-Stokes flux terms.
 */
__kernel void ns_flux_prepare_kernel(__global real* __restrict__ div_flux,
                                     __global real* __restrict__ dissipation,
                                     __global real* __restrict__ h1,
                                     __global const real* __restrict__ dudx,
                                     __global const real* __restrict__ dudy,
                                     __global const real* __restrict__ dudz,
                                     __global const real* __restrict__ dvdx,
                                     __global const real* __restrict__ dvdy,
                                     __global const real* __restrict__ dvdz,
                                     __global const real* __restrict__ dwdx,
                                     __global const real* __restrict__ dwdy,
                                     __global const real* __restrict__ dwdz,
                                     __global const real* __restrict__ mu,
                                     const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    const real div_u = dudx[i] + dvdy[i] + dwdz[i];
    const real two_thirds = 2.0 / 3.0;
    const real tau_xx = mu[i] * (2.0 * dudx[i] - two_thirds * div_u);
    const real tau_yy = mu[i] * (2.0 * dvdy[i] - two_thirds * div_u);
    const real tau_zz = mu[i] * (2.0 * dwdz[i] - two_thirds * div_u);
    const real tau_xy = mu[i] * (dudy[i] + dvdx[i]);
    const real tau_xz = mu[i] * (dudz[i] + dwdx[i]);
    const real tau_yz = mu[i] * (dvdz[i] + dwdy[i]);

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
__kernel void ns_flux_finalize_kernel(__global real* __restrict__ visc_m_x,
                                      __global real* __restrict__ visc_m_y,
                                      __global real* __restrict__ visc_m_z,
                                      __global real* __restrict__ visc_E,
                                      __global real* __restrict__ f_x,
                                      __global real* __restrict__ f_y,
                                      __global real* __restrict__ f_z,
                                      __global const real* __restrict__ opgrad_x,
                                      __global const real* __restrict__ opgrad_y,
                                      __global const real* __restrict__ opgrad_z,
                                      __global const real* __restrict__ u,
                                      __global const real* __restrict__ v,
                                      __global const real* __restrict__ w,
                                      __global const real* __restrict__ B,
                                      __global const real* __restrict__ dissipation,
                                      const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

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
__kernel void ns_flux_temperature_kernel(__global real* __restrict__ div_flux,
                                         __global real* __restrict__ h1,
                                         __global const real* __restrict__ p,
                                         __global const real* __restrict__ rho,
                                         __global const real* __restrict__ kappa,
                                         const real gamma,
                                         const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    div_flux[i] = p[i] / (rho[i] * (gamma - 1.0));
    h1[i] = kappa[i];
  }
}

#endif // __OPENCL_COMPRESSIBLE_OPS_UPDATE_KERNEL__ 
