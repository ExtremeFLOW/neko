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

#endif // __OPENCL_COMPRESSIBLE_OPS_UPDATE_KERNEL__ 
