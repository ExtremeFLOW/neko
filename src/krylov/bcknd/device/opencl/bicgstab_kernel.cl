#ifndef __KRYLOV_BICGSTAB_KERNEL_CL__
#define __KRYLOV_BICGSTAB_KERNEL_CL__
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
 * Device kernel for the BiCGStab search direction update,
 * \f$ p = r + \beta (p - \omega v) \f$
 */
__kernel void bicgstab_update_p_kernel(__global real* __restrict__ p,
                                       __global const real* __restrict__ r,
                                       __global const real* __restrict__ v,
                                       const real beta,
                                       const real omega,
                                       const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    p[i] = r[i] + beta * (p[i] - omega * v[i]);
  }
}

/**
 * Device kernel for a weighted inner product and squared norm in one pass,
 * \f$ (a^T M b, b^T M b) \f$
 *
 * The partial sums of the second value follow those of the first,
 * at @a buf_h[nb + group].
 */
__kernel
void bicgstab_product_and_norm_kernel(__global const real* __restrict__ a,
                                      __global const real* __restrict__ b,
                                      __global const real* __restrict__ mult,
                                      __global real_xp* __restrict__ buf_h,
                                      const int n,
                                      const int nb) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  __local real_xp buf_prod[256]; /* Make this nice...*/
  __local real_xp buf_nrm[256];
  real_xp tmp_prod = 0.0;
  real_xp tmp_nrm = 0.0;

  for (int i = idx; i < n; i += str) {
    const real_xp bi = (real_xp) b[i];
    const real_xp w = (real_xp) mult[i] * bi;
    tmp_prod += (real_xp) a[i] * w;
    tmp_nrm += bi * w;
  }
  buf_prod[get_local_id(0)] = tmp_prod;
  buf_nrm[get_local_id(0)] = tmp_nrm;
  barrier(CLK_LOCAL_MEM_FENCE);

  int i = (get_local_size(0)) >> 1;
  while (i != 0) {
    if (get_local_id(0) < i) {
      buf_prod[get_local_id(0)] += buf_prod[get_local_id(0) + i];
      buf_nrm[get_local_id(0)] += buf_nrm[get_local_id(0) + i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    i = i >> 1;
  }

  if (get_local_id(0) == 0) {
    buf_h[get_group_id(0)] = buf_prod[0];
    buf_h[nb + get_group_id(0)] = buf_nrm[0];
  }
}

/**
 * Device kernel for BiCGStab part 1, \f$ s = r - \alpha v \f$,
 * reducing \f$ s^T M s \f$
 */
__kernel void bicgstab_part1_kernel(__global real* __restrict__ s,
                                    __global const real* __restrict__ r,
                                    __global const real* __restrict__ v,
                                    __global const real* __restrict__ mult,
                                    __global real_xp* __restrict__ buf_h,
                                    const real alpha,
                                    const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  __local real_xp buf[256]; /* Make this nice...*/
  real_xp tmp = 0.0;

  for (int i = idx; i < n; i += str) {
    const real si = r[i] - alpha * v[i];
    s[i] = si;
    tmp += (real_xp) si * (real_xp) mult[i] * (real_xp) si;
  }
  buf[get_local_id(0)] = tmp;
  barrier(CLK_LOCAL_MEM_FENCE);

  int i = (get_local_size(0)) >> 1;
  while (i != 0) {
    if (get_local_id(0) < i) {
      buf[get_local_id(0)] += buf[get_local_id(0) + i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    i = i >> 1;
  }

  if (get_local_id(0) == 0) { buf_h[get_group_id(0)] = buf[0]; }
}

/**
 * Device kernel for BiCGStab part 2,
 * \f$ x = x + \alpha \hat{p} + \omega \hat{s} \f$ and
 * \f$ r = s - \omega t \f$, reducing \f$ (r^T M r, f^T M r) \f$
 *
 * The second value is the rho inner product of the next iteration, which
 * the recurrence would otherwise reduce separately at the top of the loop.
 */
__kernel void bicgstab_part2_kernel(__global real* __restrict__ x,
                                    __global real* __restrict__ r,
                                    __global const real* __restrict__ p_hat,
                                    __global const real* __restrict__ s_hat,
                                    __global const real* __restrict__ s,
                                    __global const real* __restrict__ t,
                                    __global const real* __restrict__ f,
                                    __global const real* __restrict__ mult,
                                    __global real_xp* __restrict__ buf_h,
                                    const real alpha,
                                    const real omega,
                                    const int n,
                                    const int nb) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  __local real_xp buf_rtr[256]; /* Make this nice...*/
  __local real_xp buf_ftr[256];
  real_xp tmp_rtr = 0.0;
  real_xp tmp_ftr = 0.0;

  for (int i = idx; i < n; i += str) {
    x[i] = x[i] + alpha * p_hat[i] + omega * s_hat[i];
    const real ri = s[i] - omega * t[i];
    r[i] = ri;
    const real_xp w = (real_xp) mult[i] * (real_xp) ri;
    tmp_rtr += (real_xp) ri * w;
    tmp_ftr += (real_xp) f[i] * w;
  }
  buf_rtr[get_local_id(0)] = tmp_rtr;
  buf_ftr[get_local_id(0)] = tmp_ftr;
  barrier(CLK_LOCAL_MEM_FENCE);

  int i = (get_local_size(0)) >> 1;
  while (i != 0) {
    if (get_local_id(0) < i) {
      buf_rtr[get_local_id(0)] += buf_rtr[get_local_id(0) + i];
      buf_ftr[get_local_id(0)] += buf_ftr[get_local_id(0) + i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    i = i >> 1;
  }

  if (get_local_id(0) == 0) {
    buf_h[get_group_id(0)] = buf_rtr[0];
    buf_h[nb + get_group_id(0)] = buf_ftr[0];
  }
}

#endif /* __KRYLOV_BICGSTAB_KERNEL_CL__ */
