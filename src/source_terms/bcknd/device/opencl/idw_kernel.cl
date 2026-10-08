#ifndef __SOURCE_TERMS_IDW_KERNEL_CL__
#define __SOURCE_TERMS_IDW_KERNEL_CL__
/*
 Copyright (c) 2024-2026, The Neko Authors
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
 * Inverse distance weight, mirrors inv_dist_weight() in direct_forcing_source_term.f90
 */
inline real inv_dist_weight(const real r, const real rmax,
                            const real p, const real eps) {
  if (r >= rmax)
    return (real) 0.0;
  return pow((rmax - r) / (rmax * r + eps), p);
}

/**
 * Gather form (atomic-free). One work-item per node of an active element;
 * each node gathers from the Lagrangian points listed for its element in the
 * CSR transpose (el_off/el_lag) and writes its slot exactly once.
 */
__kernel void idw_gather_one_sided(__global real * __restrict__ fu,
                                   __global real * __restrict__ fv,
                                   __global real * __restrict__ fw,
                                   __global const real * __restrict__ fu_ib,
                                   __global const real * __restrict__ fv_ib,
                                   __global const real * __restrict__ fw_ib,
                                   __global const real * __restrict__ fum_ib,
                                   __global const real * __restrict__ fvm_ib,
                                   __global const real * __restrict__ fwm_ib,
                                   __global const real * __restrict__ x,
                                   __global const real * __restrict__ y,
                                   __global const real * __restrict__ z,
                                   __global const real * __restrict__ ds,
                                   __global const real * __restrict__ pmsk,
                                   __global const real * __restrict__ w,
                                   __global const real * __restrict__ wm,
                                   __global const real * __restrict__ lpx,
                                   __global const real * __restrict__ lpy,
                                   __global const real * __restrict__ lpz,
                                   __global const int  * __restrict__ active_el,
                                   __global const int  * __restrict__ el_off,
                                   __global const int  * __restrict__ el_lag,
                                   const int  n_active,
                                   const int  lx3,
                                   const real dt,
                                   const real rmax,
                                   const real pwr,
                                   const real eps,
                                   const real wtol) {

  const int total = n_active * lx3;
  const int str   = get_global_size(0);

  for (int gid = get_global_id(0); gid < total; gid += str) {
    const int a   = gid / lx3;
    const int nd  = gid - a * lx3;
    const int e   = active_el[a];          /* 0-based element     */
    const int idx = e * lx3 + nd;          /* == x(j,k,l,e) flat  */

    const real xn  = x[idx];
    const real yn  = y[idx];
    const real zn  = z[idx];
    const real dsn = ds[idx];
    const int  lo  = el_off[e];
    const int  hi  = el_off[e + 1];

    real au = (real) 0.0, av = (real) 0.0, aw = (real) 0.0;

    if (pmsk[idx] > (real) 0.0) {
      for (int k = lo; k < hi; ++k) {
        const int  i  = el_lag[k];
        const real dx = xn - lpx[i];
        const real dy = yn - lpy[i];
        const real dz = zn - lpz[i];
        const real r  = sqrt(dx*dx + dy*dy + dz*dz) / dsn;
        const real idw = inv_dist_weight(r, rmax, pwr, eps);
        au -= fu_ib[i] * idw;
        av -= fv_ib[i] * idw;
        aw -= fw_ib[i] * idw;
      }
      const real wv = w[idx];
      if (fabs(wv) > wtol) {
        const real s = (real) 1.0 / (wv * dt);
        fu[idx] += au * s;                 /* unique writer -> no atomics */
        fv[idx] += av * s;
        fw[idx] += aw * s;
      }
    } else {
      for (int k = lo; k < hi; ++k) {
        const int  i  = el_lag[k];
        const real dx = xn - lpx[i];
        const real dy = yn - lpy[i];
        const real dz = zn - lpz[i];
        const real r  = sqrt(dx*dx + dy*dy + dz*dz) / dsn;
        const real idw = inv_dist_weight(r, rmax, pwr, eps);
        au -= fum_ib[i] * idw;
        av -= fvm_ib[i] * idw;
        aw -= fwm_ib[i] * idw;
      }
      const real wmv = wm[idx];
      if (fabs(wmv) > wtol) {
        const real s = (real) 1.0 / (wmv * dt);
        fu[idx] += au * s;
        fv[idx] += av * s;
        fw[idx] += aw * s;
      }
    }
  }
}

#endif // __SOURCE_TERMS_IDW_KERNEL_CL__

/**
 * Shepard / adjoint interpolation partial sums. One work-group per
 * Lagrangian point: the work-items stride over the nodes of every element
 * in the point's list (lag_off/lag_els), accumulate the eight partial sums
 * in registers (plus side: u, v, w, weight; minus side: u, v, w, weight)
 * and reduce them in local memory. Mirrors idw_interp_shepard_partials()
 * in direct_forcing_source_term.f90: with adjoint != 0 the weight of a node is
 * K * mult / (Binv * w_side), zero where |w_side| <= wtol, Binv being the
 * inverse of the assembled mass matrix (coef%Binv) and w_p / w_m the
 * assembled spread weights. The work-group
 * size must be 256.
 */
__kernel void idw_interp_partials(__global real * __restrict__ part,
                                  __global const real * __restrict__ u,
                                  __global const real * __restrict__ v,
                                  __global const real * __restrict__ w,
                                  __global const real * __restrict__ x,
                                  __global const real * __restrict__ y,
                                  __global const real * __restrict__ z,
                                  __global const real * __restrict__ ds,
                                  __global const real * __restrict__ pmsk,
                                  __global const real * __restrict__ mult,
                                  __global const real * __restrict__ Binv,
                                  __global const real * __restrict__ w_p,
                                  __global const real * __restrict__ w_m,
                                  __global const real * __restrict__ lpx,
                                  __global const real * __restrict__ lpy,
                                  __global const real * __restrict__ lpz,
                                  __global const int  * __restrict__ lag_off,
                                  __global const int  * __restrict__ lag_els,
                                  const int  n_lag,
                                  const int  lx3,
                                  const real rmax_i,
                                  const real pwr,
                                  const real eps,
                                  const real wtol,
                                  const int  adjoint) {

  __local real red[8][256];

  const int i = get_group_id(0);
  if (i >= n_lag)
    return;

  const int tid = get_local_id(0);
  const int nt  = get_local_size(0);
  const real xp = lpx[i];
  const real yp = lpy[i];
  const real zp = lpz[i];

  real acc[8];
  for (int k = 0; k < 8; ++k)
    acc[k] = (real) 0.0;

  for (int ee = lag_off[i]; ee < lag_off[i + 1]; ++ee) {
    const int e = lag_els[ee];
    for (int nd = tid; nd < lx3; nd += nt) {
      const int idx = e * lx3 + nd;
      const real dx = x[idx] - xp;
      const real dy = y[idx] - yp;
      const real dz = z[idx] - zp;
      const real r = sqrt(dx*dx + dy*dy + dz*dz) / ds[idx];
      real wgt = inv_dist_weight(r, rmax_i, pwr, eps) * mult[idx];
      if (pmsk[idx] > (real) 0.0) {
        if (adjoint) {
          const real sw = w_p[idx];
          wgt = (fabs(sw) > wtol) ? wgt / (Binv[idx] * sw) : (real) 0.0;
        }
        acc[0] += wgt * u[idx];
        acc[1] += wgt * v[idx];
        acc[2] += wgt * w[idx];
        acc[3] += wgt;
      } else {
        if (adjoint) {
          const real sw = w_m[idx];
          wgt = (fabs(sw) > wtol) ? wgt / (Binv[idx] * sw) : (real) 0.0;
        }
        acc[4] += wgt * u[idx];
        acc[5] += wgt * v[idx];
        acc[6] += wgt * w[idx];
        acc[7] += wgt;
      }
    }
  }

  for (int k = 0; k < 8; ++k)
    red[k][tid] = acc[k];
  barrier(CLK_LOCAL_MEM_FENCE);

  for (int s = nt / 2; s > 0; s >>= 1) {
    if (tid < s) {
      for (int k = 0; k < 8; ++k)
        red[k][tid] += red[k][tid + s];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }

  if (tid == 0) {
    for (int k = 0; k < 8; ++k)
      part[8 * i + k] = red[k][0];
  }
}
