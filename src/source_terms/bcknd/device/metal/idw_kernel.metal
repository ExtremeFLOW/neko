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
 * Metal compute kernel for the atomic-free IDW gather.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/**
 * Inverse distance weight, mirrors inv_dist_weight() in direct_forcing_source_term.f90
 */
inline float inv_dist_weight(const float r, const float rmax,
                             const float p, const float eps) {
  if (r >= rmax)
    return 0.0f;
  /*
   * eps is NEKO_EPS, which is precision dependent, so as r -> 0 the base
   * grows to rmax/eps: ~8.4e6 in fp32 against an overflow threshold of
   * 3.4e38, i.e. base^p is +inf for p >~ 6, and the inf/inf in the w
   * normalisation then gives NaN.  fast-math lets isinf/isnan checks be
   * folded away, so bound the base instead: 2^(127/p) is the largest base
   * with a finite p-th power.  The dp path never reaches its equivalent
   * bound, so this only saturates where fp32 would have overflowed.
   */
  const float base = (rmax - r) / (rmax * r + eps);
  return pow(fmin(base, exp2(127.0f / p)), p);
}

/**
 * Gather form (atomic-free). One thread per node of an active element; each
 * node gathers from the Lagrangian points listed for its element in the CSR
 * transpose (el_off/el_lag) and writes its slot exactly once.
 */
kernel void idw_gather_one_sided(device float * fu [[ buffer(0) ]],
                                 device float * fv [[ buffer(1) ]],
                                 device float * fw [[ buffer(2) ]],
                                 device const float * fu_ib [[ buffer(3) ]],
                                 device const float * fv_ib [[ buffer(4) ]],
                                 device const float * fw_ib [[ buffer(5) ]],
                                 device const float * fum_ib [[ buffer(6) ]],
                                 device const float * fvm_ib [[ buffer(7) ]],
                                 device const float * fwm_ib [[ buffer(8) ]],
                                 device const float * x [[ buffer(9) ]],
                                 device const float * y [[ buffer(10) ]],
                                 device const float * z [[ buffer(11) ]],
                                 device const float * ds [[ buffer(12) ]],
                                 device const float * pmsk [[ buffer(13) ]],
                                 device const float * w [[ buffer(14) ]],
                                 device const float * wm [[ buffer(15) ]],
                                 device const float * lpx [[ buffer(16) ]],
                                 device const float * lpy [[ buffer(17) ]],
                                 device const float * lpz [[ buffer(18) ]],
                                 device const int * active_el [[ buffer(19) ]],
                                 device const int * el_off [[ buffer(20) ]],
                                 device const int * el_lag [[ buffer(21) ]],
                                 constant int & n_active [[ buffer(22) ]],
                                 constant int & lx3 [[ buffer(23) ]],
                                 constant float & dt [[ buffer(24) ]],
                                 constant float & rmax [[ buffer(25) ]],
                                 constant float & pwr [[ buffer(26) ]],
                                 constant float & eps [[ buffer(27) ]],
                                 constant float & wtol [[ buffer(28) ]],
                                 uint tpig [[ thread_position_in_grid ]],
                                 uint tpg [[ threads_per_grid ]]) {

  const int total = n_active * lx3;
  const int str   = (int) tpg;

  for (int gid = (int) tpig; gid < total; gid += str) {
    const int a   = gid / lx3;
    const int nd  = gid - a * lx3;
    const int e   = active_el[a];          /* 0-based element     */
    const int idx = e * lx3 + nd;          /* == x(j,k,l,e) flat  */

    const float xn  = x[idx];
    const float yn  = y[idx];
    const float zn  = z[idx];
    const float dsn = ds[idx];
    const int   lo  = el_off[e];
    const int   hi  = el_off[e + 1];

    float au = 0.0f, av = 0.0f, aw = 0.0f;

    if (pmsk[idx] > 0.0f) {
      for (int k = lo; k < hi; ++k) {
        const int   i  = el_lag[k];
        const float dx = xn - lpx[i];
        const float dy = yn - lpy[i];
        const float dz = zn - lpz[i];
        const float r  = sqrt(dx*dx + dy*dy + dz*dz) / dsn;
        const float idw = inv_dist_weight(r, rmax, pwr, eps);
        au -= fu_ib[i] * idw;
        av -= fv_ib[i] * idw;
        aw -= fw_ib[i] * idw;
      }
      const float wv = w[idx];
      if (fabs(wv) > wtol) {
        const float s = 1.0f / (wv * dt);
        fu[idx] += au * s;                 /* unique writer -> no atomics */
        fv[idx] += av * s;
        fw[idx] += aw * s;
      }
    } else {
      for (int k = lo; k < hi; ++k) {
        const int   i  = el_lag[k];
        const float dx = xn - lpx[i];
        const float dy = yn - lpy[i];
        const float dz = zn - lpz[i];
        const float r  = sqrt(dx*dx + dy*dy + dz*dz) / dsn;
        const float idw = inv_dist_weight(r, rmax, pwr, eps);
        au -= fum_ib[i] * idw;
        av -= fvm_ib[i] * idw;
        aw -= fwm_ib[i] * idw;
      }
      const float wmv = wm[idx];
      if (fabs(wmv) > wtol) {
        const float s = 1.0f / (wmv * dt);
        fu[idx] += au * s;
        fv[idx] += av * s;
        fw[idx] += aw * s;
      }
    }
  }
}

/**
 * Shepard / adjoint interpolation partial sums. One threadgroup per
 * Lagrangian point: the threads stride over the nodes of every element in
 * the point's list (lag_off/lag_els), accumulate the eight partial sums in
 * registers (plus side: u, v, w, weight; minus side: u, v, w, weight) and
 * reduce them in threadgroup memory. Mirrors idw_interp_shepard_partials()
 * in direct_forcing_source_term.f90: with adjoint != 0 the weight of a node is
 * K * mult / (Binv * w_side), zero where |w_side| <= wtol, Binv being the
 * inverse of the assembled mass matrix (coef%Binv) and w_p / w_m the
 * assembled spread weights. The threadgroup
 * size must be 256.
 */
kernel void idw_interp_partials(device float * part [[ buffer(0) ]],
                                device const float * u [[ buffer(1) ]],
                                device const float * v [[ buffer(2) ]],
                                device const float * w [[ buffer(3) ]],
                                device const float * x [[ buffer(4) ]],
                                device const float * y [[ buffer(5) ]],
                                device const float * z [[ buffer(6) ]],
                                device const float * ds [[ buffer(7) ]],
                                device const float * pmsk [[ buffer(8) ]],
                                device const float * mult [[ buffer(9) ]],
                                device const float * Binv [[ buffer(10) ]],
                                device const float * w_p [[ buffer(11) ]],
                                device const float * w_m [[ buffer(12) ]],
                                device const float * lpx [[ buffer(13) ]],
                                device const float * lpy [[ buffer(14) ]],
                                device const float * lpz [[ buffer(15) ]],
                                device const int * lag_off [[ buffer(16) ]],
                                device const int * lag_els [[ buffer(17) ]],
                                constant int & n_lag [[ buffer(18) ]],
                                constant int & lx3 [[ buffer(19) ]],
                                constant float & rmax_i [[ buffer(20) ]],
                                constant float & pwr [[ buffer(21) ]],
                                constant float & eps [[ buffer(22) ]],
                                constant float & wtol [[ buffer(23) ]],
                                constant int & adjoint [[ buffer(24) ]],
                                uint tgid [[ threadgroup_position_in_grid ]],
                                uint tid [[ thread_position_in_threadgroup ]],
                                uint ntg [[ threads_per_threadgroup ]]) {

  threadgroup float red[8][256];

  const int i = (int) tgid;
  if (i >= n_lag)
    return;

  const int nt = (int) ntg;
  const float xp = lpx[i];
  const float yp = lpy[i];
  const float zp = lpz[i];

  float acc[8];
  for (int k = 0; k < 8; ++k)
    acc[k] = 0.0f;

  for (int ee = lag_off[i]; ee < lag_off[i + 1]; ++ee) {
    const int e = lag_els[ee];
    for (int nd = (int) tid; nd < lx3; nd += nt) {
      const int idx = e * lx3 + nd;
      const float dx = x[idx] - xp;
      const float dy = y[idx] - yp;
      const float dz = z[idx] - zp;
      const float r = sqrt(dx*dx + dy*dy + dz*dz) / ds[idx];
      float wgt = inv_dist_weight(r, rmax_i, pwr, eps) * mult[idx];
      if (pmsk[idx] > 0.0f) {
        if (adjoint) {
          const float sw = w_p[idx];
          wgt = (fabs(sw) > wtol) ? wgt / (Binv[idx] * sw) : 0.0f;
        }
        acc[0] += wgt * u[idx];
        acc[1] += wgt * v[idx];
        acc[2] += wgt * w[idx];
        acc[3] += wgt;
      } else {
        if (adjoint) {
          const float sw = w_m[idx];
          wgt = (fabs(sw) > wtol) ? wgt / (Binv[idx] * sw) : 0.0f;
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
  threadgroup_barrier(mem_flags::mem_threadgroup);

  for (int s = nt / 2; s > 0; s >>= 1) {
    if ((int) tid < s) {
      for (int k = 0; k < 8; ++k)
        red[k][tid] += red[k][tid + s];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
  }

  if (tid == 0) {
    for (int k = 0; k < 8; ++k)
      part[8 * i + k] = red[k][0];
  }
}
