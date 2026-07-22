#ifndef __SOURCE_TERMS_IDW_KERNEL_H__
#define __SOURCE_TERMS_IDW_KERNEL_H__
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
 * Inverse distance weight, mirrors inv_dist_weight() in idw_source_term.f90
 */
template< typename T >
__device__ __forceinline__ T inv_dist_weight(const T r, const T rmax,
                                              const T p, const T eps) {
  if (r >= rmax)
    return (T) 0.0;
  return pow((rmax - r) / (rmax * r + eps), p);
}

/**
 * Gather form (atomic-free). One thread per node of an active element; each
 * node gathers from the Lagrangian points listed for its element in the CSR
 * transpose (el_off/el_lag) and writes its slot exactly once.
 */
template< typename T >
__global__ void idw_gather_one_sided_kernel(T * __restrict__ fu,
                                            T * __restrict__ fv,
                                            T * __restrict__ fw,
                                            const T * __restrict__ fu_ib,
                                            const T * __restrict__ fv_ib,
                                            const T * __restrict__ fw_ib,
                                            const T * __restrict__ fum_ib,
                                            const T * __restrict__ fvm_ib,
                                            const T * __restrict__ fwm_ib,
                                            const T * __restrict__ x,
                                            const T * __restrict__ y,
                                            const T * __restrict__ z,
                                            const T * __restrict__ ds,
                                            const T * __restrict__ pmsk,
                                            const T * __restrict__ w,
                                            const T * __restrict__ wm,
                                            const T * __restrict__ lpx,
                                            const T * __restrict__ lpy,
                                            const T * __restrict__ lpz,
                                            const int * __restrict__ active_el,
                                            const int * __restrict__ el_off,
                                            const int * __restrict__ el_lag,
                                            const int n_active,
                                            const int lx3,
                                            const T dt,
                                            const T rmax,
                                            const T pwr,
                                            const T eps,
                                            const T wtol) {

  const int total = n_active * lx3;
  const int str   = blockDim.x * gridDim.x;

  for (int gid = blockIdx.x * blockDim.x + threadIdx.x; gid < total;
       gid += str) {
    const int a   = gid / lx3;
    const int nd  = gid - a * lx3;
    const int e   = active_el[a];          /* 0-based element     */
    const int idx = e * lx3 + nd;          /* == x(j,k,l,e) flat  */

    const T xn  = x[idx];
    const T yn  = y[idx];
    const T zn  = z[idx];
    const T dsn = ds[idx];
    const int lo = el_off[e];
    const int hi = el_off[e + 1];

    T au = (T) 0.0, av = (T) 0.0, aw = (T) 0.0;

    if (pmsk[idx] > (T) 0.0) {
      for (int k = lo; k < hi; ++k) {
        const int i  = el_lag[k];
        const T dx = xn - lpx[i];
        const T dy = yn - lpy[i];
        const T dz = zn - lpz[i];
        const T r  = sqrt(dx*dx + dy*dy + dz*dz) / dsn;
        const T idw = inv_dist_weight(r, rmax, pwr, eps);
        au -= fu_ib[i] * idw;
        av -= fv_ib[i] * idw;
        aw -= fw_ib[i] * idw;
      }
      const T wv = w[idx];
      if (fabs(wv) > wtol) {
        const T s = (T) 1.0 / (wv * dt);
        fu[idx] += au * s;                 /* unique writer -> no atomics */
        fv[idx] += av * s;
        fw[idx] += aw * s;
      }
    } else {
      for (int k = lo; k < hi; ++k) {
        const int i  = el_lag[k];
        const T dx = xn - lpx[i];
        const T dy = yn - lpy[i];
        const T dz = zn - lpz[i];
        const T r  = sqrt(dx*dx + dy*dy + dz*dz) / dsn;
        const T idw = inv_dist_weight(r, rmax, pwr, eps);
        au -= fum_ib[i] * idw;
        av -= fvm_ib[i] * idw;
        aw -= fwm_ib[i] * idw;
      }
      const T wmv = wm[idx];
      if (fabs(wmv) > wtol) {
        const T s = (T) 1.0 / (wmv * dt);
        fu[idx] += au * s;
        fv[idx] += av * s;
        fw[idx] += aw * s;
      }
    }
  }
}

#endif // __SOURCE_TERMS_IDW_KERNEL_H__
