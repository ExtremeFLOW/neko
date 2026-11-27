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
void compressible_res_part_visc_kernel(__global real * __restrict__ rhs_rho,
                                __global real * __restrict__ rhs_m_x,
                                __global real * __restrict__ rhs_m_y,
                                __global real * __restrict__ rhs_m_z,
                                __global real * __restrict__ rhs_E,
                                __global const real * __restrict__ visc_rho,
                                __global const real * __restrict__ visc_m_x,
                                __global const real * __restrict__ visc_m_y,
                                __global const real * __restrict__ visc_m_z,
                                __global const real * __restrict__ visc_E,
                                __global const real * __restrict__ mult,
                                __global const real * __restrict__ Binv,
                                __global real * __restrict__ h1,
                                const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  /* rhs = -(rhs * mult + Binv * visc), evaluated before the gather-scatter:
     mult and Binv are continuous, so the exchange of the sum equals the
     weighted sum of the two exchanges. */
  for (int i = idx; i < n; i += str) {
    const real mi = mult[i];
    const real Bi = Binv[i];
    rhs_rho[i] = -(rhs_rho[i] * mi + Bi * visc_rho[i]);
    rhs_m_x[i] = -(rhs_m_x[i] * mi + Bi * visc_m_x[i]);
    rhs_m_y[i] = -(rhs_m_y[i] * mi + Bi * visc_m_y[i]);
    rhs_m_z[i] = -(rhs_m_z[i] * mi + Bi * visc_m_z[i]);
    rhs_E[i] = -(rhs_E[i] * mi + Bi * visc_E[i]);
    h1[i] = 1.0;
  }
}

/*
 * Fused inviscid divergence
 *
 * One work-group per element forms the density, momentum and energy fluxes
 * from the conserved variables, stages them in local memory a few at a time
 * and differentiates them there, so no flux component round-trips through
 * global memory and the geometric factors are read once per element rather
 * than once per divergence. The momentum flux tensor is symmetric: its three
 * off-diagonal components are differentiated once and feed two divergences
 * each, 12 fluxes and 36 local derivatives per point instead of 15 and 45.
 *
 * The discrete operator is the one div() computes,
 *   jacinv * sum_i (D_r F_i drdx_i + D_s F_i dsdx_i + D_t F_i dtdx_i),
 * with the x, y and z terms accumulated in that order.
 */

/* Flux f of the 12, in staging order: the mass flux m, the six distinct
   components of the symmetric momentum flux rho u u^T + p I, and the
   energy flux (E + p) u from the velocity the stage state holds. */
inline real inviscid_div_flux(const int f, const real rho, const real m_x,
                              const real m_y, const real m_z, const real E,
                              const real p, const real u, const real v,
                              const real w) {
  switch (f) {
  case 0: return m_x;
  case 1: return m_y;
  case 2: return m_z;
  case 3: return m_x * m_x / rho + p;
  case 4: return m_x * m_y / rho;
  case 5: return m_x * m_z / rho;
  case 6: return m_y * m_y / rho + p;
  case 7: return m_y * m_z / rho;
  case 8: return m_z * m_z / rho + p;
  case 9: return (E + p) * u;
  case 10: return (E + p) * v;
  default: return (E + p) * w;
  }
}

/* Feed the local derivatives (rt, st, tt) of flux f into the divergences
   it belongs to. An off-diagonal momentum flux F_ab is the b component of
   the flux of m_a and the a component of the flux of m_b. The staging order
   keeps every divergence summed x, then y, then z, as div() does. */
inline void inviscid_div_accumulate(const int f, const real rt,
                                    const real st, const real tt,
                                    const real drdx, const real dsdx,
                                    const real dtdx, const real drdy,
                                    const real dsdy, const real dtdy,
                                    const real drdz, const real dsdz,
                                    const real dtdz, const real jacinv,
                                    real *rho, real *m_x, real *m_y,
                                    real *m_z, real *E) {
  const real tx = ((rt * drdx) + (st * dsdx) + (tt * dtdx)) * jacinv;
  const real ty = ((rt * drdy) + (st * dsdy) + (tt * dtdy)) * jacinv;
  const real tz = ((rt * drdz) + (st * dsdz) + (tt * dtdz)) * jacinv;
  switch (f) {
  case 0: *rho = tx; break;
  case 1: *rho = *rho + ty; break;
  case 2: *rho = *rho + tz; break;
  case 3: *m_x = tx; break;
  case 4: *m_x = *m_x + ty; *m_y = tx; break;
  case 5: *m_x = *m_x + tz; *m_z = tx; break;
  case 6: *m_y = *m_y + ty; break;
  case 7: *m_y = *m_y + tz; *m_z = *m_z + ty; break;
  case 8: *m_z = *m_z + tz; break;
  case 9: *E = tx; break;
  case 10: *E = *E + ty; break;
  default: *E = *E + tz; break;
  }
}

/* Points per work-item */
#define INVISCID_DIV_NCH(LX, CHUNKS)                                           \
  (((LX) * (LX) * (LX) + (CHUNKS) - 1) / (CHUNKS))

/*
 * One kernel per tuning candidate, see the host (compressible_res.c) for the
 * candidates it enumerates and the geometry each maps to. CHUNKS is the
 * work-group size, NC the number of fluxes staged per pass and R which of
 * the point state (bit 0) and the geometric factors (bit 1) of a work-item's
 * points are held in registers across the passes rather than re-read from
 * global memory, where an element's data stays cache resident.
 */
#define DEFINE_INVISCID_DIV_KERNEL(LX, CHUNKS, NC, R)                           \
__kernel void                                                                  \
compressible_res_inviscid_div_kernel_lx##LX##_ch##CHUNKS##_nc##NC##_r##R(      \
    __global real * __restrict__ rhs_rho,                                     \
    __global real * __restrict__ rhs_m_x,                                     \
    __global real * __restrict__ rhs_m_y,                                     \
    __global real * __restrict__ rhs_m_z,                                     \
    __global real * __restrict__ rhs_E,                                       \
    __global const real * __restrict__ rho,                                   \
    __global const real * __restrict__ m_x,                                   \
    __global const real * __restrict__ m_y,                                   \
    __global const real * __restrict__ m_z,                                   \
    __global const real * __restrict__ E,                                     \
    __global const real * __restrict__ p,                                     \
    __global const real * __restrict__ u,                                     \
    __global const real * __restrict__ v,                                     \
    __global const real * __restrict__ w,                                     \
    __global const real * __restrict__ drdx,                                  \
    __global const real * __restrict__ dsdx,                                  \
    __global const real * __restrict__ dtdx,                                  \
    __global const real * __restrict__ drdy,                                  \
    __global const real * __restrict__ dsdy,                                  \
    __global const real * __restrict__ dtdy,                                  \
    __global const real * __restrict__ drdz,                                  \
    __global const real * __restrict__ dsdz,                                  \
    __global const real * __restrict__ dtdz,                                  \
    __global const real * __restrict__ dx,                                    \
    __global const real * __restrict__ dy,                                    \
    __global const real * __restrict__ dz,                                    \
    __global const real * __restrict__ jacinv) {                              \
                                                                               \
  __local real shf[(NC) * (LX) * (LX) * (LX)];                                 \
  __local real shdx[(LX) * (LX)];                                              \
  __local real shdy[(LX) * (LX)];                                              \
  __local real shdz[(LX) * (LX)];                                              \
                                                                               \
  const int e = get_group_id(0);                                               \
  const int tid = get_local_id(0);                                             \
  const int ele = e * (LX) * (LX) * (LX);                                      \
                                                                               \
  for (int l = tid; l < (LX) * (LX); l += (CHUNKS)) {                          \
    shdx[l] = dx[l];                                                           \
    shdy[l] = dy[l];                                                           \
    shdz[l] = dz[l];                                                           \
  }                                                                            \
                                                                               \
  real r_rho[INVISCID_DIV_NCH(LX, CHUNKS)];                                    \
  real r_mx[INVISCID_DIV_NCH(LX, CHUNKS)];                                     \
  real r_my[INVISCID_DIV_NCH(LX, CHUNKS)];                                     \
  real r_mz[INVISCID_DIV_NCH(LX, CHUNKS)];                                     \
  real r_E[INVISCID_DIV_NCH(LX, CHUNKS)];                                      \
  real r_p[INVISCID_DIV_NCH(LX, CHUNKS)];                                      \
  real r_u[INVISCID_DIV_NCH(LX, CHUNKS)];                                      \
  real r_v[INVISCID_DIV_NCH(LX, CHUNKS)];                                      \
  real r_w[INVISCID_DIV_NCH(LX, CHUNKS)];                                      \
  real g_rx[INVISCID_DIV_NCH(LX, CHUNKS)];                                     \
  real g_sx[INVISCID_DIV_NCH(LX, CHUNKS)];                                     \
  real g_tx[INVISCID_DIV_NCH(LX, CHUNKS)];                                     \
  real g_ry[INVISCID_DIV_NCH(LX, CHUNKS)];                                     \
  real g_sy[INVISCID_DIV_NCH(LX, CHUNKS)];                                     \
  real g_ty[INVISCID_DIV_NCH(LX, CHUNKS)];                                     \
  real g_rz[INVISCID_DIV_NCH(LX, CHUNKS)];                                     \
  real g_sz[INVISCID_DIV_NCH(LX, CHUNKS)];                                     \
  real g_tz[INVISCID_DIV_NCH(LX, CHUNKS)];                                     \
  real g_j[INVISCID_DIV_NCH(LX, CHUNKS)];                                      \
  real a_rho[INVISCID_DIV_NCH(LX, CHUNKS)];                                    \
  real a_mx[INVISCID_DIV_NCH(LX, CHUNKS)];                                     \
  real a_my[INVISCID_DIV_NCH(LX, CHUNKS)];                                     \
  real a_mz[INVISCID_DIV_NCH(LX, CHUNKS)];                                     \
  real a_E[INVISCID_DIV_NCH(LX, CHUNKS)];                                      \
                                                                               \
  /* Whatever is held in registers is loaded once, the rest is re-read in   \
     every pass that needs it */                                               \
  if ((R) != 0) {                                                              \
    for (int n = 0; n < INVISCID_DIV_NCH(LX, CHUNKS); ++n) {                   \
      const int ijk = tid + n * (CHUNKS);                                      \
      if (ijk < (LX) * (LX) * (LX)) {                                          \
        const int g = ijk + ele;                                               \
        if ((R) & 1) {                                                         \
          r_rho[n] = rho[g];                                                   \
          r_mx[n] = m_x[g];                                                    \
          r_my[n] = m_y[g];                                                    \
          r_mz[n] = m_z[g];                                                    \
          r_E[n] = E[g];                                                       \
          r_p[n] = p[g];                                                       \
          r_u[n] = u[g];                                                       \
          r_v[n] = v[g];                                                       \
          r_w[n] = w[g];                                                       \
        }                                                                      \
        if ((R) & 2) {                                                         \
          g_rx[n] = drdx[g];                                                   \
          g_sx[n] = dsdx[g];                                                   \
          g_tx[n] = dtdx[g];                                                   \
          g_ry[n] = drdy[g];                                                   \
          g_sy[n] = dsdy[g];                                                   \
          g_ty[n] = dtdy[g];                                                   \
          g_rz[n] = drdz[g];                                                   \
          g_sz[n] = dsdz[g];                                                   \
          g_tz[n] = dtdz[g];                                                   \
          g_j[n] = jacinv[g];                                                  \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  }                                                                            \
                                                                               \
  for (int f = 0; f < 12; ++f) {                                               \
    const int slot = f % (NC);                                                 \
                                                                               \
    /* Stage the next NC fluxes at this work-item's points */                  \
    if (slot == 0) {                                                           \
      for (int c = 0; c < (NC); ++c) {                                         \
        for (int n = 0; n < INVISCID_DIV_NCH(LX, CHUNKS); ++n) {               \
          const int ijk = tid + n * (CHUNKS);                                  \
          if (ijk < (LX) * (LX) * (LX)) {                                      \
            const int g = ijk + ele;                                           \
            shf[c * (LX) * (LX) * (LX) + ijk] = ((R) & 1) ?                    \
              inviscid_div_flux(f + c, r_rho[n], r_mx[n], r_my[n], r_mz[n],    \
                                r_E[n], r_p[n], r_u[n], r_v[n], r_w[n]) :      \
              inviscid_div_flux(f + c, rho[g], m_x[g], m_y[g], m_z[g],         \
                                E[g], p[g], u[g], v[g], w[g]);                 \
          }                                                                    \
        }                                                                      \
      }                                                                        \
      barrier(CLK_LOCAL_MEM_FENCE);                                            \
    }                                                                          \
                                                                               \
    /* Local derivatives of flux f and its share of the divergences */         \
    for (int n = 0; n < INVISCID_DIV_NCH(LX, CHUNKS); ++n) {                   \
      const int ijk = tid + n * (CHUNKS);                                      \
      if (ijk < (LX) * (LX) * (LX)) {                                          \
        const int jk = ijk / (LX);                                             \
        const int i = ijk - jk * (LX);                                         \
        const int k = jk / (LX);                                               \
        const int j = jk - k * (LX);                                           \
        const int so = slot * (LX) * (LX) * (LX);                              \
        real rt = 0.0;                                                         \
        real st = 0.0;                                                         \
        real tt = 0.0;                                                         \
        for (int l = 0; l < (LX); ++l) {                                       \
          rt += shdx[i + l * (LX)] * shf[so + l + j * (LX) + k * (LX) * (LX)]; \
          st += shdy[j + l * (LX)] * shf[so + i + l * (LX) + k * (LX) * (LX)]; \
          tt += shdz[k + l * (LX)] * shf[so + i + j * (LX) + l * (LX) * (LX)]; \
        }                                                                      \
        if ((R) & 2) {                                                         \
          inviscid_div_accumulate(f, rt, st, tt,                               \
                                  g_rx[n], g_sx[n], g_tx[n],                   \
                                  g_ry[n], g_sy[n], g_ty[n],                   \
                                  g_rz[n], g_sz[n], g_tz[n], g_j[n],           \
                                  &a_rho[n], &a_mx[n], &a_my[n], &a_mz[n],     \
                                  &a_E[n]);                                    \
        } else {                                                               \
          const int g = ijk + ele;                                             \
          inviscid_div_accumulate(f, rt, st, tt,                               \
                                  drdx[g], dsdx[g], dtdx[g],                   \
                                  drdy[g], dsdy[g], dtdy[g],                   \
                                  drdz[g], dsdz[g], dtdz[g], jacinv[g],        \
                                  &a_rho[n], &a_mx[n], &a_my[n], &a_mz[n],     \
                                  &a_E[n]);                                    \
        }                                                                      \
      }                                                                        \
    }                                                                          \
                                                                               \
    /* The staged fluxes are consumed, the slots may be refilled */            \
    if (slot == (NC) - 1) {                                                    \
      barrier(CLK_LOCAL_MEM_FENCE);                                            \
    }                                                                          \
  }                                                                            \
                                                                               \
  for (int n = 0; n < INVISCID_DIV_NCH(LX, CHUNKS); ++n) {                     \
    const int ijk = tid + n * (CHUNKS);                                        \
    if (ijk < (LX) * (LX) * (LX)) {                                            \
      const int g = ijk + ele;                                                 \
      rhs_rho[g] = a_rho[n];                                                   \
      rhs_m_x[g] = a_mx[n];                                                    \
      rhs_m_y[g] = a_my[n];                                                    \
      rhs_m_z[g] = a_mz[n];                                                    \
      rhs_E[g] = a_E[n];                                                       \
    }                                                                          \
  }                                                                            \
}

/*
 * The distinct kernels of the tuning candidates, one line per geometry the
 * host can ask for: the work-group size covers the LX^3 points in as few
 * chunks of at most 256 work-items as needed, rounded up to 64; three fluxes
 * are staged while they fit 48 kB of local memory in double precision
 * (LX <= 12); registers are offered up to four points per work-item
 * (LX <= 10). Regenerate with the host's rules when they change.
 */
DEFINE_INVISCID_DIV_KERNEL(2, 64, 3, 0)
DEFINE_INVISCID_DIV_KERNEL(2, 64, 3, 1)
DEFINE_INVISCID_DIV_KERNEL(2, 64, 3, 2)
DEFINE_INVISCID_DIV_KERNEL(2, 64, 3, 3)
DEFINE_INVISCID_DIV_KERNEL(2, 64, 1, 0)
DEFINE_INVISCID_DIV_KERNEL(2, 64, 1, 1)
DEFINE_INVISCID_DIV_KERNEL(2, 64, 1, 2)
DEFINE_INVISCID_DIV_KERNEL(2, 64, 1, 3)
DEFINE_INVISCID_DIV_KERNEL(3, 64, 3, 0)
DEFINE_INVISCID_DIV_KERNEL(3, 64, 3, 1)
DEFINE_INVISCID_DIV_KERNEL(3, 64, 3, 2)
DEFINE_INVISCID_DIV_KERNEL(3, 64, 3, 3)
DEFINE_INVISCID_DIV_KERNEL(3, 64, 1, 0)
DEFINE_INVISCID_DIV_KERNEL(3, 64, 1, 1)
DEFINE_INVISCID_DIV_KERNEL(3, 64, 1, 2)
DEFINE_INVISCID_DIV_KERNEL(3, 64, 1, 3)
DEFINE_INVISCID_DIV_KERNEL(4, 64, 3, 0)
DEFINE_INVISCID_DIV_KERNEL(4, 64, 3, 1)
DEFINE_INVISCID_DIV_KERNEL(4, 64, 3, 2)
DEFINE_INVISCID_DIV_KERNEL(4, 64, 3, 3)
DEFINE_INVISCID_DIV_KERNEL(4, 64, 1, 0)
DEFINE_INVISCID_DIV_KERNEL(4, 64, 1, 1)
DEFINE_INVISCID_DIV_KERNEL(4, 64, 1, 2)
DEFINE_INVISCID_DIV_KERNEL(4, 64, 1, 3)
DEFINE_INVISCID_DIV_KERNEL(5, 128, 3, 0)
DEFINE_INVISCID_DIV_KERNEL(5, 128, 3, 1)
DEFINE_INVISCID_DIV_KERNEL(5, 128, 3, 2)
DEFINE_INVISCID_DIV_KERNEL(5, 128, 3, 3)
DEFINE_INVISCID_DIV_KERNEL(5, 128, 1, 0)
DEFINE_INVISCID_DIV_KERNEL(5, 128, 1, 1)
DEFINE_INVISCID_DIV_KERNEL(5, 128, 1, 2)
DEFINE_INVISCID_DIV_KERNEL(5, 128, 1, 3)
DEFINE_INVISCID_DIV_KERNEL(6, 256, 3, 0)
DEFINE_INVISCID_DIV_KERNEL(6, 256, 3, 1)
DEFINE_INVISCID_DIV_KERNEL(6, 256, 3, 2)
DEFINE_INVISCID_DIV_KERNEL(6, 256, 3, 3)
DEFINE_INVISCID_DIV_KERNEL(6, 256, 1, 0)
DEFINE_INVISCID_DIV_KERNEL(6, 256, 1, 1)
DEFINE_INVISCID_DIV_KERNEL(6, 256, 1, 2)
DEFINE_INVISCID_DIV_KERNEL(6, 256, 1, 3)
DEFINE_INVISCID_DIV_KERNEL(7, 192, 3, 0)
DEFINE_INVISCID_DIV_KERNEL(7, 192, 3, 1)
DEFINE_INVISCID_DIV_KERNEL(7, 192, 3, 2)
DEFINE_INVISCID_DIV_KERNEL(7, 192, 3, 3)
DEFINE_INVISCID_DIV_KERNEL(7, 192, 1, 0)
DEFINE_INVISCID_DIV_KERNEL(7, 192, 1, 1)
DEFINE_INVISCID_DIV_KERNEL(7, 192, 1, 2)
DEFINE_INVISCID_DIV_KERNEL(7, 192, 1, 3)
DEFINE_INVISCID_DIV_KERNEL(8, 256, 3, 0)
DEFINE_INVISCID_DIV_KERNEL(8, 256, 3, 1)
DEFINE_INVISCID_DIV_KERNEL(8, 256, 3, 2)
DEFINE_INVISCID_DIV_KERNEL(8, 256, 3, 3)
DEFINE_INVISCID_DIV_KERNEL(8, 256, 1, 0)
DEFINE_INVISCID_DIV_KERNEL(8, 256, 1, 1)
DEFINE_INVISCID_DIV_KERNEL(8, 256, 1, 2)
DEFINE_INVISCID_DIV_KERNEL(8, 256, 1, 3)
DEFINE_INVISCID_DIV_KERNEL(9, 256, 3, 0)
DEFINE_INVISCID_DIV_KERNEL(9, 256, 3, 1)
DEFINE_INVISCID_DIV_KERNEL(9, 256, 3, 2)
DEFINE_INVISCID_DIV_KERNEL(9, 256, 3, 3)
DEFINE_INVISCID_DIV_KERNEL(9, 256, 1, 0)
DEFINE_INVISCID_DIV_KERNEL(9, 256, 1, 1)
DEFINE_INVISCID_DIV_KERNEL(9, 256, 1, 2)
DEFINE_INVISCID_DIV_KERNEL(9, 256, 1, 3)
DEFINE_INVISCID_DIV_KERNEL(10, 256, 3, 0)
DEFINE_INVISCID_DIV_KERNEL(10, 256, 3, 1)
DEFINE_INVISCID_DIV_KERNEL(10, 256, 3, 2)
DEFINE_INVISCID_DIV_KERNEL(10, 256, 3, 3)
DEFINE_INVISCID_DIV_KERNEL(10, 256, 1, 0)
DEFINE_INVISCID_DIV_KERNEL(10, 256, 1, 1)
DEFINE_INVISCID_DIV_KERNEL(10, 256, 1, 2)
DEFINE_INVISCID_DIV_KERNEL(10, 256, 1, 3)
DEFINE_INVISCID_DIV_KERNEL(11, 256, 3, 0)
DEFINE_INVISCID_DIV_KERNEL(11, 256, 1, 0)
DEFINE_INVISCID_DIV_KERNEL(12, 256, 3, 0)
DEFINE_INVISCID_DIV_KERNEL(12, 256, 1, 0)
DEFINE_INVISCID_DIV_KERNEL(13, 256, 1, 0)
DEFINE_INVISCID_DIV_KERNEL(14, 256, 1, 0)
DEFINE_INVISCID_DIV_KERNEL(15, 256, 1, 0)
DEFINE_INVISCID_DIV_KERNEL(16, 256, 1, 0)

__kernel
void compressible_res_part_rk_sum_kernel(__global real * __restrict__ rho,
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
