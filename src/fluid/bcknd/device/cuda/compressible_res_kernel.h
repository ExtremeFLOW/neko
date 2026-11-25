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

#include <stdio.h>

/* Warp width the fused divergence rounds its block to, the widest block it
   launches and the static shared memory it may use per block */
#define INVISCID_DIV_WARP 32
#define INVISCID_DIV_MAX_BLOCK 512
#define INVISCID_DIV_MAX_SMEM 49152

template< typename T >
__global__ void compressible_res_part_visc_kernel(T * __restrict__ rhs_rho,
                                     T * __restrict__ rhs_m_x,
                                     T * __restrict__ rhs_m_y,
                                     T * __restrict__ rhs_m_z,
                                     T * __restrict__ rhs_E,
                                     const T * __restrict__ visc_rho,
                                     const T * __restrict__ visc_m_x,
                                     const T * __restrict__ visc_m_y,
                                     const T * __restrict__ visc_m_z,
                                     const T * __restrict__ visc_E,
                                     const T * __restrict__ mult,
                                     const T * __restrict__ Binv,
                                     T * __restrict__ h1,
                                     const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  /* rhs = -(rhs * mult + Binv * visc), evaluated before the gather-scatter:
     mult and Binv are continuous, so the exchange of the sum equals the
     weighted sum of the two exchanges. */
  for (int i = idx; i < n; i += str) {
    const T mi = mult[i];
    const T Bi = Binv[i];
    rhs_rho[i] = -(rhs_rho[i] * mi + Bi * visc_rho[i]);
    rhs_m_x[i] = -(rhs_m_x[i] * mi + Bi * visc_m_x[i]);
    rhs_m_y[i] = -(rhs_m_y[i] * mi + Bi * visc_m_y[i]);
    rhs_m_z[i] = -(rhs_m_z[i] * mi + Bi * visc_m_z[i]);
    rhs_E[i] = -(rhs_E[i] * mi + Bi * visc_E[i]);
    h1[i] = T(1);
  }
}

/*
 * Fused inviscid divergence
 *
 * One block per element forms the density, momentum and energy fluxes from
 * the conserved variables, stages them in shared memory a few at a time and
 * differentiates them there, so no flux component round-trips through
 * global memory and the geometric factors are read once per element rather
 * than once per divergence. The momentum flux tensor is symmetric: its three
 * off-diagonal components are differentiated once and feed two divergences
 * each, 12 fluxes and 36 local derivatives per point instead of 15 and 45.
 *
 * The discrete operator is the one div() computes,
 *   jacinv * sum_i (D_r F_i drdx_i + D_s F_i dsdx_i + D_t F_i dtdx_i),
 * with the x, y and z terms accumulated in that order.
 */

/*
 * Tuning candidates
 *
 * The kernel's geometry is measured rather than decided: the autotuner in
 * compressible_res.cu times INVISCID_DIV_CANDIDATES candidates on the first
 * call at each order and keeps the fastest. Candidate c combines three
 * choices:
 *
 *   chunks (bit 3)       the LX^3 points of an element are split over as few
 *                        chunks of at most INVISCID_DIV_MAX_BLOCK threads as
 *                        they need, or twice as many. A wider block covers
 *                        the element with fewer points per thread, a narrower
 *                        one lets more elements share a multiprocessor. The
 *                        narrower block is not offered where it would put
 *                        more than INVISCID_DIV_MAX_PTS points on a thread.
 *   staging (bit 2)      three fluxes staged in shared memory per pass, where
 *                        they fit next to the derivative matrices, or one: a
 *                        third of the barriers against three times the
 *                        shared memory per block.
 *   registers (bits 0-1) whether the point state (bit 0) and the geometric
 *                        factors (bit 1) of a thread's points are held in
 *                        registers across the passes or re-read from global
 *                        memory, where an element's data stays cache
 *                        resident. Registers save the re-reads and cost
 *                        occupancy, and which wins was measured to flip with
 *                        the order: on an AMD W5700X re-reading both was
 *                        1.8x faster than holding both at order 8, holding
 *                        both 1.4x faster at order 10, state alone fastest
 *                        at order 12. Neither is offered beyond
 *                        INVISCID_DIV_MAX_REG_PTS points per thread, where
 *                        the arrays would spill.
 *
 * A candidate whose choice is not offered takes the geometry of the one that
 * is, so candidates that round to the same geometry instantiate the same
 * kernel; the tuner skips the repeats, see inviscid_div_cand_dup().
 */
#define INVISCID_DIV_CANDIDATES 16
#define INVISCID_DIV_MAX_PTS 8
#define INVISCID_DIV_MAX_REG_PTS 4

/* Threads per block covering the element in nchunks chunks: the points of
   one chunk, rounded up to whole warps */
__host__ __device__ constexpr int inviscid_div_threads_of(const int lx,
                                                          const int nchunks) {
  return ((((lx * lx * lx + nchunks - 1) / nchunks)
           + INVISCID_DIV_WARP - 1) / INVISCID_DIV_WARP) * INVISCID_DIV_WARP;
}

/* Points per thread at that block */
__host__ __device__ constexpr int inviscid_div_pts_of(const int lx,
                                                      const int nchunks) {
  return (lx * lx * lx + inviscid_div_threads_of(lx, nchunks) - 1)
    / inviscid_div_threads_of(lx, nchunks);
}

/* Fewest chunks of at most INVISCID_DIV_MAX_BLOCK threads */
__host__ __device__ constexpr int inviscid_div_min_chunks(const int lx) {
  return (lx * lx * lx + INVISCID_DIV_MAX_BLOCK - 1) / INVISCID_DIV_MAX_BLOCK;
}

/* Chunks the element's points are split over for candidate c */
__host__ __device__ constexpr int inviscid_div_nchunks(const int lx,
                                                       const int c) {
  return (((c >> 3) & 1) != 0 &&
          inviscid_div_pts_of(lx, 2 * inviscid_div_min_chunks(lx))
          <= INVISCID_DIV_MAX_PTS)
    ? 2 * inviscid_div_min_chunks(lx) : inviscid_div_min_chunks(lx);
}

/* Threads per block of candidate c */
__host__ __device__ constexpr int inviscid_div_chunks(const int lx,
                                                      const int c) {
  return inviscid_div_threads_of(lx, inviscid_div_nchunks(lx, c));
}

/* Points per thread of candidate c */
__host__ __device__ constexpr int inviscid_div_pts(const int lx,
                                                   const int c) {
  return inviscid_div_pts_of(lx, inviscid_div_nchunks(lx, c));
}

/* Fluxes staged per pass, for values of size tsize bytes */
__host__ __device__ constexpr int inviscid_div_stage(const int lx,
                                                     const int c,
                                                     const int tsize) {
  return (((c >> 2) & 1) == 0 &&
          (3 * lx * lx * lx + 3 * lx * lx) * tsize <= INVISCID_DIV_MAX_SMEM)
    ? 3 : 1;
}

/* Point state held in registers? */
__host__ __device__ constexpr bool inviscid_div_regs(const int lx,
                                                     const int c) {
  return (c & 1) != 0 && inviscid_div_pts(lx, c) <= INVISCID_DIV_MAX_REG_PTS;
}

/* Geometric factors held in registers? */
__host__ __device__ constexpr bool inviscid_div_regm(const int lx,
                                                     const int c) {
  return (c & 2) != 0 && inviscid_div_pts(lx, c) <= INVISCID_DIV_MAX_REG_PTS;
}

/* Does candidate c instantiate the same kernel as a lower candidate? */
static inline bool inviscid_div_cand_dup(const int lx, const int c,
                                         const int tsize) {
  for (int d = 0; d < c; d++) {
    if (inviscid_div_chunks(lx, d) == inviscid_div_chunks(lx, c) &&
        inviscid_div_stage(lx, d, tsize) == inviscid_div_stage(lx, c, tsize) &&
        inviscid_div_regs(lx, d) == inviscid_div_regs(lx, c) &&
        inviscid_div_regm(lx, d) == inviscid_div_regm(lx, c)) {
      return true;
    }
  }
  return false;
}

/* Describe candidate c for the log */
static inline const char *inviscid_div_cand_desc(const int lx, const int c,
                                                 const int tsize, char *buf,
                                                 const int len) {
  static const char *regs[] = { "none", "state", "metrics", "both" };
  snprintf(buf, len, "%d thr, %d/pass, regs: %s",
           inviscid_div_chunks(lx, c), inviscid_div_stage(lx, c, tsize),
           regs[(inviscid_div_regs(lx, c) ? 1 : 0) +
                (inviscid_div_regm(lx, c) ? 2 : 0)]);
  return buf;
}

/* Flux f of the 12, in staging order: the mass flux m, the six distinct
   components of the symmetric momentum flux rho u u^T + p I, and the
   energy flux (E + p) u from the velocity the stage state holds. */
template< typename T >
__device__ __forceinline__ T inviscid_div_flux(const int f,
                                               const T rho, const T m_x,
                                               const T m_y, const T m_z,
                                               const T E, const T p,
                                               const T u, const T v,
                                               const T w) {
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
template< typename T >
__device__ __forceinline__ void inviscid_div_accumulate(const int f,
                                                        const T rt,
                                                        const T st,
                                                        const T tt,
                                                        const T drdx,
                                                        const T dsdx,
                                                        const T dtdx,
                                                        const T drdy,
                                                        const T dsdy,
                                                        const T dtdy,
                                                        const T drdz,
                                                        const T dsdz,
                                                        const T dtdz,
                                                        const T jacinv,
                                                        T &rho, T &m_x,
                                                        T &m_y, T &m_z,
                                                        T &E) {
  const T tx = ((rt * drdx) + (st * dsdx) + (tt * dtdx)) * jacinv;
  const T ty = ((rt * drdy) + (st * dsdy) + (tt * dtdy)) * jacinv;
  const T tz = ((rt * drdz) + (st * dsdz) + (tt * dtdz)) * jacinv;
  switch (f) {
  case 0: rho = tx; break;
  case 1: rho = rho + ty; break;
  case 2: rho = rho + tz; break;
  case 3: m_x = tx; break;
  case 4: m_x = m_x + ty; m_y = tx; break;
  case 5: m_x = m_x + tz; m_z = tx; break;
  case 6: m_y = m_y + ty; break;
  case 7: m_y = m_y + tz; m_z = m_z + ty; break;
  case 8: m_z = m_z + tz; break;
  case 9: E = tx; break;
  case 10: E = E + ty; break;
  default: E = E + tz; break;
  }
}

template< typename T, const int LX, const int CHUNKS, const int NC,
          const bool REGS, const bool REGM >
__global__ void __launch_bounds__(CHUNKS)
compressible_res_inviscid_div_kernel(T * __restrict__ rhs_rho,
                                     T * __restrict__ rhs_m_x,
                                     T * __restrict__ rhs_m_y,
                                     T * __restrict__ rhs_m_z,
                                     T * __restrict__ rhs_E,
                                     const T * __restrict__ rho,
                                     const T * __restrict__ m_x,
                                     const T * __restrict__ m_y,
                                     const T * __restrict__ m_z,
                                     const T * __restrict__ E,
                                     const T * __restrict__ p,
                                     const T * __restrict__ u,
                                     const T * __restrict__ v,
                                     const T * __restrict__ w,
                                     const T * __restrict__ drdx,
                                     const T * __restrict__ dsdx,
                                     const T * __restrict__ dtdx,
                                     const T * __restrict__ drdy,
                                     const T * __restrict__ dsdy,
                                     const T * __restrict__ dtdy,
                                     const T * __restrict__ drdz,
                                     const T * __restrict__ dsdz,
                                     const T * __restrict__ dtdz,
                                     const T * __restrict__ dx,
                                     const T * __restrict__ dy,
                                     const T * __restrict__ dz,
                                     const T * __restrict__ jacinv) {

  constexpr int NPTS = LX * LX * LX;
  /* Points per thread */
  constexpr int NCH = (NPTS + CHUNKS - 1) / CHUNKS;

  __shared__ T shf[NC * NPTS];
  __shared__ T shdx[LX * LX];
  __shared__ T shdy[LX * LX];
  __shared__ T shdz[LX * LX];

  const int e = blockIdx.x;
  const int tid = threadIdx.x;
  const int ele = e * NPTS;

  for (int l = tid; l < LX * LX; l += CHUNKS) {
    shdx[l] = dx[l];
    shdy[l] = dy[l];
    shdz[l] = dz[l];
  }

  T r_rho[NCH], r_mx[NCH], r_my[NCH], r_mz[NCH], r_E[NCH], r_p[NCH];
  T r_u[NCH], r_v[NCH], r_w[NCH];
  T g_rx[NCH], g_sx[NCH], g_tx[NCH];
  T g_ry[NCH], g_sy[NCH], g_ty[NCH];
  T g_rz[NCH], g_sz[NCH], g_tz[NCH];
  T g_j[NCH];
  T a_rho[NCH], a_mx[NCH], a_my[NCH], a_mz[NCH], a_E[NCH];

  /* Whatever is held in registers is loaded once, the rest is re-read in
     every pass that needs it: the state when a flux is staged, the
     geometric factors when a derivative is accumulated */
  if (REGS || REGM) {
#pragma unroll
    for (int n = 0; n < NCH; ++n) {
      const int ijk = tid + n * CHUNKS;
      if (ijk < NPTS) {
        const int g = ijk + ele;
        if (REGS) {
          r_rho[n] = rho[g];
          r_mx[n] = m_x[g];
          r_my[n] = m_y[g];
          r_mz[n] = m_z[g];
          r_E[n] = E[g];
          r_p[n] = p[g];
          r_u[n] = u[g];
          r_v[n] = v[g];
          r_w[n] = w[g];
        }
        if (REGM) {
          g_rx[n] = drdx[g];
          g_sx[n] = dsdx[g];
          g_tx[n] = dtdx[g];
          g_ry[n] = drdy[g];
          g_sy[n] = dsdy[g];
          g_ty[n] = dtdy[g];
          g_rz[n] = drdz[g];
          g_sz[n] = dsdz[g];
          g_tz[n] = dtdz[g];
          g_j[n] = jacinv[g];
        }
      }
    }
  }

#pragma unroll
  for (int f = 0; f < 12; ++f) {
    const int slot = f % NC;

    /* Stage the next NC fluxes at this thread's points */
    if (slot == 0) {
#pragma unroll
      for (int c = 0; c < NC; ++c) {
#pragma unroll
        for (int n = 0; n < NCH; ++n) {
          const int ijk = tid + n * CHUNKS;
          if (ijk < NPTS) {
            const int g = ijk + ele;
            shf[c * NPTS + ijk] = REGS ?
              inviscid_div_flux(f + c, r_rho[n], r_mx[n], r_my[n], r_mz[n],
                                r_E[n], r_p[n], r_u[n], r_v[n], r_w[n]) :
              inviscid_div_flux(f + c, rho[g], m_x[g], m_y[g], m_z[g],
                                E[g], p[g], u[g], v[g], w[g]);
          }
        }
      }
      __syncthreads();
    }

    /* Local derivatives of flux f and its share of the divergences */
#pragma unroll
    for (int n = 0; n < NCH; ++n) {
      const int ijk = tid + n * CHUNKS;
      if (ijk < NPTS) {
        const int jk = ijk / LX;
        const int i = ijk - jk * LX;
        const int k = jk / LX;
        const int j = jk - k * LX;
        const T *sf = shf + slot * NPTS;
        T rt = 0.0;
        T st = 0.0;
        T tt = 0.0;
#pragma unroll
        for (int l = 0; l < LX; ++l) {
          rt += shdx[i + l * LX] * sf[l + j * LX + k * LX * LX];
          st += shdy[j + l * LX] * sf[i + l * LX + k * LX * LX];
          tt += shdz[k + l * LX] * sf[i + j * LX + l * LX * LX];
        }
        if (REGM) {
          inviscid_div_accumulate(f, rt, st, tt,
                                  g_rx[n], g_sx[n], g_tx[n],
                                  g_ry[n], g_sy[n], g_ty[n],
                                  g_rz[n], g_sz[n], g_tz[n], g_j[n],
                                  a_rho[n], a_mx[n], a_my[n], a_mz[n], a_E[n]);
        } else {
          const int g = ijk + ele;
          inviscid_div_accumulate(f, rt, st, tt,
                                  drdx[g], dsdx[g], dtdx[g],
                                  drdy[g], dsdy[g], dtdy[g],
                                  drdz[g], dsdz[g], dtdz[g], jacinv[g],
                                  a_rho[n], a_mx[n], a_my[n], a_mz[n], a_E[n]);
        }
      }
    }

    /* The staged fluxes are consumed, the slots may be refilled */
    if (slot == NC - 1) {
      __syncthreads();
    }
  }

#pragma unroll
  for (int n = 0; n < NCH; ++n) {
    const int ijk = tid + n * CHUNKS;
    if (ijk < NPTS) {
      const int g = ijk + ele;
      rhs_rho[g] = a_rho[n];
      rhs_m_x[g] = a_mx[n];
      rhs_m_y[g] = a_my[n];
      rhs_m_z[g] = a_mz[n];
      rhs_E[g] = a_E[n];
    }
  }
}

template< typename T >
__global__ void compressible_res_part_rk_sum_kernel(T * __restrict__ rho,
                                    T * __restrict__ m_x,
                                    T * __restrict__ m_y,
                                    T * __restrict__ m_z,
                                    T * __restrict__ E,
                                    const T * __restrict__ k_rho_i,
                                    const T * __restrict__ k_m_x_i,
                                    const T * __restrict__ k_m_y_i,
                                    const T * __restrict__ k_m_z_i,
                                    const T * __restrict__ k_E_i,
                                    const T dt,
                                    const T c,
                                    const int n) {
  
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  
  for (int i = idx; i < n; i += str) {
    rho[i] = rho[i] + dt * c * k_rho_i[i];
    m_x[i] = m_x[i] + dt * c * k_m_x_i[i];
    m_y[i] = m_y[i] + dt * c * k_m_y_i[i];
    m_z[i] = m_z[i] + dt * c * k_m_z_i[i];
    E[i] = E[i] + dt * c * k_E_i[i];
  }
}

#endif // __FLUID_EULER_RES_KERNEL__
