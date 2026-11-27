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
 * Metal compute kernels for the compressible residual.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/**
 * Device kernel for the viscous part of the compressible residual
 */
kernel void compressible_res_part_visc_kernel(device float *rhs_rho[[ buffer(0) ]],
                                       device float *rhs_m_x[[ buffer(1) ]],
                                       device float *rhs_m_y[[ buffer(2) ]],
                                       device float *rhs_m_z[[ buffer(3) ]],
                                       device float *rhs_E[[ buffer(4) ]],
                                       device const float *visc_rho[[ buffer(5) ]],
                                       device const float *visc_m_x[[ buffer(6) ]],
                                       device const float *visc_m_y[[ buffer(7) ]],
                                       device const float *visc_m_z[[ buffer(8) ]],
                                       device const float *visc_E[[ buffer(9) ]],
                                       device const float *mult[[ buffer(10) ]],
                                       device const float *Binv[[ buffer(11) ]],
                                       device float *h1[[ buffer(12) ]],
                                       constant int &n[[ buffer(13) ]],
                                       uint idx [[ thread_position_in_grid ]]) {
  if (idx >= (uint)n) return;
  /* rhs = -(rhs * mult + Binv * visc), evaluated before the gather-scatter:
     mult and Binv are continuous, so the exchange of the sum equals the
     weighted sum of the two exchanges. */
  const float mi = mult[idx];
  const float Bi = Binv[idx];
  rhs_rho[idx] = -(rhs_rho[idx] * mi + Bi * visc_rho[idx]);
  rhs_m_x[idx] = -(rhs_m_x[idx] * mi + Bi * visc_m_x[idx]);
  rhs_m_y[idx] = -(rhs_m_y[idx] * mi + Bi * visc_m_y[idx]);
  rhs_m_z[idx] = -(rhs_m_z[idx] * mi + Bi * visc_m_z[idx]);
  rhs_E[idx] = -(rhs_E[idx] * mi + Bi * visc_E[idx]);
  h1[idx] = 1.0f;
}

/*
 * Fused inviscid divergence
 *
 * One threadgroup per element forms the density, momentum and energy fluxes
 * from the conserved variables, stages them in threadgroup memory a few at
 * a time and differentiates them there, so no flux component round-trips
 * through device memory and the geometric factors are read once per element
 * rather than once per divergence. The momentum flux tensor is symmetric:
 * its three off-diagonal components are differentiated once and feed two
 * divergences each, 12 fluxes and 36 local derivatives per point instead of
 * 15 and 45.
 *
 * The discrete operator is the one div() computes,
 *   jacinv * sum_i (D_r F_i drdx_i + D_s F_i dsdx_i + D_t F_i dtdx_i),
 * with the x, y and z terms accumulated in that order.
 */

/* One kernel is instantiated per distinct tuning candidate, see the host
   (compressible_res.m) for the candidates it enumerates: the threadgroup
   size, the fluxes staged per pass and which of the point state and the
   geometric factors a thread holds in registers across the passes. */

/* Flux f of the 12, in staging order: the mass flux m, the six distinct
   components of the symmetric momentum flux rho u u^T + p I, and the
   energy flux (E + p) u from the velocity the stage state holds. */
inline float inviscid_div_flux(const int f, const float rho, const float m_x,
                               const float m_y, const float m_z,
                               const float E, const float p, const float u,
                               const float v, const float w) {
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
inline void inviscid_div_accumulate(const int f, const float rt,
                                    const float st, const float tt,
                                    const float drdx, const float dsdx,
                                    const float dtdx, const float drdy,
                                    const float dsdy, const float dtdy,
                                    const float drdz, const float dsdz,
                                    const float dtdz, const float jacinv,
                                    thread float &rho, thread float &m_x,
                                    thread float &m_y, thread float &m_z,
                                    thread float &E) {
  const float tx = ((rt * drdx) + (st * dsdx) + (tt * dtdx)) * jacinv;
  const float ty = ((rt * drdy) + (st * dsdy) + (tt * dtdy)) * jacinv;
  const float tz = ((rt * drdz) + (st * dsdz) + (tt * dtdz)) * jacinv;
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

template <int LX, int CHUNKS, int NC, bool REGS, bool REGM>
void inviscid_div_impl(device float *rhs_rho,
                       device float *rhs_m_x,
                       device float *rhs_m_y,
                       device float *rhs_m_z,
                       device float *rhs_E,
                       device const float *rho,
                       device const float *m_x,
                       device const float *m_y,
                       device const float *m_z,
                       device const float *E,
                       device const float *p,
                       device const float *u,
                       device const float *v,
                       device const float *w,
                       device const float *drdx,
                       device const float *dsdx,
                       device const float *dtdx,
                       device const float *drdy,
                       device const float *dsdy,
                       device const float *dtdy,
                       device const float *drdz,
                       device const float *dsdz,
                       device const float *dtdz,
                       device const float *dx,
                       device const float *dy,
                       device const float *dz,
                       device const float *jacinv,
                       threadgroup float *shf,
                       threadgroup float *shdx,
                       threadgroup float *shdy,
                       threadgroup float *shdz,
                       uint e, uint tid) {

  constexpr int NPTS = LX * LX * LX;
  /* Points per thread */
  constexpr int NCH = (NPTS + CHUNKS - 1) / CHUNKS;

  const int ele = (int) e * NPTS;
  const int t = (int) tid;

  for (int l = t; l < LX * LX; l += CHUNKS) {
    shdx[l] = dx[l];
    shdy[l] = dy[l];
    shdz[l] = dz[l];
  }

  float r_rho[NCH], r_mx[NCH], r_my[NCH], r_mz[NCH], r_E[NCH], r_p[NCH];
  float r_u[NCH], r_v[NCH], r_w[NCH];
  float g_rx[NCH], g_sx[NCH], g_tx[NCH];
  float g_ry[NCH], g_sy[NCH], g_ty[NCH];
  float g_rz[NCH], g_sz[NCH], g_tz[NCH];
  float g_j[NCH];
  float a_rho[NCH], a_mx[NCH], a_my[NCH], a_mz[NCH], a_E[NCH];

  /* Whatever is held in registers is loaded once, the rest is re-read in
     every pass that needs it: the state when a flux is staged, the
     geometric factors when a derivative is accumulated */
  if (REGS || REGM) {
#pragma unroll
    for (int n = 0; n < NCH; ++n) {
      const int ijk = t + n * CHUNKS;
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

  /* Fully unrolled, so that every index into the per-thread arrays is a
     compile-time constant and they live in registers. With the loops
     kept, the Metal compiler on AMD GPUs was seen to lose the
     accumulators of all but the first point of a thread at NCH >= 3. */
#pragma unroll
  for (int f = 0; f < 12; ++f) {
    const int slot = f % NC;

    /* Stage the next NC fluxes at this thread's points */
    if (slot == 0) {
#pragma unroll
      for (int c = 0; c < NC; ++c) {
#pragma unroll
        for (int n = 0; n < NCH; ++n) {
          const int ijk = t + n * CHUNKS;
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
      threadgroup_barrier(mem_flags::mem_threadgroup);
    }

    /* Local derivatives of flux f and its share of the divergences */
#pragma unroll
    for (int n = 0; n < NCH; ++n) {
      const int ijk = t + n * CHUNKS;
      if (ijk < NPTS) {
        const int jk = ijk / LX;
        const int i = ijk - jk * LX;
        const int k = jk / LX;
        const int j = jk - k * LX;
        threadgroup const float *sf = shf + slot * NPTS;
        float rt = 0.0f;
        float st = 0.0f;
        float tt = 0.0f;
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
      threadgroup_barrier(mem_flags::mem_threadgroup);
    }
  }

#pragma unroll
  for (int n = 0; n < NCH; ++n) {
    const int ijk = t + n * CHUNKS;
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

/* Kernel for one geometry: CHUNKS threads per threadgroup, NC fluxes staged
   per pass, R the register choice (bit 0 state, bit 1 geometric factors).
   The host builds the pipeline with CHUNKS as the threadgroup limit
   (neko_metal_new_pipeline), so the compiler sizes the register use for
   exactly the threadgroup that is dispatched. */
#define INSTANTIATE_INVISCID_DIV(LX, CHUNKS, NC, R)                             \
kernel void                                                                     \
compressible_res_inviscid_div_kernel_lx##LX##_ch##CHUNKS##_nc##NC##_r##R(      \
    device float *rhs_rho [[ buffer(0) ]],                                      \
    device float *rhs_m_x [[ buffer(1) ]],                                      \
    device float *rhs_m_y [[ buffer(2) ]],                                      \
    device float *rhs_m_z [[ buffer(3) ]],                                      \
    device float *rhs_E [[ buffer(4) ]],                                        \
    device const float *rho [[ buffer(5) ]],                                    \
    device const float *m_x [[ buffer(6) ]],                                    \
    device const float *m_y [[ buffer(7) ]],                                    \
    device const float *m_z [[ buffer(8) ]],                                    \
    device const float *E [[ buffer(9) ]],                                      \
    device const float *p [[ buffer(10) ]],                                     \
    device const float *u [[ buffer(11) ]],                                     \
    device const float *v [[ buffer(12) ]],                                     \
    device const float *w [[ buffer(13) ]],                                     \
    device const float *drdx [[ buffer(14) ]],                                  \
    device const float *dsdx [[ buffer(15) ]],                                  \
    device const float *dtdx [[ buffer(16) ]],                                  \
    device const float *drdy [[ buffer(17) ]],                                  \
    device const float *dsdy [[ buffer(18) ]],                                  \
    device const float *dtdy [[ buffer(19) ]],                                  \
    device const float *drdz [[ buffer(20) ]],                                  \
    device const float *dsdz [[ buffer(21) ]],                                  \
    device const float *dtdz [[ buffer(22) ]],                                  \
    device const float *dx [[ buffer(23) ]],                                    \
    device const float *dy [[ buffer(24) ]],                                    \
    device const float *dz [[ buffer(25) ]],                                    \
    device const float *jacinv [[ buffer(26) ]],                                \
    uint e [[ threadgroup_position_in_grid ]],                                  \
    uint tid [[ thread_position_in_threadgroup ]]) {                            \
  threadgroup float shf[(NC) * LX * LX * LX];                                   \
  threadgroup float shdx[LX * LX];                                              \
  threadgroup float shdy[LX * LX];                                              \
  threadgroup float shdz[LX * LX];                                              \
  inviscid_div_impl<LX, CHUNKS, NC, ((R) & 1) != 0, ((R) & 2) != 0>(           \
      rhs_rho, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E, rho, m_x, m_y, m_z, E, p,      \
      u, v, w, drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz,            \
      dx, dy, dz, jacinv, shf, shdx, shdy, shdz, e, tid);                       \
}

/*
 * The distinct kernels of the tuning candidates, one line per geometry the
 * host can ask for: the threadgroup covers the LX^3 points in as few chunks
 * of at most 512 threads as needed, or twice as many up to eight points per
 * thread, rounded up to 32; three fluxes are staged while they fit 32 kB of
 * threadgroup memory (LX <= 13); registers are offered up to four points per
 * thread. Regenerate with the host's rules when they change.
 */
INSTANTIATE_INVISCID_DIV(2, 32, 3, 0)
INSTANTIATE_INVISCID_DIV(2, 32, 3, 1)
INSTANTIATE_INVISCID_DIV(2, 32, 3, 2)
INSTANTIATE_INVISCID_DIV(2, 32, 3, 3)
INSTANTIATE_INVISCID_DIV(2, 32, 1, 0)
INSTANTIATE_INVISCID_DIV(2, 32, 1, 1)
INSTANTIATE_INVISCID_DIV(2, 32, 1, 2)
INSTANTIATE_INVISCID_DIV(2, 32, 1, 3)
INSTANTIATE_INVISCID_DIV(3, 32, 3, 0)
INSTANTIATE_INVISCID_DIV(3, 32, 3, 1)
INSTANTIATE_INVISCID_DIV(3, 32, 3, 2)
INSTANTIATE_INVISCID_DIV(3, 32, 3, 3)
INSTANTIATE_INVISCID_DIV(3, 32, 1, 0)
INSTANTIATE_INVISCID_DIV(3, 32, 1, 1)
INSTANTIATE_INVISCID_DIV(3, 32, 1, 2)
INSTANTIATE_INVISCID_DIV(3, 32, 1, 3)
INSTANTIATE_INVISCID_DIV(4, 64, 3, 0)
INSTANTIATE_INVISCID_DIV(4, 64, 3, 1)
INSTANTIATE_INVISCID_DIV(4, 64, 3, 2)
INSTANTIATE_INVISCID_DIV(4, 64, 3, 3)
INSTANTIATE_INVISCID_DIV(4, 64, 1, 0)
INSTANTIATE_INVISCID_DIV(4, 64, 1, 1)
INSTANTIATE_INVISCID_DIV(4, 64, 1, 2)
INSTANTIATE_INVISCID_DIV(4, 64, 1, 3)
INSTANTIATE_INVISCID_DIV(4, 32, 3, 0)
INSTANTIATE_INVISCID_DIV(4, 32, 3, 1)
INSTANTIATE_INVISCID_DIV(4, 32, 3, 2)
INSTANTIATE_INVISCID_DIV(4, 32, 3, 3)
INSTANTIATE_INVISCID_DIV(4, 32, 1, 0)
INSTANTIATE_INVISCID_DIV(4, 32, 1, 1)
INSTANTIATE_INVISCID_DIV(4, 32, 1, 2)
INSTANTIATE_INVISCID_DIV(4, 32, 1, 3)
INSTANTIATE_INVISCID_DIV(5, 128, 3, 0)
INSTANTIATE_INVISCID_DIV(5, 128, 3, 1)
INSTANTIATE_INVISCID_DIV(5, 128, 3, 2)
INSTANTIATE_INVISCID_DIV(5, 128, 3, 3)
INSTANTIATE_INVISCID_DIV(5, 128, 1, 0)
INSTANTIATE_INVISCID_DIV(5, 128, 1, 1)
INSTANTIATE_INVISCID_DIV(5, 128, 1, 2)
INSTANTIATE_INVISCID_DIV(5, 128, 1, 3)
INSTANTIATE_INVISCID_DIV(5, 64, 3, 0)
INSTANTIATE_INVISCID_DIV(5, 64, 3, 1)
INSTANTIATE_INVISCID_DIV(5, 64, 3, 2)
INSTANTIATE_INVISCID_DIV(5, 64, 3, 3)
INSTANTIATE_INVISCID_DIV(5, 64, 1, 0)
INSTANTIATE_INVISCID_DIV(5, 64, 1, 1)
INSTANTIATE_INVISCID_DIV(5, 64, 1, 2)
INSTANTIATE_INVISCID_DIV(5, 64, 1, 3)
INSTANTIATE_INVISCID_DIV(6, 224, 3, 0)
INSTANTIATE_INVISCID_DIV(6, 224, 3, 1)
INSTANTIATE_INVISCID_DIV(6, 224, 3, 2)
INSTANTIATE_INVISCID_DIV(6, 224, 3, 3)
INSTANTIATE_INVISCID_DIV(6, 224, 1, 0)
INSTANTIATE_INVISCID_DIV(6, 224, 1, 1)
INSTANTIATE_INVISCID_DIV(6, 224, 1, 2)
INSTANTIATE_INVISCID_DIV(6, 224, 1, 3)
INSTANTIATE_INVISCID_DIV(6, 128, 3, 0)
INSTANTIATE_INVISCID_DIV(6, 128, 3, 1)
INSTANTIATE_INVISCID_DIV(6, 128, 3, 2)
INSTANTIATE_INVISCID_DIV(6, 128, 3, 3)
INSTANTIATE_INVISCID_DIV(6, 128, 1, 0)
INSTANTIATE_INVISCID_DIV(6, 128, 1, 1)
INSTANTIATE_INVISCID_DIV(6, 128, 1, 2)
INSTANTIATE_INVISCID_DIV(6, 128, 1, 3)
INSTANTIATE_INVISCID_DIV(7, 352, 3, 0)
INSTANTIATE_INVISCID_DIV(7, 352, 3, 1)
INSTANTIATE_INVISCID_DIV(7, 352, 3, 2)
INSTANTIATE_INVISCID_DIV(7, 352, 3, 3)
INSTANTIATE_INVISCID_DIV(7, 352, 1, 0)
INSTANTIATE_INVISCID_DIV(7, 352, 1, 1)
INSTANTIATE_INVISCID_DIV(7, 352, 1, 2)
INSTANTIATE_INVISCID_DIV(7, 352, 1, 3)
INSTANTIATE_INVISCID_DIV(7, 192, 3, 0)
INSTANTIATE_INVISCID_DIV(7, 192, 3, 1)
INSTANTIATE_INVISCID_DIV(7, 192, 3, 2)
INSTANTIATE_INVISCID_DIV(7, 192, 3, 3)
INSTANTIATE_INVISCID_DIV(7, 192, 1, 0)
INSTANTIATE_INVISCID_DIV(7, 192, 1, 1)
INSTANTIATE_INVISCID_DIV(7, 192, 1, 2)
INSTANTIATE_INVISCID_DIV(7, 192, 1, 3)
INSTANTIATE_INVISCID_DIV(8, 512, 3, 0)
INSTANTIATE_INVISCID_DIV(8, 512, 3, 1)
INSTANTIATE_INVISCID_DIV(8, 512, 3, 2)
INSTANTIATE_INVISCID_DIV(8, 512, 3, 3)
INSTANTIATE_INVISCID_DIV(8, 512, 1, 0)
INSTANTIATE_INVISCID_DIV(8, 512, 1, 1)
INSTANTIATE_INVISCID_DIV(8, 512, 1, 2)
INSTANTIATE_INVISCID_DIV(8, 512, 1, 3)
INSTANTIATE_INVISCID_DIV(8, 256, 3, 0)
INSTANTIATE_INVISCID_DIV(8, 256, 3, 1)
INSTANTIATE_INVISCID_DIV(8, 256, 3, 2)
INSTANTIATE_INVISCID_DIV(8, 256, 3, 3)
INSTANTIATE_INVISCID_DIV(8, 256, 1, 0)
INSTANTIATE_INVISCID_DIV(8, 256, 1, 1)
INSTANTIATE_INVISCID_DIV(8, 256, 1, 2)
INSTANTIATE_INVISCID_DIV(8, 256, 1, 3)
INSTANTIATE_INVISCID_DIV(9, 384, 3, 0)
INSTANTIATE_INVISCID_DIV(9, 384, 3, 1)
INSTANTIATE_INVISCID_DIV(9, 384, 3, 2)
INSTANTIATE_INVISCID_DIV(9, 384, 3, 3)
INSTANTIATE_INVISCID_DIV(9, 384, 1, 0)
INSTANTIATE_INVISCID_DIV(9, 384, 1, 1)
INSTANTIATE_INVISCID_DIV(9, 384, 1, 2)
INSTANTIATE_INVISCID_DIV(9, 384, 1, 3)
INSTANTIATE_INVISCID_DIV(9, 192, 3, 0)
INSTANTIATE_INVISCID_DIV(9, 192, 3, 1)
INSTANTIATE_INVISCID_DIV(9, 192, 3, 2)
INSTANTIATE_INVISCID_DIV(9, 192, 3, 3)
INSTANTIATE_INVISCID_DIV(9, 192, 1, 0)
INSTANTIATE_INVISCID_DIV(9, 192, 1, 1)
INSTANTIATE_INVISCID_DIV(9, 192, 1, 2)
INSTANTIATE_INVISCID_DIV(9, 192, 1, 3)
INSTANTIATE_INVISCID_DIV(10, 512, 3, 0)
INSTANTIATE_INVISCID_DIV(10, 512, 3, 1)
INSTANTIATE_INVISCID_DIV(10, 512, 3, 2)
INSTANTIATE_INVISCID_DIV(10, 512, 3, 3)
INSTANTIATE_INVISCID_DIV(10, 512, 1, 0)
INSTANTIATE_INVISCID_DIV(10, 512, 1, 1)
INSTANTIATE_INVISCID_DIV(10, 512, 1, 2)
INSTANTIATE_INVISCID_DIV(10, 512, 1, 3)
INSTANTIATE_INVISCID_DIV(10, 256, 3, 0)
INSTANTIATE_INVISCID_DIV(10, 256, 3, 1)
INSTANTIATE_INVISCID_DIV(10, 256, 3, 2)
INSTANTIATE_INVISCID_DIV(10, 256, 3, 3)
INSTANTIATE_INVISCID_DIV(10, 256, 1, 0)
INSTANTIATE_INVISCID_DIV(10, 256, 1, 1)
INSTANTIATE_INVISCID_DIV(10, 256, 1, 2)
INSTANTIATE_INVISCID_DIV(10, 256, 1, 3)
INSTANTIATE_INVISCID_DIV(11, 448, 3, 0)
INSTANTIATE_INVISCID_DIV(11, 448, 3, 1)
INSTANTIATE_INVISCID_DIV(11, 448, 3, 2)
INSTANTIATE_INVISCID_DIV(11, 448, 3, 3)
INSTANTIATE_INVISCID_DIV(11, 448, 1, 0)
INSTANTIATE_INVISCID_DIV(11, 448, 1, 1)
INSTANTIATE_INVISCID_DIV(11, 448, 1, 2)
INSTANTIATE_INVISCID_DIV(11, 448, 1, 3)
INSTANTIATE_INVISCID_DIV(11, 224, 3, 0)
INSTANTIATE_INVISCID_DIV(11, 224, 1, 0)
INSTANTIATE_INVISCID_DIV(12, 448, 3, 0)
INSTANTIATE_INVISCID_DIV(12, 448, 3, 1)
INSTANTIATE_INVISCID_DIV(12, 448, 3, 2)
INSTANTIATE_INVISCID_DIV(12, 448, 3, 3)
INSTANTIATE_INVISCID_DIV(12, 448, 1, 0)
INSTANTIATE_INVISCID_DIV(12, 448, 1, 1)
INSTANTIATE_INVISCID_DIV(12, 448, 1, 2)
INSTANTIATE_INVISCID_DIV(12, 448, 1, 3)
INSTANTIATE_INVISCID_DIV(12, 224, 3, 0)
INSTANTIATE_INVISCID_DIV(12, 224, 1, 0)
INSTANTIATE_INVISCID_DIV(13, 448, 3, 0)
INSTANTIATE_INVISCID_DIV(13, 448, 1, 0)
INSTANTIATE_INVISCID_DIV(14, 480, 1, 0)
INSTANTIATE_INVISCID_DIV(15, 512, 1, 0)
INSTANTIATE_INVISCID_DIV(16, 512, 1, 0)

/**
 * Device kernel for the Runge-Kutta stage summation
 */
kernel void compressible_res_part_rk_sum_kernel(device float *rho[[ buffer(0) ]],
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
