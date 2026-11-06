#ifndef __MATH_AX_HELM_SVV_FACTORIZED_FULL_KERNEL_H__
#define __MATH_AX_HELM_SVV_FACTORIZED_FULL_KERNEL_H__
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
 * Fused device kernel for the factorized SVV stress (full) operator
 *
 * Evaluates, in a single pass over each element, the standard stress
 * formulation built from D and h1, plus the factorized SVV stress built from
 * the complementary derivatives B = (Br, Bs, Bt) and h1_svv. This mirrors
 * ax_helm_svv_factorized_full_lx in the CPU backend: both contributions share
 * the k-plane sweep, the geometric factors are loaded once, and (au, av, aw)
 * are stored once (no read-modify-write of the output).
 *
 * The kernel is kept separate from the autotuned ax_helm kernels on purpose,
 * so that those keep their register footprint untouched. Within each k-plane
 * the standard and SVV fluxes are formed one after the other, so only the
 * three t-direction fluxes of each need to be live at the same time.
 *
 * PAD selects a shared memory pitch of LX+1, removing bank conflicts when LX
 * is a power of two.
 */

/**
 * Map reference derivatives of (u, v, w) to the reference-space fluxes of the
 * symmetric stress tensor, scaled by dj. The r- and s-fluxes are written to
 * shared memory, the t-fluxes are returned in registers.
 */
template< typename T >
__device__ __forceinline__ void
ax_helm_svv_factorized_stress_flux(const T ur, const T us, const T ut,
                                   const T vr, const T vs, const T vt,
                                   const T wr, const T ws, const T wt,
                                   const T drdx, const T drdy, const T drdz,
                                   const T dsdx, const T dsdy, const T dsdz,
                                   const T dtdx, const T dtdy, const T dtdz,
                                   const T dj,
                                   T & shur, T & shus,
                                   T & shvr, T & shvs,
                                   T & shwr, T & shws,
                                   T & rut, T & rvt, T & rwt) {

  const T u1 = ur * drdx + us * dsdx + ut * dtdx;
  const T u2 = ur * drdy + us * dsdy + ut * dtdy;
  const T u3 = ur * drdz + us * dsdz + ut * dtdz;

  const T v1 = vr * drdx + vs * dsdx + vt * dtdx;
  const T v2 = vr * drdy + vs * dsdy + vt * dtdy;
  const T v3 = vr * drdz + vs * dsdz + vt * dtdz;

  const T w1 = wr * drdx + ws * dsdx + wt * dtdx;
  const T w2 = wr * drdy + ws * dsdy + wt * dtdy;
  const T w3 = wr * drdz + ws * dsdz + wt * dtdz;

  const T s11 = dj * (u1 + u1);
  const T s22 = dj * (v2 + v2);
  const T s33 = dj * (w3 + w3);
  const T s12 = dj * (u2 + v1);
  const T s13 = dj * (u3 + w1);
  const T s23 = dj * (v3 + w2);

  shur = drdx * s11 + drdy * s12 + drdz * s13;
  shus = dsdx * s11 + dsdy * s12 + dsdz * s13;
  rut  = dtdx * s11 + dtdy * s12 + dtdz * s13;

  shvr = drdx * s12 + drdy * s22 + drdz * s23;
  shvs = dsdx * s12 + dsdy * s22 + dsdz * s23;
  rvt  = dtdx * s12 + dtdy * s22 + dtdz * s23;

  shwr = drdx * s13 + drdy * s23 + drdz * s33;
  shws = dsdx * s13 + dsdy * s23 + dsdz * s33;
  rwt  = dtdx * s13 + dtdy * s23 + dtdz * s33;
}

template< typename T, const int LX, const bool PAD >
__global__ void __launch_bounds__(LX*LX,3)
ax_helm_svv_factorized_full_kernel(T * __restrict__ au,
                                   T * __restrict__ av,
                                   T * __restrict__ aw,
                                   const T * __restrict__ u,
                                   const T * __restrict__ v,
                                   const T * __restrict__ w,
                                   const T * __restrict__ dx,
                                   const T * __restrict__ dy,
                                   const T * __restrict__ dz,
                                   const T * __restrict__ h1,
                                   const T * __restrict__ br,
                                   const T * __restrict__ bs,
                                   const T * __restrict__ bt,
                                   const T * __restrict__ h1_svv,
                                   const T * __restrict__ drdx,
                                   const T * __restrict__ drdy,
                                   const T * __restrict__ drdz,
                                   const T * __restrict__ dsdx,
                                   const T * __restrict__ dsdy,
                                   const T * __restrict__ dsdz,
                                   const T * __restrict__ dtdx,
                                   const T * __restrict__ dtdy,
                                   const T * __restrict__ dtdz,
                                   const T * __restrict__ jacinv,
                                   const T * __restrict__ weight3) {

  constexpr int P = PAD ? LX + 1 : LX;

  __shared__ T shdx[LX * P];
  __shared__ T shdy[LX * P];
  __shared__ T shdz[LX * P];
  __shared__ T shbr[LX * P];
  __shared__ T shbs[LX * P];
  __shared__ T shbt[LX * P];

  __shared__ T shu[LX * P];
  __shared__ T shv[LX * P];
  __shared__ T shw[LX * P];

  /* r-fluxes are only accessed along the fastest dimension */
  __shared__ T shur[LX * LX];
  __shared__ T shus[LX * P];
  __shared__ T shvr[LX * LX];
  __shared__ T shvs[LX * P];
  __shared__ T shwr[LX * LX];
  __shared__ T shws[LX * P];

  __shared__ T shur_svv[LX * LX];
  __shared__ T shus_svv[LX * P];
  __shared__ T shvr_svv[LX * LX];
  __shared__ T shvs_svv[LX * P];
  __shared__ T shwr_svv[LX * LX];
  __shared__ T shws_svv[LX * P];

  T ru[LX];
  T rv[LX];
  T rw[LX];

  T ruw[LX];
  T rvw[LX];
  T rww[LX];

  T rut, rvt, rwt;
  T rut_svv, rvt_svv, rwt_svv;

  const int e = blockIdx.x;
  const int j = threadIdx.y;
  const int i = threadIdx.x;
  const int ij = i + j*LX;
  const int ij_p = i + j*P;
  const int ele = e*LX*LX*LX;

  shdx[ij_p] = dx[ij];
  shdy[ij_p] = dy[ij];
  shdz[ij_p] = dz[ij];
  shbr[ij_p] = br[ij];
  shbs[ij_p] = bs[ij];
  shbt[ij_p] = bt[ij];

#pragma unroll
  for (int k = 0; k < LX; ++k) {
    ru[k] = u[ij + k*LX*LX + ele];
    ruw[k] = 0.0;

    rv[k] = v[ij + k*LX*LX + ele];
    rvw[k] = 0.0;

    rw[k] = w[ij + k*LX*LX + ele];
    rww[k] = 0.0;
  }

  __syncthreads();
#pragma unroll
  for (int k = 0; k < LX; ++k) {
    const int ijk = ij + k*LX*LX;
    const T drdx_local = drdx[ijk+ele];
    const T drdy_local = drdy[ijk+ele];
    const T drdz_local = drdz[ijk+ele];
    const T dsdx_local = dsdx[ijk+ele];
    const T dsdy_local = dsdy[ijk+ele];
    const T dsdz_local = dsdz[ijk+ele];
    const T dtdx_local = dtdx[ijk+ele];
    const T dtdy_local = dtdy[ijk+ele];
    const T dtdz_local = dtdz[ijk+ele];
    const T wj = weight3[ijk] * jacinv[ijk+ele];

    shu[ij_p] = ru[k];
    shv[ij_p] = rv[k];
    shw[ij_p] = rw[k];
    __syncthreads();

    /* Standard contribution */
    {
      T uttmp = 0.0;
      T vttmp = 0.0;
      T wttmp = 0.0;
      T urtmp = 0.0;
      T ustmp = 0.0;
      T vrtmp = 0.0;
      T vstmp = 0.0;
      T wrtmp = 0.0;
      T wstmp = 0.0;
#pragma unroll
      for (int l = 0; l < LX; l++) {
        uttmp += shdz[k+l*P] * ru[l];
        vttmp += shdz[k+l*P] * rv[l];
        wttmp += shdz[k+l*P] * rw[l];

        urtmp += shdx[i+l*P] * shu[l+j*P];
        ustmp += shdy[j+l*P] * shu[i+l*P];
        vrtmp += shdx[i+l*P] * shv[l+j*P];
        vstmp += shdy[j+l*P] * shv[i+l*P];
        wrtmp += shdx[i+l*P] * shw[l+j*P];
        wstmp += shdy[j+l*P] * shw[i+l*P];
      }

      ax_helm_svv_factorized_stress_flux<T>(
          urtmp, ustmp, uttmp, vrtmp, vstmp, vttmp, wrtmp, wstmp, wttmp,
          drdx_local, drdy_local, drdz_local,
          dsdx_local, dsdy_local, dsdz_local,
          dtdx_local, dtdy_local, dtdz_local,
          h1[ijk+ele] * wj,
          shur[ij], shus[ij_p], shvr[ij], shvs[ij_p], shwr[ij], shws[ij_p],
          rut, rvt, rwt);
    }

    /* Factorized SVV contribution */
    {
      T uttmp = 0.0;
      T vttmp = 0.0;
      T wttmp = 0.0;
      T urtmp = 0.0;
      T ustmp = 0.0;
      T vrtmp = 0.0;
      T vstmp = 0.0;
      T wrtmp = 0.0;
      T wstmp = 0.0;
#pragma unroll
      for (int l = 0; l < LX; l++) {
        uttmp += shbt[k+l*P] * ru[l];
        vttmp += shbt[k+l*P] * rv[l];
        wttmp += shbt[k+l*P] * rw[l];

        urtmp += shbr[i+l*P] * shu[l+j*P];
        ustmp += shbs[j+l*P] * shu[i+l*P];
        vrtmp += shbr[i+l*P] * shv[l+j*P];
        vstmp += shbs[j+l*P] * shv[i+l*P];
        wrtmp += shbr[i+l*P] * shw[l+j*P];
        wstmp += shbs[j+l*P] * shw[i+l*P];
      }

      ax_helm_svv_factorized_stress_flux<T>(
          urtmp, ustmp, uttmp, vrtmp, vstmp, vttmp, wrtmp, wstmp, wttmp,
          drdx_local, drdy_local, drdz_local,
          dsdx_local, dsdy_local, dsdz_local,
          dtdx_local, dtdy_local, dtdz_local,
          h1_svv[ijk+ele] * wj,
          shur_svv[ij], shus_svv[ij_p],
          shvr_svv[ij], shvs_svv[ij_p],
          shwr_svv[ij], shws_svv[ij_p],
          rut_svv, rvt_svv, rwt_svv);
    }

    __syncthreads();

    T uwijke = 0.0;
    T vwijke = 0.0;
    T wwijke = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++) {
      const T dxt = shdx[l+i*P];
      const T dyt = shdy[l+j*P];
      const T brt = shbr[l+i*P];
      const T bst = shbs[l+j*P];
      const T dzt = shdz[k+l*P];
      const T btt = shbt[k+l*P];

      uwijke += shur[l+j*LX] * dxt + shur_svv[l+j*LX] * brt
              + shus[i+l*P] * dyt + shus_svv[i+l*P] * bst;
      ruw[l] += rut * dzt + rut_svv * btt;

      vwijke += shvr[l+j*LX] * dxt + shvr_svv[l+j*LX] * brt
              + shvs[i+l*P] * dyt + shvs_svv[i+l*P] * bst;
      rvw[l] += rvt * dzt + rvt_svv * btt;

      wwijke += shwr[l+j*LX] * dxt + shwr_svv[l+j*LX] * brt
              + shws[i+l*P] * dyt + shws_svv[i+l*P] * bst;
      rww[l] += rwt * dzt + rwt_svv * btt;
    }
    ruw[k] += uwijke;
    rvw[k] += vwijke;
    rww[k] += wwijke;
  }

#pragma unroll
  for (int k = 0; k < LX; ++k) {
    au[ij + k*LX*LX + ele] = ruw[k];
    av[ij + k*LX*LX + ele] = rvw[k];
    aw[ij + k*LX*LX + ele] = rww[k];
  }
}

#endif // __MATH_AX_HELM_SVV_FACTORIZED_FULL_KERNEL_H__
