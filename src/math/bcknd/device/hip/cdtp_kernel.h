#ifndef __MATH_CDTP_KERNEL_H__
#define __MATH_CDTP_KERNEL_H__

#include "mfma_kernel.h"

/*
 Copyright (c) 2021-2026, The Neko Authors
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

#include "elem_block.h"

/**
 * Device kernel for \f$ D^T x \f$
 */

template< typename T, const int LX, const int CHUNKS >
__global__ void cdtp_kernel_1d(T * __restrict__ dtx,
                               const T * __restrict__ x,
                               const T * __restrict__ dr,
                               const T * __restrict__ ds,
                               const T * __restrict__ dt,
                               const T * __restrict__ dxt,
                               const T * __restrict__ dyt,
                               const T * __restrict__ dzt,
                               const T * __restrict__ w3) {

  __shared__ T shdxt[LX * LX];
  __shared__ T shdyt[LX * LX];
  __shared__ T shdzt[LX * LX];

  __shared__ T shtar[LX * LX * LX];
  __shared__ T shtas[LX * LX * LX];
  __shared__ T shtat[LX * LX * LX];

  const int e = blockIdx.x;
  const int iii = threadIdx.x;
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;

  if (iii < (LX * LX)) {
    shdxt[iii] = dxt[iii];
    shdyt[iii] = dyt[iii];
    shdzt[iii] = dzt[iii];
  }

  int l = iii;
  while(l < (LX * LX * LX)) {
    T wx = x[l + e * LX * LX * LX] * w3[l];

    shtar[l] = wx*dr[l + e * LX * LX * LX];
    shtas[l] = wx*ds[l + e * LX * LX * LX];
    shtat[l] = wx*dt[l + e * LX * LX * LX];

    l = l + CHUNKS;
  }

  __syncthreads();
  for (int n = 0; n < nchunks; n++) {
    const int ijk = iii + n * CHUNKS;
    const int jk = ijk / LX;
    const int i = ijk - jk * LX;
    const int k = jk / LX;
    const int j = jk - k * LX;
    if ( i < LX && j < LX && k < LX && ijk < LX*LX*LX) {
      T rtmp = 0.0;
      T stmp = 0.0;
      T ttmp = 0.0;
      for (int l = 0; l < LX; l++) {
	    rtmp += shdxt[i + l * LX] * shtar[l+j*LX+k*LX*LX];
	    stmp += shdyt[j + l * LX] * shtas[i+l*LX + k*LX*LX];
	    ttmp += shdzt[k + l * LX] * shtat[i + j*LX + l*LX*LX];
      }
      dtx[ijk + e * LX * LX * LX] = ( rtmp + stmp + ttmp );

    }
  }
}

template< typename T, const int LX, const int EB >
__global__ void NEKO_EB_BOUNDS(LX*LX*EB)
  cdtp_kernel_kstep(T * __restrict__ dtx,
                    const T * __restrict__ x,
                    const T * __restrict__ dr,
                    const T * __restrict__ ds,
                    const T * __restrict__ dt,
                    const T * __restrict__ dxt,
                    const T * __restrict__ dyt,
                    const T * __restrict__ dzt,
                    const T * __restrict__ w3,
                    const int nelv) {

  __shared__ T shdxt[LX * LX];
  __shared__ T shdyt[LX * LX];
  __shared__ T shdzt[LX * LX];

  __shared__ T shtar[EB * LX * LX];
  __shared__ T shtas[EB * LX * LX];

  static_assert(sizeof(shdxt) +
                sizeof(shdyt) +
                sizeof(shdzt) +
                sizeof(shtar) +
                sizeof(shtas)
                <= NEKO_EB_MAX_LDS,
                "kstep block exceeds the LDS budget");

  T rtar[LX];
  T rtas[LX];
  T rtat[LX];

  const int eb = (EB == 1) ? 0 : threadIdx.z;
  const int e_blk = blockIdx.x * EB + eb;
  /* Threads past the last element still have to reach the barriers in
     the k loop, so clamp their reads and drop their stores rather than
     returning early. At EB == 1 this all constant folds away */
  const bool active = (EB == 1) ? true : (e_blk < nelv);
  const int e = active ? e_blk : (nelv - 1);
  const int sh = eb * LX * LX;
  const int j = threadIdx.y;
  const int i = threadIdx.x;
  const int ij = i + j * LX;
  const int ele = e*LX*LX*LX;

  if (eb == 0) {
    shdxt[ij] = dxt[ij];
    shdyt[ij] = dyt[ij];
    shdzt[ij] = dzt[ij];
  }


#pragma unroll LX
  for (int k = 0; k < LX; ++k) {
    T wx = x[ij + k*LX*LX + ele] * w3[ij + k*LX*LX];

    rtar[k] = wx *dr[ij + k*LX*LX + ele];
    rtas[k] = wx *ds[ij + k*LX*LX + ele];
    rtat[k] = wx *dt[ij + k*LX*LX + ele];
  }

  __syncthreads();

#pragma unroll
  for (int k = 0; k < LX; ++k) {
    const int ijk = ij + k*LX*LX;
    T ttmp = 0.0;
    shtar[sh + ij] = rtar[k];
    shtas[sh + ij] = rtas[k];
#pragma unroll
    for (int l = 0; l < LX; l++) {
      ttmp += shdzt[k+l*LX] * rtat[l];
    }
    __syncthreads();

    T rtmp = 0.0;
    T stmp = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++) {
      rtmp += shdxt[i+l*LX] * shtar[sh + l+j*LX];
      stmp += shdyt[j+l*LX] * shtas[sh + i+l*LX];
    }

    if (active) {
      dtx[ijk + ele] = ( rtmp + stmp + ttmp );
    }

    __syncthreads();
  }
}



/**
 * Device kernel for D^T x on the AMD matrix cores
 *
 * The only one of these whose contractions come *after* its pointwise work
 * rather than before it, and the only one that accumulates: the three
 * weighted fields are formed first and then contracted into a single output.
 * That is axhelm's second half rather than its first.
 *
 * Because the operator is handed dxt, dyt and dzt already transposed, the
 * contraction index pattern is the ordinary forward one, so
 * mfma_contract_sel is used with TRANSPOSE = false against the staged
 * transposes rather than with TRANSPOSE = true against the forward matrices.
 * The output is zeroed with the staging and every axis accumulates into it,
 * matching the divergence half of ax_helm_mfma_elem().
 */

#if defined(__gfx90a__) || defined(__gfx942__)

template< typename T, const int LX, const int NWF >
__device__ void cdtp_mfma_elem(T * __restrict__ dtx,
                               const T * __restrict__ x,
                               const T * __restrict__ dr,
                               const T * __restrict__ ds,
                               const T * __restrict__ dt,
                               const T * __restrict__ dxt,
                               const T * __restrict__ dyt,
                               const T * __restrict__ dzt,
                               const T * __restrict__ w3,
                               const int nelv) {
  const int LX2 = LX * LX;
  const int LX3 = LX * LX * LX;

  /* NWF wavefronts per block, WPE of them cooperating on one element and the
     block covering EB elements, see the note in mfma_kernel.h. At LX = 4 the
     contraction offers one column group, so WPE is 1 and every wavefront gets
     an element of its own rather than idling. */
  enum { EB = NEKO_MFMA_EB_N(NWF, LX),
         WPE = NWF / EB };
  static_assert(WPE * EB == NWF,
                "wavefronts per block must split evenly over the elements");

  __shared__ T shdxt[LX * LX];
  __shared__ T shdyt[LX * LX];
  __shared__ T shdzt[LX * LX];
  __shared__ T shtar[EB * LX * LX * LX];   // the three weighted fields
  __shared__ T shtas[EB * LX * LX * LX];
  __shared__ T shtat[EB * LX * LX * LX];
  __shared__ T shout[EB * LX * LX * LX];   // the accumulated result

  static_assert(sizeof(shdxt) + sizeof(shdyt) + sizeof(shdzt) +
                sizeof(shtar) + sizeof(shtas) + sizeof(shtat) +
                sizeof(shout)
                <= NEKO_EB_MAX_LDS,
                "mfma block exceeds the shared memory budget");

  const int lane = threadIdx.x;          // 0..63 : lane within a wavefront
  const int wf   = threadIdx.y;          // 0..NWF-1 : which wavefront
  const int tid  = wf * 64 + lane;       // 0..NWF*64-1 : block-wide thread id
  const int nthr = NWF * 64;

  const int eb  = wf / WPE;              // which element this wavefront serves
  const int sub = wf % WPE;              // its rank among that element's waves
  const int gtid = sub * 64 + lane;      // thread id within the element group
  const int gnthr = WPE * 64;

  /* Threads past the last element still have to reach the block wide
     barriers, so clamp their reads and drop their stores rather than
     returning early. At EB == 1 the grid covers nelv exactly and this is
     constant folded away */
  const int e_blk = blockIdx.x * EB + eb;
  const bool active = (EB == 1) ? true : (e_blk < nelv);
  const int e = active ? e_blk : (nelv - 1);
  const int ele = e * LX3;
  const int sh = eb * LX3;

  /* Reference derivative matrices, already transposed, one copy per block */
  for (int p = tid; p < LX2; p += nthr) {
    shdxt[p] = dxt[p];
    shdyt[p] = dyt[p];
    shdzt[p] = dzt[p];
  }
  /* The three weighted fields, formed before anything is contracted, and the
     accumulator cleared alongside them */
  for (int p = gtid; p < LX3; p += gnthr) {
    const int gp = p + ele;
    const T wx = x[gp] * w3[p];

    shtar[sh + p] = wx * dr[gp];
    shtas[sh + p] = wx * ds[gp];
    shtat[sh + p] = wx * dt[gp];
    shout[sh + p] = 0.0;
  }

  __syncthreads();

  /* dtx = Dr^T tar + Ds^T tas + Dt^T tat, accumulated in shout */
  mfma_contract_sel<T, LX, 0, false, true, WPE>::run(shout + sh, shdxt,
                                                     shtar + sh, lane, sub);
  __syncthreads();
  mfma_contract_sel<T, LX, 1, false, true, WPE>::run(shout + sh, shdyt,
                                                     shtas + sh, lane, sub);
  __syncthreads();
  mfma_contract_sel<T, LX, 2, false, true, WPE>::run(shout + sh, shdzt,
                                                     shtat + sh, lane, sub);
  __syncthreads();

  if (active) {
    for (int p = gtid; p < LX3; p += gnthr)
      dtx[p + ele] = shout[sh + p];
  }
}

#endif // __gfx90a__ || __gfx942__

/*
 * Compile-time dispatch onto the MFMA element kernel. The launch macros in
 * opr_cdtp.hip are written for every LX the operator dispatches and for
 * whatever `real` is, so every combination has to compile; the ones the
 * strategy does not cover -- LX outside the supported range, a build without
 * a matrix-core arch -- resolve to this no-op. The autotuner never selects
 * the strategy for them, so the no-op is unreachable at runtime, see
 * mfma_lx_supported() and hip_have_mfma() in mfma_kernel.h.
 */
template< typename T, const int LX, const int NWF >
struct cdtp_mfma_dispatch {
  __device__ static void run(T *, const T *, const T *, const T *, const T *,
                             const T *, const T *, const T *, const T *,
                             const int) {}
};

#if defined(__gfx90a__) || defined(__gfx942__)

/* Keep in sync with mfma_lx_supported() in mfma_kernel.h */
#define NEKO_CDTP_MFMA_DISPATCH(TYPE, LXV)                                    \
  template< const int NWF >                                                   \
  struct cdtp_mfma_dispatch< TYPE, LXV, NWF > {                               \
    __device__ static void run(TYPE * dtx,                                    \
                               const TYPE * x,                                \
                               const TYPE * dr,                               \
                               const TYPE * ds,                               \
                               const TYPE * dt,                               \
                               const TYPE * dxt,                              \
                               const TYPE * dyt,                              \
                               const TYPE * dzt,                              \
                               const TYPE * w3,                               \
                               const int nelv) {                              \
      cdtp_mfma_elem< TYPE, LXV, NWF >(dtx, x, dr, ds, dt, dxt, dyt, dzt,     \
                                       w3, nelv);                             \
    }                                                                         \
  }

NEKO_CDTP_MFMA_DISPATCH(double, 4);
NEKO_CDTP_MFMA_DISPATCH(double, 5);
NEKO_CDTP_MFMA_DISPATCH(double, 6);
NEKO_CDTP_MFMA_DISPATCH(double, 7);
NEKO_CDTP_MFMA_DISPATCH(double, 8);
NEKO_CDTP_MFMA_DISPATCH(double, 9);
NEKO_CDTP_MFMA_DISPATCH(double, 10);
NEKO_CDTP_MFMA_DISPATCH(double, 11);
NEKO_CDTP_MFMA_DISPATCH(double, 12);

NEKO_CDTP_MFMA_DISPATCH(float, 4);
NEKO_CDTP_MFMA_DISPATCH(float, 5);
NEKO_CDTP_MFMA_DISPATCH(float, 6);
NEKO_CDTP_MFMA_DISPATCH(float, 7);
NEKO_CDTP_MFMA_DISPATCH(float, 8);
NEKO_CDTP_MFMA_DISPATCH(float, 9);
NEKO_CDTP_MFMA_DISPATCH(float, 10);
NEKO_CDTP_MFMA_DISPATCH(float, 11);
NEKO_CDTP_MFMA_DISPATCH(float, 12);

#endif // __gfx90a__ || __gfx942__

/*
 * Note the bare __launch_bounds__ rather than NEKO_EB_BOUNDS, matching
 * ax_helm_kernel_mfma: the kstep kernels ask for three waves per SIMD, and
 * the matrix core kernels were validated without that constraint.
 */
template< typename T, const int LX, const int NWF >
__global__ void __launch_bounds__(64 * NWF)
cdtp_kernel_mfma(T * __restrict__ dtx,
                 const T * __restrict__ x,
                 const T * __restrict__ dr,
                 const T * __restrict__ ds,
                 const T * __restrict__ dt,
                 const T * __restrict__ dxt,
                 const T * __restrict__ dyt,
                 const T * __restrict__ dzt,
                 const T * __restrict__ w3,
                 const int nelv) {

  cdtp_mfma_dispatch< T, LX, NWF >::run(dtx, x, dr, ds, dt, dxt, dyt, dzt,
                                        w3, nelv);
}

#endif // __MATH_CDTP_KERNEL_H__
