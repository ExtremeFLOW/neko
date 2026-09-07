#ifndef __MATH_OPGRAD_KERNEL_H__
#define __MATH_OPGRAD_KERNEL_H__

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
 * Device kernel for convective terms
 */

template< typename T, const int LX, const int CHUNKS >
__global__ void opgrad_kernel_1d(T * __restrict__ ux,
                                 T * __restrict__ uy,
                                 T * __restrict__ uz,
                                 const T * __restrict__ u,
                                 const T * __restrict__ dx,
                                 const T * __restrict__ dy,
                                 const T * __restrict__ dz,
                                 const T * __restrict__ drdx,
                                 const T * __restrict__ dsdx,
                                 const T * __restrict__ dtdx,
                                 const T * __restrict__ drdy,
                                 const T * __restrict__ dsdy,
                                 const T * __restrict__ dtdy,
                                 const T * __restrict__ drdz,
                                 const T * __restrict__ dsdz,
                                 const T * __restrict__ dtdz,
                                 const T * __restrict__ w3) {

  __shared__ T shu[LX * LX * LX];

  __shared__ T shdx[LX * LX];
  __shared__ T shdy[LX * LX];
  __shared__ T shdz[LX * LX];


  int i,j,k;

  const int e = blockIdx.x;
  const int iii = threadIdx.x;
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;

  if (iii < (LX * LX)) {
    shdx[iii] = dx[iii];
    shdy[iii] = dy[iii];
    shdz[iii] = dz[iii];
  }

  j = iii;
  while(j < (LX * LX * LX)) {
    shu[j] = u[j + e * LX * LX * LX];
    j = j + CHUNKS;
  }

  __syncthreads();

  for (int n = 0; n < nchunks; n++) {
    const int ijk = iii + n * CHUNKS;
    const int jk = ijk / LX;
    i = ijk - jk * LX;
    k = jk / LX;
    j = jk - k * LX;
    if ( i < LX && j < LX && k < LX ) {
      T rtmp = 0.0;
      T stmp = 0.0;
      T ttmp = 0.0;
      for (int l = 0; l < LX; l++) {
	rtmp += shdx[i + l * LX] * shu[l + j * LX + k * LX * LX];
	stmp += shdy[j + l * LX] * shu[i + l * LX + k * LX * LX];
	ttmp += shdz[k + l * LX] * shu[i + j * LX + l * LX * LX];
      }

      ux[ijk + e * LX * LX * LX] = w3[ijk]
	* (drdx[ijk + e * LX * LX * LX] * rtmp
	   + dsdx[ijk + e * LX * LX * LX] * stmp
	   + dtdx[ijk + e * LX * LX * LX] * ttmp);

      uy[ijk + e * LX * LX * LX] = w3[ijk]
	* (drdy[ijk + e * LX * LX * LX] * rtmp
	   + dsdy[ijk + e * LX * LX * LX] * stmp
	   + dtdy[ijk + e * LX * LX * LX] * ttmp);

      uz[ijk + e * LX * LX * LX] = w3[ijk]
	* (drdz[ijk + e * LX * LX * LX] * rtmp
	   + dsdz[ijk + e * LX * LX * LX] * stmp
	   + dtdz[ijk + e * LX * LX * LX] * ttmp);

    }
  }

}

template< typename T, const int LX, const int EB >
__global__ void NEKO_EB_BOUNDS(LX*LX*EB)
opgrad_kernel_kstep(T * __restrict__ ux,
                    T * __restrict__ uy,
                    T * __restrict__ uz,
                    const T * __restrict__ u,
                    const T * __restrict__ dx,
                    const T * __restrict__ dy,
                    const T * __restrict__ dz,
                    const T * __restrict__ drdx,
                    const T * __restrict__ dsdx,
                    const T * __restrict__ dtdx,
                    const T * __restrict__ drdy,
                    const T * __restrict__ dsdy,
                    const T * __restrict__ dtdy,
                    const T * __restrict__ drdz,
                    const T * __restrict__ dsdz,
                    const T * __restrict__ dtdz,
                    const T * __restrict__ w3,
                    const int nelv) {

  /* One slice per element in the block */
  __shared__ T shu[EB * LX * LX];

  /* Element independent, one copy per block */
  __shared__ T shdx[LX * LX];
  __shared__ T shdy[LX * LX];
  __shared__ T shdz[LX * LX];

  static_assert(sizeof(shu) + sizeof(shdx) + sizeof(shdy) + sizeof(shdz)
                <= NEKO_EB_MAX_LDS,
                "kstep block exceeds the LDS budget");

  /* Threads past the last element still have to reach the barriers in the k
     loop, so clamp their reads and drop their stores rather than returning
     early. At EB == 1 all of this is constant folded away */
  const int eb = (EB == 1) ? 0 : threadIdx.z;
  const int e_blk = blockIdx.x * EB + eb;
  const bool active = (EB == 1) ? true : (e_blk < nelv);
  const int e = active ? e_blk : (nelv - 1);
  const int j = threadIdx.y;
  const int i = threadIdx.x;
  const int ij = i + j * LX;
  const int sh = eb * LX * LX;
  const int ele = e*LX*LX*LX;

  if (eb == 0) {
    shdx[ij] = dx[ij];
    shdy[ij] = dy[ij];
    shdz[ij] = dz[ij];
  }

  T ru[LX];

#pragma unroll LX
  for (int k = 0; k < LX; ++k) {
    ru[k] = u[ij + k*LX*LX + ele];
  }

  __syncthreads();

  #pragma unroll
  for (int k = 0; k < LX; ++k) {
    const int ijk = ij + k*LX*LX;
    const T W3 = w3[ijk];
    T ttmp = 0.0;
    shu[sh + ij] = ru[k];
#pragma unroll
    for (int l = 0; l < LX; l++) {
      ttmp += shdz[k+l*LX] * ru[l];
    }
    __syncthreads();

    T rtmp = 0.0;
    T stmp = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++) {
      rtmp += shdx[i+l*LX] * shu[sh + l+j*LX];
      stmp += shdy[j+l*LX] * shu[sh + i+l*LX];
    }

    if (active) {
      ux[ijk + ele] = W3 * (drdx[ijk + ele] * rtmp
                            + dsdx[ijk + ele] * stmp
                            + dtdx[ijk + ele] * ttmp);

      uy[ijk + ele] = W3 * (drdy[ijk + ele] * rtmp
                            + dsdy[ijk + ele] * stmp
                            + dtdy[ijk + ele] * ttmp);

      uz[ijk + ele] = W3 * (drdz[ijk + ele] * rtmp
                            + dsdz[ijk + ele] * stmp
                            + dtdz[ijk + ele] * ttmp);
    }
    __syncthreads();
  }
}



/**
 * Device kernel for the weak gradient on the AMD matrix cores
 *
 * The same staging and the same three contractions as dudxyz_mfma_elem(),
 * with a wider tail: three outputs formed pointwise from nine metric cubes
 * instead of one from three. See mfma_kernel.h for the tiles, the lane
 * layouts, the wavefronts-per-element scheme and the arch and LX bounds.
 *
 * w3 is a quadrature weight shared by every element, so it is indexed by the
 * point within the element rather than by the global offset. Unlike the CUDA
 * counterpart there is no padded cube here, so that index is simply p.
 */

#if defined(__gfx90a__) || defined(__gfx942__)

template< typename T, const int LX, const int NWF >
__device__ void opgrad_mfma_elem(T * __restrict__ ux,
                                 T * __restrict__ uy,
                                 T * __restrict__ uz,
                                 const T * __restrict__ u,
                                 const T * __restrict__ dx,
                                 const T * __restrict__ dy,
                                 const T * __restrict__ dz,
                                 const T * __restrict__ drdx,
                                 const T * __restrict__ dsdx,
                                 const T * __restrict__ dtdx,
                                 const T * __restrict__ drdy,
                                 const T * __restrict__ dsdy,
                                 const T * __restrict__ dtdy,
                                 const T * __restrict__ drdz,
                                 const T * __restrict__ dsdz,
                                 const T * __restrict__ dtdz,
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

  __shared__ T shdx[LX * LX];
  __shared__ T shdy[LX * LX];
  __shared__ T shdz[LX * LX];
  __shared__ T shu[EB * LX * LX * LX];   // the staged field
  __shared__ T shr[EB * LX * LX * LX];   // d/dr
  __shared__ T shs[EB * LX * LX * LX];   // d/ds
  __shared__ T sht[EB * LX * LX * LX];   // d/dt

  static_assert(sizeof(shdx) + sizeof(shdy) + sizeof(shdz) +
                sizeof(shu) + sizeof(shr) + sizeof(shs) + sizeof(sht)
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

  /* Reference derivative matrices, one copy shared by every element */
  for (int p = tid; p < LX2; p += nthr) {
    shdx[p] = dx[p];
    shdy[p] = dy[p];
    shdz[p] = dz[p];
  }
  /* Element-local field, staged by the wavefronts that own it */
  for (int p = gtid; p < LX3; p += gnthr)
    shu[sh + p] = u[p + ele];

  __syncthreads();

  /* Reference space derivatives ur, us, ut, striped across the cooperating
     wavefronts by mfma_contract_sel */
  mfma_contract_sel<T, LX, 0, false, false, WPE>::run(shr + sh, shdx,
                                                      shu + sh, lane, sub);
  mfma_contract_sel<T, LX, 1, false, false, WPE>::run(shs + sh, shdy,
                                                      shu + sh, lane, sub);
  mfma_contract_sel<T, LX, 2, false, false, WPE>::run(sht + sh, shdz,
                                                      shu + sh, lane, sub);

  __syncthreads();

  if (active) {
    for (int p = gtid; p < LX3; p += gnthr) {
      const int gp = p + ele;
      const T W3 = w3[p];
      const T rr = shr[sh + p], ss = shs[sh + p], tt = sht[sh + p];

      ux[gp] = W3 * (drdx[gp] * rr + dsdx[gp] * ss + dtdx[gp] * tt);
      uy[gp] = W3 * (drdy[gp] * rr + dsdy[gp] * ss + dtdy[gp] * tt);
      uz[gp] = W3 * (drdz[gp] * rr + dsdz[gp] * ss + dtdz[gp] * tt);
    }
  }
}

#endif // __gfx90a__ || __gfx942__

/*
 * Compile-time dispatch onto the MFMA element kernel. The launch macros in
 * opr_opgrad.hip are written for every LX the operator dispatches and for
 * whatever `real` is, so every combination has to compile; the ones the
 * strategy does not cover -- LX outside the supported range, a build without
 * a matrix-core arch -- resolve to this no-op. The autotuner never selects
 * the strategy for them, so the no-op is unreachable at runtime, see
 * mfma_lx_supported() and hip_have_mfma() in mfma_kernel.h.
 */
template< typename T, const int LX, const int NWF >
struct opgrad_mfma_dispatch {
  __device__ static void run(T *, T *, T *, const T *, const T *, const T *,
                             const T *, const T *, const T *, const T *,
                             const T *, const T *, const T *, const T *,
                             const T *, const T *, const T *, const int) {}
};

#if defined(__gfx90a__) || defined(__gfx942__)

/* Keep in sync with mfma_lx_supported() in mfma_kernel.h */
#define NEKO_OPGRAD_MFMA_DISPATCH(TYPE, LXV)                                  \
  template< const int NWF >                                                   \
  struct opgrad_mfma_dispatch< TYPE, LXV, NWF > {                             \
    __device__ static void run(TYPE * ux,                                     \
                               TYPE * uy,                                     \
                               TYPE * uz,                                     \
                               const TYPE * u,                                \
                               const TYPE * dx,                               \
                               const TYPE * dy,                               \
                               const TYPE * dz,                               \
                               const TYPE * drdx,                             \
                               const TYPE * dsdx,                             \
                               const TYPE * dtdx,                             \
                               const TYPE * drdy,                             \
                               const TYPE * dsdy,                             \
                               const TYPE * dtdy,                             \
                               const TYPE * drdz,                             \
                               const TYPE * dsdz,                             \
                               const TYPE * dtdz,                             \
                               const TYPE * w3,                               \
                               const int nelv) {                              \
      opgrad_mfma_elem< TYPE, LXV, NWF >(ux, uy, uz, u, dx, dy, dz, drdx,     \
                                         dsdx, dtdx, drdy, dsdy, dtdy,        \
                                         drdz, dsdz, dtdz, w3, nelv);         \
    }                                                                         \
  }

NEKO_OPGRAD_MFMA_DISPATCH(double, 4);
NEKO_OPGRAD_MFMA_DISPATCH(double, 5);
NEKO_OPGRAD_MFMA_DISPATCH(double, 6);
NEKO_OPGRAD_MFMA_DISPATCH(double, 7);
NEKO_OPGRAD_MFMA_DISPATCH(double, 8);
NEKO_OPGRAD_MFMA_DISPATCH(double, 9);
NEKO_OPGRAD_MFMA_DISPATCH(double, 10);
NEKO_OPGRAD_MFMA_DISPATCH(double, 11);
NEKO_OPGRAD_MFMA_DISPATCH(double, 12);

NEKO_OPGRAD_MFMA_DISPATCH(float, 4);
NEKO_OPGRAD_MFMA_DISPATCH(float, 5);
NEKO_OPGRAD_MFMA_DISPATCH(float, 6);
NEKO_OPGRAD_MFMA_DISPATCH(float, 7);
NEKO_OPGRAD_MFMA_DISPATCH(float, 8);
NEKO_OPGRAD_MFMA_DISPATCH(float, 9);
NEKO_OPGRAD_MFMA_DISPATCH(float, 10);
NEKO_OPGRAD_MFMA_DISPATCH(float, 11);
NEKO_OPGRAD_MFMA_DISPATCH(float, 12);

#endif // __gfx90a__ || __gfx942__

/*
 * Note the bare __launch_bounds__ rather than NEKO_EB_BOUNDS, matching
 * ax_helm_kernel_mfma: the kstep kernels ask for three waves per SIMD, and
 * the matrix core kernels were validated without that constraint.
 */
template< typename T, const int LX, const int NWF >
__global__ void __launch_bounds__(64 * NWF)
opgrad_kernel_mfma(T * __restrict__ ux,
                   T * __restrict__ uy,
                   T * __restrict__ uz,
                   const T * __restrict__ u,
                   const T * __restrict__ dx,
                   const T * __restrict__ dy,
                   const T * __restrict__ dz,
                   const T * __restrict__ drdx,
                   const T * __restrict__ dsdx,
                   const T * __restrict__ dtdx,
                   const T * __restrict__ drdy,
                   const T * __restrict__ dsdy,
                   const T * __restrict__ dtdy,
                   const T * __restrict__ drdz,
                   const T * __restrict__ dsdz,
                   const T * __restrict__ dtdz,
                   const T * __restrict__ w3,
                   const int nelv) {

  opgrad_mfma_dispatch< T, LX, NWF >::run(ux, uy, uz, u, dx, dy, dz, drdx,
                                          dsdx, dtdx, drdy, dsdy, dtdy, drdz,
                                          dsdz, dtdz, w3, nelv);
}

#endif // __MATH_OPGRAD_KERNEL_H__
