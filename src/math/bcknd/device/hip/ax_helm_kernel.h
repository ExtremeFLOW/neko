#ifndef __MATH_AX_HELM_KERNEL_H__
#define __MATH_AX_HELM_KERNEL_H__

/*
 Copyright (c) 2021-2024, The Neko Authors
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
#include "mfma_kernel.h"

/**
 * Device kernels for Ax helm
 */

template< typename T, const int LX, const int CHUNKS >
__global__ void ax_helm_kernel_1d(T * __restrict__ w,
                                  const T * __restrict__ u,
                                  const T * __restrict__ dx,
                                  const T * __restrict__ dy,
                                  const T * __restrict__ dz,
                                  const T * __restrict__ dxt,
                                  const T * __restrict__ dyt,
                                  const T * __restrict__ dzt,
                                  const T * __restrict__ h1,
                                  const T * __restrict__ g11,
                                  const T * __restrict__ g22,
                                  const T * __restrict__ g33,
                                  const T * __restrict__ g12,
                                  const T * __restrict__ g13,
                                  const T * __restrict__ g23) {

  __shared__ T shdx[LX*LX];
  __shared__ T shdy[LX*LX];
  __shared__ T shdzt[LX*LX];

  __shared__ T shdxt[LX*LX];
  __shared__ T shdyt[LX*LX];
  __shared__ T shdz[LX*LX];

  __shared__ T shu[LX*LX*LX];
  __shared__ T shur[LX*LX*LX];
  __shared__ T shus[LX*LX*LX];
  __shared__ T shut[LX*LX*LX];

  const int e = blockIdx.x;
  const int iii = threadIdx.x;
  const int nchunks = (LX * LX * LX - 1)/CHUNKS + 1;

  if (iii<LX*LX) {
    shdx[iii] = dx[iii];
    shdy[iii] = dy[iii];
    shdz[iii] = dz[iii];
  }

  {
    int i = iii;
    while (i < LX * LX * LX){
      shu[i] = u[i+e*LX*LX*LX];
      i = i + CHUNKS;
    }
  }

  __syncthreads();

  if (iii<LX*LX){
    shdxt[iii] = dxt[iii];
    shdyt[iii] = dyt[iii];
    shdzt[iii] = dzt[iii];
  }

  for (int n=0; n<nchunks; n++){
    const int ijk = iii+n*CHUNKS;
    const int jk = ijk/LX;
    const int i = ijk-jk*LX;
    const int k = jk/LX;
    const int j = jk-k*LX;
    if (i<LX && j<LX && k<LX && ijk < LX*LX*LX){
      const T G00 = g11[ijk+e*LX*LX*LX];
      const T G11 = g22[ijk+e*LX*LX*LX];
      const T G22 = g33[ijk+e*LX*LX*LX];
      const T G01 = g12[ijk+e*LX*LX*LX];
      const T G02 = g13[ijk+e*LX*LX*LX];
      const T G12 = g23[ijk+e*LX*LX*LX];
      const T H1 = h1[ijk+e*LX*LX*LX];
      T rtmp = 0.0;
      T stmp = 0.0;
      T ttmp = 0.0;
#pragma unroll
      for (int l = 0; l<LX; l++){
        rtmp = rtmp + shdx[i+l*LX] * shu[l+j*LX+k*LX*LX];
        stmp = stmp + shdy[j+l*LX] * shu[i+l*LX+k*LX*LX];
        ttmp = ttmp + shdz[k+l*LX] * shu[i+j*LX+l*LX*LX];
      }
      shur[ijk] = H1 * (G00 * rtmp + G01 * stmp + G02 * ttmp);
      shus[ijk] = H1 * (G01 * rtmp + G11 * stmp + G12 * ttmp);
      shut[ijk] = H1 * (G02 * rtmp + G12 * stmp + G22 * ttmp);
    }
  }

  __syncthreads();

  for (int n=0; n<nchunks; n++){
    const int ijk = iii+n*CHUNKS;
    const int jk = ijk/LX;
    const int k = jk/LX;
    const int j = jk-k*LX;
    const int i = ijk-jk*LX;
    if (i<LX && j<LX && k<LX && ijk <LX*LX*LX){
      T wijke = 0.0;
#pragma unroll
      for (int l = 0; l<LX; l++){
        wijke = wijke
              + shdxt[i+l*LX] * shur[l+j*LX+k*LX*LX]
              + shdyt[j+l*LX] * shus[i+l*LX+k*LX*LX]
              + shdzt[k+l*LX] * shut[i+j*LX+l*LX*LX];
      }
      w[ijk+e*LX*LX*LX] = wijke;
    }
  }
}

template< typename T, const int LX, const int EB >
__global__ void NEKO_EB_BOUNDS(LX*LX*EB)
ax_helm_kernel_kstep(T * __restrict__ w,
                     const T * __restrict__ u,
                     const T * __restrict__ dx,
                     const T * __restrict__ dy,
                     const T * __restrict__ dz,
                     const T * __restrict__ h1,
                     const T * __restrict__ g11,
                     const T * __restrict__ g22,
                     const T * __restrict__ g33,
                     const T * __restrict__ g12,
                     const T * __restrict__ g13,
                     const T * __restrict__ g23,
                     const int nelv) {

  /* Element independent, one copy per block */
  __shared__ T shdx[LX * LX];
  __shared__ T shdy[LX * LX];
  __shared__ T shdz[LX * LX];

  /* One slice per element in the block */
  __shared__ T shu[EB * LX * LX];
  __shared__ T shur[EB * LX * LX];
  __shared__ T shus[EB * LX * LX];

  static_assert(sizeof(shdx) +
                sizeof(shdy) +
                sizeof(shdz) +
                sizeof(shu) +
                sizeof(shur) +
                sizeof(shus)
                <= NEKO_EB_MAX_LDS,
                "kstep block exceeds the LDS budget");

  T ru[LX];
  T rw[LX];
  T rut;

  /* Threads past the last element still have to reach the block wide
     barriers below, so clamp their reads and drop their stores rather
     than returning early. At EB == 1 the grid covers nelv exactly, so all
     of this bookkeeping is constant folded and the kernel is exactly what
     it was before blocking */
  const int eb = (EB == 1) ? 0 : threadIdx.z;
  const int e_blk = blockIdx.x * EB + eb;
  const bool active = (EB == 1) ? true : (e_blk < nelv);
  const int e = active ? e_blk : (nelv - 1);
  const int j = threadIdx.y;
  const int i = threadIdx.x;
  const int ij = i + j*LX;
  const int sh = eb*LX*LX;
  const int ele = e*LX*LX*LX;

  if (eb == 0) {
    shdx[ij] = dx[ij];
    shdy[ij] = dy[ij];
    shdz[ij] = dz[ij];
  }

#pragma unroll
  for(int k = 0; k < LX; ++k){
    ru[k] = u[ij + k*LX*LX + ele];
    rw[k] = 0.0;
  }


  __syncthreads();
#pragma unroll
  for (int k = 0; k < LX; ++k){
    const int ijk = ij + k*LX*LX;
    const T G00 = g11[ijk+ele];
    const T G11 = g22[ijk+ele];
    const T G22 = g33[ijk+ele];
    const T G01 = g12[ijk+ele];
    const T G02 = g13[ijk+ele];
    const T G12 = g23[ijk+ele];
    const T H1  = h1[ijk+ele];
    T ttmp = 0.0;
    shu[sh + ij] = ru[k];
#pragma unroll
    for (int l = 0; l < LX; l++){
      ttmp += shdz[k+l*LX] * ru[l];
    }
    __syncthreads();

    T rtmp = 0.0;
    T stmp = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++){
      rtmp += shdx[i+l*LX] * shu[sh + l+j*LX];
      stmp += shdy[j+l*LX] * shu[sh + i+l*LX];
    }
    shur[sh + ij] = H1
                  * (G00 * rtmp
                     + G01 * stmp
                     + G02 * ttmp);
    shus[sh + ij] = H1
                  * (G01 * rtmp
                     + G11 * stmp
                     + G12 * ttmp);
    rut           = H1
                  * (G02 * rtmp
                     + G12 * stmp
                     + G22 * ttmp);

    __syncthreads();

    T wijke = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++){
      wijke += shur[sh + l+j*LX] * shdx[l+i*LX];
      rw[l] += rut * shdz[k+l*LX];
      wijke += shus[sh + i+l*LX] * shdy[l + j*LX];
    }
    rw[k] += wijke;
  }
  if (active) {
#pragma unroll
    for (int k = 0; k < LX; ++k){
      w[ij + k*LX*LX + ele] = rw[k];
    }
  }
}

/**
 * Device kernel for axhelm with padding in shared memory to
 * remove bank conflicts when LX is a power of 2
 */

template< typename T, const int LX, const int EB >
__global__ void NEKO_EB_BOUNDS(LX*LX*EB)
ax_helm_kernel_kstep_padded(T * __restrict__ w,
                            const T * __restrict__ u,
                            const T * __restrict__ dx,
                            const T * __restrict__ dy,
                            const T * __restrict__ dz,
                            const T * __restrict__ h1,
                            const T * __restrict__ g11,
                            const T * __restrict__ g22,
                            const T * __restrict__ g33,
                            const T * __restrict__ g12,
                            const T * __restrict__ g13,
                            const T * __restrict__ g23,
                            const int nelv) {

  /* Element independent, one copy per block */
  __shared__ T shdx[LX * (LX+1)];
  __shared__ T shdy[LX * (LX+1)];
  __shared__ T shdz[LX * (LX+1)];

  /* One slice per element in the block */
  __shared__ T shu[EB * LX * (LX+1)];
  __shared__ T shur[EB * LX * LX];  // only accessed using fastest dimension
  __shared__ T shus[EB * LX * (LX+1)];

  static_assert(sizeof(shdx) +
                sizeof(shdy) +
                sizeof(shdz) +
                sizeof(shu) +
                sizeof(shur) +
                sizeof(shus)
                <= NEKO_EB_MAX_LDS,
                "kstep block exceeds the LDS budget");

  T ru[LX];
  T rw[LX];
  T rut;

  /* At EB == 1 the grid covers nelv exactly, so the blocking bookkeeping is
     constant folded away and the kernel is exactly what it was before */
  const int eb = (EB == 1) ? 0 : threadIdx.z;
  const int e_blk = blockIdx.x * EB + eb;
  const bool active = (EB == 1) ? true : (e_blk < nelv);
  const int e = active ? e_blk : (nelv - 1);
  const int j = threadIdx.y;
  const int i = threadIdx.x;
  const int ij = i + j*LX;
  const int ij_p = i + j*(LX+1);
  const int sh = eb*LX*LX;
  const int sh_p = eb*LX*(LX+1);
  const int ele = e*LX*LX*LX;

  if (eb == 0) {
    shdx[ij_p] = dx[ij];
    shdy[ij_p] = dy[ij];
    shdz[ij_p] = dz[ij];
  }

#pragma unroll
  for(int k = 0; k < LX; ++k){
    ru[k] = u[ij + k*LX*LX + ele];
    rw[k] = 0.0;
  }


  __syncthreads();
#pragma unroll
  for (int k = 0; k < LX; ++k){
    const int ijk = ij + k*LX*LX;
    const T G00 = g11[ijk+ele];
    const T G11 = g22[ijk+ele];
    const T G22 = g33[ijk+ele];
    const T G01 = g12[ijk+ele];
    const T G02 = g13[ijk+ele];
    const T G12 = g23[ijk+ele];
    const T H1  = h1[ijk+ele];
    T ttmp = 0.0;
    shu[sh_p + ij_p] = ru[k];
#pragma unroll
    for (int l = 0; l < LX; l++){
      ttmp += shdz[k+l*(LX+1)] * ru[l];
    }
    __syncthreads();

    T rtmp = 0.0;
    T stmp = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++){
      rtmp += shdx[i+l*(LX+1)] * shu[sh_p + l+j*(LX+1)];
      stmp += shdy[j+l*(LX+1)] * shu[sh_p + i+l*(LX+1)];
    }
    shur[sh + ij]      = H1
                       * (G00 * rtmp
                          + G01 * stmp
                          + G02 * ttmp);
    shus[sh_p + ij_p]  = H1
                       * (G01 * rtmp
                          + G11 * stmp
                          + G12 * ttmp);
    rut                = H1
                       * (G02 * rtmp
                          + G12 * stmp
                          + G22 * ttmp);

    __syncthreads();

    T wijke = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++){
      wijke += shur[sh + l+j*LX] * shdx[l+i*(LX+1)];
      rw[l] += rut * shdz[k+l*(LX+1)];
      wijke += shus[sh_p + i+l*(LX+1)] * shdy[l + j*(LX+1)];
    }
    rw[k] += wijke;
  }
  if (active) {
#pragma unroll
    for (int k = 0; k < LX; ++k){
      w[ij + k*LX*LX + ele] = rw[k];
    }
  }
}


/**
 * Matrix-core (MFMA) device kernel for axhelm.
 *
 * Additional autotuner strategy that maps the six spectral-element tensor
 * contractions (3 gradient + 3 divergence) onto the AMD matrix cores via the
 * precision-dispatched mfma_contract_sel primitive -- double precision uses the
 * batched 4x4x4 matrix core (full M-utilisation), single precision the 16x16x4
 * tile (see mfma_kernel.h for the precision traits, lane layout and supported
 * (precision, LX) set).  AX_HELM_MFMA_NWF cooperating wavefronts process one
 * element, sharing the LX^3 field staged in LDS; decoupling the block shape
 * from LX lets one implementation cover 4 <= LX <= 12 rather than only LX = 8.
 * Unsupported (T, LX) instantiate to a no-op; the autotuner only launches this
 * kernel for the supported set (see mfma_lx_supported() and hip_have_mfma() in
 * mfma_kernel.h), so the no-op is never reached at runtime.
 */
#if defined(__gfx90a__) || defined(__gfx942__)

/* Matrix-core axhelm for one element of order LX-1 (T = float or double).
 *
 * NWF cooperating wavefronts (blockDim = (64, NWF, 1)) process one element,
 * sharing a single staged cube in LDS.  The matrix-core column tiles are
 * striped across the wavefronts (wf = threadIdx.y) via mfma_contract_sel, while
 * the staging, geometry and write-back passes parallelise over all NWF*64
 * threads (tid).  The three accumulating divergence contractions stay separated
 * by __syncthreads and each wavefront owns disjoint output columns, so the
 * accumulation is race-free.  NWF = 1 reproduces the single-wavefront kernel.
 *
 * mfma_contract_sel routes double precision through the batched 4x4x4 matrix
 * core (full M-utilisation) and single precision through the 16x16x4 tile. */
template< typename T, const int LX, const int NWF >
__device__ void ax_helm_mfma_elem(T * __restrict__ w,
                                  const T * __restrict__ u,
                                  const T * __restrict__ dx,
                                  const T * __restrict__ dy,
                                  const T * __restrict__ dz,
                                  const T * __restrict__ h1,
                                  const T * __restrict__ g11,
                                  const T * __restrict__ g22,
                                  const T * __restrict__ g33,
                                  const T * __restrict__ g12,
                                  const T * __restrict__ g13,
                                  const T * __restrict__ g23,
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
  __shared__ T shu[EB * LX * LX * LX];   // input u, later reused as output w
  __shared__ T shr[EB * LX * LX * LX];   // d/dr -> Sr
  __shared__ T shs[EB * LX * LX * LX];   // d/ds -> Ss
  __shared__ T sht[EB * LX * LX * LX];   // d/dt -> St

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

  /* Gradient: ur, us, ut in canonical i + LX*j + LX*LX*k layout. */
  mfma_contract_sel<T, LX, 0, false, false, WPE>::run(shr + sh, shdx,
                                                      shu + sh, lane, sub);
  mfma_contract_sel<T, LX, 1, false, false, WPE>::run(shs + sh, shdy,
                                                      shu + sh, lane, sub);
  mfma_contract_sel<T, LX, 2, false, false, WPE>::run(sht + sh, shdz,
                                                      shu + sh, lane, sub);

  __syncthreads();

  /* Geometry (pointwise): (ur,us,ut) -> (Sr,Ss,St), reusing shr/shs/sht. */
  for (int p = gtid; p < LX3; p += gnthr) {
    const int gp = p + ele;
    const T G00 = g11[gp], G11 = g22[gp], G22 = g33[gp];
    const T G01 = g12[gp], G02 = g13[gp], G12 = g23[gp];
    const T H1  = h1[gp];
    const T rr = shr[sh + p], ss = shs[sh + p], tt = sht[sh + p];
    shr[sh + p] = H1 * (G00 * rr + G01 * ss + G02 * tt);
    shs[sh + p] = H1 * (G01 * rr + G11 * ss + G12 * tt);
    sht[sh + p] = H1 * (G02 * rr + G12 * ss + G22 * tt);
  }

  __syncthreads();

  /* Divergence: w = Dr^T Sr + Ds^T Ss + Dt^T St, accumulated in shu (= w). */
  for (int p = gtid; p < LX3; p += gnthr)
    shu[sh + p] = 0.0;
  __syncthreads();

  mfma_contract_sel<T, LX, 0, true, true, WPE>::run(shu + sh, shdx,
                                                    shr + sh, lane, sub);
  __syncthreads();
  mfma_contract_sel<T, LX, 1, true, true, WPE>::run(shu + sh, shdy,
                                                    shs + sh, lane, sub);
  __syncthreads();
  mfma_contract_sel<T, LX, 2, true, true, WPE>::run(shu + sh, shdz,
                                                    sht + sh, lane, sub);
  __syncthreads();

  if (active) {
    for (int p = gtid; p < LX3; p += gnthr)
      w[p + ele] = shu[sh + p];
  }
}
#endif // __gfx90a__ || __gfx942__

/*
 * Compile-time dispatch onto the MFMA element kernel. The launch macros in
 * ax_helm.hip are written for every LX the operator dispatches and for
 * whatever `real` is, so every combination has to compile; the ones the
 * strategy does not cover -- LX outside the supported range, a build without
 * a matrix-core arch -- resolve to this no-op. The autotuner never selects
 * the strategy for them, so the no-op is unreachable at runtime, see
 * mfma_lx_supported() and hip_have_mfma() in mfma_kernel.h.
 */
template< typename T, const int LX, const int NWF >
struct ax_helm_mfma_dispatch {
  __device__ static void run(T *, const T *, const T *, const T *, const T *,
                             const T *, const T *, const T *, const T *,
                             const T *, const T *, const T *, const int) {}
};

#if defined(__gfx90a__) || defined(__gfx942__)

/* Keep in sync with mfma_lx_supported() in mfma_kernel.h */
#define NEKO_AX_HELM_MFMA_DISPATCH(TYPE, LXV)                                  \
  template< const int NWF >                                                    \
  struct ax_helm_mfma_dispatch< TYPE, LXV, NWF > {                             \
    __device__ static void run(TYPE *w, const TYPE *u,                         \
                               const TYPE *dx, const TYPE *dy,                 \
                               const TYPE *dz, const TYPE *h1,                 \
                               const TYPE *g11, const TYPE *g22,               \
                               const TYPE *g33, const TYPE *g12,               \
                               const TYPE *g13, const TYPE *g23,               \
                               const int nelv) {                               \
      ax_helm_mfma_elem< TYPE, LXV, NWF >(w, u, dx, dy, dz, h1,                \
                                          g11, g22, g33, g12, g13, g23,        \
                                          nelv);                               \
    }                                                                          \
  }

NEKO_AX_HELM_MFMA_DISPATCH(double, 4);
NEKO_AX_HELM_MFMA_DISPATCH(double, 5);
NEKO_AX_HELM_MFMA_DISPATCH(double, 6);
NEKO_AX_HELM_MFMA_DISPATCH(double, 7);
NEKO_AX_HELM_MFMA_DISPATCH(double, 8);
NEKO_AX_HELM_MFMA_DISPATCH(double, 9);
NEKO_AX_HELM_MFMA_DISPATCH(double, 10);
NEKO_AX_HELM_MFMA_DISPATCH(double, 11);
NEKO_AX_HELM_MFMA_DISPATCH(double, 12);
NEKO_AX_HELM_MFMA_DISPATCH(float, 4);
NEKO_AX_HELM_MFMA_DISPATCH(float, 5);
NEKO_AX_HELM_MFMA_DISPATCH(float, 6);
NEKO_AX_HELM_MFMA_DISPATCH(float, 7);
NEKO_AX_HELM_MFMA_DISPATCH(float, 8);
NEKO_AX_HELM_MFMA_DISPATCH(float, 9);
NEKO_AX_HELM_MFMA_DISPATCH(float, 10);
NEKO_AX_HELM_MFMA_DISPATCH(float, 11);
NEKO_AX_HELM_MFMA_DISPATCH(float, 12);

#endif // __gfx90a__ || __gfx942__

/*
 * Note the bare __launch_bounds__ rather than NEKO_EB_BOUNDS: the kstep
 * kernels ask for three waves per SIMD, and this kernel was validated on
 * gfx90a/gfx942 without that constraint. Keep it byte-identical to the
 * configuration that was confirmed on hardware.
 */
template< typename T, const int LX, const int NWF >
__global__ void __launch_bounds__(64 * NWF)
ax_helm_kernel_mfma(T * __restrict__ w,
                    const T * __restrict__ u,
                    const T * __restrict__ dx,
                    const T * __restrict__ dy,
                    const T * __restrict__ dz,
                    const T * __restrict__ h1,
                    const T * __restrict__ g11,
                    const T * __restrict__ g22,
                    const T * __restrict__ g33,
                    const T * __restrict__ g12,
                    const T * __restrict__ g13,
                    const T * __restrict__ g23,
                    const int nelv) {

  ax_helm_mfma_dispatch< T, LX, NWF >::run(w, u, dx, dy, dz, h1,
                                           g11, g22, g33, g12, g13, g23, nelv);
}

/*
 * Vector versions
 */

template< typename T, const int LX, const int EB >
__global__ void NEKO_EB_BOUNDS(LX*LX*EB)
ax_helm_kernel_vector_kstep(T * __restrict__ au,
                            T * __restrict__ av,
                            T * __restrict__ aw,
                            const T * __restrict__ u,
                            const T * __restrict__ v,
                            const T * __restrict__ w,
                            const T * __restrict__ dx,
                            const T * __restrict__ dy,
                            const T * __restrict__ dz,
                            const T * __restrict__ h1,
                            const T * __restrict__ g11,
                            const T * __restrict__ g22,
                            const T * __restrict__ g33,
                            const T * __restrict__ g12,
                            const T * __restrict__ g13,
                            const T * __restrict__ g23,
                            const int nelv) {

  /* Element independent, one copy per block */
  __shared__ T shdx[LX * LX];
  __shared__ T shdy[LX * LX];
  __shared__ T shdz[LX * LX];

  /* One slice per element in the block */
  __shared__ T shu[EB * LX * LX];
  __shared__ T shur[EB * LX * LX];
  __shared__ T shus[EB * LX * LX];

  __shared__ T shv[EB * LX * LX];
  __shared__ T shvr[EB * LX * LX];
  __shared__ T shvs[EB * LX * LX];

  __shared__ T shw[EB * LX * LX];
  __shared__ T shwr[EB * LX * LX];
  __shared__ T shws[EB * LX * LX];

  static_assert(sizeof(shdx) +
                sizeof(shdy) +
                sizeof(shdz) +
                sizeof(shu) +
                sizeof(shur) +
                sizeof(shus) +
                sizeof(shv) +
                sizeof(shvr) +
                sizeof(shvs) +
                sizeof(shw) +
                sizeof(shwr) +
                sizeof(shws)
                <= NEKO_EB_MAX_LDS,
                "kstep block exceeds the LDS budget");

  T ru[LX];
  T rv[LX];
  T rw[LX];

  T ruw[LX];
  T rvw[LX];
  T rww[LX];

  T rut;
  T rvt;
  T rwt;

  /* At EB == 1 the grid covers nelv exactly, so the blocking bookkeeping is
     constant folded away and the kernel is exactly what it was before */
  const int eb = (EB == 1) ? 0 : threadIdx.z;
  const int e_blk = blockIdx.x * EB + eb;
  const bool active = (EB == 1) ? true : (e_blk < nelv);
  const int e = active ? e_blk : (nelv - 1);
  const int j = threadIdx.y;
  const int i = threadIdx.x;
  const int ij = i + j*LX;
  const int sh = eb*LX*LX;
  const int ele = e*LX*LX*LX;

  if (eb == 0) {
    shdx[ij] = dx[ij];
    shdy[ij] = dy[ij];
    shdz[ij] = dz[ij];
  }

#pragma unroll
  for(int k = 0; k < LX; ++k){
    ru[k] = u[ij + k*LX*LX + ele];
    ruw[k] = 0.0;

    rv[k] = v[ij + k*LX*LX + ele];
    rvw[k] = 0.0;

    rw[k] = w[ij + k*LX*LX + ele];
    rww[k] = 0.0;
  }


  __syncthreads();
#pragma unroll
  for (int k = 0; k < LX; ++k){
    const int ijk = ij + k*LX*LX;
    const T G00 = g11[ijk+ele];
    const T G11 = g22[ijk+ele];
    const T G22 = g33[ijk+ele];
    const T G01 = g12[ijk+ele];
    const T G02 = g13[ijk+ele];
    const T G12 = g23[ijk+ele];
    const T H1  = h1[ijk+ele];
    T uttmp = 0.0;
    T vttmp = 0.0;
    T wttmp = 0.0;
    shu[sh + ij] = ru[k];
    shv[sh + ij] = rv[k];
    shw[sh + ij] = rw[k];
#pragma unroll
    for (int l = 0; l < LX; l++){
      uttmp += shdz[k+l*LX] * ru[l];
      vttmp += shdz[k+l*LX] * rv[l];
      wttmp += shdz[k+l*LX] * rw[l];
    }
    __syncthreads();

    T urtmp = 0.0;
    T ustmp = 0.0;

    T vrtmp = 0.0;
    T vstmp = 0.0;

    T wrtmp = 0.0;
    T wstmp = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++){
      urtmp += shdx[i+l*LX] * shu[sh + l+j*LX];
      ustmp += shdy[j+l*LX] * shu[sh + i+l*LX];

      vrtmp += shdx[i+l*LX] * shv[sh + l+j*LX];
      vstmp += shdy[j+l*LX] * shv[sh + i+l*LX];

      wrtmp += shdx[i+l*LX] * shw[sh + l+j*LX];
      wstmp += shdy[j+l*LX] * shw[sh + i+l*LX];
    }

    shur[sh + ij] = H1
                  * (G00 * urtmp
                     + G01 * ustmp
                     + G02 * uttmp);
    shus[sh + ij] = H1
                  * (G01 * urtmp
                     + G11 * ustmp
                     + G12 * uttmp);
    rut           = H1
                  * (G02 * urtmp
                     + G12 * ustmp
                     + G22 * uttmp);

    shvr[sh + ij] = H1
                  * (G00 * vrtmp
                     + G01 * vstmp
                     + G02 * vttmp);
    shvs[sh + ij] = H1
                  * (G01 * vrtmp
                     + G11 * vstmp
                     + G12 * vttmp);
    rvt           = H1
                  * (G02 * vrtmp
                     + G12 * vstmp
                     + G22 * vttmp);

    shwr[sh + ij] = H1
                  * (G00 * wrtmp
                     + G01 * wstmp
                     + G02 * wttmp);
    shws[sh + ij] = H1
                  * (G01 * wrtmp
                     + G11 * wstmp
                     + G12 * wttmp);
    rwt           = H1
                  * (G02 * wrtmp
                     + G12 * wstmp
                     + G22 * wttmp);

    __syncthreads();

    T uwijke = 0.0;
    T vwijke = 0.0;
    T wwijke = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++){
      uwijke += shur[sh + l+j*LX] * shdx[l+i*LX];
      ruw[l] += rut * shdz[k+l*LX];
      uwijke += shus[sh + i+l*LX] * shdy[l + j*LX];

      vwijke += shvr[sh + l+j*LX] * shdx[l+i*LX];
      rvw[l] += rvt * shdz[k+l*LX];
      vwijke += shvs[sh + i+l*LX] * shdy[l + j*LX];

      wwijke += shwr[sh + l+j*LX] * shdx[l+i*LX];
      rww[l] += rwt * shdz[k+l*LX];
      wwijke += shws[sh + i+l*LX] * shdy[l + j*LX];
    }
    ruw[k] += uwijke;
    rvw[k] += vwijke;
    rww[k] += wwijke;
  }
  if (active) {
#pragma unroll
    for (int k = 0; k < LX; ++k){
      au[ij + k*LX*LX + ele] = ruw[k];
      av[ij + k*LX*LX + ele] = rvw[k];
      aw[ij + k*LX*LX + ele] = rww[k];
    }
  }
}

template< typename T, const int LX, const int EB >
__global__ void NEKO_EB_BOUNDS(LX*LX*EB)
ax_helm_kernel_vector_kstep_padded(T * __restrict__ au,
                                   T * __restrict__ av,
                                   T * __restrict__ aw,
                                   const T * __restrict__ u,
                                   const T * __restrict__ v,
                                   const T * __restrict__ w,
                                   const T * __restrict__ dx,
                                   const T * __restrict__ dy,
                                   const T * __restrict__ dz,
                                   const T * __restrict__ h1,
                                   const T * __restrict__ g11,
                                   const T * __restrict__ g22,
                                   const T * __restrict__ g33,
                                   const T * __restrict__ g12,
                                   const T * __restrict__ g13,
                                   const T * __restrict__ g23,
                                   const int nelv) {

  /* Element independent, one copy per block */
  __shared__ T shdx[LX * (LX+1)];
  __shared__ T shdy[LX * (LX+1)];
  __shared__ T shdz[LX * (LX+1)];

  /* One slice per element in the block */
  __shared__ T shu[EB * LX * (LX+1)];
  __shared__ T shur[EB * LX * LX];
  __shared__ T shus[EB * LX * (LX+1)];

  __shared__ T shv[EB * LX * (LX+1)];
  __shared__ T shvr[EB * LX * LX];
  __shared__ T shvs[EB * LX * (LX+1)];

  __shared__ T shw[EB * LX * (LX+1)];
  __shared__ T shwr[EB * LX * LX];
  __shared__ T shws[EB * LX * (LX+1)];

  static_assert(sizeof(shdx) +
                sizeof(shdy) +
                sizeof(shdz) +
                sizeof(shu) +
                sizeof(shur) +
                sizeof(shus) +
                sizeof(shv) +
                sizeof(shvr) +
                sizeof(shvs) +
                sizeof(shw) +
                sizeof(shwr) +
                sizeof(shws)
                <= NEKO_EB_MAX_LDS,
                "kstep block exceeds the LDS budget");

  T ru[LX];
  T rv[LX];
  T rw[LX];

  T ruw[LX];
  T rvw[LX];
  T rww[LX];

  T rut;
  T rvt;
  T rwt;

  /* At EB == 1 the grid covers nelv exactly, so the blocking bookkeeping is
     constant folded away and the kernel is exactly what it was before */
  const int eb = (EB == 1) ? 0 : threadIdx.z;
  const int e_blk = blockIdx.x * EB + eb;
  const bool active = (EB == 1) ? true : (e_blk < nelv);
  const int e = active ? e_blk : (nelv - 1);
  const int j = threadIdx.y;
  const int i = threadIdx.x;
  const int ij = i + j*LX;
  const int ij_p = i + j*(LX+1);
  const int sh = eb*LX*LX;
  const int sh_p = eb*LX*(LX+1);
  const int ele = e*LX*LX*LX;

  if (eb == 0) {
    shdx[ij_p] = dx[ij];
    shdy[ij_p] = dy[ij];
    shdz[ij_p] = dz[ij];
  }

#pragma unroll
  for(int k = 0; k < LX; ++k){
    ru[k] = u[ij + k*LX*LX + ele];
    ruw[k] = 0.0;

    rv[k] = v[ij + k*LX*LX + ele];
    rvw[k] = 0.0;

    rw[k] = w[ij + k*LX*LX + ele];
    rww[k] = 0.0;
  }


  __syncthreads();
#pragma unroll
  for (int k = 0; k < LX; ++k){
    const int ijk = ij + k*LX*LX;
    const T G00 = g11[ijk+ele];
    const T G11 = g22[ijk+ele];
    const T G22 = g33[ijk+ele];
    const T G01 = g12[ijk+ele];
    const T G02 = g13[ijk+ele];
    const T G12 = g23[ijk+ele];
    const T H1  = h1[ijk+ele];
    T uttmp = 0.0;
    T vttmp = 0.0;
    T wttmp = 0.0;
    shu[sh_p + ij_p] = ru[k];
    shv[sh_p + ij_p] = rv[k];
    shw[sh_p + ij_p] = rw[k];
#pragma unroll
    for (int l = 0; l < LX; l++){
      uttmp += shdz[k+l*(LX+1)] * ru[l];
      vttmp += shdz[k+l*(LX+1)] * rv[l];
      wttmp += shdz[k+l*(LX+1)] * rw[l];
    }
    __syncthreads();

    T urtmp = 0.0;
    T ustmp = 0.0;

    T vrtmp = 0.0;
    T vstmp = 0.0;

    T wrtmp = 0.0;
    T wstmp = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++){
      urtmp += shdx[i+l*(LX+1)] * shu[sh_p + l+j*(LX+1)];
      ustmp += shdy[j+l*(LX+1)] * shu[sh_p + i+l*(LX+1)];

      vrtmp += shdx[i+l*(LX+1)] * shv[sh_p + l+j*(LX+1)];
      vstmp += shdy[j+l*(LX+1)] * shv[sh_p + i+l*(LX+1)];

      wrtmp += shdx[i+l*(LX+1)] * shw[sh_p + l+j*(LX+1)];
      wstmp += shdy[j+l*(LX+1)] * shw[sh_p + i+l*(LX+1)];
    }

    shur[sh + ij]     = H1
                      * (G00 * urtmp
                         + G01 * ustmp
                         + G02 * uttmp);
    shus[sh_p + ij_p] = H1
                      * (G01 * urtmp
                         + G11 * ustmp
                         + G12 * uttmp);
    rut               = H1
                      * (G02 * urtmp
                         + G12 * ustmp
                         + G22 * uttmp);

    shvr[sh + ij]     = H1
                      * (G00 * vrtmp
                         + G01 * vstmp
                         + G02 * vttmp);
    shvs[sh_p + ij_p] = H1
                      * (G01 * vrtmp
                         + G11 * vstmp
                         + G12 * vttmp);
    rvt               = H1
                      * (G02 * vrtmp
                         + G12 * vstmp
                         + G22 * vttmp);

    shwr[sh + ij]     = H1
                      * (G00 * wrtmp
                         + G01 * wstmp
                         + G02 * wttmp);
    shws[sh_p + ij_p] = H1
                      * (G01 * wrtmp
                         + G11 * wstmp
                         + G12 * wttmp);
    rwt               = H1
                      * (G02 * wrtmp
                         + G12 * wstmp
                         + G22 * wttmp);

    __syncthreads();

    T uwijke = 0.0;
    T vwijke = 0.0;
    T wwijke = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++){
      uwijke += shur[sh + l+j*LX] * shdx[l+i*(LX+1)];
      ruw[l] += rut * shdz[k+l*(LX+1)];
      uwijke += shus[sh_p + i+l*(LX+1)] * shdy[l + j*(LX+1)];

      vwijke += shvr[sh + l+j*LX] * shdx[l+i*(LX+1)];
      rvw[l] += rvt * shdz[k+l*(LX+1)];
      vwijke += shvs[sh_p + i+l*(LX+1)] * shdy[l + j*(LX+1)];

      wwijke += shwr[sh + l+j*LX] * shdx[l+i*(LX+1)];
      rww[l] += rwt * shdz[k+l*(LX+1)];
      wwijke += shws[sh_p + i+l*(LX+1)] * shdy[l + j*(LX+1)];
    }
    ruw[k] += uwijke;
    rvw[k] += vwijke;
    rww[k] += wwijke;
  }
  if (active) {
#pragma unroll
    for (int k = 0; k < LX; ++k){
      au[ij + k*LX*LX + ele] = ruw[k];
      av[ij + k*LX*LX + ele] = rvw[k];
      aw[ij + k*LX*LX + ele] = rww[k];
    }
  }
}

template< typename T >
__global__ void ax_helm_kernel_vector_part2(T * __restrict__ au,
                                            T * __restrict__ av,
                                            T * __restrict__ aw,
                                            const T * __restrict__ u,
                                            const T * __restrict__ v,
                                            const T * __restrict__ w,
                                            const T * __restrict__ h2,
                                            const T * __restrict__ B,
                                            const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    au[i] = au[i] + h2[i] * B[i] * u[i];
    av[i] = av[i] + h2[i] * B[i] * v[i];
    aw[i] = aw[i] + h2[i] * B[i] * w[i];
  }

}
#endif // __MATH_AX_HELM_KERNEL_H__
