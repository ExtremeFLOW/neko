#ifndef __MATH_AX_HELM_KERNEL_H__
#define __MATH_AX_HELM_KERNEL_H__
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
#include "dmma_kernel.h"
#include "dmma_tma_kernel.h"

/*
 * A note on elements per block for the vector kstep kernels.
 *
 * These hold three times the register blocked state of the scalar ones -- six
 * T[LX] arrays rather than two -- and measure at 254-255 registers on sm_90
 * for every lx from 8 up, spilling at 10, 12, 14 and 16. Blocking is therefore
 * not expected to buy much: it does not change registers per thread, only
 * threads per block, so at 255 registers the occupancy is the same either way
 * and all it saves is loading the derivative matrices once per block instead
 * of once per element.
 *
 * They are swept anyway, over the same elem_block<> candidates as the scalar
 * kernels. "Not expected to buy much" is a prediction, and pinning it at build
 * time is a prediction the tuner can never check. elem_block<>'s thread clamp
 * keeps every candidate inside the shared memory budget at every lx -- the
 * widest case, lx = 11 at four elements per block, comes to 37 kB against the
 * 48 kB cap -- so nothing here needs a bound of its own.
 */

/**
 * Device kernel for axhelm
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
                <= NEKO_EB_MAX_SMEM,
                "kstep block exceeds the shared memory budget");

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
                <= NEKO_EB_MAX_SMEM,
                "kstep block exceeds the shared memory budget");

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
 * Device kernel for axhelm on the fp64 tensor cores
 *
 * One element per block, NW warps, and the whole element resident in shared
 * memory as four DMMA_P^3 cubes: the staged input (reused as the output), and
 * the three reference derivatives. Unlike the kstep variants, which stream a
 * k plane at a time and keep the k contraction in registers, all six
 * contractions here are full D * U GEMMs handed to dmma_contract(), with the
 * geometric factors applied pointwise in between. See dmma_kernel.h for the
 * padded staging, the per axis matrix views and the arch and LX bounds.
 */

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 800) && (__CUDA_ARCH__ < 1000)

template< const int LX, const int NW >
__device__ __forceinline__
void ax_helm_dmma_elem(double * __restrict__ w,
                       const double * __restrict__ u,
                       const double * __restrict__ dx,
                       const double * __restrict__ dy,
                       const double * __restrict__ dz,
                       const double * __restrict__ h1,
                       const double * __restrict__ g11,
                       const double * __restrict__ g22,
                       const double * __restrict__ g33,
                       const double * __restrict__ g12,
                       const double * __restrict__ g13,
                       const double * __restrict__ g23,
                       const int nelv) {

  /* Element independent, one copy per block. Block diagonal when more than
     one element is packed: one LX x LX copy of D per sub-cube */
  __shared__ __align__(16) double shdx[DMMA_MAT];
  __shared__ __align__(16) double shdy[DMMA_MAT];
  __shared__ __align__(16) double shdz[DMMA_MAT];

  /* The pack, padded to DMMA_P^3. shu carries the input and then the
     output, shr, shs and sht the reference derivatives */
  __shared__ __align__(16) double shu[DMMA_CUBE];
  __shared__ __align__(16) double shr[DMMA_CUBE];
  __shared__ __align__(16) double shs[DMMA_CUBE];
  __shared__ __align__(16) double sht[DMMA_CUBE];

  static_assert(sizeof(shdx) +
                sizeof(shdy) +
                sizeof(shdz) +
                sizeof(shu) +
                sizeof(shr) +
                sizeof(shs) +
                sizeof(sht)
                <= NEKO_EB_MAX_SMEM,
                "dmma block exceeds the shared memory budget");

  /* Elements per sub-cube axis, and per cube, see the note in dmma_kernel.h.
     At PPA == 1 the addressing is not left to constant fold -- dmma_pack is
     specialised so the tail clamp and the guarded store are never emitted at
     all, see the note there */
  enum { PPA = (DMMA_P % LX == 0) ? (DMMA_P / LX) : 1,
         PACK = PPA * PPA * PPA,
         LX3 = LX * LX * LX,
         NP = PACK * LX3 };

  typedef dmma_pack< LX, PPA > pack;

  const int nthrds = 32 * NW;
  const int tid = threadIdx.x;
  const int wf = tid >> 5;
  const int ebase = pack::ebase();

  /* The padding has to be finite, see the note above. At LX == DMMA_P with
     one element packed there is none and this is folded away */
  if (PACK * LX3 < DMMA_CUBE) {
    for (int p = tid; p < DMMA_CUBE; p += nthrds) {
      shu[p] = 0.0;
    }
  }
  if (LX < DMMA_P) {
    for (int p = tid; p < DMMA_MAT; p += nthrds) {
      shdx[p] = 0.0;
      shdy[p] = 0.0;
      shdz[p] = 0.0;
    }
  }
  if (LX < DMMA_P) {
    __syncthreads();
  }

  /* One copy of D per sub-cube, down the diagonal */
  for (int p = tid; p < LX * LX; p += nthrds) {
    const int i = p % LX;
    const int l = p / LX;
#pragma unroll
    for (int b = 0; b < PPA; b++) {
      const int m = (b * LX + i) + DMMA_P * (b * LX + l);
      shdx[m] = dx[p];
      shdy[m] = dy[p];
      shdz[m] = dz[p];
    }
  }

  for (int p = tid; p < NP; p += nthrds) {
    const dmma_idx x = pack::map(p, ebase, nelv);

    shu[x.c] = u[x.g];
  }

  __syncthreads();

  dmma_contract< 0, false, false, NW >(shr, shdx, shu, wf);
  dmma_contract< 1, false, false, NW >(shs, shdy, shu, wf);
  dmma_contract< 2, false, false, NW >(sht, shdz, shu, wf);

  __syncthreads();

  for (int p = tid; p < NP; p += nthrds) {
    const dmma_idx x = pack::map(p, ebase, nelv);
    const int c = x.c;
    const int gp = x.g;

    const double G00 = g11[gp];
    const double G11 = g22[gp];
    const double G22 = g33[gp];
    const double G01 = g12[gp];
    const double G02 = g13[gp];
    const double G12 = g23[gp];
    const double H1  = h1[gp];

    const double rtmp = shr[c];
    const double stmp = shs[c];
    const double ttmp = sht[c];

    shr[c] = H1
           * (G00 * rtmp
              + G01 * stmp
              + G02 * ttmp);
    shs[c] = H1
           * (G01 * rtmp
              + G11 * stmp
              + G12 * ttmp);
    sht[c] = H1
           * (G02 * rtmp
              + G12 * stmp
              + G22 * ttmp);
  }

  __syncthreads();

  /* The result overwrites the staged input; the first contraction writes
     every tile of the cube, so nothing has to be cleared first. The barriers
     are needed because the j slab tiles of axis 2 span what every warp wrote
     along axis 0 and 1 */
  dmma_contract< 0, true, false, NW >(shu, shdx, shr, wf);
  __syncthreads();
  dmma_contract< 1, true, true, NW >(shu, shdy, shs, wf);
  __syncthreads();
  dmma_contract< 2, true, true, NW >(shu, shdz, sht, wf);
  __syncthreads();

  for (int p = tid; p < NP; p += nthrds) {
    const dmma_idx x = pack::map(p, ebase, nelv);

    if (x.live) {
      w[x.g] = shu[x.c];
    }
  }
}

#endif // __CUDA_ARCH__ in [800, 1000)

/*
 * Compile-time dispatch onto the DMMA element kernel. The launch macros in
 * ax_helm.cu are written for every LX the operator dispatches and for whatever
 * `real` is, so every combination has to compile; the ones the strategy does
 * not cover -- single precision, LX outside the supported range, a build
 * without an fp64 tensor core arch -- resolve to this no-op. The autotuner
 * never selects the strategy for them, so the no-op is unreachable at runtime,
 * see dmma_lx_supported() and cuda_have_dmma() in dmma_kernel.h.
 */
template< typename T, const int LX, const int NW >
struct ax_helm_dmma_dispatch {
  __device__ static void run(T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const int) { }
};

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 800) && (__CUDA_ARCH__ < 1000)

/* Keep in sync with dmma_lx_supported() in dmma_kernel.h */
#define NEKO_AX_HELM_DMMA_DISPATCH(LXV)                                        \
  template< const int NW >                                                     \
  struct ax_helm_dmma_dispatch< double, LXV, NW > {                            \
    __device__ static void run(double * __restrict__ w,                        \
                               const double * __restrict__ u,                  \
                               const double * __restrict__ dx,                 \
                               const double * __restrict__ dy,                 \
                               const double * __restrict__ dz,                 \
                               const double * __restrict__ h1,                 \
                               const double * __restrict__ g11,                \
                               const double * __restrict__ g22,                \
                               const double * __restrict__ g33,                \
                               const double * __restrict__ g12,                \
                               const double * __restrict__ g13,                \
                               const double * __restrict__ g23,                \
                               const int nelv) {                               \
      ax_helm_dmma_elem< LXV, NW >(w, u, dx, dy, dz, h1,                       \
                                   g11, g22, g33, g12, g13, g23, nelv);        \
    }                                                                          \
  }

NEKO_AX_HELM_DMMA_DISPATCH(2);
NEKO_AX_HELM_DMMA_DISPATCH(3);
NEKO_AX_HELM_DMMA_DISPATCH(4);
NEKO_AX_HELM_DMMA_DISPATCH(5);
NEKO_AX_HELM_DMMA_DISPATCH(6);
NEKO_AX_HELM_DMMA_DISPATCH(7);
NEKO_AX_HELM_DMMA_DISPATCH(8);

#endif // __CUDA_ARCH__ in [800, 1000)

template< typename T, const int LX, const int NW >
__global__ void NEKO_EB_BOUNDS(32 * NW)
ax_helm_kernel_dmma(T * __restrict__ w,
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

  ax_helm_dmma_dispatch< T, LX, NW >::run(w, u, dx, dy, dz, h1,
                                          g11, g22, g33, g12, g13, g23, nelv);
}

/**
 * Device kernel for axhelm on the fp64 tensor cores, with the element staged
 * by the TMA engine
 *
 * Same six contractions and the same padded cube layout as
 * ax_helm_dmma_elem(), and at lx == DMMA_P the padding is empty and the cube
 * offset is the flat point index, so the arithmetic is identical -- the only
 * thing that differs is how the eleven cubes get in and out of shared memory.
 *
 * Here u arrives as one bulk copy on its own mbarrier and the seven geometric
 * factors as seven more on a second, and only the first is waited on before
 * the contractions start. The factors land in shared memory across the three
 * contractions and are waited on at the pointwise step that consumes them,
 * which is 78% of the read traffic moved out from between two barriers and
 * put underneath the tensor cores. The result leaves as a single bulk store.
 *
 * See dmma_tma_kernel.h for the primitives, the sm_90 and toolkit guards, the
 * lx == DMMA_P bound and the occupancy trade this buys the overlap with.
 */

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 900) &&                        \
    (__CUDA_ARCH__ < 1000) && NEKO_TMA_TOOLKIT

template< const int LX, const int NW >
__device__ __forceinline__
void ax_helm_dmma_tma_elem(double * __restrict__ w,
                           const double * __restrict__ u,
                           const double * __restrict__ dx,
                           const double * __restrict__ dy,
                           const double * __restrict__ dz,
                           const double * __restrict__ h1,
                           const double * __restrict__ g11,
                           const double * __restrict__ g22,
                           const double * __restrict__ g33,
                           const double * __restrict__ g12,
                           const double * __restrict__ g13,
                           const double * __restrict__ g23) {

  /* A bulk copy needs 16 byte alignment at both ends; the cubes are given the
     128 the tensor variants would want anyway, which costs nothing here since
     every one of them is a whole number of 128 byte lines */
  __shared__ __align__(128) double shdx[DMMA_MAT];
  __shared__ __align__(128) double shdy[DMMA_MAT];
  __shared__ __align__(128) double shdz[DMMA_MAT];

  /* shu carries the input and then the output, shr, shs and sht the reference
     derivatives */
  __shared__ __align__(128) double shu[DMMA_CUBE];
  __shared__ __align__(128) double shr[DMMA_CUBE];
  __shared__ __align__(128) double shs[DMMA_CUBE];
  __shared__ __align__(128) double sht[DMMA_CUBE];

  /* h1, g11, g22, g33, g12, g13, g23 */
  __shared__ __align__(128) double shg[DMMA_NG][DMMA_CUBE];

  __shared__ __align__(8) unsigned long long bar_u;
  __shared__ __align__(8) unsigned long long bar_g;

  static_assert(sizeof(shdx) +
                sizeof(shdy) +
                sizeof(shdz) +
                sizeof(shu) +
                sizeof(shr) +
                sizeof(shs) +
                sizeof(sht) +
                sizeof(shg) +
                sizeof(bar_u) +
                sizeof(bar_g)
                <= NEKO_EB_MAX_SMEM,
                "dmma tma block exceeds the shared memory budget");

  /* Only lx == DMMA_P stages as a contiguous run of bytes, which is what a
     bulk copy moves; see the scope note in dmma_tma_kernel.h. Keep in step
     with dmma_tma_lx_supported() */
  static_assert(LX == DMMA_P,
                "the dmma tma variant stages whole cubes only");

  enum { CUBE_BYTES = DMMA_CUBE * (int) sizeof(double) };

  const int nthrds = 32 * NW;
  const int tid = threadIdx.x;
  const int wf = tid >> 5;
  const int ebase = blockIdx.x * DMMA_CUBE;

  if (tid == 0) {
    tma_barrier_init(&bar_u, 1);
    tma_barrier_init(&bar_g, 1);
  }

  /* At LX == DMMA_P the derivative matrix fills the staged one exactly, and
     there is no padding anywhere to zero */
  for (int p = tid; p < DMMA_MAT; p += nthrds) {
    shdx[p] = dx[p];
    shdy[p] = dy[p];
    shdz[p] = dz[p];
  }

  /* Both barriers initialised and both derivative matrices staged before
     anyone waits on the one or contracts with the other */
  __syncthreads();

  if (tid == 0) {
    tma_expect(&bar_u, CUBE_BYTES);
    tma_load(shu, u + ebase, CUBE_BYTES, &bar_u);

    tma_expect(&bar_g, DMMA_NG * CUBE_BYTES);
    tma_load(shg[0], h1  + ebase, CUBE_BYTES, &bar_g);
    tma_load(shg[1], g11 + ebase, CUBE_BYTES, &bar_g);
    tma_load(shg[2], g22 + ebase, CUBE_BYTES, &bar_g);
    tma_load(shg[3], g33 + ebase, CUBE_BYTES, &bar_g);
    tma_load(shg[4], g12 + ebase, CUBE_BYTES, &bar_g);
    tma_load(shg[5], g13 + ebase, CUBE_BYTES, &bar_g);
    tma_load(shg[6], g23 + ebase, CUBE_BYTES, &bar_g);
  }

  /* u only. The seven factor cubes are still arriving */
  tma_wait(&bar_u, 0);

  dmma_contract< 0, false, false, NW >(shr, shdx, shu, wf);
  dmma_contract< 1, false, false, NW >(shs, shdy, shu, wf);
  dmma_contract< 2, false, false, NW >(sht, shdz, shu, wf);

  __syncthreads();

  tma_wait(&bar_g, 0);

  for (int p = tid; p < DMMA_CUBE; p += nthrds) {
    const double H1  = shg[0][p];
    const double G00 = shg[1][p];
    const double G11 = shg[2][p];
    const double G22 = shg[3][p];
    const double G01 = shg[4][p];
    const double G02 = shg[5][p];
    const double G12 = shg[6][p];

    const double rtmp = shr[p];
    const double stmp = shs[p];
    const double ttmp = sht[p];

    shr[p] = H1
           * (G00 * rtmp
              + G01 * stmp
              + G02 * ttmp);
    shs[p] = H1
           * (G01 * rtmp
              + G11 * stmp
              + G12 * ttmp);
    sht[p] = H1
           * (G02 * rtmp
              + G12 * stmp
              + G22 * ttmp);
  }

  __syncthreads();

  /* The result overwrites the staged input, exactly as in the scalar dmma
     kernel; see the note there on why the barriers between these are needed */
  dmma_contract< 0, true, false, NW >(shu, shdx, shr, wf);
  __syncthreads();
  dmma_contract< 1, true, true, NW >(shu, shdy, shs, wf);
  __syncthreads();
  dmma_contract< 2, true, true, NW >(shu, shdz, sht, wf);

  /* The contractions wrote shu through the generic proxy and the bulk store
     reads it through the async one, so the fence is needed on top of the
     barrier; see tma_fence_shared() */
  tma_fence_shared();
  __syncthreads();

  if (tid == 0) {
    tma_store(w + ebase, shu, CUBE_BYTES);
    tma_store_wait();
  }

  /* The store reads shu asynchronously, and shared memory lives only as long
     as the block does: nothing may retire until it has been read out */
  __syncthreads();
}

#endif // __CUDA_ARCH__ == sm_90 with a CUDA 12 toolkit

/*
 * Compile-time dispatch onto the TMA staged DMMA element kernel, see the note
 * on ax_helm_dmma_dispatch above. The no-op covers everything the strategy
 * does not: single precision, any lx but DMMA_P, a build without sm_90, and a
 * toolkit older than CUDA 12.
 */
template< typename T, const int LX, const int NW >
struct ax_helm_dmma_tma_dispatch {
  __device__ static void run(T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__) { }
};

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 900) &&                        \
    (__CUDA_ARCH__ < 1000) && NEKO_TMA_TOOLKIT

/* Keep in sync with dmma_tma_lx_supported() in dmma_tma_kernel.h */
#define NEKO_AX_HELM_DMMA_TMA_DISPATCH(LXV)                                    \
  template< const int NW >                                                     \
  struct ax_helm_dmma_tma_dispatch< double, LXV, NW > {                        \
    __device__ static void run(double * __restrict__ w,                        \
                               const double * __restrict__ u,                  \
                               const double * __restrict__ dx,                 \
                               const double * __restrict__ dy,                 \
                               const double * __restrict__ dz,                 \
                               const double * __restrict__ h1,                 \
                               const double * __restrict__ g11,                \
                               const double * __restrict__ g22,                \
                               const double * __restrict__ g33,                \
                               const double * __restrict__ g12,                \
                               const double * __restrict__ g13,                \
                               const double * __restrict__ g23) {              \
      ax_helm_dmma_tma_elem< LXV, NW >(w, u, dx, dy, dz, h1,                   \
                                       g11, g22, g33, g12, g13, g23);          \
    }                                                                          \
  }

NEKO_AX_HELM_DMMA_TMA_DISPATCH(8);

#endif // __CUDA_ARCH__ == sm_90 with a CUDA 12 toolkit

template< typename T, const int LX, const int NW >
__global__ void NEKO_EB_BOUNDS(32 * NW)
ax_helm_kernel_dmma_tma(T * __restrict__ w,
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
                        const T * __restrict__ g23) {

  ax_helm_dmma_tma_dispatch< T, LX, NW >::run(w, u, dx, dy, dz, h1,
                                              g11, g22, g33, g12, g13, g23);
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
                <= NEKO_EB_MAX_SMEM,
                "kstep block exceeds the shared memory budget");

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
                <= NEKO_EB_MAX_SMEM,
                "kstep block exceeds the shared memory budget");

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


/**
 * Device kernel for the vector axhelm on the fp64 tensor cores
 *
 * The three components share one set of geometric factors, which is the whole
 * point of the vector operator: reading them once and applying them to u, v
 * and w moves 13 arrays per element instead of the 27 that three scalar calls
 * would. Twelve staged cubes will not fit in a block's shared memory, so the
 * components are run one after another through the same four cubes and the
 * factors are hoisted into registers instead -- PPT points per thread, held
 * across all three passes. Everything else is the scalar kernel; see
 * dmma_kernel.h for the padded staging and the per axis matrix views.
 *
 * Note the padded entries of a staged cube only have to be *finite*, not zero.
 * A padded row of the derivative matrix is zero and kills the contribution of
 * a padded row of the cube, but 0 * NaN would not, so the cube is zeroed once
 * before its first use and then left to carry whatever the contractions put
 * there. That is why the component loop below re-stages only the LX^3 real
 * points and never clears the padding again.
 *
 * Measured on GH200 at lx = 8 this loses to the kstep variant: holding the
 * factors costs enough occupancy (~60 registers of metrics at nw = 4, against
 * the scalar kernel's handful) to give back more than the access pattern wins.
 * It is kept for the low order end, where PPT falls to 1 and the register cost
 * with it, while the shared memory footprint stays flat. The autotuner decides.
 */

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 800) && (__CUDA_ARCH__ < 1000)

template< const int LX, const int NW >
__device__ __forceinline__
void ax_helm_dmma_vector_elem(double * __restrict__ au,
                              double * __restrict__ av,
                              double * __restrict__ aw,
                              const double * __restrict__ u,
                              const double * __restrict__ v,
                              const double * __restrict__ w,
                              const double * __restrict__ dx,
                              const double * __restrict__ dy,
                              const double * __restrict__ dz,
                              const double * __restrict__ h1,
                              const double * __restrict__ g11,
                              const double * __restrict__ g22,
                              const double * __restrict__ g33,
                              const double * __restrict__ g12,
                              const double * __restrict__ g13,
                              const double * __restrict__ g23) {

  /* Element independent, one copy per block */
  __shared__ __align__(16) double shdx[DMMA_MAT];
  __shared__ __align__(16) double shdy[DMMA_MAT];
  __shared__ __align__(16) double shdz[DMMA_MAT];

  /* One component at a time: shc carries it in and the result out, shr, shs
     and sht its reference derivatives */
  __shared__ __align__(16) double shc[DMMA_CUBE];
  __shared__ __align__(16) double shr[DMMA_CUBE];
  __shared__ __align__(16) double shs[DMMA_CUBE];
  __shared__ __align__(16) double sht[DMMA_CUBE];

  static_assert(sizeof(shdx) +
                sizeof(shdy) +
                sizeof(shdz) +
                sizeof(shc) +
                sizeof(shr) +
                sizeof(shs) +
                sizeof(sht)
                <= NEKO_EB_MAX_SMEM,
                "dmma vector block exceeds the shared memory budget");

  enum { NTHRDS = 32 * NW,
         LX3 = LX * LX * LX,
         PPT = (LX3 + NTHRDS - 1) / NTHRDS };

  const int tid = threadIdx.x;
  const int wf = tid >> 5;
  const int ele = blockIdx.x * LX3;

  /* The geometric factors, read once and reused by all three components */
  double rG00[PPT], rG11[PPT], rG22[PPT];
  double rG01[PPT], rG02[PPT], rG12[PPT];
  double rH1[PPT];
  int rc[PPT];

  /* The padding only has to be finite, see the note above. At LX == DMMA_P
     there is none and this is folded away */
  if (LX < DMMA_P) {
    for (int p = tid; p < DMMA_MAT; p += NTHRDS) {
      shdx[p] = 0.0;
      shdy[p] = 0.0;
      shdz[p] = 0.0;
    }
    for (int p = tid; p < DMMA_CUBE; p += NTHRDS) {
      shc[p] = 0.0;
    }
    __syncthreads();
  }

  for (int p = tid; p < LX * LX; p += NTHRDS) {
    const int i = p % LX;
    const int l = p / LX;
    const int m = i + DMMA_P * l;
    shdx[m] = dx[p];
    shdy[m] = dy[p];
    shdz[m] = dz[p];
  }

#pragma unroll
  for (int q = 0; q < PPT; q++) {
    const int p = tid + q * NTHRDS;

    if (p < LX3) {
      const int i = p % LX;
      const int jk = p / LX;
      const int j = jk % LX;
      const int k = jk / LX;
      rc[q]   = i + DMMA_SI * j + DMMA_SJ * k;
      rG00[q] = g11[p + ele];
      rG11[q] = g22[p + ele];
      rG22[q] = g33[p + ele];
      rG01[q] = g12[p + ele];
      rG02[q] = g13[p + ele];
      rG12[q] = g23[p + ele];
      rH1[q]  = h1[p + ele];
    } else {
      rc[q]   = 0;
      rG00[q] = 0.0;
      rG11[q] = 0.0;
      rG22[q] = 0.0;
      rG01[q] = 0.0;
      rG02[q] = 0.0;
      rG12[q] = 0.0;
      rH1[q]  = 0.0;
    }
  }

  __syncthreads();

  const double * const cin[3]  = { u, v, w };
  double * const       cout[3] = { au, av, aw };

#pragma unroll
  for (int c = 0; c < 3; c++) {

    for (int p = tid; p < LX3; p += NTHRDS) {
      const int i = p % LX;
      const int jk = p / LX;
      const int j = jk % LX;
      const int k = jk / LX;
      shc[i + DMMA_SI * j + DMMA_SJ * k] = cin[c][p + ele];
    }

    __syncthreads();

    dmma_contract< 0, false, false, NW >(shr, shdx, shc, wf);
    dmma_contract< 1, false, false, NW >(shs, shdy, shc, wf);
    dmma_contract< 2, false, false, NW >(sht, shdz, shc, wf);

    __syncthreads();

#pragma unroll
    for (int q = 0; q < PPT; q++) {
      const int p = tid + q * NTHRDS;

      if (p < LX3) {
        const int idx = rc[q];
        const double rtmp = shr[idx];
        const double stmp = shs[idx];
        const double ttmp = sht[idx];

        shr[idx] = rH1[q]
                 * (rG00[q] * rtmp
                    + rG01[q] * stmp
                    + rG02[q] * ttmp);
        shs[idx] = rH1[q]
                 * (rG01[q] * rtmp
                    + rG11[q] * stmp
                    + rG12[q] * ttmp);
        sht[idx] = rH1[q]
                 * (rG02[q] * rtmp
                    + rG12[q] * stmp
                    + rG22[q] * ttmp);
      }
    }

    __syncthreads();

    dmma_contract< 0, true, false, NW >(shc, shdx, shr, wf);
    __syncthreads();
    dmma_contract< 1, true, true, NW >(shc, shdy, shs, wf);
    __syncthreads();
    dmma_contract< 2, true, true, NW >(shc, shdz, sht, wf);
    __syncthreads();

    for (int p = tid; p < LX3; p += NTHRDS) {
      const int i = p % LX;
      const int jk = p / LX;
      const int j = jk % LX;
      const int k = jk / LX;
      cout[c][p + ele] = shc[i + DMMA_SI * j + DMMA_SJ * k];
    }

    /* shc is restaged by the next component */
    __syncthreads();
  }
}

#endif // __CUDA_ARCH__ in [800, 1000)

/*
 * Compile-time dispatch onto the vector DMMA element kernel, see the note on
 * ax_helm_dmma_dispatch above.
 */
template< typename T, const int LX, const int NW >
struct ax_helm_dmma_vector_dispatch {
  __device__ static void run(T * __restrict__,
                             T * __restrict__,
                             T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__) { }
};

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 800) && (__CUDA_ARCH__ < 1000)

/* Keep in sync with dmma_vector_lx_supported() in dmma_kernel.h, which is
   4 <= LX <= 8 and not the 2 <= LX <= 8 of the packed scalar kernel */
#define NEKO_AX_HELM_DMMA_VECTOR_DISPATCH(LXV)                                 \
  template< const int NW >                                                     \
  struct ax_helm_dmma_vector_dispatch< double, LXV, NW > {                     \
    __device__ static void run(double * __restrict__ au,                       \
                               double * __restrict__ av,                       \
                               double * __restrict__ aw,                       \
                               const double * __restrict__ u,                  \
                               const double * __restrict__ v,                  \
                               const double * __restrict__ w,                  \
                               const double * __restrict__ dx,                 \
                               const double * __restrict__ dy,                 \
                               const double * __restrict__ dz,                 \
                               const double * __restrict__ h1,                 \
                               const double * __restrict__ g11,                \
                               const double * __restrict__ g22,                \
                               const double * __restrict__ g33,                \
                               const double * __restrict__ g12,                \
                               const double * __restrict__ g13,                \
                               const double * __restrict__ g23) {              \
      ax_helm_dmma_vector_elem< LXV, NW >(au, av, aw, u, v, w, dx, dy, dz, h1, \
                                          g11, g22, g33, g12, g13, g23);       \
    }                                                                          \
  }

NEKO_AX_HELM_DMMA_VECTOR_DISPATCH(4);
NEKO_AX_HELM_DMMA_VECTOR_DISPATCH(5);
NEKO_AX_HELM_DMMA_VECTOR_DISPATCH(6);
NEKO_AX_HELM_DMMA_VECTOR_DISPATCH(7);
NEKO_AX_HELM_DMMA_VECTOR_DISPATCH(8);

#endif // __CUDA_ARCH__ in [800, 1000)

template< typename T, const int LX, const int NW >
__global__ void NEKO_EB_BOUNDS(32 * NW)
ax_helm_kernel_dmma_vector(T * __restrict__ au,
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
                           const T * __restrict__ g23) {

  ax_helm_dmma_vector_dispatch< T, LX, NW >::run(au, av, aw, u, v, w,
                                                 dx, dy, dz, h1,
                                                 g11, g22, g33,
                                                 g12, g13, g23);
}

/**
 * Device kernel for the vector axhelm on the fp64 tensor cores, with the
 * element staged by the TMA engine
 *
 * The same trade as ax_helm_dmma_tma_elem() makes for the scalar operator,
 * applied to the thing that was holding the vector one back. The register
 * variant above reads the seven geometric factors once and keeps them in
 * registers across all three components, which is the right access pattern
 * and the wrong place to put it: ~60 registers at nw = 4, measured on GH200
 * as a loss against the kstep variant even though it moves 13 arrays per
 * element rather than 27. Here they are read once into shared memory instead,
 * by the TMA engine, and the register file is left alone.
 *
 * Everything else follows the scalar TMA kernel. The factor copies are issued
 * before the component loop and waited on only at the first pointwise step,
 * so they arrive underneath the staging and the first three contractions of
 * component 0; after that they are simply resident for components 1 and 2.
 * Each component is staged by one bulk copy and stored by another.
 *
 * The block footprint is the scalar kernel's exactly -- four working cubes
 * plus seven factor cubes -- because the components run one at a time through
 * the same four. That leaves no room for a second component cube, so a
 * component's staging is not overlapped with the previous component's
 * contractions; doing that needs either the opt-in dynamic shared memory path
 * or h1 moved back into registers to free a cube. Neither is done here.
 *
 * See dmma_tma_kernel.h for the primitives, the sm_90 and toolkit guards and
 * the lx == DMMA_P bound.
 */

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 900) &&                        \
    (__CUDA_ARCH__ < 1000) && NEKO_TMA_TOOLKIT

template< const int LX, const int NW >
__device__ __forceinline__
void ax_helm_dmma_tma_vector_elem(double * __restrict__ au,
                                  double * __restrict__ av,
                                  double * __restrict__ aw,
                                  const double * __restrict__ u,
                                  const double * __restrict__ v,
                                  const double * __restrict__ w,
                                  const double * __restrict__ dx,
                                  const double * __restrict__ dy,
                                  const double * __restrict__ dz,
                                  const double * __restrict__ h1,
                                  const double * __restrict__ g11,
                                  const double * __restrict__ g22,
                                  const double * __restrict__ g33,
                                  const double * __restrict__ g12,
                                  const double * __restrict__ g13,
                                  const double * __restrict__ g23) {

  /* Element independent, one copy per block */
  __shared__ __align__(128) double shdx[DMMA_MAT];
  __shared__ __align__(128) double shdy[DMMA_MAT];
  __shared__ __align__(128) double shdz[DMMA_MAT];

  /* One component at a time: shc carries it in and the result out, shr, shs
     and sht its reference derivatives */
  __shared__ __align__(128) double shc[DMMA_CUBE];
  __shared__ __align__(128) double shr[DMMA_CUBE];
  __shared__ __align__(128) double shs[DMMA_CUBE];
  __shared__ __align__(128) double sht[DMMA_CUBE];

  /* h1, g11, g22, g33, g12, g13, g23, shared by all three components */
  __shared__ __align__(128) double shg[DMMA_NG][DMMA_CUBE];

  __shared__ __align__(8) unsigned long long bar_c;
  __shared__ __align__(8) unsigned long long bar_g;

  static_assert(sizeof(shdx) +
                sizeof(shdy) +
                sizeof(shdz) +
                sizeof(shc) +
                sizeof(shr) +
                sizeof(shs) +
                sizeof(sht) +
                sizeof(shg) +
                sizeof(bar_c) +
                sizeof(bar_g)
                <= NEKO_EB_MAX_SMEM,
                "dmma tma vector block exceeds the shared memory budget");

  /* Keep in step with dmma_tma_vector_lx_supported() */
  static_assert(LX == DMMA_P,
                "the dmma tma vector variant stages whole cubes only");

  enum { CUBE_BYTES = DMMA_CUBE * (int) sizeof(double) };

  const int nthrds = 32 * NW;
  const int tid = threadIdx.x;
  const int wf = tid >> 5;
  const int ele = blockIdx.x * DMMA_CUBE;

  if (tid == 0) {
    tma_barrier_init(&bar_c, 1);
    tma_barrier_init(&bar_g, 1);
  }

  /* At LX == DMMA_P the derivative matrix fills the staged one exactly, and
     there is no padding anywhere to zero */
  for (int p = tid; p < DMMA_MAT; p += nthrds) {
    shdx[p] = dx[p];
    shdy[p] = dy[p];
    shdz[p] = dz[p];
  }

  /* Both barriers initialised and both derivative matrices staged before
     anyone waits on the one or contracts with the other */
  __syncthreads();

  /* Read once, outside the component loop, and waited on inside it at the
     first step that needs them */
  if (tid == 0) {
    tma_expect(&bar_g, DMMA_NG * CUBE_BYTES);
    tma_load(shg[0], h1  + ele, CUBE_BYTES, &bar_g);
    tma_load(shg[1], g11 + ele, CUBE_BYTES, &bar_g);
    tma_load(shg[2], g22 + ele, CUBE_BYTES, &bar_g);
    tma_load(shg[3], g33 + ele, CUBE_BYTES, &bar_g);
    tma_load(shg[4], g12 + ele, CUBE_BYTES, &bar_g);
    tma_load(shg[5], g13 + ele, CUBE_BYTES, &bar_g);
    tma_load(shg[6], g23 + ele, CUBE_BYTES, &bar_g);
  }

  const double * const cin[3]  = { u, v, w };
  double * const       cout[3] = { au, av, aw };

#pragma unroll
  for (int c = 0; c < 3; c++) {

    if (tid == 0) {
      tma_expect(&bar_c, CUBE_BYTES);
      tma_load(shc, cin[c] + ele, CUBE_BYTES, &bar_c);
    }

    /* bar_c is re-armed by every completion, so the parity of the phase each
       component's copy completes alternates. The __syncthreads() at the foot
       of the loop is what keeps the next component's arrive from racing a
       thread that has not yet observed this one */
    tma_wait(&bar_c, c & 1);

    dmma_contract< 0, false, false, NW >(shr, shdx, shc, wf);
    dmma_contract< 1, false, false, NW >(shs, shdy, shc, wf);
    dmma_contract< 2, false, false, NW >(sht, shdz, shc, wf);

    __syncthreads();

    /* Resident from here on; only the first component ever waits */
    if (c == 0) {
      tma_wait(&bar_g, 0);
    }

    for (int p = tid; p < DMMA_CUBE; p += nthrds) {
      const double H1  = shg[0][p];
      const double G00 = shg[1][p];
      const double G11 = shg[2][p];
      const double G22 = shg[3][p];
      const double G01 = shg[4][p];
      const double G02 = shg[5][p];
      const double G12 = shg[6][p];

      const double rtmp = shr[p];
      const double stmp = shs[p];
      const double ttmp = sht[p];

      shr[p] = H1
             * (G00 * rtmp
                + G01 * stmp
                + G02 * ttmp);
      shs[p] = H1
             * (G01 * rtmp
                + G11 * stmp
                + G12 * ttmp);
      sht[p] = H1
             * (G02 * rtmp
                + G12 * stmp
                + G22 * ttmp);
    }

    __syncthreads();

    dmma_contract< 0, true, false, NW >(shc, shdx, shr, wf);
    __syncthreads();
    dmma_contract< 1, true, true, NW >(shc, shdy, shs, wf);
    __syncthreads();
    dmma_contract< 2, true, true, NW >(shc, shdz, sht, wf);

    /* The contractions wrote shc through the generic proxy and the bulk store
       reads it through the async one; see tma_fence_shared() */
    tma_fence_shared();
    __syncthreads();

    if (tid == 0) {
      tma_store(cout[c] + ele, shc, CUBE_BYTES);
      tma_store_wait();
    }

    /* shc is restaged by the next component, and nothing may retire while the
       store is still reading it */
    __syncthreads();
  }
}

#endif // __CUDA_ARCH__ == sm_90 with a CUDA 12 toolkit

/*
 * Compile-time dispatch onto the TMA staged vector DMMA element kernel, see
 * the note on ax_helm_dmma_tma_dispatch above.
 */
template< typename T, const int LX, const int NW >
struct ax_helm_dmma_tma_vector_dispatch {
  __device__ static void run(T * __restrict__,
                             T * __restrict__,
                             T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__) { }
};

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 900) &&                        \
    (__CUDA_ARCH__ < 1000) && NEKO_TMA_TOOLKIT

/* Keep in sync with dmma_tma_vector_lx_supported() in dmma_tma_kernel.h */
#define NEKO_AX_HELM_DMMA_TMA_VECTOR_DISPATCH(LXV)                             \
  template< const int NW >                                                     \
  struct ax_helm_dmma_tma_vector_dispatch< double, LXV, NW > {                 \
    __device__ static void run(double * __restrict__ au,                       \
                               double * __restrict__ av,                       \
                               double * __restrict__ aw,                       \
                               const double * __restrict__ u,                  \
                               const double * __restrict__ v,                  \
                               const double * __restrict__ w,                  \
                               const double * __restrict__ dx,                 \
                               const double * __restrict__ dy,                 \
                               const double * __restrict__ dz,                 \
                               const double * __restrict__ h1,                 \
                               const double * __restrict__ g11,                \
                               const double * __restrict__ g22,                \
                               const double * __restrict__ g33,                \
                               const double * __restrict__ g12,                \
                               const double * __restrict__ g13,                \
                               const double * __restrict__ g23) {              \
      ax_helm_dmma_tma_vector_elem< LXV, NW >(au, av, aw, u, v, w,             \
                                              dx, dy, dz, h1,                  \
                                              g11, g22, g33,                   \
                                              g12, g13, g23);                  \
    }                                                                          \
  }

NEKO_AX_HELM_DMMA_TMA_VECTOR_DISPATCH(8);

#endif // __CUDA_ARCH__ == sm_90 with a CUDA 12 toolkit

template< typename T, const int LX, const int NW >
__global__ void NEKO_EB_BOUNDS(32 * NW)
ax_helm_kernel_dmma_tma_vector(T * __restrict__ au,
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
                               const T * __restrict__ g23) {

  ax_helm_dmma_tma_vector_dispatch< T, LX, NW >::run(au, av, aw, u, v, w,
                                                     dx, dy, dz, h1,
                                                     g11, g22, g33,
                                                     g12, g13, g23);
}


/**
 * Device kernel for the vector axhelm on the fp64 tensor cores, with the whole
 * element staged by the TMA engine in one batch
 *
 * The component-at-a-time variant above measured 4.9% behind the register
 * hoisting DMMA kernel at nw = 4 on GH200 while running *more* blocks per SM
 * than it (five against four, ptxas confirmed) -- so residency is not what
 * separates them. What separates them is how many bulk copies are outstanding.
 * The scalar TMA kernel issues eight at once, covering every byte the element
 * reads, and wins at 20 warps per SM against 48; the component-at-a-time
 * vector kernel issues seven once and then, for components 1 and 2, a lone
 * 4 kB copy that one thread issues and every other thread then waits on. A
 * 4 kB copy is small, its fixed cost only disappears when several are in
 * flight, and paying it alone three times an element is the deficit.
 *
 * So this variant issues all ten of an element's input copies at entry -- the
 * three components and the seven geometric factors -- and stores the three
 * results without waiting on each other. Each component keeps its own cube, in
 * and then out, which is what removes both serialisations at once: nothing has
 * to be restaged, so no load waits on a store.
 *
 * The cost is thirteen cubes rather than four, 54800 B, past the 48 kB a block
 * gets for free -- hence dynamic shared memory and the one-time opt-in in
 * ax_helm_dmma_tma_batch_optin() below. That lands at four blocks per SM at
 * nw = 4, which is exactly what the DMMA kernel it is competing with already
 * runs at, so the batching is bought with no occupancy at all.
 *
 * Only the first component's cube is waited on before work starts; the other
 * nine copies are waited on at the first pointwise step, and components 1 and
 * 2 never wait at all. See dmma_tma_kernel.h for the primitives, the struct
 * that fixes the layout, the sm_90 and toolkit guards and the lx bound.
 */

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 900) &&                        \
    (__CUDA_ARCH__ < 1000) && NEKO_TMA_TOOLKIT

template< const int LX, const int NW >
__device__ __forceinline__
void ax_helm_dmma_tma_batch_elem(double * __restrict__ au,
                                 double * __restrict__ av,
                                 double * __restrict__ aw,
                                 const double * __restrict__ u,
                                 const double * __restrict__ v,
                                 const double * __restrict__ w,
                                 const double * __restrict__ dx,
                                 const double * __restrict__ dy,
                                 const double * __restrict__ dz,
                                 const double * __restrict__ h1,
                                 const double * __restrict__ g11,
                                 const double * __restrict__ g22,
                                 const double * __restrict__ g33,
                                 const double * __restrict__ g12,
                                 const double * __restrict__ g13,
                                 const double * __restrict__ g23) {

  /* Keep in step with dmma_tma_batch_lx_supported() */
  static_assert(LX == DMMA_P,
                "the dmma tma batch variant stages whole cubes only");

  extern __shared__ __align__(128) char neko_ax_dyn_smem[];
  dmma_tma_batch_smem &sm =
    *reinterpret_cast< dmma_tma_batch_smem * >(neko_ax_dyn_smem);

  enum { CUBE_BYTES = DMMA_CUBE * (int) sizeof(double),
         /* v, w and the seven factors -- everything but the first cube */
         REST_BYTES = (DMMA_NG + 2) * CUBE_BYTES };

  const int nthrds = 32 * NW;
  const int tid = threadIdx.x;
  const int wf = tid >> 5;
  const int ele = blockIdx.x * DMMA_CUBE;

  if (tid == 0) {
    tma_barrier_init(&sm.bar_a, 1);
    tma_barrier_init(&sm.bar_b, 1);
  }

  /* At LX == DMMA_P the derivative matrix fills the staged one exactly, and
     there is no padding anywhere to zero */
  for (int p = tid; p < DMMA_MAT; p += nthrds) {
    sm.dx[p] = dx[p];
    sm.dy[p] = dy[p];
    sm.dz[p] = dz[p];
  }

  /* Both barriers initialised and both derivative matrices staged before
     anyone waits on the one or contracts with the other */
  __syncthreads();

  /* The whole element, ten copies, one batch. u leads because it is the only
     one anything waits for before the contractions start */
  if (tid == 0) {
    tma_expect(&sm.bar_a, CUBE_BYTES);
    tma_load(sm.c[0], u + ele, CUBE_BYTES, &sm.bar_a);

    tma_expect(&sm.bar_b, REST_BYTES);
    tma_load(sm.c[1], v + ele, CUBE_BYTES, &sm.bar_b);
    tma_load(sm.c[2], w + ele, CUBE_BYTES, &sm.bar_b);
    tma_load(sm.g[0], h1  + ele, CUBE_BYTES, &sm.bar_b);
    tma_load(sm.g[1], g11 + ele, CUBE_BYTES, &sm.bar_b);
    tma_load(sm.g[2], g22 + ele, CUBE_BYTES, &sm.bar_b);
    tma_load(sm.g[3], g33 + ele, CUBE_BYTES, &sm.bar_b);
    tma_load(sm.g[4], g12 + ele, CUBE_BYTES, &sm.bar_b);
    tma_load(sm.g[5], g13 + ele, CUBE_BYTES, &sm.bar_b);
    tma_load(sm.g[6], g23 + ele, CUBE_BYTES, &sm.bar_b);
  }

  /* The first component only. The other nine copies are still arriving */
  tma_wait(&sm.bar_a, 0);

  double * const cout[3] = { au, av, aw };

#pragma unroll
  for (int c = 0; c < 3; c++) {

    dmma_contract< 0, false, false, NW >(sm.r, sm.dx, sm.c[c], wf);
    dmma_contract< 1, false, false, NW >(sm.s, sm.dy, sm.c[c], wf);
    dmma_contract< 2, false, false, NW >(sm.t, sm.dz, sm.c[c], wf);

    __syncthreads();

    /* Everything else has landed by here; components 1 and 2 never wait */
    if (c == 0) {
      tma_wait(&sm.bar_b, 0);
    }

    for (int p = tid; p < DMMA_CUBE; p += nthrds) {
      const double H1  = sm.g[0][p];
      const double G00 = sm.g[1][p];
      const double G11 = sm.g[2][p];
      const double G22 = sm.g[3][p];
      const double G01 = sm.g[4][p];
      const double G02 = sm.g[5][p];
      const double G12 = sm.g[6][p];

      const double rtmp = sm.r[p];
      const double stmp = sm.s[p];
      const double ttmp = sm.t[p];

      sm.r[p] = H1
              * (G00 * rtmp
                 + G01 * stmp
                 + G02 * ttmp);
      sm.s[p] = H1
              * (G01 * rtmp
                 + G11 * stmp
                 + G12 * ttmp);
      sm.t[p] = H1
              * (G02 * rtmp
                 + G12 * stmp
                 + G22 * ttmp);
    }

    __syncthreads();

    /* The result overwrites the component's own staged input, which nothing
       reads again */
    dmma_contract< 0, true, false, NW >(sm.c[c], sm.dx, sm.r, wf);
    __syncthreads();
    dmma_contract< 1, true, true, NW >(sm.c[c], sm.dy, sm.s, wf);
    __syncthreads();
    dmma_contract< 2, true, true, NW >(sm.c[c], sm.dz, sm.t, wf);

    /* The contractions wrote sm.c[c] through the generic proxy and the bulk
       store reads it through the async one; see tma_fence_shared(). The
       barrier that follows also separates this component's last reads of
       r, s and t from the next one's writes to them */
    tma_fence_shared();
    __syncthreads();

    /* Issued and left uncommitted: the next component reads its own cube, so
       nothing here waits on this store. All three are committed together
       below */
    if (tid == 0) {
      tma_store(cout[c] + ele, sm.c[c], CUBE_BYTES);
    }
  }

  /* Commit the three stores as one group and wait for them to have read their
     cubes out; shared memory lives only as long as the block does */
  if (tid == 0) {
    tma_store_wait();
  }

  __syncthreads();
}

#endif // __CUDA_ARCH__ == sm_90 with a CUDA 12 toolkit

/*
 * Compile-time dispatch onto the batched TMA element kernel, see the note on
 * ax_helm_dmma_tma_dispatch above.
 */
template< typename T, const int LX, const int NW >
struct ax_helm_dmma_tma_batch_dispatch {
  __device__ static void run(T * __restrict__,
                             T * __restrict__,
                             T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__,
                             const T * __restrict__) { }
};

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 900) &&                        \
    (__CUDA_ARCH__ < 1000) && NEKO_TMA_TOOLKIT

/* Keep in sync with dmma_tma_batch_lx_supported() in dmma_tma_kernel.h */
#define NEKO_AX_HELM_DMMA_TMA_BATCH_DISPATCH(LXV)                              \
  template< const int NW >                                                     \
  struct ax_helm_dmma_tma_batch_dispatch< double, LXV, NW > {                  \
    __device__ static void run(double * __restrict__ au,                       \
                               double * __restrict__ av,                       \
                               double * __restrict__ aw,                       \
                               const double * __restrict__ u,                  \
                               const double * __restrict__ v,                  \
                               const double * __restrict__ w,                  \
                               const double * __restrict__ dx,                 \
                               const double * __restrict__ dy,                 \
                               const double * __restrict__ dz,                 \
                               const double * __restrict__ h1,                 \
                               const double * __restrict__ g11,                \
                               const double * __restrict__ g22,                \
                               const double * __restrict__ g33,                \
                               const double * __restrict__ g12,                \
                               const double * __restrict__ g13,                \
                               const double * __restrict__ g23) {              \
      ax_helm_dmma_tma_batch_elem< LXV, NW >(au, av, aw, u, v, w,              \
                                             dx, dy, dz, h1,                   \
                                             g11, g22, g33,                    \
                                             g12, g13, g23);                   \
    }                                                                          \
  }

NEKO_AX_HELM_DMMA_TMA_BATCH_DISPATCH(8);

#endif // __CUDA_ARCH__ == sm_90 with a CUDA 12 toolkit

template< typename T, const int LX, const int NW >
__global__ void NEKO_EB_BOUNDS(32 * NW)
ax_helm_kernel_dmma_tma_batch(T * __restrict__ au,
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
                              const T * __restrict__ g23) {

  ax_helm_dmma_tma_batch_dispatch< T, LX, NW >::run(au, av, aw, u, v, w,
                                                    dx, dy, dz, h1,
                                                    g11, g22, g33,
                                                    g12, g13, g23);
}

/*
 * Opt into the batched variant's dynamic allocation, once per specialisation.
 *
 * A block gets 48 kB of shared memory without asking; anything past that has
 * to be requested per kernel, and the carveout moved with it or the extra is
 * granted at the expense of nothing. Both are properties of the function, not
 * of the launch, so this is a function local static -- one flag per
 * <T, LX, NW> -- and the launch macro calls it ahead of every launch, where
 * after the first it is a single predictable branch.
 *
 * Returns false if the device refuses, which the tuner has already ruled out
 * via cuda_have_tma_batch(); the error is cleared rather than left to surface
 * against an unrelated CUDA_CHECK later.
 */
template< typename T, const int LX, const int NW >
static inline bool ax_helm_dmma_tma_batch_optin()
{
  static int state = -1;

  if (state < 0) {
    const void * const fn =
      (const void *) ax_helm_kernel_dmma_tma_batch< T, LX, NW >;
    cudaError_t err =
      cudaFuncSetAttribute(fn, cudaFuncAttributeMaxDynamicSharedMemorySize,
                           NEKO_DMMA_TMA_BATCH_SMEM);

    if (err == cudaSuccess) {
      err = cudaFuncSetAttribute(fn,
                                 cudaFuncAttributePreferredSharedMemoryCarveout,
                                 cudaSharedmemCarveoutMaxShared);
    }
    state = (err == cudaSuccess) ? 1 : 0;
    if (state == 0) {
      cudaGetLastError();
    }
  }
  return state == 1;
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
