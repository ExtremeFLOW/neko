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
