#ifndef __MATH_DUDXYZ_KERNEL_H__
#define __MATH_DUDXYZ_KERNEL_H__

#include "elem_block.h"
#include "dmma_kernel.h"
#include "dmma_tma_kernel.h"
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

/**
 * Device kernel for derivative 
 */

template< typename T, const int LX, const int CHUNKS >
__global__ void dudxyz_kernel_1d(T * __restrict__ du,
                                 const T * __restrict__ u,
                                 const T * __restrict__ dr,
                                 const T * __restrict__ ds,
                                 const T * __restrict__ dt,
                                 const T * __restrict__ dx,
                                 const T * __restrict__ dy,
                                 const T * __restrict__ dz,
                                 const T * __restrict__ jacinv) { 
  
  __shared__ T shu[LX * LX * LX];
  __shared__ T shdr[LX * LX * LX];
  __shared__ T shds[LX * LX * LX];
  __shared__ T shdt[LX * LX * LX];

  __shared__ T shdx[LX * LX];
  __shared__ T shdy[LX * LX];
  __shared__ T shdz[LX * LX];
  
  __shared__ T shjacinv[LX * LX * LX];
  
  const int e = blockIdx.x;
  const int iii = threadIdx.x;
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;

  if (iii < (LX * LX)) {
    shdx[iii] = dx[iii];
    shdy[iii] = dy[iii];
    shdz[iii] = dz[iii];
  }

  int l = iii;
  while(l < (LX * LX * LX)) {
    shu[l] = u[l + e * LX * LX * LX];
    shdr[l] = dr[l + e * LX * LX * LX];
    shds[l] = ds[l + e * LX * LX * LX];
    shdt[l] = dt[l + e * LX * LX * LX];
    shjacinv[l] = jacinv[l + e * LX * LX * LX];
    l = l + CHUNKS;
  }

  __syncthreads();
  
  for (int n = 0; n < nchunks; n++) {
    const int ijk = iii + n * CHUNKS;
    const int jk = ijk / LX;
    const int i = ijk - jk * LX;
    const int k = jk / LX;
    const int j = jk - k * LX;
    if ( i < LX && j < LX && k < LX) {
      T rtmp = 0.0;
      T stmp = 0.0;
      T ttmp = 0.0;
      for (int l = 0; l < LX; l++) {
	rtmp += shdx[i + l * LX] * shu[l + j * LX + k * LX * LX];
	stmp += shdy[j + l * LX] * shu[i + l * LX + k * LX * LX];
	ttmp += shdz[k + l * LX] * shu[i + j * LX + l * LX * LX];
      }
      du[ijk + e * LX * LX * LX] = ((rtmp * shdr[ijk])
				    + (stmp * shds[ijk])
				    + (ttmp * shdt[ijk]))
	                           * shjacinv[ijk];
      
    }
  }
}

template< typename T, const int LX, const int EB >
__global__ void NEKO_EB_BOUNDS(LX*LX*EB)
  dudxyz_kernel_kstep(T * __restrict__ du,
                      const T * __restrict__ u,
                      const T * __restrict__ dr,
                      const T * __restrict__ ds,
                      const T * __restrict__ dt,
                      const T * __restrict__ dx,
                      const T * __restrict__ dy,
                      const T * __restrict__ dz,
                      const T * __restrict__ jacinv,
                    const int nelv) { 
  
  __shared__ T shu[EB * LX * LX];

  __shared__ T shdx[LX * LX];
  __shared__ T shdy[LX * LX];
  __shared__ T shdz[LX * LX];

  static_assert(sizeof(shu) +
                sizeof(shdx) +
                sizeof(shdy) +
                sizeof(shdz)
                <= NEKO_EB_MAX_SMEM,
                "kstep block exceeds the shared memory budget");

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
    shdx[ij] = dx[ij];
    shdy[ij] = dy[ij];
    shdz[ij] = dz[ij];
  }
  
  T ru[LX];
  T rdr[LX];
  T rds[LX];
  T rdt[LX];  
  T rjacinv[LX];

  #pragma unroll LX
  for (int k = 0; k < LX; ++k) {
    ru[k] = u[ij + k*LX*LX + ele];
    rdr[k] = dr[ij + k*LX*LX + ele];
    rds[k] = ds[ij + k*LX*LX + ele];
    rdt[k] = dt[ij + k*LX*LX + ele];
    rjacinv[k] = jacinv[ij + k*LX*LX + ele];
  }
    
  __syncthreads();

  #pragma unroll
  for (int k = 0; k < LX; ++k) {
    const int ijk = ij + k*LX*LX;
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
      du[ijk + ele] = rjacinv[k] * ((rtmp * rdr[k])
                                    + (stmp * rds[k])
                                    + (ttmp * rdt[k]));
    }
    __syncthreads();
  }
}


/**
 * Device kernel for the derivative on the fp64 tensor cores
 *
 * One element per block -- or a packed cube of them at low order -- NW warps,
 * and the element staged into shared memory as four DMMA_P^3 cubes: the
 * staged input and the three reference derivatives. This is the first half of
 * ax_helm_dmma_elem() and nothing else: the same three D * U GEMMs handed to
 * dmma_contract(), but the geometric factors turn the reference derivatives
 * into the physical one pointwise and it goes straight back out to global
 * memory, where axhelm contracts a second time with D^T. See dmma_kernel.h
 * for the padded staging, the per axis matrix views and the arch and LX
 * bounds.
 *
 * Half the contractions for two thirds of the bytes is the reason to expect
 * anything here. What won axhelm was not the tensor cores -- it reached 82%
 * of the achievable streaming bandwidth on GH200 at 7.8% of the fp64 tensor
 * peak -- but the access pattern the staging buys, and the cost of that
 * staging is paid per contraction while the memory time it hides behind is
 * paid per byte. At lx = 8 this operator has 8192 bytes per contraction
 * against axhelm's 6144, so the trade is 33% better here, not worse.
 */

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 800) && (__CUDA_ARCH__ < 1000)

template< const int LX, const int NW >
__device__ __forceinline__
void dudxyz_dmma_elem(double * __restrict__ du,
                      const double * __restrict__ u,
                      const double * __restrict__ dr,
                      const double * __restrict__ ds,
                      const double * __restrict__ dt,
                      const double * __restrict__ dx,
                      const double * __restrict__ dy,
                      const double * __restrict__ dz,
                      const double * __restrict__ jacinv,
                      const int nelv) {

  /* Element independent, one copy per block. Block diagonal when more than
     one element is packed: one LX x LX copy of D per sub-cube */
  __shared__ __align__(16) double shdx[DMMA_MAT];
  __shared__ __align__(16) double shdy[DMMA_MAT];
  __shared__ __align__(16) double shdz[DMMA_MAT];

  /* The pack, padded to DMMA_P^3. shu carries the input, shr, shs and sht
     the reference derivatives */
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

  /* The padding has to be finite, see the note in dmma_kernel.h. At
     LX == DMMA_P with one element packed there is none and this is folded
     away */
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

  /* No second set of contractions and so no restaging: the metrics are
     applied on the way out. A dead tail slot reads nothing at all here,
     unlike the axhelm kernel where the clamped read still has a shared write
     to feed; at PACK == 1 the guard is a compile time true and folds */
  for (int p = tid; p < NP; p += nthrds) {
    const dmma_idx x = pack::map(p, ebase, nelv);

    if (x.live) {
      const int c = x.c;
      const int gp = x.g;

      du[gp] = jacinv[gp] * ((shr[c] * dr[gp])
                             + (shs[c] * ds[gp])
                             + (sht[c] * dt[gp]));
    }
  }
}

#endif // __CUDA_ARCH__ in [800, 1000)

/*
 * Compile-time dispatch onto the DMMA element kernel. The launch macros in
 * opr_dudxyz.cu are written for every LX the operator dispatches and for
 * whatever `real` is, so every combination has to compile; the ones the
 * strategy does not cover -- single precision, LX outside the supported
 * range, a build without an fp64 tensor core arch -- resolve to this no-op.
 * The autotuner never selects the strategy for them, so the no-op is
 * unreachable at runtime, see dmma_lx_supported() and cuda_have_dmma() in
 * dmma_kernel.h.
 */
template< typename T, const int LX, const int NW >
struct dudxyz_dmma_dispatch {
  __device__ static void run(T * __restrict__,
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
#define NEKO_DUDXYZ_DMMA_DISPATCH(LXV)                                         \
  template< const int NW >                                                     \
  struct dudxyz_dmma_dispatch< double, LXV, NW > {                             \
    __device__ static void run(double * __restrict__ du,                       \
                               const double * __restrict__ u,                  \
                               const double * __restrict__ dr,                 \
                               const double * __restrict__ ds,                 \
                               const double * __restrict__ dt,                 \
                               const double * __restrict__ dx,                 \
                               const double * __restrict__ dy,                 \
                               const double * __restrict__ dz,                 \
                               const double * __restrict__ jacinv,             \
                               const int nelv) {                               \
      dudxyz_dmma_elem< LXV, NW >(du, u, dr, ds, dt,                           \
                                  dx, dy, dz, jacinv, nelv);                   \
    }                                                                          \
  }

NEKO_DUDXYZ_DMMA_DISPATCH(2);
NEKO_DUDXYZ_DMMA_DISPATCH(3);
NEKO_DUDXYZ_DMMA_DISPATCH(4);
NEKO_DUDXYZ_DMMA_DISPATCH(5);
NEKO_DUDXYZ_DMMA_DISPATCH(6);
NEKO_DUDXYZ_DMMA_DISPATCH(7);
NEKO_DUDXYZ_DMMA_DISPATCH(8);

#endif // __CUDA_ARCH__ in [800, 1000)

template< typename T, const int LX, const int NW >
__global__ void NEKO_EB_BOUNDS(32 * NW)
dudxyz_kernel_dmma(T * __restrict__ du,
                   const T * __restrict__ u,
                   const T * __restrict__ dr,
                   const T * __restrict__ ds,
                   const T * __restrict__ dt,
                   const T * __restrict__ dx,
                   const T * __restrict__ dy,
                   const T * __restrict__ dz,
                   const T * __restrict__ jacinv,
                   const int nelv) {

  dudxyz_dmma_dispatch< T, LX, NW >::run(du, u, dr, ds, dt,
                                         dx, dy, dz, jacinv, nelv);
}

/**
 * Device kernel for the derivative on the fp64 tensor cores, with the element
 * staged by the TMA engine
 *
 * Same three contractions and the same cube layout as dudxyz_dmma_elem(), and
 * at lx == DMMA_P the padding is empty and the cube offset is the flat point
 * index, so the arithmetic is identical -- only how the six cubes get in and
 * out of shared memory differs.
 *
 * u arrives as one bulk copy on its own mbarrier and the four factors dr, ds,
 * dt and jacinv as four more on a second, waited on only at the pointwise
 * step, so they are in flight underneath the three contractions. That is 80%
 * of the read traffic moved out from between two barriers, the same share the
 * axhelm variant moves. Five copies issued together is a smaller batch than
 * that kernel's eight but the same shape -- every one of the element's reads
 * in one issue -- which is what separates the scalar axhelm variant that won
 * from the component-at-a-time vector one that did not.
 *
 * The result is formed in the staged input, dead once the contractions are
 * done, and leaves as a single bulk store. Nothing is contracted back, so the
 * block is 33.8 kB against the axhelm variant's 45.5 -- six resident blocks
 * per SM rather than five.
 *
 * See dmma_tma_kernel.h for the primitives, the sm_90 and toolkit guards and
 * the lx == DMMA_P bound.
 */

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 900) &&                        \
    (__CUDA_ARCH__ < 1000) && NEKO_TMA_TOOLKIT

template< const int LX, const int NW >
__device__ __forceinline__
void dudxyz_dmma_tma_elem(double * __restrict__ du,
                          const double * __restrict__ u,
                          const double * __restrict__ dr,
                          const double * __restrict__ ds,
                          const double * __restrict__ dt,
                          const double * __restrict__ dx,
                          const double * __restrict__ dy,
                          const double * __restrict__ dz,
                          const double * __restrict__ jacinv) {

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

  /* dr, ds, dt and jacinv */
  __shared__ __align__(128) double shgr[DMMA_CUBE];
  __shared__ __align__(128) double shgs[DMMA_CUBE];
  __shared__ __align__(128) double shgt[DMMA_CUBE];
  __shared__ __align__(128) double shgj[DMMA_CUBE];

  __shared__ __align__(8) unsigned long long bar_u;
  __shared__ __align__(8) unsigned long long bar_g;

  static_assert(sizeof(shdx) +
                sizeof(shdy) +
                sizeof(shdz) +
                sizeof(shu) +
                sizeof(shr) +
                sizeof(shs) +
                sizeof(sht) +
                sizeof(shgr) +
                sizeof(shgs) +
                sizeof(shgt) +
                sizeof(shgj) +
                sizeof(bar_u) +
                sizeof(bar_g)
                <= NEKO_EB_MAX_SMEM,
                "dudxyz dmma tma block exceeds the shared memory budget");

  /* Only lx == DMMA_P stages as a contiguous run of bytes, which is what a
     bulk copy moves; see the scope note in dmma_tma_kernel.h. Keep in step
     with dmma_tma_dudxyz_lx_supported() */
  static_assert(LX == DMMA_P,
                "the dmma tma variant stages whole cubes only");

  enum { CUBE_BYTES = DMMA_CUBE * (int) sizeof(double),
         NG = 4 };

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

  /* Both barriers initialised and every derivative matrix staged before
     anyone waits on the one or contracts with the other */
  __syncthreads();

  if (tid == 0) {
    tma_expect(&bar_u, CUBE_BYTES);
    tma_load(shu, u + ebase, CUBE_BYTES, &bar_u);

    tma_expect(&bar_g, NG * CUBE_BYTES);
    tma_load(shgr, dr     + ebase, CUBE_BYTES, &bar_g);
    tma_load(shgs, ds     + ebase, CUBE_BYTES, &bar_g);
    tma_load(shgt, dt     + ebase, CUBE_BYTES, &bar_g);
    tma_load(shgj, jacinv + ebase, CUBE_BYTES, &bar_g);
  }

  /* u only. The four factor cubes are still arriving */
  tma_wait(&bar_u, 0);

  dmma_contract< 0, false, false, NW >(shr, shdx, shu, wf);
  dmma_contract< 1, false, false, NW >(shs, shdy, shu, wf);
  dmma_contract< 2, false, false, NW >(sht, shdz, shu, wf);

  __syncthreads();

  tma_wait(&bar_g, 0);

  /* The staged input is dead once the contractions have been read out of it,
     which the barrier above guarantees, so the result is formed in it and
     leaves as one bulk store rather than DMMA_CUBE scalar ones */
  for (int p = tid; p < DMMA_CUBE; p += nthrds) {
    shu[p] = shgj[p] * ((shr[p] * shgr[p])
                        + (shs[p] * shgs[p])
                        + (sht[p] * shgt[p]));
  }

  /* The pointwise pass wrote shu through the generic proxy and the bulk store
     reads it through the async one, so the fence is needed on top of the
     barrier; see tma_fence_shared() */
  tma_fence_shared();
  __syncthreads();

  if (tid == 0) {
    tma_store(du + ebase, shu, CUBE_BYTES);
    tma_store_wait();
  }

  /* The store reads shu asynchronously, and shared memory lives only as long
     as the block does: nothing may retire until it has been read out */
  __syncthreads();
}

#endif // __CUDA_ARCH__ == sm_90 with a CUDA 12 toolkit

/*
 * Compile-time dispatch onto the TMA staged DMMA element kernel, see the note
 * on dudxyz_dmma_dispatch above. The no-op covers everything the strategy
 * does not: single precision, any lx but DMMA_P, a build without sm_90, and a
 * toolkit older than CUDA 12.
 */
template< typename T, const int LX, const int NW >
struct dudxyz_dmma_tma_dispatch {
  __device__ static void run(T * __restrict__,
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

/* Keep in sync with dmma_tma_dudxyz_lx_supported() in dmma_tma_kernel.h */
#define NEKO_DUDXYZ_DMMA_TMA_DISPATCH(LXV)                                     \
  template< const int NW >                                                     \
  struct dudxyz_dmma_tma_dispatch< double, LXV, NW > {                         \
    __device__ static void run(double * __restrict__ du,                       \
                               const double * __restrict__ u,                  \
                               const double * __restrict__ dr,                 \
                               const double * __restrict__ ds,                 \
                               const double * __restrict__ dt,                 \
                               const double * __restrict__ dx,                 \
                               const double * __restrict__ dy,                 \
                               const double * __restrict__ dz,                 \
                               const double * __restrict__ jacinv) {           \
      dudxyz_dmma_tma_elem< LXV, NW >(du, u, dr, ds, dt,                       \
                                      dx, dy, dz, jacinv);                     \
    }                                                                          \
  }

NEKO_DUDXYZ_DMMA_TMA_DISPATCH(8);

#endif // __CUDA_ARCH__ == sm_90 with a CUDA 12 toolkit

template< typename T, const int LX, const int NW >
__global__ void NEKO_EB_BOUNDS(32 * NW)
dudxyz_kernel_dmma_tma(T * __restrict__ du,
                       const T * __restrict__ u,
                       const T * __restrict__ dr,
                       const T * __restrict__ ds,
                       const T * __restrict__ dt,
                       const T * __restrict__ dx,
                       const T * __restrict__ dy,
                       const T * __restrict__ dz,
                       const T * __restrict__ jacinv) {

  dudxyz_dmma_tma_dispatch< T, LX, NW >::run(du, u, dr, ds, dt,
                                             dx, dy, dz, jacinv);
}

#endif // __MATH_DUDXYZ_KERNEL_H__
