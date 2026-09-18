#ifndef __MATH_CDTP_KERNEL_H__
#define __MATH_CDTP_KERNEL_H__

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
                <= NEKO_EB_MAX_SMEM,
                "kstep block exceeds the shared memory budget");

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
 * Device kernel for D^T x on the fp64 tensor cores
 *
 * The only operator here whose contractions come *after* its pointwise work
 * rather than before it, and the only one that accumulates: the three weighted
 * fields tar, tas and tat are formed first and then contracted with the
 * transposed derivative matrices into one output. That is axhelm's phase 2
 * rather than its phase 1, and it is what makes cdtp a different test of the
 * strategy from dudxyz, opgrad and conv1 rather than a fourth copy of one.
 *
 * Because the operator is handed dxt, dyt and dzt already transposed, the
 * contraction index pattern is the ordinary forward one -- dxt[i + l*LX]
 * against tar[l,j,k] -- so dmma_contract is used with TRANSPOSE = false and
 * the staged transposes, not with TRANSPOSE = true and the forward matrices.
 * The accumulation is ACCUM: axis 0 initialises every tile of the output and
 * axes 1 and 2 add into it.
 *
 * Four cubes, 17920 B, the same as every other DMMA variant here. w3 is a
 * quadrature weight shared by every element, so it is indexed by the point
 * within the element (dmma_idx::l) rather than the global offset, as in
 * opgrad_dmma_elem().
 *
 * Expectation, recorded before measuring: cdtp carries 6827 bytes per
 * contraction at lx = 8, *below* dudxyz's 8192 and conv1's 20480, both of
 * which measured as ties against their incumbents on GH200 because they were
 * already at the bandwidth roof. This one should tie as well. It is offered
 * because the shape is untested, not because the ranking predicts a win.
 */

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 800) && (__CUDA_ARCH__ < 1000)

template< const int LX, const int NW >
__device__ __forceinline__
void cdtp_dmma_elem(double * __restrict__ dtx,
                    const double * __restrict__ x,
                    const double * __restrict__ dr,
                    const double * __restrict__ ds,
                    const double * __restrict__ dt,
                    const double * __restrict__ dxt,
                    const double * __restrict__ dyt,
                    const double * __restrict__ dzt,
                    const double * __restrict__ w3,
                    const int nelv) {

  /* Element independent, one copy per block. Block diagonal when more than
     one element is packed: one LX x LX copy of D^T per sub-cube */
  __shared__ __align__(16) double shdxt[DMMA_MAT];
  __shared__ __align__(16) double shdyt[DMMA_MAT];
  __shared__ __align__(16) double shdzt[DMMA_MAT];

  /* The pack, padded to DMMA_P^3: the three weighted fields and the
     accumulated result */
  __shared__ __align__(16) double shtar[DMMA_CUBE];
  __shared__ __align__(16) double shtas[DMMA_CUBE];
  __shared__ __align__(16) double shtat[DMMA_CUBE];
  __shared__ __align__(16) double shout[DMMA_CUBE];

  static_assert(sizeof(shdxt) +
                sizeof(shdyt) +
                sizeof(shdzt) +
                sizeof(shtar) +
                sizeof(shtas) +
                sizeof(shtat) +
                sizeof(shout)
                <= NEKO_EB_MAX_SMEM,
                "dmma block exceeds the shared memory budget");

  /* Elements per sub-cube axis, and per cube, see the note in dmma_kernel.h */
  enum { PPA = (DMMA_P % LX == 0) ? (DMMA_P / LX) : 1,
         PACK = PPA * PPA * PPA,
         LX3 = LX * LX * LX,
         NP = PACK * LX3 };

  typedef dmma_pack< LX, PPA > pack;

  const int nthrds = 32 * NW;
  const int tid = threadIdx.x;
  const int wf = tid >> 5;
  const int ebase = pack::ebase();

  /* The padding has to be finite, see the note in dmma_kernel.h. All three
     inputs are contracted, so all three are zero filled -- unlike the phase 1
     kernels, which stage a single field. shout needs none: the axis 0
     contraction initialises every tile of it */
  if (PACK * LX3 < DMMA_CUBE) {
    for (int p = tid; p < DMMA_CUBE; p += nthrds) {
      shtar[p] = 0.0;
      shtas[p] = 0.0;
      shtat[p] = 0.0;
    }
  }
  if (LX < DMMA_P) {
    for (int p = tid; p < DMMA_MAT; p += nthrds) {
      shdxt[p] = 0.0;
      shdyt[p] = 0.0;
      shdzt[p] = 0.0;
    }
  }
  if (LX < DMMA_P) {
    __syncthreads();
  }

  /* One copy of D^T per sub-cube, down the diagonal */
  for (int p = tid; p < LX * LX; p += nthrds) {
    const int i = p % LX;
    const int l = p / LX;
#pragma unroll
    for (int b = 0; b < PPA; b++) {
      const int m = (b * LX + i) + DMMA_P * (b * LX + l);
      shdxt[m] = dxt[p];
      shdyt[m] = dyt[p];
      shdzt[m] = dzt[p];
    }
  }

  for (int p = tid; p < NP; p += nthrds) {
    const dmma_idx idx = pack::map(p, ebase, nelv);
    const int c = idx.c;
    const int gp = idx.g;
    const double wx = x[gp] * w3[idx.l];

    shtar[c] = wx * dr[gp];
    shtas[c] = wx * ds[gp];
    shtat[c] = wx * dt[gp];
  }

  __syncthreads();

  /* Axis 0 initialises the output, the other two add into it. The barriers
     are needed because the j slab tiles of axis 2 span what every warp wrote
     along axes 0 and 1; see the same sequence in ax_helm_dmma_elem() */
  dmma_contract< 0, false, false, NW >(shout, shdxt, shtar, wf);
  __syncthreads();
  dmma_contract< 1, false, true, NW >(shout, shdyt, shtas, wf);
  __syncthreads();
  dmma_contract< 2, false, true, NW >(shout, shdzt, shtat, wf);
  __syncthreads();

  for (int p = tid; p < NP; p += nthrds) {
    const dmma_idx idx = pack::map(p, ebase, nelv);

    if (idx.live) {
      dtx[idx.g] = shout[idx.c];
    }
  }
}

#endif // __CUDA_ARCH__ in [800, 1000)

/*
 * Compile-time dispatch onto the DMMA element kernel. The launch macros in
 * opr_cdtp.cu are written for every LX the operator dispatches and for
 * whatever `real` is, so every combination has to compile; the ones the
 * strategy does not cover -- single precision, LX outside the supported
 * range, a build without an fp64 tensor core arch -- resolve to this no-op.
 * The autotuner never selects the strategy for them, see dmma_lx_supported()
 * and cuda_have_dmma() in dmma_kernel.h.
 */
template< typename T, const int LX, const int NW >
struct cdtp_dmma_dispatch {
  __device__ static void run(T * __restrict__,
                             const T * __restrict__, const T * __restrict__,
                             const T * __restrict__, const T * __restrict__,
                             const T * __restrict__, const T * __restrict__,
                             const T * __restrict__, const T * __restrict__,
                             const int) { }
};

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 800) && (__CUDA_ARCH__ < 1000)

/* Keep in sync with dmma_lx_supported() in dmma_kernel.h */
#define NEKO_CDTP_DMMA_DISPATCH(LXV)                                           \
  template< const int NW >                                                     \
  struct cdtp_dmma_dispatch< double, LXV, NW > {                               \
    __device__ static void run(double * __restrict__ dtx,                      \
                               const double * __restrict__ x,                  \
                               const double * __restrict__ dr,                 \
                               const double * __restrict__ ds,                 \
                               const double * __restrict__ dt,                 \
                               const double * __restrict__ dxt,                \
                               const double * __restrict__ dyt,                \
                               const double * __restrict__ dzt,                \
                               const double * __restrict__ w3,                 \
                               const int nelv) {                               \
      cdtp_dmma_elem< LXV, NW >(dtx, x, dr, ds, dt,                            \
                                dxt, dyt, dzt, w3, nelv);                      \
    }                                                                          \
  }

NEKO_CDTP_DMMA_DISPATCH(2);
NEKO_CDTP_DMMA_DISPATCH(3);
NEKO_CDTP_DMMA_DISPATCH(4);
NEKO_CDTP_DMMA_DISPATCH(5);
NEKO_CDTP_DMMA_DISPATCH(6);
NEKO_CDTP_DMMA_DISPATCH(7);
NEKO_CDTP_DMMA_DISPATCH(8);

#endif // __CUDA_ARCH__ in [800, 1000)

template< typename T, const int LX, const int NW >
__global__ void NEKO_EB_BOUNDS(32 * NW)
cdtp_kernel_dmma(T * __restrict__ dtx,
                 const T * __restrict__ x,
                 const T * __restrict__ dr,
                 const T * __restrict__ ds,
                 const T * __restrict__ dt,
                 const T * __restrict__ dxt,
                 const T * __restrict__ dyt,
                 const T * __restrict__ dzt,
                 const T * __restrict__ w3,
                 const int nelv) {

  cdtp_dmma_dispatch< T, LX, NW >::run(dtx, x, dr, ds, dt,
                                       dxt, dyt, dzt, w3, nelv);
}

/**
 * Device kernel for D^T x on the fp64 tensor cores, with the element staged by
 * the TMA engine
 *
 * Same three accumulating contractions and the same cube layout as
 * cdtp_dmma_elem(); at lx == DMMA_P the padding is empty and the cube offset
 * is the flat point index, so the arithmetic is identical.
 *
 * **The overlap has to be arranged differently here, and this is the one
 * interesting thing about the variant.** Every other TMA kernel in the tree
 * contracts a field that arrives first and consumes its factors afterwards, so
 * the factor copies fly underneath three contractions for free. cdtp needs x
 * *and* a factor before it can form anything at all, and needs all three
 * factors before the last contraction -- a single wait on everything would
 * leave the copies overlapping nothing.
 *
 * So the four copies are all issued at entry but waited in two stages: x and
 * dr on the first barrier, ds and dt on the second, which is waited only after
 * the axis 0 contraction has been issued. Half the factor traffic flies
 * underneath a contraction; the other half cannot, by the shape of the
 * operator.
 *
 * That costs one cube over the phase 1 kernels -- the weighted fields
 * overwrite their own staged factors, but x has to outlive all three of them,
 * so the result needs its own -- for **22032 B and 10 blocks per SM**, the
 * smallest TMA block and the best occupancy of any of these variants, and the
 * only one that stays under the 48 kB a block gets without opting in.
 */

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 900) &&                        \
    (__CUDA_ARCH__ < 1000) && NEKO_TMA_TOOLKIT

template< const int LX, const int NW >
__device__ __forceinline__
void cdtp_dmma_tma_elem(double * __restrict__ dtx,
                        const double * __restrict__ x,
                        const double * __restrict__ dr,
                        const double * __restrict__ ds,
                        const double * __restrict__ dt,
                        const double * __restrict__ dxt,
                        const double * __restrict__ dyt,
                        const double * __restrict__ dzt,
                        const double * __restrict__ w3) {

  /* A bulk copy needs 16 byte alignment at both ends; the cubes are given the
     128 the tensor variants would want anyway */
  __shared__ __align__(128) double shdxt[DMMA_MAT];
  __shared__ __align__(128) double shdyt[DMMA_MAT];
  __shared__ __align__(128) double shdzt[DMMA_MAT];

  /* x, and the three factors which become the weighted fields in place */
  __shared__ __align__(128) double shx[DMMA_CUBE];
  __shared__ __align__(128) double shtar[DMMA_CUBE];
  __shared__ __align__(128) double shtas[DMMA_CUBE];
  __shared__ __align__(128) double shtat[DMMA_CUBE];
  __shared__ __align__(128) double shout[DMMA_CUBE];

  __shared__ __align__(8) unsigned long long bar_a;
  __shared__ __align__(8) unsigned long long bar_b;

  static_assert(sizeof(shdxt) +
                sizeof(shdyt) +
                sizeof(shdzt) +
                sizeof(shx) +
                sizeof(shtar) +
                sizeof(shtas) +
                sizeof(shtat) +
                sizeof(shout) +
                sizeof(bar_a) +
                sizeof(bar_b)
                <= NEKO_EB_MAX_SMEM,
                "cdtp dmma tma block exceeds the shared memory budget");

  /* Only lx == DMMA_P stages as a contiguous run of bytes, which is what a
     bulk copy moves; see the scope note in dmma_tma_kernel.h. Keep in step
     with dmma_tma_cdtp_lx_supported() */
  static_assert(LX == DMMA_P,
                "the dmma tma variant stages whole cubes only");

  enum { CUBE_BYTES = DMMA_CUBE * (int) sizeof(double) };

  const int nthrds = 32 * NW;
  const int tid = threadIdx.x;
  const int wf = tid >> 5;
  const int ebase = blockIdx.x * DMMA_CUBE;

  if (tid == 0) {
    tma_barrier_init(&bar_a, 1);
    tma_barrier_init(&bar_b, 1);
  }

  /* At LX == DMMA_P the derivative matrix fills the staged one exactly, and
     there is no padding anywhere to zero */
  for (int p = tid; p < DMMA_MAT; p += nthrds) {
    shdxt[p] = dxt[p];
    shdyt[p] = dyt[p];
    shdzt[p] = dzt[p];
  }

  /* Both barriers initialised and every derivative matrix staged before
     anyone waits on the one or contracts with the other */
  __syncthreads();

  /* All four copies issued together; only the first two are waited on now */
  if (tid == 0) {
    tma_expect(&bar_a, 2 * CUBE_BYTES);
    tma_load(shx,   x  + ebase, CUBE_BYTES, &bar_a);
    tma_load(shtar, dr + ebase, CUBE_BYTES, &bar_a);

    tma_expect(&bar_b, 2 * CUBE_BYTES);
    tma_load(shtas, ds + ebase, CUBE_BYTES, &bar_b);
    tma_load(shtat, dt + ebase, CUBE_BYTES, &bar_b);
  }

  tma_wait(&bar_a, 0);

  /* At LX == DMMA_P the cube offset, the point within the element and the
     loop counter are all the same index, so w3 is read at p. The weighted
     field overwrites the staged factor it was formed from */
  for (int p = tid; p < DMMA_CUBE; p += nthrds) {
    shtar[p] *= shx[p] * w3[p];
  }

  __syncthreads();

  dmma_contract< 0, false, false, NW >(shout, shdxt, shtar, wf);

  /* ds and dt have been in flight underneath that contraction */
  tma_wait(&bar_b, 0);

  for (int p = tid; p < DMMA_CUBE; p += nthrds) {
    const double wx = shx[p] * w3[p];

    shtas[p] *= wx;
    shtat[p] *= wx;
  }

  __syncthreads();

  dmma_contract< 1, false, true, NW >(shout, shdyt, shtas, wf);
  __syncthreads();
  dmma_contract< 2, false, true, NW >(shout, shdzt, shtat, wf);

  /* The contractions wrote shout through the generic proxy and the bulk store
     reads it through the async one, so the fence is needed on top of the
     barrier; see tma_fence_shared() */
  tma_fence_shared();
  __syncthreads();

  if (tid == 0) {
    tma_store(dtx + ebase, shout, CUBE_BYTES);
    tma_store_wait();
  }

  /* The store reads shout asynchronously, and shared memory lives only as long
     as the block does: nothing may retire until it has been read out */
  __syncthreads();
}

#endif // __CUDA_ARCH__ == sm_90 with a CUDA 12 toolkit

/*
 * Compile-time dispatch onto the TMA staged DMMA element kernel, see the note
 * on cdtp_dmma_dispatch above. The no-op covers everything the strategy does
 * not: single precision, any lx but DMMA_P, a build without sm_90, and a
 * toolkit older than CUDA 12.
 */
template< typename T, const int LX, const int NW >
struct cdtp_dmma_tma_dispatch {
  __device__ static void run(T * __restrict__,
                             const T * __restrict__, const T * __restrict__,
                             const T * __restrict__, const T * __restrict__,
                             const T * __restrict__, const T * __restrict__,
                             const T * __restrict__, const T * __restrict__) { }
};

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 900) &&                        \
    (__CUDA_ARCH__ < 1000) && NEKO_TMA_TOOLKIT

/* Keep in sync with dmma_tma_cdtp_lx_supported() in dmma_tma_kernel.h */
#define NEKO_CDTP_DMMA_TMA_DISPATCH(LXV)                                       \
  template< const int NW >                                                     \
  struct cdtp_dmma_tma_dispatch< double, LXV, NW > {                           \
    __device__ static void run(double * __restrict__ dtx,                      \
                               const double * __restrict__ x,                  \
                               const double * __restrict__ dr,                 \
                               const double * __restrict__ ds,                 \
                               const double * __restrict__ dt,                 \
                               const double * __restrict__ dxt,                \
                               const double * __restrict__ dyt,                \
                               const double * __restrict__ dzt,                \
                               const double * __restrict__ w3) {               \
      cdtp_dmma_tma_elem< LXV, NW >(dtx, x, dr, ds, dt,                        \
                                    dxt, dyt, dzt, w3);                        \
    }                                                                          \
  }

NEKO_CDTP_DMMA_TMA_DISPATCH(8);

#endif // __CUDA_ARCH__ == sm_90 with a CUDA 12 toolkit

template< typename T, const int LX, const int NW >
__global__ void NEKO_EB_BOUNDS(32 * NW)
cdtp_kernel_dmma_tma(T * __restrict__ dtx,
                     const T * __restrict__ x,
                     const T * __restrict__ dr,
                     const T * __restrict__ ds,
                     const T * __restrict__ dt,
                     const T * __restrict__ dxt,
                     const T * __restrict__ dyt,
                     const T * __restrict__ dzt,
                     const T * __restrict__ w3) {

  cdtp_dmma_tma_dispatch< T, LX, NW >::run(dtx, x, dr, ds, dt,
                                           dxt, dyt, dzt, w3);
}

#endif // __MATH_CDTP_KERNEL_H__
