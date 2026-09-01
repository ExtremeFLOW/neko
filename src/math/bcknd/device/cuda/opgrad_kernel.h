#ifndef __MATH_OPGRAD_KERNEL_H__
#define __MATH_OPGRAD_KERNEL_H__

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
                <= NEKO_EB_MAX_SMEM,
                "kstep block exceeds the shared memory budget");

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
 * Device kernel for the weak gradient on the fp64 tensor cores
 *
 * Structurally this is dudxyz_dmma_elem() with a wider tail: the same four
 * DMMA_P^3 cubes, the same stage-u-then-contract-three-axes, and then three
 * outputs formed pointwise from nine streamed metric cubes instead of one from
 * three. See dmma_kernel.h for the padded staging, the per axis matrix views
 * and the arch and LX bounds.
 *
 * This is the best placed scalar operator in the tree for the strategy. What
 * won axhelm was the access pattern rather than the tensor cores -- 82% of
 * achievable bandwidth at 7.8% of fp64 tensor peak -- and the cost of the
 * staging is paid per contraction while the memory time hiding it is paid per
 * byte. At lx = 8 opgrad carries 17749 bytes per contraction against axhelm's
 * 6144, so there is 2.9x as much memory time here to hide the same staging
 * behind. It is also in the default pnpn path, three calls per timestep in the
 * dealiased advection plus the pressure gradient, where dudxyz is not.
 *
 * w3 is the quadrature weight and is shared by every element, so it is indexed
 * by the point within the element rather than the global offset -- dmma_idx::l
 * rather than ::g. That is the only thing here the derivative kernel does not
 * already do.
 */

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 800) && (__CUDA_ARCH__ < 1000)

template< const int LX, const int NW >
__device__ __forceinline__
void opgrad_dmma_elem(double * __restrict__ ux,
                      double * __restrict__ uy,
                      double * __restrict__ uz,
                      const double * __restrict__ u,
                      const double * __restrict__ dx,
                      const double * __restrict__ dy,
                      const double * __restrict__ dz,
                      const double * __restrict__ drdx,
                      const double * __restrict__ dsdx,
                      const double * __restrict__ dtdx,
                      const double * __restrict__ drdy,
                      const double * __restrict__ dsdy,
                      const double * __restrict__ dtdy,
                      const double * __restrict__ drdz,
                      const double * __restrict__ dsdz,
                      const double * __restrict__ dtdz,
                      const double * __restrict__ w3,
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

  /* One staging, three outputs: the nine metrics stream from global exactly
     as the seven geometric factors do in the scalar axhelm dmma kernel */
  for (int p = tid; p < NP; p += nthrds) {
    const dmma_idx x = pack::map(p, ebase, nelv);

    if (x.live) {
      const int c = x.c;
      const int gp = x.g;
      const double W3 = w3[x.l];
      const double rtmp = shr[c];
      const double stmp = shs[c];
      const double ttmp = sht[c];

      ux[gp] = W3 * (drdx[gp] * rtmp
                     + dsdx[gp] * stmp
                     + dtdx[gp] * ttmp);

      uy[gp] = W3 * (drdy[gp] * rtmp
                     + dsdy[gp] * stmp
                     + dtdy[gp] * ttmp);

      uz[gp] = W3 * (drdz[gp] * rtmp
                     + dsdz[gp] * stmp
                     + dtdz[gp] * ttmp);
    }
  }
}

#endif // __CUDA_ARCH__ in [800, 1000)

/*
 * Compile-time dispatch onto the DMMA element kernel. The launch macros in
 * opr_opgrad.cu are written for every LX the operator dispatches and for
 * whatever `real` is, so every combination has to compile; the ones the
 * strategy does not cover -- single precision, LX outside the supported
 * range, a build without an fp64 tensor core arch -- resolve to this no-op.
 * The autotuner never selects the strategy for them, see dmma_lx_supported()
 * and cuda_have_dmma() in dmma_kernel.h.
 */
template< typename T, const int LX, const int NW >
struct opgrad_dmma_dispatch {
  __device__ static void run(T * __restrict__, T * __restrict__,
                             T * __restrict__,
                             const T * __restrict__, const T * __restrict__,
                             const T * __restrict__, const T * __restrict__,
                             const T * __restrict__, const T * __restrict__,
                             const T * __restrict__, const T * __restrict__,
                             const T * __restrict__, const T * __restrict__,
                             const T * __restrict__, const T * __restrict__,
                             const T * __restrict__, const T * __restrict__,
                             const int) { }
};

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 800) && (__CUDA_ARCH__ < 1000)

/* Keep in sync with dmma_lx_supported() in dmma_kernel.h */
#define NEKO_OPGRAD_DMMA_DISPATCH(LXV)                                         \
  template< const int NW >                                                     \
  struct opgrad_dmma_dispatch< double, LXV, NW > {                             \
    __device__ static void run(double * __restrict__ ux,                       \
                               double * __restrict__ uy,                       \
                               double * __restrict__ uz,                       \
                               const double * __restrict__ u,                  \
                               const double * __restrict__ dx,                 \
                               const double * __restrict__ dy,                 \
                               const double * __restrict__ dz,                 \
                               const double * __restrict__ drdx,               \
                               const double * __restrict__ dsdx,               \
                               const double * __restrict__ dtdx,               \
                               const double * __restrict__ drdy,               \
                               const double * __restrict__ dsdy,               \
                               const double * __restrict__ dtdy,               \
                               const double * __restrict__ drdz,               \
                               const double * __restrict__ dsdz,               \
                               const double * __restrict__ dtdz,               \
                               const double * __restrict__ w3,                 \
                               const int nelv) {                               \
      opgrad_dmma_elem< LXV, NW >(ux, uy, uz, u, dx, dy, dz,                   \
                                  drdx, dsdx, dtdx, drdy, dsdy, dtdy,          \
                                  drdz, dsdz, dtdz, w3, nelv);                 \
    }                                                                          \
  }

NEKO_OPGRAD_DMMA_DISPATCH(2);
NEKO_OPGRAD_DMMA_DISPATCH(3);
NEKO_OPGRAD_DMMA_DISPATCH(4);
NEKO_OPGRAD_DMMA_DISPATCH(5);
NEKO_OPGRAD_DMMA_DISPATCH(6);
NEKO_OPGRAD_DMMA_DISPATCH(7);
NEKO_OPGRAD_DMMA_DISPATCH(8);

#endif // __CUDA_ARCH__ in [800, 1000)

template< typename T, const int LX, const int NW >
__global__ void NEKO_EB_BOUNDS(32 * NW)
opgrad_kernel_dmma(T * __restrict__ ux,
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

  opgrad_dmma_dispatch< T, LX, NW >::run(ux, uy, uz, u, dx, dy, dz,
                                         drdx, dsdx, dtdx, drdy, dsdy, dtdy,
                                         drdz, dsdz, dtdz, w3, nelv);
}

/**
 * Device kernel for the weak gradient on the fp64 tensor cores, with the
 * element staged by the TMA engine
 *
 * Same three contractions and the same cube layout as opgrad_dmma_elem(); at
 * lx == DMMA_P the padding is empty and the cube offset is the flat point
 * index, so the arithmetic is identical and only the staging differs.
 *
 * All ten input copies -- u and the nine metric cubes, 77% of the traffic --
 * are issued together at kernel entry, u on its own mbarrier and waited
 * immediately, the metrics on a second and waited only at the pointwise step
 * that consumes them. Batching every read into one issue is what separated
 * the scalar axhelm TMA kernel that won from the component-at-a-time vector
 * one that did not: a lone 4 kB bulk copy never amortises the engine's fixed
 * latency.
 *
 * Thirteen cubes is 54800 B, past the 48 kB a block gets for free, so the
 * allocation is dynamic and its layout is opgrad_tma_smem in
 * dmma_tma_kernel.h. That is byte for byte the batched axhelm block, so it
 * runs at the same four blocks per SM and reuses the same device gate.
 *
 * The three outputs leave as ordinary coalesced stores, see the note on the
 * struct: only u's cube is free by then, and the stores are 19% of the
 * traffic against the loads' 77%.
 */

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 900) &&                        \
    (__CUDA_ARCH__ < 1000) && NEKO_TMA_TOOLKIT

template< const int LX, const int NW >
__device__ __forceinline__
void opgrad_dmma_tma_elem(double * __restrict__ ux,
                          double * __restrict__ uy,
                          double * __restrict__ uz,
                          const double * __restrict__ u,
                          const double * __restrict__ dx,
                          const double * __restrict__ dy,
                          const double * __restrict__ dz,
                          const double * __restrict__ drdx,
                          const double * __restrict__ dsdx,
                          const double * __restrict__ dtdx,
                          const double * __restrict__ drdy,
                          const double * __restrict__ dsdy,
                          const double * __restrict__ dtdy,
                          const double * __restrict__ drdz,
                          const double * __restrict__ dsdz,
                          const double * __restrict__ dtdz,
                          const double * __restrict__ w3) {

  /* Only lx == DMMA_P stages as a contiguous run of bytes, which is what a
     bulk copy moves; see the scope note in dmma_tma_kernel.h. Keep in step
     with dmma_tma_opgrad_lx_supported() */
  static_assert(LX == DMMA_P,
                "the dmma tma variant stages whole cubes only");

  extern __shared__ __align__(128) char neko_opgrad_dyn_smem[];
  opgrad_tma_smem &sm =
    *reinterpret_cast< opgrad_tma_smem * >(neko_opgrad_dyn_smem);

  enum { CUBE_BYTES = DMMA_CUBE * (int) sizeof(double) };

  const int nthrds = 32 * NW;
  const int tid = threadIdx.x;
  const int wf = tid >> 5;
  const int ebase = blockIdx.x * DMMA_CUBE;

  if (tid == 0) {
    tma_barrier_init(&sm.bar_u, 1);
    tma_barrier_init(&sm.bar_g, 1);
  }

  /* At LX == DMMA_P the derivative matrix fills the staged one exactly, and
     there is no padding anywhere to zero */
  for (int p = tid; p < DMMA_MAT; p += nthrds) {
    sm.dx[p] = dx[p];
    sm.dy[p] = dy[p];
    sm.dz[p] = dz[p];
  }

  /* Both barriers initialised and every derivative matrix staged before
     anyone waits on the one or contracts with the other */
  __syncthreads();

  if (tid == 0) {
    tma_expect(&sm.bar_u, CUBE_BYTES);
    tma_load(sm.u, u + ebase, CUBE_BYTES, &sm.bar_u);

    tma_expect(&sm.bar_g, DMMA_NG_OPGRAD * CUBE_BYTES);
    tma_load(sm.g[0], drdx + ebase, CUBE_BYTES, &sm.bar_g);
    tma_load(sm.g[1], dsdx + ebase, CUBE_BYTES, &sm.bar_g);
    tma_load(sm.g[2], dtdx + ebase, CUBE_BYTES, &sm.bar_g);
    tma_load(sm.g[3], drdy + ebase, CUBE_BYTES, &sm.bar_g);
    tma_load(sm.g[4], dsdy + ebase, CUBE_BYTES, &sm.bar_g);
    tma_load(sm.g[5], dtdy + ebase, CUBE_BYTES, &sm.bar_g);
    tma_load(sm.g[6], drdz + ebase, CUBE_BYTES, &sm.bar_g);
    tma_load(sm.g[7], dsdz + ebase, CUBE_BYTES, &sm.bar_g);
    tma_load(sm.g[8], dtdz + ebase, CUBE_BYTES, &sm.bar_g);
  }

  /* u only. The nine metric cubes are still arriving */
  tma_wait(&sm.bar_u, 0);

  dmma_contract< 0, false, false, NW >(sm.r, sm.dx, sm.u, wf);
  dmma_contract< 1, false, false, NW >(sm.s, sm.dy, sm.u, wf);
  dmma_contract< 2, false, false, NW >(sm.t, sm.dz, sm.u, wf);

  __syncthreads();

  tma_wait(&sm.bar_g, 0);

  /* At LX == DMMA_P the cube offset, the point within the element and the
     loop counter are all the same index, so w3 is read at p */
  for (int p = tid; p < DMMA_CUBE; p += nthrds) {
    const int gp = ebase + p;
    const double W3 = w3[p];
    const double rtmp = sm.r[p];
    const double stmp = sm.s[p];
    const double ttmp = sm.t[p];

    ux[gp] = W3 * (sm.g[0][p] * rtmp
                   + sm.g[1][p] * stmp
                   + sm.g[2][p] * ttmp);

    uy[gp] = W3 * (sm.g[3][p] * rtmp
                   + sm.g[4][p] * stmp
                   + sm.g[5][p] * ttmp);

    uz[gp] = W3 * (sm.g[6][p] * rtmp
                   + sm.g[7][p] * stmp
                   + sm.g[8][p] * ttmp);
  }
}

#endif // __CUDA_ARCH__ == sm_90 with a CUDA 12 toolkit

/*
 * Compile-time dispatch onto the TMA staged DMMA element kernel, see the note
 * on opgrad_dmma_dispatch above. The no-op covers everything the strategy does
 * not: single precision, any lx but DMMA_P, a build without sm_90, and a
 * toolkit older than CUDA 12.
 */
template< typename T, const int LX, const int NW >
struct opgrad_dmma_tma_dispatch {
  __device__ static void run(T * __restrict__, T * __restrict__,
                             T * __restrict__,
                             const T * __restrict__, const T * __restrict__,
                             const T * __restrict__, const T * __restrict__,
                             const T * __restrict__, const T * __restrict__,
                             const T * __restrict__, const T * __restrict__,
                             const T * __restrict__, const T * __restrict__,
                             const T * __restrict__, const T * __restrict__,
                             const T * __restrict__, const T * __restrict__) { }
};

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 900) &&                        \
    (__CUDA_ARCH__ < 1000) && NEKO_TMA_TOOLKIT

/* Keep in sync with dmma_tma_opgrad_lx_supported() in dmma_tma_kernel.h */
#define NEKO_OPGRAD_DMMA_TMA_DISPATCH(LXV)                                     \
  template< const int NW >                                                     \
  struct opgrad_dmma_tma_dispatch< double, LXV, NW > {                         \
    __device__ static void run(double * __restrict__ ux,                       \
                               double * __restrict__ uy,                       \
                               double * __restrict__ uz,                       \
                               const double * __restrict__ u,                  \
                               const double * __restrict__ dx,                 \
                               const double * __restrict__ dy,                 \
                               const double * __restrict__ dz,                 \
                               const double * __restrict__ drdx,               \
                               const double * __restrict__ dsdx,               \
                               const double * __restrict__ dtdx,               \
                               const double * __restrict__ drdy,               \
                               const double * __restrict__ dsdy,               \
                               const double * __restrict__ dtdy,               \
                               const double * __restrict__ drdz,               \
                               const double * __restrict__ dsdz,               \
                               const double * __restrict__ dtdz,               \
                               const double * __restrict__ w3) {               \
      opgrad_dmma_tma_elem< LXV, NW >(ux, uy, uz, u, dx, dy, dz,               \
                                      drdx, dsdx, dtdx, drdy, dsdy, dtdy,      \
                                      drdz, dsdz, dtdz, w3);                   \
    }                                                                          \
  }

NEKO_OPGRAD_DMMA_TMA_DISPATCH(8);

#endif // __CUDA_ARCH__ == sm_90 with a CUDA 12 toolkit

template< typename T, const int LX, const int NW >
__global__ void NEKO_EB_BOUNDS(32 * NW)
opgrad_kernel_dmma_tma(T * __restrict__ ux,
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

  opgrad_dmma_tma_dispatch< T, LX, NW >::run(ux, uy, uz, u, dx, dy, dz,
                                             drdx, dsdx, dtdx,
                                             drdy, dsdy, dtdy,
                                             drdz, dsdz, dtdz, w3);
}

/*
 * Opt into the TMA variant's dynamic allocation, once per specialisation.
 *
 * A block gets 48 kB of shared memory without asking; anything past that has
 * to be requested per kernel, and the carveout moved with it. Both are
 * properties of the function rather than of the launch, so this is a function
 * local static -- one flag per <T, LX, NW> -- and the launch macro calls it
 * ahead of every launch, where after the first it is a predictable branch.
 * Mirrors ax_helm_dmma_tma_batch_optin().
 *
 * Returns false if the device refuses, which the tuner has already ruled out
 * via cuda_have_tma_opgrad(); the error is cleared rather than left to surface
 * against an unrelated CUDA_CHECK later.
 */
template< typename T, const int LX, const int NW >
static inline bool opgrad_dmma_tma_optin()
{
  static int state = -1;

  if (state < 0) {
    const void * const fn =
      (const void *) opgrad_kernel_dmma_tma< T, LX, NW >;
    cudaError_t err =
      cudaFuncSetAttribute(fn, cudaFuncAttributeMaxDynamicSharedMemorySize,
                           NEKO_OPGRAD_TMA_SMEM);

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

#endif // __MATH_OPGRAD_KERNEL_H__
