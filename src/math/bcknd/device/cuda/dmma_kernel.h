#ifndef __MATH_DMMA_KERNEL_H__
#define __MATH_DMMA_KERNEL_H__
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
 * Double precision matrix-core (DMMA) primitives for the spectral element
 * tensor contractions, the CUDA counterpart of the AMD MFMA strategies.
 *
 * Maps a reference derivative matrix * field contraction onto the fp64 tensor
 * cores, which are exposed as a single fixed 8x8x4 WMMA tile lowering to
 * mma.sync.aligned.m8n8k4.f64. Each contraction is a D * U GEMM with M = LX,
 * N = LX^2, K = LX, and rather than mask the partial M, N and K tiles that a
 * general LX would leave, every cube is staged into shared memory padded to
 * DMMA_P^3 and every derivative matrix to DMMA_P^2, with the padding zero
 * filled. The GEMM is then always 8 x 64 x 8 and every tile is full:
 *
 *   - zero rows of D contribute nothing to a padded output row,
 *   - zero rows of the staged cube contribute nothing to any output,
 *   - padded output entries are never read back out of shared memory,
 *
 * so no lane ever has to test an index. At LX = 8 the padding is empty and
 * the staging is a straight copy, which is the case the strategy is aimed at.
 *
 * That fixed tile is also what bounds the strategy: LX > DMMA_P would need
 * the cubes padded to 16^3, i.e. 4 * 16^3 * 8 = 128 kB of shared memory for
 * the four cubes, well past the 48 kB a block gets without the opt-in dynamic
 * path. Supported for double precision and 2 <= LX <= 8 only -- the lower
 * bound is 2 rather than 4 because of the packing below, see
 * dmma_lx_supported(); note that fp32 has no tensor core equivalent, only
 * TF32 with an 11 bit significand, which is not usable for the operator
 * inside a Krylov solve, so a single precision build has no DMMA strategy.
 *
 * The fp64 tensor cores themselves only exist on the data centre parts:
 * sm_80 (A100) and sm_90 (H100 / GH200), where they run at twice the fp64
 * FMA rate. Consumer sm_86 / sm_89 assemble the instruction but have no fp64
 * tensor hardware at all, and Blackwell moves the fp64 rates again, so both
 * the compile time guard and cuda_have_dmma() below allow-list [sm_80, sm_90]
 * rather than testing for "at least sm_80". Widen both together if a later
 * part turns out to be worth measuring.
 *
 * The result is bit-for-bit an fp64 computation -- mma.m8n8k4.f64 is IEEE
 * double throughout, unlike the reduced precision tensor tiles -- but the
 * summation order differs from the scalar variants, so the Ax output differs
 * in the last bits, exactly as it already does between the 1d and kstep
 * variants.
 *
 * Cube and matrix layout. A staged cube is
 *
 *   cube[i + DMMA_SI * j + DMMA_SJ * k],  i, j, k < DMMA_P
 *
 * and a staged derivative matrix is dmat[m + DMMA_P * l] = D(m, l), matching
 * the dx[i + l*LX] = D(i,l) convention of the scalar kernels. Contracting
 * axis AXIS turns the cube into an M x N matrix whose (m, n) element sits at a
 * constant stride, which is all load_matrix_sync() needs:
 *
 *   AXIS = 0: m = i, n = (j,k), element at m + DMMA_SI * n -> column major,
 *             ldm = DMMA_SI, one tile per group of 8 (j,k) pairs
 *   AXIS = 1: m = j, n = i, element at n + DMMA_SI * m + DMMA_SJ * k -> row
 *             major, ldm = DMMA_SI, one tile per k slab
 *   AXIS = 2: m = k, n = i, element at n + DMMA_SI * j + DMMA_SJ * m -> row
 *             major, ldm = DMMA_SJ, one tile per j slab
 *
 * Each view has DMMA_P tiles and the same view describes the input and the
 * output of a contraction, so one traits struct covers both.
 */

#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime.h>
#include <mma.h>
#include <device/device_config.h>

enum {
  DMMA_P = 8,                    /* Tile M and N, and the padded cube extent */
  DMMA_KS = 4,                   /* Tile K */
  DMMA_SI = DMMA_P,              /* Padded cube stride between j */
  DMMA_SJ = DMMA_P * DMMA_P,     /* Padded cube stride between k */
  DMMA_CUBE = DMMA_P * DMMA_P * DMMA_P,
  DMMA_MAT = DMMA_P * DMMA_P,
  DMMA_NTILES = DMMA_P,          /* Tiles (or slabs) per contraction */
  DMMA_KSTEPS = DMMA_P / DMMA_KS
};

/*
 * Warps per block for the DMMA kernels. Candidate C selects 2^(C+1) warps,
 * which stripe the DMMA_NTILES tiles of every contraction among themselves;
 * at the top candidate each warp owns exactly one tile. Picked between at
 * runtime by the operator's autotuner, like the elements per block candidates
 * of the kstep variants.
 */
#define NEKO_DMMA_CANDIDATES 3
#define NEKO_DMMA_NW(C) (1 << ((C) + 1))
#define NEKO_DMMA_NTHRDS(C) dim3(32 * NEKO_DMMA_NW(C), 1, 1)

/*
 * Elements packed into one padded cube.
 *
 * The staging pads every cube to DMMA_P^3 so that the fixed 8x8x4 tile is
 * always full, which at LX = DMMA_P is exact and below it is waste: at LX = 4
 * only 64 of 512 points are real, so seven eighths of every contraction is
 * spent on zeros. Where LX divides DMMA_P that waste can be turned back into
 * work by packing DMMA_P/LX elements along each axis -- one per sub-cube --
 * and making the staged derivative matrix block diagonal, one LX x LX copy of
 * D per sub-cube. Because D is then block diagonal, the contraction along any
 * axis couples only indices within the same sub-cube, so the single padded
 * contraction computes every packed element independently and correctly.
 *
 * LX = 8 packs one element and is exactly what it was; LX = 4 packs eight,
 * LX = 2 packs sixty-four. LX = 3, 5, 6 and 7 do not divide DMMA_P and keep
 * one element with the old padding waste. This matters because p-multigrid
 * smooths at LX = 4 and 2, so low order Ax is hot rather than incidental --
 * and note the packing costs no extra shared memory at all, the cube is
 * DMMA_P^3 either way.
 */
#define NEKO_DMMA_PPA(LX) ((DMMA_P % (LX) == 0) ? (DMMA_P / (LX)) : 1)
#define NEKO_DMMA_PACK(LX)                                                    \
  (NEKO_DMMA_PPA(LX) * NEKO_DMMA_PPA(LX) * NEKO_DMMA_PPA(LX))
#define NEKO_DMMA_NBLCKS(NELV, LX)                                            \
  dim3(((NELV) + NEKO_DMMA_PACK(LX) - 1) / NEKO_DMMA_PACK(LX), 1, 1)

/* Forced candidate, used when NEKO_AUTOTUNE pins the DMMA variant */
static int neko_dmma_env()
{
  const char *v = getenv("NEKO_DMMA_NW");
  int c = (v != NULL) ? atoi(v) : 0;

  if (c < 0 || c >= NEKO_DMMA_CANDIDATES) {
    c = 0;
  }
  return c;
}

/* Report every measured DMMA candidate, see NEKO_TUNE_LOG in
   elem_block_tune.h. The label is padded to 13 as a whole rather than by
   hand counted field widths, so it stays column aligned with the 1D, KSTEP
   and Chose lines whatever the warp and element counts render as */
#define NEKO_TUNE_LOG_DMMA(LX, T3)                                            \
  do {                                                                        \
    for (int c = 0; c < NEKO_DMMA_CANDIDATES; c++) {                          \
      if ((T3)[c] >= NEKO_TUNE_INIT) { continue; }                            \
      char lbl_[16];                                                          \
      sprintf(lbl_, "DMMA  %dw %de",                                          \
              NEKO_DMMA_NW(c), NEKO_DMMA_PACK(LX));                           \
      sprintf(neko_log_buf, "%-13s: %9.2f us/call", lbl_,                     \
              NEKO_TUNE_US((T3)[c], iters));                                  \
      log_message(neko_log_buf);                                              \
    }                                                                         \
  } while (0)

/*
 * Was any architecture with fp64 tensor cores among the ones this translation
 * unit was compiled for? Unlike HIP, a binary built for an older arch still
 * runs on a newer device by JIT compiling its PTX, and the kernel body below
 * is guarded on __CUDA_ARCH__ -- so a build for, say, sm_70 running on an
 * H100 would JIT the *no-op* body and silently return a zero Ax. Checking the
 * compiled arch list keeps the strategy from ever being offered in that case.
 *
 * __CUDA_ARCH_LIST__ needs CUDA >= 11.5; on an older toolkit the strategy is
 * off unless forced with -DNEKO_DMMA_ARCH_COMPILED=1.
 */
#ifndef NEKO_DMMA_ARCH_COMPILED
#if defined(__CUDA_ARCH_LIST__)
static inline bool dmma_arch_compiled()
{
  const int arch[] = { __CUDA_ARCH_LIST__ };

  for (int i = 0; i < (int) (sizeof(arch)/sizeof(arch[0])); i++) {
    if (arch[i] >= 800 && arch[i] < 1000) {
      return true;
    }
  }
  return false;
}
#else
static inline bool dmma_arch_compiled()
{
  return false;
}
#endif
#else
static inline bool dmma_arch_compiled()
{
  return (NEKO_DMMA_ARCH_COMPILED != 0);
}
#endif

/**
 * Returns true if the current device exposes the fp64 tensor cores used by the
 * DMMA strategies and this build can actually reach them. Result is cached
 * after the first query.
 */
static inline bool cuda_have_dmma()
{
  static int cached = -1;

  if (cached < 0) {
    int dev = 0;
    cudaDeviceProp prop;

    if (!dmma_arch_compiled()) {
      cached = 0;
    } else if (cudaGetDevice(&dev) == cudaSuccess &&
               cudaGetDeviceProperties(&prop, dev) == cudaSuccess) {
      cached = ((prop.major == 8 && prop.minor == 0) ||
                (prop.major == 9)) ? 1 : 0;
    } else {
      cached = 0;
    }
  }
  return cached == 1;
}

/**
 * Compile-time predicate for the LX values that the DMMA strategy supports.
 * Double precision and 2 <= LX <= DMMA_P, see the note on the padded staging
 * above. The lower bound is 2 rather than 4 because packing makes LX = 2
 * useful rather than absurd -- sixty-four elements to a cube, and it is a
 * p-multigrid level. LX = 3 is supported but packs one element and wastes
 * 95% of every contraction; the tuner will reject it.
 *
 * MUST match the dispatch specialisations in each operator's kernel
 * header.
 *
 * The double precision requirement is deliberate and is NOT an oversight to
 * be brought in line with the HIP side, which does offer both precisions.
 * That difference is hardware: AMD's v_mfma_f32_16x16x4f32 is true IEEE fp32,
 * so an sp build there loses nothing, while NVIDIA has no fp32 tensor path at
 * all -- only TF32, whose 11 bit significand (unit roundoff 4.9e-4, against
 * fp32's 6e-8) is not usable for the operator inside a Krylov solve. The
 * accuracy preserving alternative, 3xTF32 splitting, costs three MMAs to buy
 * back precision in a kernel that is bandwidth bound with flops to spare.
 *
 * The empirical argument is stronger than the roofline one. AMD ran the
 * experiment NVIDIA cannot: with a true fp32 matrix core and therefore no
 * accuracy compromise whatsoever, MFMA still measured 4.5% BEHIND a plain 1d
 * kernel at lx = 8, and the whole fp32 speedup over fp64 (2.05x) came from
 * moving half the bytes at identical achieved bandwidth. A path that loses
 * where it is free will not win where it is lossy.
 */
template< const int LX >
static inline bool dmma_lx_supported()
{
  return (sizeof(real) == 8) && (LX >= 2) && (LX <= DMMA_P);
}

/**
 * The same predicate for the vector operator, whose supported range is not the
 * same and must not be assumed to be.
 *
 * The vector kernel stages one element per block and does not pack, so the
 * lower bound stays at 4: below it the padding waste that packing removes for
 * the scalar kernel is still there, and there is nothing to gain. That is a
 * performance argument, but the predicate is a correctness one -- the tuner
 * launches whatever it is told is supported, and an LX the vector dispatch
 * does not specialise resolves to the no-op primary template, which times as
 * free, wins the comparison and leaves stale values in au/av/aw. Lowering
 * dmma_lx_supported() to 2 for the packed scalar kernel is exactly how that
 * happened.
 *
 * MUST match NEKO_AX_HELM_DMMA_VECTOR_DISPATCH in ax_helm_kernel.h.
 */
template< const int LX >
static inline bool dmma_vector_lx_supported()
{
  return (sizeof(real) == 8) && (LX >= 4) && (LX <= DMMA_P);
}

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 800) && (__CUDA_ARCH__ < 1000)

/*
 * Index map from the flat staging loop counter p to the padded cube offset and
 * the global element offset, specialised on whether the cube holds more than
 * one element.
 *
 * The two cases have to *generate* different code, not merely fold to it. With
 * PACK > 1 the grid is ceil(nelv/PACK) blocks, so the last block can own slots
 * past nelv: those are clamped on the way in -- the cube has to be finite, see
 * the note above -- and dropped on the way out. With PACK == 1 the grid is
 * exactly nelv blocks and neither is needed, but leaving that to the optimiser
 * does NOT work: nvcc cannot prove blockIdx.x < nelv, so the select survives on
 * every staging and metric load where the unpacked kernel had a single hoisted
 * blockIdx.x * LX3, and the guarded store keeps its branch. Hence the explicit
 * PPA == 1 specialisation below, which restores that addressing verbatim and
 * makes `live` a compile time true.
 */
struct dmma_idx {
  int c;                     /* offset into the padded cube */
  int g;                     /* offset into the global element arrays */
  bool live;                 /* false only for a padded tail slot, PACK > 1 */
};

template< const int LX, const int PPA >
struct dmma_pack {
  enum { PACK = PPA * PPA * PPA,
         LX3 = LX * LX * LX,
         NP = PACK * LX3 };

  /* Loop invariant part of the addressing, hoisted by the caller */
  __device__ __forceinline__ static int ebase() {
    return blockIdx.x * PACK;
  }

  __device__ __forceinline__ static dmma_idx map(const int p, const int ebase,
                                                 const int nelv) {
    const int q = p / LX3;
    const int r = p - q * LX3;
    const int i = r % LX;
    const int jk = r / LX;
    const int j = jk % LX;
    const int k = jk / LX;
    const int qa = q % PPA;
    const int qb = (q / PPA) % PPA;
    const int qc = q / (PPA * PPA);
    const int eq = ebase + q;
    dmma_idx x;

    x.c = (qa * LX + i) + DMMA_SI * (qb * LX + j) + DMMA_SJ * (qc * LX + k);
    x.g = r + (eq < nelv ? eq : nelv - 1) * LX3;
    x.live = (eq < nelv);
    return x;
  }
};

/*
 * One element per cube: the grid covers nelv exactly, so there is no tail to
 * clamp, no store to predicate, and ebase() is the element's base offset
 * outright rather than an element index. This is the addressing the kernel had
 * before packing existed, and at LX == DMMA_P the cube offset reduces to p.
 */
template< const int LX >
struct dmma_pack< LX, 1 > {
  enum { PACK = 1,
         LX3 = LX * LX * LX,
         NP = LX3 };

  __device__ __forceinline__ static int ebase() {
    return blockIdx.x * LX3;
  }

  __device__ __forceinline__ static dmma_idx map(const int p, const int ebase,
                                                 const int) {
    const int i = p % LX;
    const int jk = p / LX;
    const int j = jk % LX;
    const int k = jk / LX;
    dmma_idx x;

    x.c = i + DMMA_SI * j + DMMA_SJ * k;
    x.g = ebase + p;
    x.live = true;
    return x;
  }
};

/*
 * How the staged cube is read as an M x N matrix when contracting axis AXIS,
 * see the layout note above. COL_MAJOR is carried alongside the layout tag
 * because the pointer offset of a K step depends on it: the contraction index
 * is the row index of B, so it advances by one element in a column major view
 * and by ldm in a row major one.
 */
template< const int AXIS >
struct dmma_view;

template< >
struct dmma_view< 0 > {
  typedef nvcuda::wmma::col_major layout;
  enum { LDM = DMMA_SI, COL_MAJOR = 1 };
  __device__ __forceinline__ static int base(const int t) {
    return DMMA_SI * DMMA_P * t;
  }
};

template< >
struct dmma_view< 1 > {
  typedef nvcuda::wmma::row_major layout;
  enum { LDM = DMMA_SI, COL_MAJOR = 0 };
  __device__ __forceinline__ static int base(const int t) {
    return DMMA_SJ * t;
  }
};

template< >
struct dmma_view< 2 > {
  typedef nvcuda::wmma::row_major layout;
  enum { LDM = DMMA_SJ, COL_MAJOR = 0 };
  __device__ __forceinline__ static int base(const int t) {
    return DMMA_SI * t;
  }
};

/*
 * The A operand is the staged derivative matrix, D(m,l) = dmat[m + DMMA_P*l].
 * Untransposed that is a column major 8x8 with ldm = DMMA_P; transposed it is
 * the same storage read row major. Here the contraction index is the column
 * index of A, so the K step offset is the mirror of the B one above.
 */
template< const bool TRANSPOSE >
struct dmma_amat;

template< >
struct dmma_amat< false > {
  typedef nvcuda::wmma::col_major layout;
  __device__ __forceinline__ static int koff(const int l0) {
    return DMMA_P * l0;
  }
};

template< >
struct dmma_amat< true > {
  typedef nvcuda::wmma::row_major layout;
  __device__ __forceinline__ static int koff(const int l0) {
    return l0;
  }
};

/*
 * NW cooperating warps (wf = threadIdx.x / 32) contract the staged derivative
 * matrix 'dmat' with the staged cube 'in' along axis AXIS into the staged cube
 * 'out':
 *
 *   out[idx(m,n)] (+)= sum_l  D(m,l) * in[idx(l,n)]   (TRANSPOSE = false)
 *   out[idx(m,n)] (+)= sum_l  D(l,m) * in[idx(l,n)]   (TRANSPOSE = true)
 *
 * ACCUM selects += over =. Warp wf takes tiles wf, wf+NW, ..., so every lane
 * of a warp follows the same tile sequence and the warp wide WMMA operations
 * are never reached divergently.
 */
template< const int AXIS, const bool TRANSPOSE, const bool ACCUM,
          const int NW >
__device__ __forceinline__
void dmma_contract(double * __restrict__ out,
                   const double * __restrict__ dmat,
                   const double * __restrict__ in,
                   const int wf)
{
  namespace wmma = nvcuda::wmma;
  typedef dmma_view< AXIS > view;
  typedef dmma_amat< TRANSPOSE > amat;

  const wmma::layout_t mem_layout =
    view::COL_MAJOR ? wmma::mem_col_major : wmma::mem_row_major;

  /* A is the same for every tile of the contraction, so it is loaded once
     rather than once per tile */
  wmma::fragment< wmma::matrix_a, DMMA_P, DMMA_P, DMMA_KS, double,
                  typename amat::layout > a[DMMA_KSTEPS];
#pragma unroll
  for (int ks = 0; ks < DMMA_KSTEPS; ks++) {
    wmma::load_matrix_sync(a[ks], dmat + amat::koff(ks * DMMA_KS), DMMA_P);
  }

  const int npass = (DMMA_NTILES + NW - 1) / NW;

#pragma unroll
  for (int p = 0; p < npass; p++) {
    const int t = wf + p * NW;

    if (t < DMMA_NTILES) {
      wmma::fragment< wmma::accumulator, DMMA_P, DMMA_P, DMMA_KS, double > acc;

      if (ACCUM) {
        wmma::load_matrix_sync(acc, out + view::base(t), view::LDM, mem_layout);
      } else {
        wmma::fill_fragment(acc, 0.0);
      }

#pragma unroll
      for (int ks = 0; ks < DMMA_KSTEPS; ks++) {
        const int l0 = ks * DMMA_KS;
        const int koff = view::COL_MAJOR ? l0 : (l0 * (int) view::LDM);
        wmma::fragment< wmma::matrix_b, DMMA_P, DMMA_P, DMMA_KS, double,
                        typename view::layout > b;

        wmma::load_matrix_sync(b, in + view::base(t) + koff, view::LDM);
        wmma::mma_sync(acc, a[ks], b, acc);
      }

      wmma::store_matrix_sync(out + view::base(t), acc, view::LDM, mem_layout);
    }
  }
}

#endif // __CUDA_ARCH__ in [800, 1000)

#endif // __MATH_DMMA_KERNEL_H__
