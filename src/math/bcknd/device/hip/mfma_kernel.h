#ifndef __MATH_MFMA_KERNEL_H__
#define __MATH_MFMA_KERNEL_H__
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
 * Shared matrix-core (MFMA) primitives for the spectral-element tensor
 * contractions, used by the Ax-helm autotuner strategy.
 *
 * Recovered onto develop 2026-08-22 from the feature/mfma branch, scoped down
 * to Ax-helm: the other operators it covered (dudxyz, cdtp, conv1,
 * convect_scalar) are the four worst placed in the arithmetic intensity
 * ranking and are not worth carrying without a measurement. The contraction
 * primitives below are unchanged from the hardware-validated versions -- see
 * the layout note on mfma_contract_4x4 for why they must not be re-derived
 * from host simulation alone.
 *
 * Maps a reference derivative-matrix * field contraction onto the AMD matrix
 * cores on CDNA2 (MI250X / gfx90a) and CDNA3 (MI300A / MI300X / gfx942).  Double
 * precision defaults to the batched 4x4x4 tile (v_mfma_f64_4x4x4f64, full
 * M-utilisation for M = LX < 16; see mfma_contract_4x4), single precision uses
 * the 16x16x4 tile (v_mfma_f32_16x16x4f32; no f32 4x4x4 equivalent).  Double
 * precision can be forced back onto 16x16x4 with -DMFMA_F64_USE_16X16.  One
 * thread block is launched as a single 64-lane wavefront (blockDim = (64,1,1),
 * or (64, NWF, 1) for the multi-wavefront Ax-helm) and processes one element;
 * each contraction is a D * U GEMM with M = LX, N = LX^2, K = LX, the partial
 * tiles masked off.
 *
 * Supported for single and double precision and 4 <= LX <= 12; the upper
 * bound is set by the LDS needed to keep the cubes resident (operators stage
 * up to 4*LX^3 + 3*LX^2 elements, ~57 KB of f64 at LX = 12).
 *
 * 16x16x4 register/lane layout (used by the f32 tile and the f64 16x16x4
 * fallback; the 4x4x4 layout is documented at mfma_contract_4x4).
 *
 * Verification status, since this tile came from the same source that got the
 * 4x4x4 block selector and contraction index the wrong way round: CONFIRMED
 * on gfx90a 2026-08-22, both accumulator packings. The f32 sweep of mfma_probe
 * covers the A and B layouts (shared between precisions) and the i = 4*g + r
 * packing at ~1e-7, fp32 epsilon; a second run with -DMFMA_F64_USE_16X16
 * covers the i = g + 4*r packing at ~1e-16. A layout error reads as O(1) here,
 * so neither pass rests on a loose tolerance.
 *
 * Note the f64 results are bit-identical between this tile and the 4x4x4 one.
 * That is expected rather than suspicious: both decompose K into chunks of
 * four with one MFMA per chunk accumulated in sequence, so the summation order
 * is the same. It is not the signature of dead code, which is what identical
 * results across NWF would be -- NWF changes which wavefront takes which
 * column group and so must perturb the schedule.  D = A*B + C,
 * wave of 64 lanes, g = lane/16 in 0..3, c = lane%16 in 0..15.  A and B share
 * the same layout for both precisions; only the accumulator packing of the
 * 16x16 result differs:
 *   A[i][k] : lane holds A[i = c][k = g]
 *   B[k][j] : lane holds B[k = g][j = c]
 *   D[i][j] : lane holds D[j = c] in accumulator slot r = 0..3, with row
 *             i = g + 4*r for f64 (rows spread with stride 4) and
 *             i = 4*g + r for f32 (four contiguous rows per lane group)
 */

#include <stdlib.h>
#include <string.h>
#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>

/*
 * Reports whether the device code really was compiled for a matrix core
 * architecture, i.e. whether the __gfx90a__ / __gfx942__ guard below was true
 * in the device pass.
 *
 * This is not the same question as "does the device have matrix cores". The
 * contraction primitives and their call sites are all guarded on those
 * macros, so a build whose offload arch does not include the running device's
 * -- or a code object selected from a fat binary built for something else --
 * turns the whole strategy into a silent no-op: the kernel launches, writes
 * nothing, and leaves stale values in the output. That fails as bad results
 * rather than as an error, which is the worst way for it to fail. Checking it
 * from the device removes the guesswork.
 */
/* static, not merely file scope by convention: this header is included by
   every operator that offers an MFMA strategy -- ax_helm, dudxyz, opgrad,
   conv1 and cdtp -- and a non-template __global__ with external linkage is
   then defined once per translation unit, which the linker rejects as a
   multiple definition. Internal linkage gives each unit its own copy, which
   is what the rest of this header already relies on. */
static __global__ void hip_mfma_arch_probe(int * flag) {
#if defined(__gfx90a__) || defined(__gfx942__)
  *flag = 1;
#else
  *flag = 0;
#endif
}

/**
 * Returns true if the current device exposes the matrix cores used by the
 * MFMA strategies (currently gfx90a / MI250X and gfx942 / MI300A / MI300X)
 * *and* this build actually compiled the device code for them.  The f64
 * 4x4x4 / 16x16x4 and f32 16x16x4 instructions are all available on these
 * arches, so a single check gates either precision.  Result is cached after
 * the first query.
 */
static inline bool hip_have_mfma() {
  static int cached = -1;
  if (cached < 0) {
    int dev = 0;
    hipDeviceProp_t prop;
    cached = 0;
    if (hipGetDevice(&dev) == hipSuccess &&
        hipGetDeviceProperties(&prop, dev) == hipSuccess &&
        (strstr(prop.gcnArchName, "gfx90a") != NULL ||
         strstr(prop.gcnArchName, "gfx942") != NULL)) {
      int *d_flag = NULL;
      int flag = 0;
      if (hipMalloc(&d_flag, sizeof(int)) == hipSuccess) {
        if (hipMemcpy(d_flag, &flag, sizeof(int),
                      hipMemcpyHostToDevice) == hipSuccess) {
          hipLaunchKernelGGL(hip_mfma_arch_probe, dim3(1), dim3(1), 0, 0,
                             d_flag);
          if (hipGetLastError() == hipSuccess &&
              hipMemcpy(&flag, d_flag, sizeof(int),
                        hipMemcpyDeviceToHost) == hipSuccess) {
            cached = flag;
          }
        }
        /* Unlike the queries above, a failure here is not "the strategy is
           unavailable" -- the pointer came from a hipMalloc that succeeded, so
           a bad free means the context is broken. Checked rather than folded
           into cached, and checked rather than discarded: hipFree is
           nodiscard */
        HIP_CHECK(hipFree(d_flag));
      }
    }
  }
  return cached == 1;
}

/**
 * Compile-time predicate for the LX values that the MFMA strategy supports.
 * Single or double precision and 4 <= LX <= 12 (bounded above by the LDS
 * needed to keep one element resident).  MUST match the dispatch
 * specialisations in each operator's kernel header.
 */
/*
 * Both precisions are offered, and both are verified against a reference on
 * gfx90a (mfma_probe, 144 configurations, 2026-08-22). They do not share a
 * code path: f64 goes through the batched 4x4x4 tile, f32 has no 4x4x4
 * instruction and uses the 16x16x4 one. If a future part disagrees, excluding
 * a precision here is a one-line change.
 */
template < const int LX >
static inline bool mfma_lx_supported() {
  return (sizeof(real) == 8 || sizeof(real) == 4) && (LX >= 4) && (LX <= 12);
}

/*
 * Wavefronts per block for the MFMA kernels. Candidate C selects 2^C
 * wavefronts, i.e. 1, 2, 4 or 8.
 *
 * Two things are being traded. The contraction stripes its N column groups
 * across wavefronts, and there are only NGROUPS = ceil(LX^2/16) of them --
 * one at LX = 4, four at LX = 8 -- so past that count the extra wavefronts
 * idle through the matrix core work. The staging and pointwise loops, on the
 * other hand, keep scaling: at LX = 8, eight wavefronts is 512 threads for
 * 512 points, one each, which is why the branch this came from defaulted to
 * eight. The measured LX = 8 curve was still improving at four, hence the
 * fourth candidate.
 */
#define NEKO_MFMA_CANDIDATES 4
#define NEKO_MFMA_NWF(C) (1 << (C))
#define NEKO_MFMA_NTHRDS(C) dim3(64, NEKO_MFMA_NWF(C), 1)

/*
 * Column groups per contraction -- the wavefront-parallel work one element
 * offers. Both tiles group N the same way, 16 columns at a time, so this is
 * ceil(LX^2/16) either way: 1 at LX = 4, 4 at LX = 8, 9 at LX = 12.
 *
 * A wavefront beyond that count has no matrix core work left on that element,
 * which is what the LX = 4 single precision sweep measured: 20.5 / 23.6 /
 * 34.0 / 60.2 us as NWF went 1, 2, 4, 8, monotonically worse, against
 * 218 / 157 / 136.4 / 136.4 at LX = 8 where four groups exist. Rather than
 * cap the sweep, the surplus wavefronts are given their own element: NWF is
 * read as wavefronts per block, the block covers EB elements and WPE =
 * NWF/EB wavefronts cooperate on each. At LX = 4 with eight wavefronts that
 * is eight elements, one each, with nothing idle; at LX = 12 it is one
 * element and eight cooperating wavefronts, exactly as before. This matters
 * because p-multigrid smooths at LX = 4 and 2, so low order Ax is hot rather
 * than incidental.
 *
 * EB is the driving quantity and WPE follows from it, not the other way
 * round: the block is partitioned into EB equal groups, so EB has to divide
 * NWF exactly or the leftover wavefronts address an element the block does
 * not own -- past the end of the shared staging arrays, and past the end of
 * global storage in the last block. EB = NWF/NGROUPS is therefore rounded
 * down to a power of two, which divides NWF for every candidate since NWF is
 * itself 2^C. WPE may then exceed NGROUPS -- at LX = 6, NGROUPS = 3 and four
 * wavefronts give EB = 1, WPE = 4 -- which is harmless: mfma_contract_4x4()
 * strides the groups with `ng = wf + gp * NWF` under `ng < NGROUPS`, so a
 * wavefront without a group of its own simply issues no matrix core work,
 * while still taking its share of the staging and pointwise passes. The
 * alternative, capping WPE at NGROUPS and letting EB absorb the remainder,
 * would grow the shared footprint (LX = 10 with eight wavefronts would want
 * 66 kB) for no gain.
 */
#define NEKO_MFMA_NGROUPS(LX) (((LX) * (LX) + 15) / 16)
/* Elements per block: surplus wavefronts, rounded down to a power of two so
   that WPE * EB == NWF exactly */
#define NEKO_MFMA_EB_N(NWF, LX)                                               \
  ((NWF) / NEKO_MFMA_NGROUPS(LX) >= 8 ? 8 :                                   \
   (NWF) / NEKO_MFMA_NGROUPS(LX) >= 4 ? 4 :                                   \
   (NWF) / NEKO_MFMA_NGROUPS(LX) >= 2 ? 2 : 1)
#define NEKO_MFMA_EB(LX, C) NEKO_MFMA_EB_N(NEKO_MFMA_NWF(C), LX)
/* Wavefronts cooperating on one element */
#define NEKO_MFMA_WPE(LX, C) (NEKO_MFMA_NWF(C) / NEKO_MFMA_EB(LX, C))
#define NEKO_MFMA_NBLCKS(NELV, LX, C)                                         \
  dim3(((NELV) + NEKO_MFMA_EB(LX, C) - 1) / NEKO_MFMA_EB(LX, C), 1, 1)

/*
 * Whether the autotuner sweeps the MFMA strategy, on by default wherever the
 * hardware and the polynomial order allow it.
 *
 * It was briefly off while the matrix core contraction was known broken -- the
 * lane layout had the block selector and the contraction index interchanged,
 * which a "the solver converges" check had failed to catch for a long time.
 * The layout was measured on gfx90a, corrected, and mfma_contract_4x4 now
 * reproduces a CPU reference to ~1e-16 over every supported order, axis,
 * transpose/accumulate mode and wavefront count, so there is no reason to
 * withhold it from the sweep. Kept as an off switch in the shape of
 * NEKO_EB_TUNE, for A/B work.
 */
static int neko_mfma_sweep()
{
  const char *v = getenv("NEKO_MFMA_TUNE");

  if (v != NULL) {
    return (atoi(v) != 0);
  }
  return 1;
}

/* Forced candidate, used when NEKO_AUTOTUNE pins the MFMA variant */
static int neko_mfma_env()
{
  const char *v = getenv("NEKO_MFMA_NWF");
  int c = (v != NULL) ? atoi(v) : 0;

  if (c < 0 || c >= NEKO_MFMA_CANDIDATES) {
    c = 0;
  }
  return c;
}

/* Report every measured MFMA candidate, see NEKO_TUNE_LOG in
   elem_block_tune.h */
#define NEKO_TUNE_LOG_MFMA(LX, T3)                                            \
  do {                                                                        \
    for (int c = 0; c < NEKO_MFMA_CANDIDATES; c++) {                          \
      if ((T3)[c] >= NEKO_TUNE_INIT) { continue; }                            \
      sprintf(neko_log_buf, "MFMA  %dwf %-2de: %9.2f us/call",                \
              NEKO_MFMA_NWF(c), NEKO_MFMA_EB(LX, c),                          \
              NEKO_TUNE_US((T3)[c], iters));                                  \
      log_message(neko_log_buf);                                              \
    }                                                                         \
  } while (0)

#if defined(__gfx90a__) || defined(__gfx942__)

/* 4-wide accumulators for the gfx90a / gfx942 matrix cores. */
typedef double mfma_f64x4 __attribute__((ext_vector_type(4)));
typedef float  mfma_f32x4 __attribute__((ext_vector_type(4)));

/*
 * Per-precision matrix-core traits: the 4-wide accumulator type, the MFMA
 * builtin, and the accumulator-slot -> output-row mapping (see the layout
 * note above; f64 spreads rows with stride 4, f32 packs four contiguous rows).
 */
template< typename T >
struct mfma_traits;

template< >
struct mfma_traits< double > {
  typedef mfma_f64x4 acc_t;
  __device__ __forceinline__ static acc_t mma(double a, double b, acc_t c) {
    return __builtin_amdgcn_mfma_f64_16x16x4f64(a, b, c, 0, 0, 0);
  }
  __device__ __forceinline__ static int out_row(const int g, const int r) {
    return g + 4 * r;
  }
};

template< >
struct mfma_traits< float > {
  typedef mfma_f32x4 acc_t;
  __device__ __forceinline__ static acc_t mma(float a, float b, acc_t c) {
    return __builtin_amdgcn_mfma_f32_16x16x4f32(a, b, c, 0, 0, 0);
  }
  __device__ __forceinline__ static int out_row(const int g, const int r) {
    return 4 * g + r;
  }
};

/*
 * Linearised index into an LX^3 cube stored as i + LX*j + LX*LX*k, where the
 * coordinate 'p' lies on contraction axis AXIS and 'n' enumerates the two
 * remaining axes as n = a + LX*b.
 */
template< const int LX, const int AXIS >
__device__ __forceinline__ int mfma_cube_idx(const int p, const int n) {
  const int a = n % LX;
  const int b = n / LX;
  if (AXIS == 0) return p + LX * a + LX * LX * b;   // contract i; n = (j,k)
  if (AXIS == 1) return a + LX * p + LX * LX * b;   // contract j; n = (i,k)
  return a + LX * b + LX * LX * p;                  // contract k; n = (i,j)
}

/*
 * One wavefront contracts the reference derivative matrix 'dmat' (LX x LX,
 * stored column-major: D(row,col) = dmat[row + LX*col]) with the cube 'in'
 * along axis AXIS, writing the cube 'out':
 *
 *   out[idx(m,n)] (+)= sum_l  D(m,l) * in[idx(l,n)]   (TRANSPOSE = false)
 *   out[idx(m,n)] (+)= sum_l  D(l,m) * in[idx(l,n)]   (TRANSPOSE = true)
 *
 * GEMM dimensions M = LX, N = LX*LX, K = LX, tiled over 16x16x4 MFMA tiles
 * with the partial M/N/K tiles masked.  ACCUM selects += over = .
 */
/*
 * NWF cooperating wavefronts (wf = 0..NWF-1) stripe the NTILES N-tiles among
 * themselves -- wavefront wf handles nt = wf, wf+NWF, ...  NWF = 1, wf = 0
 * (the defaults) reproduces the single-wavefront contraction.
 */
template< typename T, const int LX, const int AXIS,
          const bool TRANSPOSE, const bool ACCUM, const int NWF = 1 >
__device__ __forceinline__
void mfma_contract(T * __restrict__ out,
                   const T * __restrict__ dmat,
                   const T * __restrict__ in,
                   const int lane, const int wf = 0) {
  typedef mfma_traits<T> mma_t;
  const int g = lane >> 4;             // lane / 16 -> 0..3
  const int c = lane & 15;             // lane % 16 -> 0..15
  const int NTILES = (LX * LX + 15) / 16;
  const int KSTEPS = (LX + 3) / 4;
  const int NPASS  = (NTILES + NWF - 1) / NWF;  // N-tiles handled by this wave

#pragma unroll
  for (int p = 0; p < NPASS; ++p) {
    const int nt = wf + p * NWF;                // this wavefront's N-tile
    if (nt < NTILES) {
      const int n = nt * 16 + c;                // free (column) index
      typename mma_t::acc_t acc = {0, 0, 0, 0};
#pragma unroll
      for (int ks = 0; ks < KSTEPS; ++ks) {
        const int l = ks * 4 + g;               // contraction index
        T a = 0;
        if (c < LX && l < LX)
          a = TRANSPOSE ? dmat[l + c * LX]      // D(l,c) = D^T(c,l)
                        : dmat[c + l * LX];     // D(c,l)
        T b = 0;
        if (l < LX && n < LX * LX)
          b = in[mfma_cube_idx<LX, AXIS>(l, n)];
        acc = mma_t::mma(a, b, acc);
      }
      const T dvals[4] = { acc[0], acc[1], acc[2], acc[3] };
#pragma unroll
      for (int r = 0; r < 4; ++r) {
        const int m = mma_t::out_row(g, r);     // output coordinate on AXIS
        if (m < LX && n < LX * LX) {
          const int idx = mfma_cube_idx<LX, AXIS>(m, n);
          if (ACCUM) out[idx] += dvals[r];
          else       out[idx]  = dvals[r];
        }
      }
    }
  }
}

/*
 * Double-precision batched matrix-core contraction using v_mfma_f64_4x4x4f64.
 *
 * Same contract as mfma_contract() -- out(m,n) (+)= sum_l D(m,l) in(l,n), with
 * D(l,m) when TRANSPOSE -- but tiled with 4x4x4 MFMA tiles and the
 * instruction's four blocks assigned to four consecutive 4-column N-subtiles.
 * GEMM M = LX, N = LX*LX, K = LX, tiled as MT = ceil(LX/4) M-tiles,
 * KSTEPS = ceil(LX/4) K-steps and NGROUPS = ceil(LX^2/16) column groups (each
 * group = 4 blocks x 4 columns); partial M/N/K masked with zeros.
 *
 * Versus the 16x16x4 tile, the 4-wide M granularity fills M = LX < 16 exactly
 * (LX = 8 runs the matrix core at 100% M-utilisation instead of 50%), at the
 * cost of 4x as many, 1/4-sized MFMA issues that feed the same f64 matrix
 * pipeline.  Double precision only: the f32 4x4 instruction has K = 1, so there
 * is no single-precision counterpart -- single precision keeps the 16x16x4 tile
 * (see mfma_contract_sel below).
 *
 * Wave of 64 lanes = 4 blocks x 16 lanes.  v_mfma_f64_4x4x4f64 lane layout,
 * with lo = lane%4, gemm = (lane/4)%4 and kq = lane/16:
 *   A[i][k] : lane holds A[i = lo ][k = kq]
 *   B[k][j] : lane holds B[k = kq ][j = lo]
 *   D[i][j] : lane holds D[i = kq ][j = lo]   (scalar f64 accumulator)
 * and the four independent 4x4x4 blocks are selected by 'gemm', not by
 * lane/16.
 *
 * This was MEASURED on gfx90a, not assumed: a one-hot read-out of which lane
 * receives which element (128 launches, no assumptions) produced exactly this
 * mapping.  The previous version had 'gemm' and 'kq' interchanged -- it used
 * lane/16 as the block selector and (lane/4)%4 as the contraction index -- and
 * every one of 72 probe configurations disagreed with a CPU reference.
 *
 * It had been believed validated because the Ax-helm fluid solver converged
 * with it.  It does: a wrong but symmetric operator makes CG converge happily
 * to the solution of a different system.  Convergence is not correctness, and
 * only a diff against a reference settles a layout.  Any future correction
 * stays localised to the four index expressions (m_a, l, n, m_d).
 *
 * NWF cooperating wavefronts (wf = 0..NWF-1) stripe the NGROUPS column groups
 * among themselves -- wavefront wf handles ng = wf, wf+NWF, ...  NWF = 1,
 * wf = 0 (the defaults) reproduces the single-wavefront contraction.
 */
template< const int LX, const int AXIS,
          const bool TRANSPOSE, const bool ACCUM, const int NWF = 1 >
__device__ __forceinline__
void mfma_contract_4x4(double * __restrict__ out,
                       const double * __restrict__ dmat,
                       const double * __restrict__ in,
                       const int lane, const int wf = 0) {
  const int lo   = lane & 3;           // 0..3 : A row, B column, D column
  const int gemm = (lane >> 2) & 3;    // 0..3 : which of the four 4x4x4 blocks
  const int kq   = lane >> 4;          // 0..3 : contraction index within a step
  const int MT      = (LX + 3) / 4;
  const int KSTEPS  = (LX + 3) / 4;
  const int NGROUPS = (LX * LX + 15) / 16;
  const int NPASS   = (NGROUPS + NWF - 1) / NWF; // groups handled by this wave

#pragma unroll
  for (int mt = 0; mt < MT; ++mt) {
    const int m_a = mt * 4 + lo;       // A row (input lane layout: i = lo)
    const int m_d = mt * 4 + kq;       // D row (output lane layout: i = kq)
#pragma unroll
    for (int gp = 0; gp < NPASS; ++gp) {
      const int ng = wf + gp * NWF;    // this wavefront's column group
      if (ng < NGROUPS) {
        /* the four blocks take four consecutive 4-column N-subtiles */
        const int n = ng * 16 + gemm * 4 + lo;  // column (N) for B in / D out
        double acc = 0.0;
#pragma unroll
        for (int ks = 0; ks < KSTEPS; ++ks) {
          const int l = ks * 4 + kq;   // contraction index (k = kq for A and B)
          double a = 0.0;
          if (m_a < LX && l < LX)
            a = TRANSPOSE ? dmat[l + m_a * LX]  // D(l,m) = D^T(m,l)
                          : dmat[m_a + l * LX]; // D(m,l)
          double b = 0.0;
          if (l < LX && n < LX * LX)
            b = in[mfma_cube_idx<LX, AXIS>(l, n)];
          acc = __builtin_amdgcn_mfma_f64_4x4x4f64(a, b, acc, 0, 0, 0);
        }
        if (m_d < LX && n < LX * LX) {
          const int idx = mfma_cube_idx<LX, AXIS>(m_d, n);
          if (ACCUM) out[idx] += acc;
          else       out[idx]  = acc;
        }
      }
    }
  }
}

/*
 * Precision-dispatched tensor contraction: double precision uses the batched
 * 4x4x4 matrix core (full M-utilisation for M = LX < 16), single precision the
 * 16x16x4 tile (no f32 4x4x4 equivalent).  Lets one call site cover both.
 */
template< typename T, const int LX, const int AXIS,
          const bool TRANSPOSE, const bool ACCUM, const int NWF = 1 >
struct mfma_contract_sel {
  __device__ __forceinline__
  static void run(T * __restrict__ out, const T * __restrict__ dmat,
                  const T * __restrict__ in, const int lane,
                  const int wf = 0) {
    mfma_contract<T, LX, AXIS, TRANSPOSE, ACCUM, NWF>(out, dmat, in, lane, wf);
  }
};

template< const int LX, const int AXIS, const bool TRANSPOSE, const bool ACCUM,
          const int NWF >
struct mfma_contract_sel<double, LX, AXIS, TRANSPOSE, ACCUM, NWF> {
  __device__ __forceinline__
  static void run(double * __restrict__ out, const double * __restrict__ dmat,
                  const double * __restrict__ in, const int lane,
                  const int wf = 0) {
#ifdef MFMA_F64_USE_16X16
    /*
     * Opt-out fallback: the original 16x16x4 tile (M-utilisation 50% at LX=8,
     * 75% at LX=12).  Kept as an escape hatch behind -DMFMA_F64_USE_16X16; the
     * batched 4x4x4 path below is the hardware-validated default.
     */
    mfma_contract<double, LX, AXIS, TRANSPOSE, ACCUM, NWF>(out, dmat, in,
                                                           lane, wf);
#else
    /*
     * Default f64 path: batched 4x4x4 matrix-core tile, full M-utilisation.
     * Hardware-validated -- the ax_helm fluid solver converges on gfx90a/
     * gfx942.  Still multi-wavefront via NWF.
     */
    mfma_contract_4x4<LX, AXIS, TRANSPOSE, ACCUM, NWF>(out, dmat, in, lane, wf);
#endif
  }
};

#endif // __gfx90a__ || __gfx942__
#endif // __MATH_MFMA_KERNEL_H__
