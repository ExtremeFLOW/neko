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
 * four cubes plus three LX^2 derivative matrices, 63360 B of f64 at LX = 12
 * with the padded Ax-helm layout, 58752 B without it).
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
/* NEKO_EB_MAX_LDS, which the elements-per-block ladder below is clamped by */
#include "elem_block.h"

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
/*
 * The candidate space is two dimensional: the wavefronts per block above, and
 * which matrix core tile the contraction is issued on. Candidate C encodes
 * both -- NWF = 2^(C mod 4) and TILE = C / 4 -- so one selector still names a
 * whole geometry and candidates 0..3 mean exactly what they always did.
 *
 * TILE 0 is the precision's default tile: the batched 4x4x4 in double
 * precision, the 16x16x4 in single, where there is no 4x4x4 counterpart
 * (its f32 sibling has K = 1). TILE 1 is the 16x16x4 tile in both, so in
 * single precision the two are the same code and only TILE 0 is offered --
 * see mfma_tile_offered().
 *
 * WHY THE TILE IS MEASURED RATHER THAN CHOSEN. It used to be a compile-time
 * #ifdef, on the argument that the 4x4x4 tile fills M = LX < 16 exactly where
 * the 16x16x4 one wastes half its rows at LX = 8. That argument tacitly
 * assumes the two instructions run at the same rate, and they do not. From
 * AMD's own matrix instruction calculator, identically on gfx90a and gfx942:
 *
 *   v_mfma_f64_16x16x4f64  2048 flop / 32 cyc = 256 flop/CU/cycle, 0 wait
 *   v_mfma_f64_4x4x4f64     512 flop / 16 cyc = 128 flop/CU/cycle, 4 waits
 *
 * The small tile is rate-halved -- 128 flop/CU/cycle is the f64 *vector* rate,
 * so it surrenders the whole reason to use a matrix core in double precision
 * -- and a chain of them into one accumulator, which is exactly what
 * mfma_contract_4x4's K loop is, owes four cycles per step where the large
 * tile owes none. Counting issues per element per contraction (4x4x4:
 * ceil(LX/4) M-tiles x groups x ceil(LX/4) K-steps at 16 cycles; 16x16x4:
 * groups x K-steps at 32) the M-utilisation gain is exactly cancelled at
 * LX = 5..8 and reversed beyond it: 4x4x4 wins 2x at LX = 4, loses ~11% at
 * LX = 5..8 once the accumulate waits are counted, and loses 1.5x at
 * LX = 9..12. It also re-reads the B operand once per M-tile, so it costs 2x
 * the operand traffic at LX = 5..8 and 3x at LX = 9..12 on a kernel measured
 * at 69-72% of the memory roof.
 *
 * So the analytical answer is "4x4x4 at LX = 4, 16x16x4 above it", which is
 * not what the old default did -- and it is an analytical answer, of the same
 * kind as the one it replaces. Both tiles are hardware verified against a CPU
 * reference and produce bit-identical f64 results, so trying both costs
 * nothing but tuning time. Hence: measure it.
 */
/*
 * THE 16-WAVEFRONT RUNG. The ladder used to stop at eight wavefronts, 512
 * threads, and that ceiling is what the strategy was losing on at LX = 8
 * rather than anything to do with matrix cores. The 1D kernel that beats it
 * there stages MORE LDS -- 6*LX^2 + 4*LX^3 against this kernel's
 * 3*LX^2 + 4*LX^3 -- but runs CHUNKS threads on one element, so its 1024- and
 * 512-thread candidates reach 8.00 and 6.00 wavefronts per SIMD where this
 * ladder tops out at 3.00. A block shape the tuner cannot express is a block
 * shape it cannot choose, so the rung is added rather than argued about.
 *
 * 1024 threads is the CDNA workgroup maximum, so this is the last rung there
 * is. It is only reachable where the LDS budget allows it, which is what the
 * clamp in NEKO_MFMA_EB_N below is for: at LX = 8 in double precision the
 * unclamped EB = NWF/NGROUPS = 4 would ask for 4 elements x 4 cubes and blow
 * the 64 kB workgroup limit -- as a compile error in the kernel's
 * static_assert, not as a bad launch, but a compile error all the same.
 */
#define NEKO_MFMA_NWF_CANDIDATES 5
#define NEKO_MFMA_TILE_CANDIDATES 2
#define NEKO_MFMA_CANDIDATES                                                  \
  (NEKO_MFMA_NWF_CANDIDATES * NEKO_MFMA_TILE_CANDIDATES)
#define NEKO_MFMA_NWF(C) (1 << ((C) % NEKO_MFMA_NWF_CANDIDATES))
#define NEKO_MFMA_TILE(C) ((C) / NEKO_MFMA_NWF_CANDIDATES)
#define NEKO_MFMA_NTHRDS(C) dim3(64, NEKO_MFMA_NWF(C), 1)

/*
 * LDS PADDING OF THE STAGED CUBE.
 *
 * The cube is stored i + SJ*j + SK*k with SJ = LX + PAD and SK = SJ*LX. PAD
 * is not a free parameter: CDNA's LDS is 32 banks of one dword, a b32 access
 * is serviced 32 lanes at a time and a b64 access 16 lanes at a time (16
 * lanes x 2 dwords = 32 dwords), and distinct addresses landing on one bank
 * inside such a group serialise. A contraction along AXIS = 0 walks its
 * free index n = (j,k) across the lanes, so consecutive lanes are SJ elements
 * apart, and the whole question is what SJ does modulo the bank count:
 *
 *   SJ odd  -> f64 lanes cover 16 distinct bank pairs, f32 lanes 32 distinct
 *              banks. Conflict free.
 *   SJ even -> the stride and the bank count share a factor and the group
 *              collapses onto a few banks. At LX = 8 in double precision that
 *              is an 8-WAY conflict on both the B read and the D store.
 *
 * So padding an even LX by one is what removes the conflict, and padding an
 * odd LX would introduce one -- the opposite of the usual "always pad by one"
 * reflex. AXIS = 2 is conflict free either way (its index is n + SK*p, unit
 * stride across lanes), AXIS = 1 sits in between.
 *
 * MODELLED, NOT MEASURED: a bank simulator over every lane of every LDS
 * access to the CUBE -- the six contractions' operand reads and result
 * stores, and the three linear passes -- for every LX in 4..12, both
 * precisions, both tiles. The derivative matrix is padded by the same policy
 * and is counted separately, see NEKO_MFMA_SD_N. Cycles per element, best
 * wavefront shape, B-hoist on, padded-slot linear passes:
 *
 *   LX = 8  f64  4x4x4 : 1696 -> 836   (-51%)
 *   LX = 8  f32  16x16 :  816 -> 642   (-21%)
 *   LX = 12 f64  4x4x4 : 3996 -> 2808  (-30%)
 *   LX = 4  f64  4x4x4 :  180 -> 136   (-24%)
 *   LX = 4  f32  16x16 :  102 -> 108   (+6%, the one regression)
 *   LX = 6, 10         : within 1% either way
 *
 * Hence the policy below: pad where LX is a multiple of four, which is where
 * the conflict is worst and the win is large, and leave LX = 4 in single
 * precision alone. That last exception is not cosmetic -- LX = 4 f32 is the
 * only configuration where this strategy has ever been measured to win
 * anything (1.9% over 1D on gfx90a), so it does not get handed a 6%
 * regression on a model. LX odd is never padded; LX = 6 and 10 gain nothing
 * because 6 and 10 are already only 2-way, and padding them costs a bigger
 * cube.
 *
 * -DNEKO_MFMA_PAD=0 forces the old unpadded layout and -DNEKO_MFMA_PAD=1 pads
 * every even LX, for an A/B against the model. Odd LX is never padded whatever
 * the setting, because there it is known to make things worse.
 *
 * Padding is a property of the layout, not of the lane mapping: which lane
 * handles which (m, n) is untouched, so it is exactly the kind of change a
 * host reference catches -- unlike a lane-layout change, where an internally
 * consistent error cancels in write-back (see mfma_contract_4x4).
 */
#ifdef NEKO_MFMA_PAD
#define NEKO_MFMA_PAD_N(LX, SZ) (((LX) % 2 == 0) ? (NEKO_MFMA_PAD) : 0)
#else
#define NEKO_MFMA_PAD_N(LX, SZ)                                               \
  ((((LX) % 4 == 0) && !((LX) == 4 && (SZ) == 4)) ? 1 : 0)
#endif
/* j-stride, k-stride and total slots of one staged cube */
#define NEKO_MFMA_SJ_N(LX, SZ) ((LX) + NEKO_MFMA_PAD_N(LX, SZ))
#define NEKO_MFMA_SK_N(LX, SZ) (NEKO_MFMA_SJ_N(LX, SZ) * (LX))
#define NEKO_MFMA_CUBE_N(LX, SZ) (NEKO_MFMA_SK_N(LX, SZ) * (LX))
/*
 * Column stride and size of a staged reference derivative matrix, which takes
 * the same padding for the same reason.
 *
 * D(row, col) = dmat[row + SD*col]. The three divergence contractions read it
 * transposed, D(l, m), so there the LANE index m walks the column stride --
 * 16 lanes SD elements apart in the 16x16x4 tile, four in the 4x4x4 one -- and
 * an even stride puts them on the same banks exactly as it does in the cube.
 * It is the same policy rather than a second one because it is win or neutral
 * at every order the policy pads, modelled the same way:
 *
 *   LX = 8  f64 : 4x4x4 576 -> 384 cycles, 16x16 480 -> 192
 *   LX = 12 f64 : 16x16 1296 -> 648
 *   LX = 4, and both precisions at every unpadded order : unchanged
 *
 * A rule of "pad every even order" would also take LX = 10 f64 from 714 to
 * 504, but costs LX = 10 f32 a 42% increase, which is the same reason the cube
 * policy stops at multiples of four.
 */
#define NEKO_MFMA_SD_N(LX, SZ) NEKO_MFMA_SJ_N(LX, SZ)
#define NEKO_MFMA_DMAT_N(LX, SZ) (NEKO_MFMA_SD_N(LX, SZ) * (LX))

/*
 * Whether candidate C's tile is a distinct thing to measure in this build.
 *
 * TILE 1 is the 16x16x4 tile, which is what TILE 0 already resolves to in
 * single precision and under -DMFMA_F64_USE_16X16, so offering it there would
 * time the same kernel twice. The launch macros instantiate both regardless --
 * every (precision, LX, candidate) has to compile -- and this only decides
 * what the sweep and the env pin will run.
 */
static inline bool mfma_tile_offered(const int c) {
  if (NEKO_MFMA_TILE(c) == 0) {
    return true;
  }
#ifdef MFMA_F64_USE_16X16
  return false;
#else
  return sizeof(real) == 8;
#endif
}

/* Which tile candidate C names, for the tuner log */
static inline const char *mfma_tile_name(const int c) {
#ifdef MFMA_F64_USE_16X16
  (void) c;
  return "16x16";
#else
  if (NEKO_MFMA_TILE(c) != 0) {
    return "16x16";
  }
  return (sizeof(real) == 8) ? "4x4x4" : "16x16";
#endif
}

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

/*
 * LDS a block of EB elements occupies: three reference derivative matrices
 * shared by the whole block, plus four staged cubes per element. Kept here
 * rather than in the kernels so that the EB ladder below and the kernels'
 * own static_assert cannot disagree about what fits.
 */
#define NEKO_MFMA_LDS_N(EB, LX, SZ)                                           \
  ((3 * NEKO_MFMA_DMAT_N(LX, SZ) + 4 * (EB) * NEKO_MFMA_CUBE_N(LX, SZ)) * (SZ))
#define NEKO_MFMA_LDS_FITS(EB, LX, SZ)                                        \
  (NEKO_MFMA_LDS_N(EB, LX, SZ) <= NEKO_EB_MAX_LDS)

/* Elements per block: surplus wavefronts, rounded down to a power of two so
   that WPE * EB == NWF exactly, and clamped to what the 64 kB workgroup LDS
   limit allows. The LDS clamp is not defensive tidiness: at 16 wavefronts the
   unclamped ladder asks for EB = 4 at LX = 8 and EB = 2 at LX = 10 and 11 in
   double precision, all three of which overflow the limit and stop the build
   in the kernels' static_assert. */
#define NEKO_MFMA_EB_N(NWF, LX, SZ)                                           \
  ((NWF) / NEKO_MFMA_NGROUPS(LX) >= 16 && NEKO_MFMA_LDS_FITS(16, LX, SZ) ? 16 :\
   (NWF) / NEKO_MFMA_NGROUPS(LX) >= 8 && NEKO_MFMA_LDS_FITS(8, LX, SZ) ? 8 :  \
   (NWF) / NEKO_MFMA_NGROUPS(LX) >= 4 && NEKO_MFMA_LDS_FITS(4, LX, SZ) ? 4 :  \
   (NWF) / NEKO_MFMA_NGROUPS(LX) >= 2 && NEKO_MFMA_LDS_FITS(2, LX, SZ) ? 2 : 1)
#define NEKO_MFMA_EB(LX, C)                                                   \
  NEKO_MFMA_EB_N(NEKO_MFMA_NWF(C), LX, sizeof(real))
/* Wavefronts cooperating on one element */
#define NEKO_MFMA_WPE(LX, C) (NEKO_MFMA_NWF(C) / NEKO_MFMA_EB(LX, C))

/*
 * Slots per thread: the CUBE_N slots of one staged cube shared out over the
 * WPE * 64 threads that serve it. Slots rather than points, because the
 * staging, pointwise and write-back passes walk the LDS cube at unit stride
 * and skip the pad slots -- see mfma_slot_point(). Where the cube is not
 * padded a slot is a point and this is the old points-per-thread count.
 *
 * Only the vector operator needs this on the host side, but it belongs here
 * with the rest of the geometry so that the device enum and the host launcher
 * cannot drift apart -- which is how the WPE/EB split came to address past the
 * end of its staging arrays.
 */
#define NEKO_MFMA_SPT_N(WPE, LX, SZ)                                          \
  ((NEKO_MFMA_CUBE_N(LX, SZ) + (WPE) * 64 - 1) / ((WPE) * 64))
#define NEKO_MFMA_SPT(LX, C)                                                  \
  NEKO_MFMA_SPT_N(NEKO_MFMA_WPE(LX, C), LX, sizeof(real))

/*
 * Whether the seven geometric factors of each of a thread's points are held in
 * registers, or read from global memory where they are used. The cost is
 * 7 * SPT * (sizeof(T)/4) VGPRs, budgeted at a quarter of the 256 a lane
 * addresses -- which leaves the four cube base pointers, the two index bases,
 * the contraction accumulator and the seven pointwise temporaries room to
 * live alongside it.
 *
 * It buys two different things in the two Ax-helm kernels, and the budget is
 * the same for both. In the vector kernel (ax_helm_mfma_vector_elem) it is
 * reuse: the factors are read once for the three components instead of once
 * each, which is the reason that operator exists at all. In the scalar kernel
 * (ax_helm_mfma_elem) there is no reuse to be had -- it is latency: the loads
 * are issued before the gradient contractions instead of behind the barrier
 * that follows them, so they are in flight across the matrix core work rather
 * than exposed after it.
 *
 * It follows from LX and the wavefront count rather than being a candidate of
 * its own, so the two modes cannot be compared directly -- which is why the
 * tuner reports which one each candidate ran, see NEKO_TUNE_LOG_MFMA_VEC.
 * Without that the sweep reads as a matrix core result when the step between
 * two candidates is really a change in memory traffic.
 */
#ifndef NEKO_MFMA_VECTOR_GREG_VGPRS
#define NEKO_MFMA_VECTOR_GREG_VGPRS 64
#endif
#define NEKO_MFMA_VECTOR_GREG_N(SPT, SZ)                                      \
  ((7 * (SPT) * ((SZ) / 4)) <= NEKO_MFMA_VECTOR_GREG_VGPRS)
#define NEKO_MFMA_VECTOR_GREG(LX, C)                                          \
  NEKO_MFMA_VECTOR_GREG_N(NEKO_MFMA_SPT(LX, C), sizeof(real))

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

/*
 * Candidate pinned by NEKO_MFMA_NWF and NEKO_MFMA_TILE, or -1 to leave it to
 * the sweep, see neko_eb_pin() in elem_block_tune.h.
 *
 * The two dimensions stay separate variables -- NEKO_MFMA_NWF keeps its old
 * range and meaning, NEKO_MFMA_TILE adds the tile -- rather than one index
 * into the combined space, so that a pin written before the tile existed
 * still selects what it used to. Setting either one pins the candidate;
 * leaving both unset leaves the whole geometry to the sweep. A tile this
 * build does not offer falls back to the default one at the requested
 * wavefront count.
 */
static int neko_mfma_pin()
{
  const char *v = getenv("NEKO_MFMA_NWF");
  const char *t = getenv("NEKO_MFMA_TILE");
  int nwf, tile, c;

  if (v == NULL && t == NULL) {
    return -1;
  }

  nwf = (v != NULL) ? atoi(v) : 0;
  tile = (t != NULL) ? atoi(t) : 0;
  if (nwf < 0 || nwf >= NEKO_MFMA_NWF_CANDIDATES) {
    nwf = 0;
  }
  if (tile < 0 || tile >= NEKO_MFMA_TILE_CANDIDATES) {
    tile = 0;
  }
  c = tile * NEKO_MFMA_NWF_CANDIDATES + nwf;
  if (!mfma_tile_offered(c)) {
    c = nwf;
  }
  return c;
}

/*
 * The same pin restricted to the wavefront dimension, for the operators that
 * offer only that one.
 *
 * The tile dimension is carried by the Helmholtz operator alone, scalar and
 * vector: it is the operator with the headroom to be worth the doubled
 * instantiation count, and the one every measurement so far is about. The
 * gradient-type operators keep the four wavefront candidates and their
 * contractions stay on the default tile. This exists so that a
 * NEKO_MFMA_TILE=1 pin is ignored there *explicitly* rather than by falling
 * through a switch onto some other candidate.
 */
static int neko_mfma_nwf_pin()
{
  const char *v = getenv("NEKO_MFMA_NWF");
  int c;

  if (v == NULL) {
    return -1;
  }

  c = atoi(v);
  if (c < 0 || c >= NEKO_MFMA_NWF_CANDIDATES) {
    c = 0;
  }
  return c;
}

/*
 * Candidates the sweep runs: the wavefront counts on the default tile, and
 * the same again on the 16x16x4 tile wherever this build offers that as a
 * distinct thing to measure, see mfma_tile_offered(). The offered set is
 * contiguous from 0, so it is a count rather than a skip and NEKO_TUNE_FOR()
 * in elem_block_tune.h can take it directly.
 */
static inline int neko_mfma_candidates()
{
  return mfma_tile_offered(NEKO_MFMA_NWF_CANDIDATES) ?
    NEKO_MFMA_CANDIDATES : NEKO_MFMA_NWF_CANDIDATES;
}

/* Report every measured MFMA candidate, see NEKO_TUNE_LOG in
   elem_block_tune.h. The tile is named because it is the dimension whose
   default was wrong for eight of the nine supported orders */
#define NEKO_TUNE_LOG_MFMA(LX, T3)                                            \
  do {                                                                        \
    for (int c = 0; c < NEKO_MFMA_CANDIDATES; c++) {                          \
      if ((T3)[c] >= NEKO_TUNE_INIT) { continue; }                            \
      sprintf(neko_log_buf, "MFMA  %s %2dwf %-2de %-5s: %9.2f us/call",      \
              mfma_tile_name(c), NEKO_MFMA_NWF(c), NEKO_MFMA_EB(LX, c),       \
              NEKO_MFMA_PAD_N(LX, sizeof(real)) ? "pad" : "plain",            \
              NEKO_TUNE_US((T3)[c], iters));                                  \
      log_message(neko_log_buf);                                              \
    }                                                                         \
  } while (0)

/*
 * The same for the vector operator, plus the register mode each candidate
 * ran in -- 'reg' where the geometric factors stay in registers across the
 * three components, 'glob' where they are re-read from global memory for each
 * one, see NEKO_MFMA_VECTOR_GREG. It is reported because it is derived from
 * the wavefront count rather than swept, so two neighbouring candidates can
 * differ in memory traffic as well as in block shape, and a step between them
 * would otherwise be read as a matrix core effect.
 */
#define NEKO_TUNE_LOG_MFMA_VEC(LX, T3)                                        \
  do {                                                                        \
    for (int c = 0; c < NEKO_MFMA_CANDIDATES; c++) {                          \
      if ((T3)[c] >= NEKO_TUNE_INIT) { continue; }                            \
      sprintf(neko_log_buf, "MFMA  %s %2dwf %-2de %-5s %-4s: %9.2f us/call",  \
              mfma_tile_name(c), NEKO_MFMA_NWF(c), NEKO_MFMA_EB(LX, c),       \
              NEKO_MFMA_PAD_N(LX, sizeof(real)) ? "pad" : "plain",            \
              NEKO_MFMA_VECTOR_GREG(LX, c) ? "reg" : "glob",                  \
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
 * Storage layout of one staged cube: i + SJ*j + SK*k over SIZE slots, with the
 * j-stride padded to an odd number of elements where that removes the LDS bank
 * conflicts -- see the padding note above for the model and the policy.
 *
 * PAD is carried as a template parameter rather than read from the macro so
 * that the operators which stage a plain LX^3 cube (dudxyz, opgrad, conv1,
 * cdtp) keep the layout they were written for: the default is 0 and only
 * Ax-helm opts in.
 */
template< const int LX, const int PAD >
struct mfma_cube {
  enum { SJ = LX + PAD,
         SK = (LX + PAD) * LX,
         SIZE = (LX + PAD) * LX * LX,
         /* column stride of a staged derivative matrix, padded with the cube
            and for the same reason -- see the note on NEKO_MFMA_SD_N */
         SD = LX + PAD,
         DSIZE = (LX + PAD) * LX };
};

/*
 * Offset of D(row, col) in a staged derivative matrix. The caller's copy is
 * LX x LX with a column stride of LX; this is where it lives once staged.
 */
template< const int LX, const int PAD >
__device__ __forceinline__ int mfma_dmat_idx(const int row, const int col) {
  return row + mfma_cube<LX, PAD>::SD * col;
}

/*
 * Linearised index into a staged cube, where the coordinate 'p' lies on
 * contraction axis AXIS and 'n' enumerates the two remaining axes as
 * n = a + LX*b.
 */
template< const int LX, const int AXIS, const int PAD = 0 >
__device__ __forceinline__ int mfma_cube_idx(const int p, const int n) {
  typedef mfma_cube<LX, PAD> C;
  const int a = n % LX;
  const int b = n / LX;
  if (AXIS == 0) return p + C::SJ * a + C::SK * b;  // contract i; n = (j,k)
  if (AXIS == 1) return a + C::SJ * p + C::SK * b;  // contract j; n = (i,k)
  return a + C::SJ * b + C::SK * p;                 // contract k; n = (i,j)
}

/*
 * The point of the element that LDS slot 's' holds, as an offset into the
 * caller's LX^3 global storage, or -1 for a pad slot that holds nothing.
 *
 * The staging, pointwise and write-back passes walk SLOTS at unit stride and
 * translate here, rather than walking points and translating the other way.
 * Both directions are correct; this one is the one that keeps those passes
 * conflict free. Walking points would leave the LDS side striding across the
 * padding -- a 2-way conflict on every pass, which the bank model prices at
 * more than the contraction gains back at LX = 6 and 10 -- while walking slots
 * leaves the GLOBAL side with a hole every SJ lanes, which costs nothing: the
 * lanes either side of the hole still touch the same cache lines, the hole
 * lane simply takes no part.
 *
 * The division is by a compile-time constant and is meant to be hoisted: each
 * kernel decodes its own slots once into a small register array and reuses it
 * across all three passes.
 */
template< const int LX, const int PAD >
__device__ __forceinline__ int mfma_slot_point(const int s) {
  typedef mfma_cube<LX, PAD> C;
  if (PAD == 0) {
    return s;
  }
  const int k = s / C::SK;
  const int rem = s - k * C::SK;
  const int j = rem / C::SJ;
  const int i = rem - j * C::SJ;
  /* j < LX holds by construction (rem < SK = SJ*LX), i and k do not */
  return (i < LX && k < LX) ? (i + LX * j + LX * LX * k) : -1;
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
          const bool TRANSPOSE, const bool ACCUM, const int NWF = 1,
          const int PAD = 0 >
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
          a = TRANSPOSE ? dmat[mfma_dmat_idx<LX, PAD>(l, c)]  // D(l,c)
                        : dmat[mfma_dmat_idx<LX, PAD>(c, l)]; // D(c,l)
        T b = 0;
        if (l < LX && n < LX * LX)
          b = in[mfma_cube_idx<LX, AXIS, PAD>(l, n)];
        acc = mma_t::mma(a, b, acc);
      }
      const T dvals[4] = { acc[0], acc[1], acc[2], acc[3] };
#pragma unroll
      for (int r = 0; r < 4; ++r) {
        const int m = mma_t::out_row(g, r);     // output coordinate on AXIS
        if (m < LX && n < LX * LX) {
          const int idx = mfma_cube_idx<LX, AXIS, PAD>(m, n);
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
          const bool TRANSPOSE, const bool ACCUM, const int NWF = 1,
          const int PAD = 0 >
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

  /*
   * M-tile innermost, with one accumulator per M-tile live across the K loop.
   *
   * The obvious nest is M-tile outermost, one accumulator, and that is what
   * this was. It re-reads the whole B operand once per M-tile -- B depends on
   * (group, K-step) and not on the M-tile at all -- which costs 2x the LDS
   * operand traffic at LX = 5..8 and 3x at LX = 9..12, on a kernel that
   * measures at 69-72% of the memory roof. The bank model puts the B read at
   * 1408 of the 2400 LDS cycles an LX = 8 f64 element spends, and hoisting it
   * out takes that to 704.
   *
   * It also breaks up the accumulate chain. AMD documents a four cycle wait on
   * a V_MFMA_4x4x4_F64 -> V_MFMA_4x4x4_F64 SrcC dependency, and the single
   * accumulator made the K loop exactly that chain; consecutive issues now
   * write different accumulators wherever MT > 1, which is every order above
   * four, so the waits overlap instead of adding up.
   *
   * Each accumulator still sums over ks in the same order it did, so the
   * results are bit-identical to the previous nest, not merely equivalent.
   */
#pragma unroll
  for (int gp = 0; gp < NPASS; ++gp) {
    const int ng = wf + gp * NWF;      // this wavefront's column group
    if (ng < NGROUPS) {
      /* the four blocks take four consecutive 4-column N-subtiles */
      const int n = ng * 16 + gemm * 4 + lo;  // column (N) for B in / D out
      double acc[MT];
#pragma unroll
      for (int mt = 0; mt < MT; ++mt)
        acc[mt] = 0.0;
#pragma unroll
      for (int ks = 0; ks < KSTEPS; ++ks) {
        const int l = ks * 4 + kq;     // contraction index (k = kq for A and B)
        double b = 0.0;
        if (l < LX && n < LX * LX)
          b = in[mfma_cube_idx<LX, AXIS, PAD>(l, n)];
#pragma unroll
        for (int mt = 0; mt < MT; ++mt) {
          const int m_a = mt * 4 + lo; // A row (input lane layout: i = lo)
          double a = 0.0;
          if (m_a < LX && l < LX)
            a = TRANSPOSE ? dmat[mfma_dmat_idx<LX, PAD>(l, m_a)]  // D(l,m)
                          : dmat[mfma_dmat_idx<LX, PAD>(m_a, l)]; // D(m,l)
          acc[mt] = __builtin_amdgcn_mfma_f64_4x4x4f64(a, b, acc[mt], 0, 0, 0);
        }
      }
#pragma unroll
      for (int mt = 0; mt < MT; ++mt) {
        const int m_d = mt * 4 + kq;   // D row (output lane layout: i = kq)
        if (m_d < LX && n < LX * LX) {
          const int idx = mfma_cube_idx<LX, AXIS, PAD>(m_d, n);
          if (ACCUM) out[idx] += acc[mt];
          else       out[idx]  = acc[mt];
        }
      }
    }
  }
}

/*
 * Precision- and tile-dispatched tensor contraction, so that one call site
 * covers both precisions and both matrix core tiles.
 *
 * TILE 0 is the precision's default: the batched 4x4x4 tile in double
 * precision (full M-utilisation for M = LX < 16, at half the FLOP rate --
 * see the candidate space note above for why that trade is worth measuring),
 * the 16x16x4 tile in single, which has no 4x4x4 counterpart. TILE 1 is the
 * 16x16x4 tile in both, so in single precision the two are the same code.
 *
 * -DMFMA_F64_USE_16X16 makes TILE 0 the 16x16x4 tile in double precision as
 * well, which collapses the tile dimension; the sweep then offers TILE 0
 * alone, see mfma_tile_offered(). It is retained as a way to build without
 * the 4x4x4 path at all, not as the way to choose between them -- that is now
 * the autotuner's job.
 *
 * Both tiles are hardware verified against a CPU reference on gfx90a
 * (mfma_probe, 144 configurations per tile) and give bit-identical f64
 * results, so which one runs is a performance question only.
 */
template< typename T, const int LX, const int AXIS,
          const bool TRANSPOSE, const bool ACCUM, const int NWF = 1,
          const int TILE = 0, const int PAD = 0 >
struct mfma_contract_sel {
  __device__ __forceinline__
  static void run(T * __restrict__ out, const T * __restrict__ dmat,
                  const T * __restrict__ in, const int lane,
                  const int wf = 0) {
    mfma_contract<T, LX, AXIS, TRANSPOSE, ACCUM, NWF, PAD>(out, dmat, in,
                                                           lane, wf);
  }
};

/* Double precision, TILE 0: the batched 4x4x4 tile, or the 16x16x4 one if the
   build asked for it */
template< const int LX, const int AXIS, const bool TRANSPOSE, const bool ACCUM,
          const int NWF, const int PAD >
struct mfma_contract_sel<double, LX, AXIS, TRANSPOSE, ACCUM, NWF, 0, PAD> {
  __device__ __forceinline__
  static void run(double * __restrict__ out, const double * __restrict__ dmat,
                  const double * __restrict__ in, const int lane,
                  const int wf = 0) {
#ifdef MFMA_F64_USE_16X16
    mfma_contract<double, LX, AXIS, TRANSPOSE, ACCUM, NWF, PAD>(out, dmat, in,
                                                                lane, wf);
#else
    mfma_contract_4x4<LX, AXIS, TRANSPOSE, ACCUM, NWF, PAD>(out, dmat, in,
                                                            lane, wf);
#endif
  }
};

/* Double precision, TILE 1: the 16x16x4 tile, M-utilisation 50% at LX = 8 and
   75% at LX = 12, but at twice the FLOP rate of the batched tile, a free
   accumulate chain and a third of its operand traffic at high order */
template< const int LX, const int AXIS, const bool TRANSPOSE, const bool ACCUM,
          const int NWF, const int PAD >
struct mfma_contract_sel<double, LX, AXIS, TRANSPOSE, ACCUM, NWF, 1, PAD> {
  __device__ __forceinline__
  static void run(double * __restrict__ out, const double * __restrict__ dmat,
                  const double * __restrict__ in, const int lane,
                  const int wf = 0) {
    mfma_contract<double, LX, AXIS, TRANSPOSE, ACCUM, NWF, PAD>(out, dmat, in,
                                                                lane, wf);
  }
};

#endif // __gfx90a__ || __gfx942__
#endif // __MATH_MFMA_KERNEL_H__
