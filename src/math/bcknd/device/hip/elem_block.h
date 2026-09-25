#ifndef __MATH_ELEM_BLOCK_H__
#define __MATH_ELEM_BLOCK_H__
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

#include "wave.h"

/**
 * Elements per thread block for the SEM operator kstep kernels (HIP)
 *
 * A kstep block is one (LX,LX) thread plane per element. Against a wavefront
 * that masks off a lot of lanes for every LX whose square is not a multiple
 * of the wave: at the 64 lanes of CDNA, 25% useful at LX = 4, 39% at LX = 5,
 * 63% at LX = 9.
 * Stacking EB elements along threadIdx.z packs the block back up.
 *
 * MEASURED TO LOSE on gfx90a and gfx942 for ax_helm: the blocked
 * specialisations roughly double VGPR usage and spill to scratch
 * (kstep_padded<double,8,4>: 168 VGPRs, 2596 B/lane scratch, against 110/0
 * unblocked), and the cause of the register blow-up is still not understood.
 * A single precision run at lx = 8 puts numbers on it: 130.7 us/call at one
 * element per block against 1212 us at four and 1213 at eight, a factor of
 * 9.3.
 *
 * The sweep is nevertheless on by default here as of 2026-08-22, as it is on
 * CUDA. Encoding the verdict in a default meant the tuner could never
 * re-measure it, and the ranking does invert on NVIDIA, where blocking
 * measures a 1.6x win. The tuner rejects the blocked variants on its own; the
 * price is tuning time, and at 9.3x slower those two candidates are most of
 * it -- roughly 0.7 s per polynomial order at the default sampling.
 *
 * Candidate 0 is one element per block. Candidate C > 0 fits as many whole
 * elements as it can into 2^(C+1) wavefronts.
 *
 * NEKO_EB_WAVE is deliberately NEKO_WAVE_SIZE_UNIFORM and not NEKO_WAVE_SIZE:
 * it sizes a kernel template argument and the matching launch geometry, which
 * the host and device passes have to agree on, and only the command line form
 * of the width is visible to both. On RDNA, pass -DNEKO_WAVE_SIZE=32 to size
 * the candidates for a 32 lane wave; getting it wrong costs occupancy, not
 * correctness, and this sweep defaults off on AMD regardless.
 */
#ifndef NEKO_EB_WAVE
#define NEKO_EB_WAVE NEKO_WAVE_SIZE_UNIFORM
#endif

/* LDS per CU, and the per workgroup maximum, on CDNA */
#ifndef NEKO_EB_MAX_LDS
#define NEKO_EB_MAX_LDS 65536
#endif

#define NEKO_EB_CANDIDATES 3

template< int LX, int C >
struct elem_block {
  static const int waves = 2 << C;
  static const int fit = (NEKO_EB_WAVE * waves)/(LX * LX);
  static const int value = (C == 0) ? 1 : (fit > 0 ? fit : 1);
};

/*
 * Launch bounds.
 *
 * HIP's second argument is MIN_WARPS_PER_EXECUTION_UNIT -- minimum wavefronts
 * resident per SIMD -- and not CUDA's min blocks per SM. The first argument
 * must track the real block size: launching more threads than declared is a
 * launch failure, not a slowdown.
 *
 * The second argument defaults to 1, i.e. no constraint. It carried 3
 * historically, which holds the compiler to a 512/3 VGPR budget on every
 * kstep launch. The kstep kernels keep ru[LX] and rw[LX] -- 16 doubles, 32
 * VGPRs -- live across a fully unrolled k loop that also issues seven metric
 * loads per step, and the vector variants keep six such arrays, ~96 VGPRs,
 * before any addressing. That is the shape that spills to scratch under a
 * cap.
 *
 * Measured on MI250X, lx = 8, 13824 elements, on two builds differing only
 * in this setting: kstep eb=1 runs 1317 us with the cap at 3 and 456 us
 * without it, eb=4 4824 against 2254. Over the same pair the 1d kernel,
 * which declares no launch bounds at all, moves 449.82 -> 455.53, and the
 * MFMA kernels, which carry a bare __launch_bounds__ with no second
 * argument, move less than 4%. Only the kernels this macro touches moved.
 *
 * Note what the fix does not buy: kstep eb=1 uncapped only draws level with
 * the 1d kernel (456 against 450), because both are at the DRAM roof. It
 * matters for the vector operator, which has no 1d variant at all, and for
 * not paying for spilled candidates on every tune.
 *
 * The CUDA header records the same failure mode from its own second argument
 * -- "squeezed registers to 80 at lx >= 14 and made cdtp spill 736 bytes per
 * thread" -- and dropped it; see cuda/elem_block.h.
 *
 * Set -DNEKO_EB_MIN_WAVES_EU=3 to restore the old behaviour, and check
 * ScratchSize with -Rpass-analysis=kernel-resource-usage when changing it.
 */
#ifndef NEKO_EB_MIN_WAVES_EU
#define NEKO_EB_MIN_WAVES_EU 1
#endif

/* A minimum of one wave per SIMD is no constraint at all, so emit the
   single argument form rather than a bound that only looks like one. */
#if NEKO_EB_MIN_WAVES_EU > 1
#define NEKO_EB_BOUNDS(NT) __launch_bounds__((NT), NEKO_EB_MIN_WAVES_EU)
#else
#define NEKO_EB_BOUNDS(NT) __launch_bounds__((NT))
#endif

/**
 * Barrier for the kstep k loop
 *
 * A kstep block is dim3(LX, LX, EB), and every LDS buffer the k loop touches
 * -- shu, shur, shus, and the nine of the vector variant -- is indexed from
 * eb*LX*LX. One private slice per element, nothing shared between slices.
 * The only cross slice data is shdx/shdy/shdz, written by eb == 0 *before*
 * the loop, and the __syncthreads() guarding that stays where it is.
 *
 * So when a slice cannot straddle a wavefront, the two barriers inside the
 * loop are synchronising waves that share nothing, and a wave scoped barrier
 * is enough. The alignment test is exact rather than conservative: the linear
 * thread id of slice z runs [z*LX*LX, (z+1)*LX*LX), so slices land on
 * wavefront boundaries exactly when LX*LX divides the wave. On wave64 that
 * admits LX = 2, 4 and 8.
 *
 * Note LX = 4, where one wave holds four whole slices: that is also where the
 * packing has something to gain, 16 of 64 lanes busy unblocked. This is not a
 * knob that is only legal where it is useless -- unlike readfirstlane on the
 * element base, whose validity condition picks out exactly the orders with no
 * upside.
 *
 * OFF BY DEFAULT, and MEASURED TO BUY NOTHING at lx = 8. On MI250X, 13824
 * elements: kstep eb=4 runs 2254, 2280 and 2240 us over three builds, one
 * of which has this on. If the 2*LX workgroup barriers were what the
 * blocked variants were paying for, removing them should have collapsed
 * eb=4 back toward eb=1 -- a factor of five. It did not move outside the
 * noise.
 *
 * Calibrate against that noise before reading anything into a single pair.
 * Over three nominally identical runs the spread is 1.7% on kstep eb=1,
 * 1.8% on eb=4, but 8.1% on the MFMA kernels and 8.4% on kstep eb=8 --
 * and the MFMA kernels contain no barrier at all, so 8% swings there are
 * the machine, not the code. Anything under ~10% on the low occupancy
 * kernels is not a result.
 *
 * So the residual eb cost is NOT barrier scope -- look to the divergent
 * element base (nelv enters the address, defeating the gfx9 SADDR form)
 * and to the hardware co-residency floor at EB = 8, where 8 wave64s on 4
 * SIMDs cap occupancy at 2 waves/SIMD regardless of any launch bound.
 *
 * Before trusting that refutation, confirm the mechanism actually engaged:
 * the s_barrier count in the ISA must drop between a build with this off
 * and one with it on. Equal counts mean the barrier was never removed and
 * the experiment said nothing.
 *
 * Kept because it is off, cheap and documented, and because lx = 4 is
 * untested -- though the prior is now weak, since lx = 4 blocking is also
 * 4 wavefronts, the same shape that showed nothing at lx = 8.
 *
 * It weakens a synchronisation scope, so it fails as wrong results rather
 * than as an error, and wrong results inside a Krylov solve present as poor
 * convergence rather than as anything that looks like a bug. The benchmark
 * reports no correctness signal at all. Validate by pinning two strategies
 * and requiring bit identical output. Never by checking that the solver
 * still converges.
 *
 * Ordering is spelled out rather than left to an inline asm memory clobber.
 * A clobber constrains the compiler and says nothing to the memory model, and
 * __builtin_amdgcn_wave_barrier() is itself only a scheduling barrier -- it
 * emits no instruction. The release/acquire pair around it is what states the
 * requirement, the same shape the backend lowers a workgroup barrier into,
 * at wavefront scope instead.
 */
#ifndef NEKO_EB_WAVE_BARRIER
#define NEKO_EB_WAVE_BARRIER 0
#endif

/*
 * Is a kstep slice wavefront contained, so that a wave barrier covers it?
 *
 * Deliberately NOT keyed on NEKO_EB_WAVE. That macro is the command line
 * width, and it sizes kernel template arguments and launch geometry, where
 * the host and device passes have to agree on one number; getting it wrong
 * there costs occupancy, not correctness. Here it would cost correctness, so
 * this asks the target instead of being told.
 *
 * __builtin_amdgcn_wavefrontsize() is the supported spelling.
 * __AMDGCN_WAVEFRONT_SIZE__ is deprecated -- clang warns that compile-time
 * constant access to the wavefront size is going away -- so there is no
 * static_assert to be had, and wanting one was the wrong instinct: under a
 * fixed --offload-arch the builtin folds to a constant and the branch below
 * disappears, and under a multi target fatbin it correctly resolves per
 * target, which a build time check could not do at all.
 */
__device__ __host__ constexpr bool neko_eb_slice_in_wave(unsigned wave, int lx)
{
  return ((unsigned) (lx * lx) <= wave) &&
         ((wave % (unsigned) (lx * lx)) == 0);
}

template< int LX >
__device__ __forceinline__ void neko_eb_kstep_barrier()
{
#if NEKO_EB_WAVE_BARRIER && defined(__HIP_DEVICE_COMPILE__)
  /* Uniform across the whole block, so the __syncthreads() arm is never
     reached by only part of it */
  if (neko_eb_slice_in_wave(__builtin_amdgcn_wavefrontsize(), LX)) {
    __builtin_amdgcn_fence(__ATOMIC_RELEASE, "wavefront");
    __builtin_amdgcn_wave_barrier();
    __builtin_amdgcn_fence(__ATOMIC_ACQUIRE, "wavefront");
  } else {
    __syncthreads();
  }
#else
  __syncthreads();
#endif
}

/*
 * Launch geometry helpers. NELV is the element count; the grid is sized so
 * the tail block is partially filled and clamps, rather than the kernel
 * returning early -- these kernels have __syncthreads() in the k loop, so an
 * early return would leave the barrier unmatched.
 */
#define NEKO_EB(LX, C) (elem_block<LX, C>::value)
#define NEKO_EB_NTHRDS(LX, C) dim3((LX), (LX), NEKO_EB(LX, C))
#define NEKO_EB_NBLCKS(NELV, LX, C)                                           \
  dim3(((NELV) + NEKO_EB(LX, C) - 1)/NEKO_EB(LX, C), 1, 1)
#define NEKO_EB_SEL(LX, SEL)                                                  \
  ((SEL) == 0 ? NEKO_EB(LX, 0) : (SEL) == 1 ? NEKO_EB(LX, 1) : NEKO_EB(LX, 2))

/**
 * Chunk size for the 1d kernels
 *
 * CHUNKS is both the thread block size and the stride over the LX^3 points of
 * an element. At the historical 1024 that is badly mismatched at low order:
 * only 64 of 1024 threads do work at LX = 4, 125 at LX = 5, 512 at LX = 8.
 * Shared memory is sized by LX rather than by CHUNKS, so a smaller block also
 * raises the number of elements resident per SM at no extra cost.
 *
 * CONSTRAINT: the 1d kernels stage the derivative matrices with a single
 * `if (iii < LX*LX)` guard, so a block smaller than one matrix would leave
 * part of it unwritten. Candidates below LX*LX are therefore rejected
 * outright and fall back to 1024, which always satisfies it for LX <= 16.
 */
#define NEKO_CHUNKS_CANDIDATES 4

/*
 * Candidate C selects 1024 >> C, i.e. 1024, 512, 256, 128. Candidate 0 is the
 * historical value, so it is always available as an A/B baseline.
 *
 * A heuristic is deliberately avoided here. The best chunk depends on how the
 * shared memory footprint (which grows as LX^3) caps the resident block count
 * on the target device, and that differs between architectures -- 164 kB per
 * SM on A100, 228 kB on H100, 64 kB of LDS per CU on CDNA. Picking by thread
 * utilisation alone gets it wrong from LX = 6 upward, because once the block
 * count is capped by shared memory a smaller block simply wastes thread slots.
 * So all four are instantiated and the tuner measures them.
 */
template< int LX, int C >
struct chunk_block {
  static const int raw = 1024 >> C;
  static const int value = (raw >= LX * LX) ? raw : 1024;
};

#define NEKO_CHUNKS(LX, C) (chunk_block<LX, C>::value)
#define NEKO_CHUNKS_NTHRDS(LX, C) dim3(NEKO_CHUNKS(LX, C), 1, 1)
#define NEKO_CHUNKS_SEL(LX, SEL)                                              \
  ((SEL) == 0 ? NEKO_CHUNKS(LX, 0) : (SEL) == 1 ? NEKO_CHUNKS(LX, 1) :        \
   (SEL) == 2 ? NEKO_CHUNKS(LX, 2) : NEKO_CHUNKS(LX, 3))

#endif // __MATH_ELEM_BLOCK_H__
