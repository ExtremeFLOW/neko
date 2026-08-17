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

/**
 * Elements per thread block for the SEM operator kstep kernels (CUDA)
 *
 * A kstep block is one (LX,LX) thread plane per element, which is half a warp
 * at LX = 4 and two warps at LX = 8. Stacking EB elements along threadIdx.z
 * packs the block up and lets the element independent derivative matrices be
 * loaded once per block rather than once per element.
 *
 * Whether that pays depends on the kernel's register footprint and on the
 * architecture, so the variants are instantiated for the candidates below and
 * picked between at runtime by each operator's autotuner. Measured on GH200
 * for ax_helm at lx = 4: 13.55 -> 8.34 us/call at EB = 4.
 *
 * Candidate C selects EB = 2^C, clamped so a block stays within
 * NEKO_EB_MAX_BLOCK_THREADS. Clamped candidates collapse onto the same
 * specialisation, which just means the tuner times it twice.
 */
#ifndef NEKO_EB_MAX_BLOCK_THREADS
#define NEKO_EB_MAX_BLOCK_THREADS 512
#endif

/* Static shared memory limit per block without the opt-in dynamic path */
#ifndef NEKO_EB_MAX_SMEM
#define NEKO_EB_MAX_SMEM 49152
#endif

#define NEKO_EB_CANDIDATES 3

template< int LX, int C >
struct elem_block {
  static const int raw = 1 << C;
  static const int value =
    (LX * LX * raw <= NEKO_EB_MAX_BLOCK_THREADS) ? raw : 1;
};

/*
 * Launch bounds.
 *
 * Max threads only, deliberately. CUDA's second argument is
 * minBlocksPerMultiprocessor, and setting it caps registers at
 * 65536/(minBlocks * LX*LX*EB) -- which under blocking gets tight enough to
 * force spilling (at lx = 8, EB = 4 with minBlocks = 3 that is 85 registers,
 * against opgrad's measured 152). The form below can only ever relax the
 * cap, never tighten it.
 *
 * Note the operator kernels historically carried __launch_bounds__(LX*LX, 3),
 * which is the HIP form where the second argument means minimum waves per
 * execution unit -- something quite different. On CUDA that squeezed
 * registers to 80 at lx >= 14 and made cdtp spill 736 bytes per thread.
 *
 * Defining NEKO_EB_THREADS_SM pins an occupancy target instead: ptxas is held
 * to 65536/NEKO_EB_THREADS_SM registers per thread and spills to meet it.
 */
#ifdef NEKO_EB_THREADS_SM
#define NEKO_EB_BOUNDS(NT)                                                    \
  __launch_bounds__((NT), ((NEKO_EB_THREADS_SM)/(NT) > 0 ?                    \
                           (NEKO_EB_THREADS_SM)/(NT) : 1))
#else
#define NEKO_EB_BOUNDS(NT) __launch_bounds__((NT))
#endif

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
