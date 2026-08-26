#ifndef __MATH_HIP_WAVE_H__
#define __MATH_HIP_WAVE_H__
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
 * Wavefront width (HIP)
 *
 * GCN and CDNA run 64 lane wavefronts, RDNA runs 32, and gfx10 and later can
 * be built either way with -mwavefrontsize32/-mwavefrontsize64. The warp level
 * shuffle reductions have to know the width at compile time, so that the
 * ladder stays the same fully unrolled instruction sequence it has always been
 * on CDNA rather than becoming a runtime loop.
 *
 * `warpSize` is not usable for that: in device code it expands to
 * __builtin_amdgcn_wavefrontsize(), which is only folded to a literal in the
 * backend, long after the preprocessor and the template instantiation that
 * would need it.
 *
 * IMPORTANT: NEKO_WAVE_SIZE is resolved per compilation *pass*, and hipcc
 * compiles every .hip twice. None of the macros it keys off exist in the host
 * pass, so the host falls back to 64 whatever the device is. It may therefore
 * only be used inside device code. Anything that has to agree between the two
 * passes -- launch geometry, kernel template arguments, __launch_bounds__ --
 * must use NEKO_WAVE_SIZE_UNIFORM instead, or the host will launch a kernel
 * instantiation the device pass never emitted.
 */
#ifdef NEKO_WAVE_SIZE

/* Set on the command line, hence identical in both passes */
#define NEKO_WAVE_SIZE_UNIFORM NEKO_WAVE_SIZE

#else

/*
 * An explicit wavefront macro is tried first: it is the only form that tracks
 * an explicit -mwavefrontsize flag. Both spellings have to be checked. The
 * bare __AMDGCN_WAVEFRONT_SIZE is the original and is what every ROCm up to
 * the 6.x series defines; __AMDGCN_WAVEFRONT_SIZE__ was added later, and ROCm
 * 6.3 deprecated both without providing a compile time replacement. They are
 * expanded here rather than at the point of use, so that a deprecation warning
 * costs one diagnostic per translation unit instead of one per reduction.
 *
 * The gfx family macros are the fallback for a toolchain that defines neither.
 * gfx6 through gfx9 (GCN, Vega, CDNA) are 64 lanes; gfx10 and later (RDNA)
 * default to 32.
 *
 * A target that matches nothing at all is a build error rather than a silent
 * 64: on RDNA that would not fail, because __shfl_down past the end of a
 * wavefront returns the lane's own value, so the reduction would quietly
 * double every partial sum instead.
 */
#  if defined(__AMDGCN_WAVEFRONT_SIZE__)
#    if __AMDGCN_WAVEFRONT_SIZE__ == 32
#      define NEKO_WAVE_SIZE 32
#    else
#      define NEKO_WAVE_SIZE 64
#    endif
#  elif defined(__AMDGCN_WAVEFRONT_SIZE)
#    if __AMDGCN_WAVEFRONT_SIZE == 32
#      define NEKO_WAVE_SIZE 32
#    else
#      define NEKO_WAVE_SIZE 64
#    endif
#  elif defined(__GFX10__) || defined(__GFX11__) || defined(__GFX12__)
#    define NEKO_WAVE_SIZE 32
#  elif defined(__GFX6__) || defined(__GFX7__) || defined(__GFX8__) ||         \
        defined(__GFX9__)
#    define NEKO_WAVE_SIZE 64
#  elif defined(__HIP_DEVICE_COMPILE__)
#    error "Unknown AMD wavefront width, rebuild with -DNEKO_WAVE_SIZE=32 or 64"
#  else
#    define NEKO_WAVE_SIZE 64 /* host pass, see the note above */
#  endif

/*
 * The pass independent width. Auto detection cannot be used here, so this
 * stays at the CDNA 64 unless the width is pinned on the command line. It only
 * feeds tuning heuristics (elem_block.h), never correctness, and the blocked
 * candidates it sizes are measured off by default on this backend anyway.
 */
#define NEKO_WAVE_SIZE_UNIFORM 64

#endif /* NEKO_WAVE_SIZE */

#if (NEKO_WAVE_SIZE != 32) && (NEKO_WAVE_SIZE != 64)
#error "NEKO_WAVE_SIZE must be 32 or 64"
#endif

/*
 * The two stage block reductions stage one partial per wavefront in a flat
 * __shared__ T [64]. Every reducing kernel is launched at the 1024 thread
 * maximum, so that is 16 entries at wave64 and 32 at wave32: the existing 64
 * covers both widths and is left alone.
 */

#endif // __MATH_HIP_WAVE_H__
