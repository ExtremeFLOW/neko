#ifndef __MATH_ELEM_BLOCK_TUNE_H__
#define __MATH_ELEM_BLOCK_TUNE_H__
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
 * Shared autotuner scaffolding for the SEM operator kernels
 *
 * Each operator's tune() picks between the 1d variant and the kstep variant,
 * and -- when sweeping is enabled -- between the elements per block
 * candidates in elem_block.h. The sampling is interleaved and min reduced:
 * timing the variants one after another in a fixed order lets clock drift
 * bias the comparison by candidate position, which no amount of extra
 * iterations removes.
 *
 * NEKO_AUTOTUNE takes a formulation out of the tuner's hands, but only the
 * formulation: the geometry candidates of the one it names are still measured
 * and reported, and it takes that formulation's own variable (NEKO_EB,
 * NEKO_CHUNKS, ...) to fix the geometry too. Pinning the formulation used to
 * imply candidate 0, which silently answered a different question than the one
 * an A/B run of two formulations is asking. See NEKO_TUNE_FOR() below.
 *
 * A tune function using these macros is expected to have `start`, `stop` and
 * `stream` in scope --- and `iters` for the reporting macros --- plus a
 * CASE_1D(LX) macro and a kstep launch macro taking (LX, C).
 */

#include <stdlib.h>

/* Sweep the elements per block candidates?
   On by default, as on CUDA. It was off here because the blocked variants
   measured a loss on MI250X and MI300A, where they spill -- but baking that
   in also stopped the tuner from ever re-testing it, and the kstep family has
   since measured badly enough on gfx90a (788 us/call at lx = 8 against 268
   for the 1d variant) that its geometry is worth measuring per case rather
   than assuming. Where the blocked variants still spill the tuner simply
   rejects them; the cost is tuning time, not run time.

   Note NEKO_TUNE_ROUNDS and NEKO_TUNE_ITERS below are not specific to
   this sweep: they control the sampling of every candidate the tuner
   times, chunk sizes included. */
#ifndef NEKO_EB_SWEEP_DEFAULT
#define NEKO_EB_SWEEP_DEFAULT 1
#endif

#define NEKO_TUNE_INIT 1.0e30f
#define NEKO_TUNE_WARMUP 20

static int neko_eb_sweep()
{
  const char *v = getenv("NEKO_EB_TUNE");

  if (v != NULL) {
    return (atoi(v) != 0);
  }
  return NEKO_EB_SWEEP_DEFAULT;
}

/*
 * Elements per block candidate pinned by NEKO_EB, or -1 to leave it to the
 * sweep.
 *
 * NEKO_AUTOTUNE selects the formulation and nothing more -- the geometry
 * inside it is still measured unless it is pinned here -- so "unset" has to
 * be distinguishable from candidate 0, which is a candidate like any other.
 * Hence the -1 rather than a plain default of 0. Out of range values clamp
 * to 0 rather than releasing the pin, as they always have.
 */
static int neko_eb_pin()
{
  const char *v = getenv("NEKO_EB");
  int c;

  if (v == NULL) {
    return -1;
  }

  c = atoi(v);
  if (c < 0 || c >= NEKO_EB_CANDIDATES) {
    c = 0;
  }
  return c;
}

/* Chunk candidate pinned by NEKO_CHUNKS, or -1 to sweep, see neko_eb_pin() */
static int neko_chunks_pin()
{
  const char *v = getenv("NEKO_CHUNKS");
  int c;

  if (v == NULL) {
    return -1;
  }

  c = atoi(v);
  if (c < 0 || c >= NEKO_CHUNKS_CANDIDATES) {
    c = 0;
  }
  return c;
}

static int neko_tune_rounds()
{
  const char *v = getenv("NEKO_TUNE_ROUNDS");
  int n = (v != NULL) ? atoi(v) : 3;

  return (n < 1) ? 1 : n;
}

static int neko_tune_iters()
{
  const char *v = getenv("NEKO_TUNE_ITERS");
  int n = (v != NULL) ? atoi(v) : 100;

  return (n < 1) ? 1 : n;
}

/*
 * Loop over the candidates of one formulation, binding C to each in turn.
 *
 * ON gates the formulation itself: hardware support, and -- when NEKO_AUTOTUNE
 * is set -- whether it is the formulation that was asked for. PIN is the
 * candidate that formulation's own variable forces, or -1 to measure all N of
 * them. So NEKO_AUTOTUNE narrows the search to one kernel family without
 * deciding the geometry within it, which is still swept and reported.
 *
 * An out of play formulation yields an empty range, leaving its candidates at
 * NEKO_TUNE_INIT so it loses the comparison at the end rather than having to
 * be excluded from it.
 */
#define NEKO_TUNE_FOR(C, ON, PIN, N)                                          \
  for (int C = (((PIN) >= 0) ? (PIN) : 0),                                    \
       C##_last_ = ((ON) ? (((PIN) >= 0) ? (PIN) + 1 : (N)) : 0);             \
       C < C##_last_; C++)

/* One timed round of LAUNCH at candidate C, min reduced into T[C] */
#define NEKO_TUNE_TIME(T, LAUNCH, LX, C, ITERS)                           \
  do {                                                                        \
    float t_;                                                                 \
    HIP_CHECK(hipEventRecord(start, stream));                                           \
    for (int i_ = 0; i_ < (ITERS); i_++) { LAUNCH(LX, C); }                   \
    HIP_CHECK(hipEventRecord(stop, stream));                                            \
    HIP_CHECK(hipEventSynchronize(stop));                                               \
    HIP_CHECK(hipEventElapsedTime(&t_, start, stop));                                   \
    if (t_ < (T)[C]) { (T)[C] = t_; }                                         \
  } while (0)

/*
 * Elapsed time (ms, over ITERS launches) as microseconds per call.
 *
 * A division rather than a constant because ITERS is NEKO_TUNE_ITERS, which
 * is settable: the factor of 10 this used to carry is only correct at the
 * default of 100 and silently rescales every reported time otherwise. The
 * ranking never moves, every candidate sharing the divisor, but the number
 * the log prints does.
 */
#define NEKO_TUNE_US(T, ITERS) ((T) * 1000.0 / (double) (ITERS))

#define NEKO_TUNE_BEST(T, BEST, N)                                              \
  do {                                                                        \
    for (int c = 1; c < (N); c++) {                                           \
      if ((T)[c] < (T)[BEST]) { BEST = c; }                                   \
    }                                                                         \
  } while (0)

/* Report every measured candidate of both sweeps, not just the winner */
#define NEKO_TUNE_LOG(LX, T1, T2)                                               \
  do {                                                                        \
    for (int c = 0; c < NEKO_CHUNKS_CANDIDATES; c++) {                        \
      if ((T1)[c] >= NEKO_TUNE_INIT) { continue; }                         \
      sprintf(neko_log_buf, "1D    ch=%-4d: %9.2f us/call",                   \
              NEKO_CHUNKS_SEL(LX, c), NEKO_TUNE_US((T1)[c], iters));          \
      log_message(neko_log_buf);                                              \
    }                                                                         \
    for (int c = 0; c < NEKO_EB_CANDIDATES; c++) {                            \
      if ((T2)[c] >= NEKO_TUNE_INIT) { continue; }                         \
      sprintf(neko_log_buf, "KSTEP eb=%-4d: %9.2f us/call",                   \
              NEKO_EB_SEL(LX, c), NEKO_TUNE_US((T2)[c], iters));              \
      log_message(neko_log_buf);                                              \
    }                                                                         \
  } while (0)

#endif // __MATH_ELEM_BLOCK_TUNE_H__
