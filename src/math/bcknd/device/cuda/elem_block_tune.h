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
 * A tune function using these macros is expected to have `start`, `stop` and
 * `stream` in scope, plus a CASE_1D(LX) macro and a kstep launch macro taking
 * (LX, C).
 */

#include <stdlib.h>

/* Sweep the elements per block candidates?
   Measured a 1.6x win on GH200 and
   a loss on MI250X/MI300A, hence the per backend default

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

/* Forced candidate, used when NEKO_AUTOTUNE pins the kstep variant */
static int neko_eb_env()
{
  const char *v = getenv("NEKO_EB");
  int c = (v != NULL) ? atoi(v) : 0;

  if (c < 0 || c >= NEKO_EB_CANDIDATES) {
    c = 0;
  }
  return c;
}

/* Forced chunk candidate, used when NEKO_AUTOTUNE pins the 1d variant */
static int neko_chunks_env()
{
  const char *v = getenv("NEKO_CHUNKS");
  int c = (v != NULL) ? atoi(v) : 0;

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

/* One timed round of LAUNCH at candidate C, min reduced into T[C] */
#define NEKO_TUNE_TIME(T, LAUNCH, LX, C, ITERS)                           \
  do {                                                                        \
    float t_;                                                                 \
    cudaEventRecord(start, stream);                                           \
    for (int i_ = 0; i_ < (ITERS); i_++) { LAUNCH(LX, C); }                   \
    cudaEventRecord(stop, stream);                                            \
    cudaEventSynchronize(stop);                                               \
    cudaEventElapsedTime(&t_, start, stop);                                   \
    if (t_ < (T)[C]) { (T)[C] = t_; }                                         \
  } while (0)

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
              NEKO_CHUNKS_SEL(LX, c), (T1)[c] * 10.0);                        \
      log_message(neko_log_buf);                                              \
    }                                                                         \
    for (int c = 0; c < NEKO_EB_CANDIDATES; c++) {                            \
      if ((T2)[c] >= NEKO_TUNE_INIT) { continue; }                         \
      sprintf(neko_log_buf, "KSTEP eb=%-4d: %9.2f us/call",                   \
              NEKO_EB_SEL(LX, c), (T2)[c] * 10.0);                            \
      log_message(neko_log_buf);                                              \
    }                                                                         \
  } while (0)

#endif // __MATH_ELEM_BLOCK_TUNE_H__
