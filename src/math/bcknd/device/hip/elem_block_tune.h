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
#include <stdio.h>
#include <math.h>
#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>

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

/* Report every measured candidate of one sweep, not just the winner. Split
   per formulation because the vector operator has no 1d variant and so wants
   the kstep half on its own */
#define NEKO_TUNE_LOG_1D(LX, T1)                                              \
  do {                                                                        \
    for (int c = 0; c < NEKO_CHUNKS_CANDIDATES; c++) {                        \
      if ((T1)[c] >= NEKO_TUNE_INIT) { continue; }                            \
      sprintf(neko_log_buf, "1D    ch=%-4d: %9.2f us/call",                   \
              NEKO_CHUNKS_SEL(LX, c), NEKO_TUNE_US((T1)[c], iters));          \
      log_message(neko_log_buf);                                              \
    }                                                                         \
  } while (0)

#define NEKO_TUNE_LOG_KSTEP(LX, T2)                                           \
  do {                                                                        \
    for (int c = 0; c < NEKO_EB_CANDIDATES; c++) {                            \
      if ((T2)[c] >= NEKO_TUNE_INIT) { continue; }                            \
      sprintf(neko_log_buf, "KSTEP eb=%-4d: %9.2f us/call",                   \
              NEKO_EB_SEL(LX, c), NEKO_TUNE_US((T2)[c], iters));              \
      log_message(neko_log_buf);                                              \
    }                                                                         \
  } while (0)

/* Both sweeps, for the operators that have a 1d variant */
#define NEKO_TUNE_LOG(LX, T1, T2)                                             \
  do {                                                                        \
    NEKO_TUNE_LOG_1D(LX, T1);                                                 \
    NEKO_TUNE_LOG_KSTEP(LX, T2);                                              \
  } while (0)


/*
 * TUNE-TIME VERIFICATION
 *
 * Every candidate the sweep is about to time is first run once on a synthetic
 * input and compared with the reference formulation's result on the same
 * input: the 1d kernel's first candidate where there is one, the kstep
 * kernel's where there is not. A candidate whose max-norm difference exceeds
 * NEKO_TUNE_VERIFY_TOL times the reference's max-norm is reported as REJECTED
 * and excluded from the comparison, so it can be measured but cannot win.
 *
 * Why it exists: the matrix core contraction once ran for a long time with
 * its lane layout wrong -- a wrong but symmetric operator makes CG converge
 * happily to the solution of a different system (see mfma_contract_4x4 in
 * mfma_kernel.h), so "the solver converges" never caught it. The single
 * precision 4x4x1 tile added since carries a layout that is documented by
 * AMD's calculator but has not been read out on hardware, and this is what
 * keeps it honest: a wrong register index costs that candidate at tuning,
 * loudly, not the solve.
 *
 * Why a synthetic input: the tuner runs on the first call of an operator,
 * which can be on a zero field (a zero initial guess, a zero velocity), and
 * every layout reproduces zero. The fields the operator differentiates are
 * swapped for a deterministic non-zero pattern in (0.5, 1.5) for the duration
 * of the pass -- the tune function hands over the addresses of its input
 * pointer variables, which its launch macros name, so the swap is invisible
 * to them -- while geometry, metrics and weights stay the case's own. The
 * reference max-norm is logged with every rejection, so an inconclusive pass
 * is visible as such.
 *
 * The tolerance is a layout detector, not an accuracy test: the variants sum
 * the same terms in a different order, a few ulp per point, where a layout
 * error reads as O(1). NaN or Inf anywhere in a candidate's output counts as
 * an infinite difference, so it cannot slip through a comparison that is
 * false for NaN.
 *
 * Cost: one launch per candidate plus a reduction over the output, once per
 * (operator, order), and nin + nout arrays of n reals for the duration of the
 * pass; a failed allocation skips the pass and says so in the log, it never
 * stops the tune. NEKO_TUNE_VERIFY=0 turns it off, NEKO_TUNE_VERIFY_TOL
 * overrides the tolerance.
 */
#ifndef NEKO_TUNE_VERIFY_DEFAULT
#define NEKO_TUNE_VERIFY_DEFAULT 1
#endif
#define NEKO_TUNE_VERIFY_BLOCKS 256
#define NEKO_TUNE_VERIFY_THREADS 256
#define NEKO_TUNE_VERIFY_MAXIN 8
/* Upper bound on any formulation's candidate count, for the reject flags */
#define NEKO_TUNE_MAXCAND 16

static int neko_tune_verify_on()
{
  const char *v = getenv("NEKO_TUNE_VERIFY");

  if (v != NULL) {
    return (atoi(v) != 0);
  }
  return NEKO_TUNE_VERIFY_DEFAULT;
}

static double neko_tune_verify_tol()
{
  const char *v = getenv("NEKO_TUNE_VERIFY_TOL");

  if (v != NULL && atof(v) > 0.0) {
    return atof(v);
  }
  return (sizeof(real) == 4) ? 1.0e-4 : 1.0e-10;
}

/* Deterministic non-zero pattern in (0.5, 1.5), different per input slot */
static __global__ void neko_tune_fill_kernel(real * __restrict__ x,
                                             const int n,
                                             const unsigned int seed) {
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n;
       i += gridDim.x * blockDim.x) {
    unsigned int h = (unsigned int) i * 2654435761u + seed * 40503u;
    h ^= h >> 15;
    h *= 2246822519u;
    h ^= h >> 13;
    x[i] = (real) 0.5 + (real) (h & 0xffffffu) / (real) 0x1000000u;
  }
}

/* Per block maxima of |a - b| and |b|, NaN or Inf in a counted as infinite */
static __global__ void neko_tune_diff_kernel(const real * __restrict__ a,
                                             const real * __restrict__ b,
                                             const int n,
                                             double * __restrict__ part) {
  __shared__ double sd[NEKO_TUNE_VERIFY_THREADS];
  __shared__ double sr[NEKO_TUNE_VERIFY_THREADS];
  double md = 0.0, mr = 0.0;

  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n;
       i += gridDim.x * blockDim.x) {
    const double x = (double) a[i];
    const double y = (double) b[i];
    const double d = fabs(x - y);
    const double r = fabs(y);
    if (!(d < HUGE_VAL)) {
      md = HUGE_VAL;
    } else if (d > md) {
      md = d;
    }
    if (r > mr) {
      mr = r;
    }
  }
  sd[threadIdx.x] = md;
  sr[threadIdx.x] = mr;
  __syncthreads();
  for (int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s) {
      if (sd[threadIdx.x + s] > sd[threadIdx.x]) {
        sd[threadIdx.x] = sd[threadIdx.x + s];
      }
      if (sr[threadIdx.x + s] > sr[threadIdx.x]) {
        sr[threadIdx.x] = sr[threadIdx.x + s];
      }
    }
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    part[2 * blockIdx.x] = sd[0];
    part[2 * blockIdx.x + 1] = sr[0];
  }
}

typedef struct {
  real *buf;                        /* nin synthetic inputs, nout references */
  double *part;                     /* per block maxima */
  real **outs;                      /* the tune function's output arrays */
  void ***ins;                      /* addresses of its input pointers */
  void *saved[NEKO_TUNE_VERIFY_MAXIN];
  int nin, nout, n;
  double tol;
} neko_tune_verify_t;

/* Reject flags of one formulation's candidates, see NEKO_TUNE_VERIFY_CAND */
typedef struct {
  int rej[NEKO_TUNE_MAXCAND];
} neko_tune_rej_t;

/*
 * Set the pass up: allocate, fill the synthetic inputs and point the tune
 * function's input pointer variables at them. Returns 0, with nothing
 * swapped, if the pass is off or cannot be set up.
 */
static int neko_tune_verify_begin(neko_tune_verify_t *v, real **outs,
                                  const int nout, void ***ins, const int nin,
                                  const int n, hipStream_t stream)
{
  const size_t bytes = (size_t) (nin + nout) * (size_t) n * sizeof(real);

  if (!neko_tune_verify_on() || n <= 0 || nin > NEKO_TUNE_VERIFY_MAXIN) {
    return 0;
  }
  v->outs = outs;
  v->ins = ins;
  v->nin = nin;
  v->nout = nout;
  v->n = n;
  v->tol = neko_tune_verify_tol();
  v->buf = NULL;
  v->part = NULL;
  if (hipMalloc((void **) &v->buf, bytes) != hipSuccess) {
    /* Clear the sticky error: the tuner checks hipGetLastError() after
       every launch and would otherwise report this one there */
    (void) hipGetLastError();
    return 0;
  }
  if (hipMalloc((void **) &v->part,
                2 * NEKO_TUNE_VERIFY_BLOCKS * sizeof(double)) != hipSuccess) {
    (void) hipGetLastError();
    HIP_CHECK(hipFree(v->buf));
    return 0;
  }
  for (int i = 0; i < nin; i++) {
    real *x = v->buf + (size_t) i * n;
    hipLaunchKernelGGL(neko_tune_fill_kernel, dim3(NEKO_TUNE_VERIFY_BLOCKS),
                       dim3(NEKO_TUNE_VERIFY_THREADS), 0, stream,
                       x, n, (unsigned int) (i + 1));
    HIP_CHECK(hipGetLastError());
    v->saved[i] = *ins[i];
    *ins[i] = (void *) x;
  }
  return 1;
}

/* Snapshot the outputs of the reference launch the caller has just made */
static void neko_tune_verify_ref(neko_tune_verify_t *v, hipStream_t stream)
{
  for (int j = 0; j < v->nout; j++) {
    /* hipMemcpyDefault rather than DeviceToDevice: the runtime then
       resolves the kind from the pointers, which keeps this correct if
       the operator's arrays are ever host-mapped on an APU */
    HIP_CHECK(hipMemcpyAsync(v->buf + (size_t) (v->nin + j) * v->n,
                             v->outs[j], (size_t) v->n * sizeof(real),
                             hipMemcpyDefault, stream));
  }
  HIP_CHECK(hipStreamSynchronize(stream));
}

/*
 * Compare the outputs of the launch the caller has just made with the
 * snapshot. Returns 1 if every output is within tolerance, and leaves the
 * worst max-norm difference and the reference max-norm behind for the log
 * either way. Written so that an infinite difference fails.
 */
static int neko_tune_verify_check(neko_tune_verify_t *v, hipStream_t stream,
                                  double *maxdiff, double *maxref)
{
  double part[2 * NEKO_TUNE_VERIFY_BLOCKS];

  *maxdiff = 0.0;
  *maxref = 0.0;
  for (int j = 0; j < v->nout; j++) {
    const real *ref = v->buf + (size_t) (v->nin + j) * v->n;
    hipLaunchKernelGGL(neko_tune_diff_kernel, dim3(NEKO_TUNE_VERIFY_BLOCKS),
                       dim3(NEKO_TUNE_VERIFY_THREADS), 0, stream,
                       (const real *) v->outs[j], ref, v->n, v->part);
    HIP_CHECK(hipGetLastError());
    HIP_CHECK(hipMemcpyAsync(part, v->part, sizeof(part),
                             hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipStreamSynchronize(stream));
    for (int b = 0; b < NEKO_TUNE_VERIFY_BLOCKS; b++) {
      if (part[2 * b] > *maxdiff) {
        *maxdiff = part[2 * b];
      }
      if (part[2 * b + 1] > *maxref) {
        *maxref = part[2 * b + 1];
      }
    }
  }
  return (*maxdiff <= v->tol * (*maxref));
}

/* Restore the input pointers and release the buffers */
static void neko_tune_verify_end(neko_tune_verify_t *v, hipStream_t stream)
{
  HIP_CHECK(hipStreamSynchronize(stream));
  for (int i = 0; i < v->nin; i++) {
    *v->ins[i] = v->saved[i];
  }
  HIP_CHECK(hipFree(v->buf));
  HIP_CHECK(hipFree(v->part));
}

/* Launch candidate C of LAUNCH, compare it with the reference and record the
   verdict in REJ; a rejection is logged with the two norms behind it. Expects
   `stream` and `neko_log_buf` in scope like the reporting macros above */
#define NEKO_TUNE_VERIFY_CAND(VF, REJ, LAUNCH, LX, C, NAME)                   \
  do {                                                                        \
    double md_, mr_;                                                          \
    LAUNCH(LX, C);                                                            \
    (REJ).rej[C] = !neko_tune_verify_check(&(VF), stream, &md_, &mr_);        \
    if ((REJ).rej[C]) {                                                       \
      sprintf(neko_log_buf,                                                   \
              "%-5s cand %2d REJECTED: |diff| %.1e, |ref| %.1e",              \
              NAME, C, md_, mr_);                                             \
      log_message(neko_log_buf);                                              \
    }                                                                         \
  } while (0)

/* Take the rejected candidates of one formulation out of the comparison */
#define NEKO_TUNE_APPLY_REJ(T, REJ, N)                                        \
  do {                                                                        \
    for (int c = 0; c < (N); c++) {                                           \
      if ((REJ).rej[c]) { (T)[c] = NEKO_TUNE_INIT; }                          \
    }                                                                         \
  } while (0)

#endif // __MATH_ELEM_BLOCK_TUNE_H__
