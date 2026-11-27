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
 * Metal host-side dispatch for the compressible residual.
 *
 * @note Apple GPUs do not support FP64. This backend operates in FP32.
 */

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <device/device_config.h>
#include <device/metal/check.h>
#include <common/neko_log.h>

/* Defined in device/metal/metal.m */
extern id<MTLDevice> neko_metal_device(void);
extern id<MTLLibrary> neko_metal_library(void);

/* Cached pipeline states */
static id<MTLComputePipelineState> pso_visc = nil;
/* Fused inviscid divergence, indexed by LX and tuning candidate */
#define INVISCID_DIV_CANDIDATES 16
static id<MTLComputePipelineState>
  pso_inviscid_div[17][INVISCID_DIV_CANDIDATES] = { { nil } };
/* Candidate chosen by the tuner plus one, zero until it has run */
static int autotune_inviscid_div[17] = { 0 };
static id<MTLComputePipelineState> pso_rk_sum = nil;

/**
 * Create a compute pipeline state for the named kernel, compiled to support
 * max_threads threads per threadgroup, see neko_metal_new_pipeline.
 */
static id<MTLComputePipelineState>
get_compressible_res_pipeline(const char *name, const NSUInteger max_threads) {
  id<MTLLibrary> lib = neko_metal_library();
  NSString *nsName = [NSString stringWithUTF8String:name];
  id<MTLFunction> func = [lib newFunctionWithName:nsName];
  if (func == nil) {
    fprintf(stderr, "Metal: kernel '%s' not found in metallib\n", name);
    exit(EXIT_FAILURE);
  }
  NSError *error = nil;
  id<MTLComputePipelineState> pso =
    neko_metal_new_pipeline(func, max_threads, &error);
  METAL_CHECK(error);
  return pso;
}

/**
 * Create a compute command encoder on the global command queue.
 */
static id<MTLComputeCommandEncoder>
compressible_res_encoder(id<MTLCommandBuffer> *cmdBuf) {
  id<MTLCommandQueue> queue =
    (__bridge id<MTLCommandQueue>)glb_cmd_queue;
  *cmdBuf = [queue commandBuffer];
  return [*cmdBuf computeCommandEncoder];
}

/**
 * Dispatch one thread per point and wait for completion.
 */
static void compressible_res_dispatch(id<MTLCommandBuffer> cmdBuf,
                               id<MTLComputeCommandEncoder> enc, int n) {
  const NSUInteger nthrds = 256;
  MTLSize groupSize = MTLSizeMake(nthrds, 1, 1);
  MTLSize numGroups = MTLSizeMake(((NSUInteger)n + nthrds - 1) / nthrds, 1, 1);

  [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
  METAL_CHECK_CB(cmdBuf);
}

/**
 * Viscous part of the compressible residual on the Metal GPU.
 */
void compressible_res_part_visc_metal(void *rhs_rho, void *rhs_m_x,
                               void *rhs_m_y, void *rhs_m_z, void *rhs_E,
                               void *visc_rho, void *visc_m_x, void *visc_m_y,
                               void *visc_m_z, void *visc_E,
                               void *mult, void *Binv, void *h1, int *n) {

  if (*n <= 0)
    return;

  if (pso_visc == nil)
    pso_visc =
      get_compressible_res_pipeline("compressible_res_part_visc_kernel",
                                    NEKO_METAL_MAX_THREADS);

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = compressible_res_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_visc];

  [enc setBuffer:(__bridge id<MTLBuffer>)rhs_rho  offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)rhs_m_x  offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)rhs_m_y  offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)rhs_m_z  offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)rhs_E    offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)visc_rho offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)visc_m_x offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)visc_m_y offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)visc_m_z offset:0 atIndex:8];
  [enc setBuffer:(__bridge id<MTLBuffer>)visc_E   offset:0 atIndex:9];
  [enc setBuffer:(__bridge id<MTLBuffer>)mult     offset:0 atIndex:10];
  [enc setBuffer:(__bridge id<MTLBuffer>)Binv     offset:0 atIndex:11];
  [enc setBuffer:(__bridge id<MTLBuffer>)h1       offset:0 atIndex:12];
  [enc setBytes:n length:sizeof(int) atIndex:13];

  compressible_res_dispatch(cmdBuf, enc, *n);
}

/*
 * Tuning candidates of the fused divergence
 *
 * The kernel's geometry is measured rather than decided: the tuner below
 * times the candidates on the first call at each order and keeps the
 * fastest. Candidate c combines three choices:
 *
 *   chunks (bit 3)       the LX^3 points of an element are split over as few
 *                        chunks of at most INVISCID_DIV_MAX_BLOCK threads as
 *                        they need, or twice as many. The narrower
 *                        threadgroup is not offered where it would put more
 *                        than INVISCID_DIV_MAX_PTS points on a thread.
 *   staging (bit 2)      three fluxes staged in threadgroup memory per pass,
 *                        where they fit next to the derivative matrices, or
 *                        one.
 *   registers (bits 0-1) whether the point state (bit 0) and the geometric
 *                        factors (bit 1) of a thread's points are held in
 *                        registers across the passes or re-read from device
 *                        memory. Neither is offered beyond
 *                        INVISCID_DIV_MAX_REG_PTS points per thread. On a
 *                        W5700X the winner flipped with the order: neither
 *                        at order 8, both at order 10, the state alone at
 *                        order 12, by 1.4-1.8x over the worst.
 *
 * A candidate whose choice is not offered takes the geometry of the one that
 * is, so candidates that round to the same geometry name the same kernel;
 * the tuner skips the repeats. The kernel names are generated from these
 * rules, see the instantiation list in compressible_res_kernel.metal.
 */
#define INVISCID_DIV_WARP 32
#define INVISCID_DIV_MAX_BLOCK 512
#define INVISCID_DIV_MAX_SMEM 32768
#define INVISCID_DIV_MAX_PTS 8
#define INVISCID_DIV_MAX_REG_PTS 4
#define INVISCID_DIV_TUNE_WARMUP 20

/* Threads per threadgroup covering the element in nchunks chunks: the
   points of one chunk, rounded up to whole SIMD groups */
static int inviscid_div_threads_of(const int lx, const int nchunks) {
  return ((((lx * lx * lx + nchunks - 1) / nchunks)
           + INVISCID_DIV_WARP - 1) / INVISCID_DIV_WARP) * INVISCID_DIV_WARP;
}

/* Points per thread at that threadgroup */
static int inviscid_div_pts_of(const int lx, const int nchunks) {
  return (lx * lx * lx + inviscid_div_threads_of(lx, nchunks) - 1)
    / inviscid_div_threads_of(lx, nchunks);
}

/* Chunks the element's points are split over for candidate c */
static int inviscid_div_nchunks(const int lx, const int c) {
  const int min_chunks = (lx * lx * lx + INVISCID_DIV_MAX_BLOCK - 1)
    / INVISCID_DIV_MAX_BLOCK;
  return (((c >> 3) & 1) != 0 &&
          inviscid_div_pts_of(lx, 2 * min_chunks) <= INVISCID_DIV_MAX_PTS)
    ? 2 * min_chunks : min_chunks;
}

/* Threads per threadgroup of candidate c */
static int inviscid_div_chunks(const int lx, const int c) {
  return inviscid_div_threads_of(lx, inviscid_div_nchunks(lx, c));
}

/* Fluxes staged per pass */
static int inviscid_div_stage(const int lx, const int c) {
  return (((c >> 2) & 1) == 0 &&
          (3 * lx * lx * lx + 3 * lx * lx) * (int) sizeof(float)
          <= INVISCID_DIV_MAX_SMEM) ? 3 : 1;
}

/* Register choice: bit 0 the state, bit 1 the geometric factors */
static int inviscid_div_reg(const int lx, const int c) {
  return (inviscid_div_pts_of(lx, inviscid_div_nchunks(lx, c))
          <= INVISCID_DIV_MAX_REG_PTS) ? (c & 3) : 0;
}

/* Does candidate c name the same kernel as a lower candidate? */
static int inviscid_div_cand_dup(const int lx, const int c) {
  for (int d = 0; d < c; d++) {
    if (inviscid_div_chunks(lx, d) == inviscid_div_chunks(lx, c) &&
        inviscid_div_stage(lx, d) == inviscid_div_stage(lx, c) &&
        inviscid_div_reg(lx, d) == inviscid_div_reg(lx, c)) {
      return 1;
    }
  }
  return 0;
}

/* Describe candidate c for the log */
static const char *inviscid_div_cand_desc(const int lx, const int c,
                                          char *buf, const int len) {
  static const char *regs[] = { "none", "state", "metrics", "both" };
  snprintf(buf, len, "%d thr, %d/pass, regs: %s", inviscid_div_chunks(lx, c),
           inviscid_div_stage(lx, c), regs[inviscid_div_reg(lx, c)]);
  return buf;
}

/* Candidate pinned by NEKO_INVISCID_DIV, or -1 to sweep */
static int neko_inviscid_div_pin(void) {
  const char *v = getenv("NEKO_INVISCID_DIV");
  int c;

  if (v == NULL) {
    return -1;
  }

  c = atoi(v);
  if (c < 0 || c >= INVISCID_DIV_CANDIDATES) {
    c = 0;
  }
  return c;
}

/* Sampling of the sweep and the verification pass, as for the CUDA and HIP
   operator tuners */
static int neko_tune_rounds(void) {
  const char *v = getenv("NEKO_TUNE_ROUNDS");
  int n = (v != NULL) ? atoi(v) : 3;

  return (n < 1) ? 1 : n;
}

static int neko_tune_iters(void) {
  const char *v = getenv("NEKO_TUNE_ITERS");
  int n = (v != NULL) ? atoi(v) : 100;

  return (n < 1) ? 1 : n;
}

static int neko_tune_verify_on(void) {
  const char *v = getenv("NEKO_TUNE_VERIFY");

  return (v != NULL) ? (atoi(v) != 0) : 1;
}

/* Pipeline of candidate c at order lx, created on first use */
static id<MTLComputePipelineState> inviscid_div_pipeline(const int lx,
                                                        const int c) {
  if (pso_inviscid_div[lx][c] == nil) {
    char name[96];
    snprintf(name, sizeof(name),
             "compressible_res_inviscid_div_kernel_lx%d_ch%d_nc%d_r%d",
             lx, inviscid_div_chunks(lx, c), inviscid_div_stage(lx, c),
             inviscid_div_reg(lx, c));
    pso_inviscid_div[lx][c] =
      get_compressible_res_pipeline(name, inviscid_div_chunks(lx, c));
  }
  return pso_inviscid_div[lx][c];
}

/* Encode iters dispatches of candidate c, one threadgroup per element */
static void inviscid_div_encode(id<MTLComputeCommandEncoder> enc,
                                const int lx, const int c, void **bufs,
                                const int nel, const int iters) {
  [enc setComputePipelineState:inviscid_div_pipeline(lx, c)];
  for (int i = 0; i < 27; i++) {
    [enc setBuffer:(__bridge id<MTLBuffer>)bufs[i] offset:0 atIndex:i];
  }

  MTLSize groupSize = MTLSizeMake(inviscid_div_chunks(lx, c), 1, 1);
  MTLSize numGroups = MTLSizeMake((NSUInteger) nel, 1, 1);

  for (int i = 0; i < iters; i++) {
    [enc dispatchThreadgroups:numGroups threadsPerThreadgroup:groupSize];
  }
}

/* Run iters dispatches of candidate c and return the GPU time per
   dispatch in seconds, HUGE_VAL if the command buffer failed */
static double inviscid_div_run(const int lx, const int c, void **bufs,
                               const int nel, const int iters) {
  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = compressible_res_encoder(&cmdBuf);

  inviscid_div_encode(enc, lx, c, bufs, nel, iters);
  [enc endEncoding];
  [cmdBuf commit];
  [cmdBuf waitUntilCompleted];
  if (cmdBuf.error != nil) {
    return HUGE_VAL;
  }
  return (cmdBuf.GPUEndTime - cmdBuf.GPUStartTime) / (double) iters;
}

/*
 * Time every distinct candidate at this order and return the fastest:
 * interleaved rounds, min reduced, after a warm-up of every candidate, as
 * the CUDA and HIP operator tuners sample. Before anything is timed every
 * candidate is run on a synthetic state and compared with candidate 0: the
 * variants compute the same sums in the same order, so a difference beyond
 * single precision rounding is a kernel the compiler got wrong, which is
 * reported and kept from winning rather than left to the solve. (Without
 * the unrolling the kernels now carry, the Metal compiler on AMD GPUs was
 * seen to lose the accumulators of all but the first point of a thread.)
 * NEKO_INVISCID_DIV pins a candidate, NEKO_TUNE_ROUNDS and NEKO_TUNE_ITERS
 * set the sampling, NEKO_TUNE_VERIFY=0 skips the verification.
 *
 * The tuner leaves the output of the chosen candidate in place: it stands
 * in for the evaluation it was called for.
 */
static int tune_inviscid_div(const int lx, void **bufs, const int nel) {
  char neko_log_buf[80];
  char desc[48];
  double time[INVISCID_DIV_CANDIDATES];
  int rej[INVISCID_DIV_CANDIDATES] = { 0 };
  const int rounds = neko_tune_rounds();
  const int iters = neko_tune_iters();
  const int pin = neko_inviscid_div_pin();
  const int n = nel * lx * lx * lx;
  int best = 0;

  sprintf(neko_log_buf, "Autotune inviscid div (lx: %d)", lx);
  log_section(neko_log_buf);

  if (pin >= 0) {
    inviscid_div_run(lx, pin, bufs, nel, 1);
    sprintf(neko_log_buf, "Set by env   : %d (%s)", pin,
            inviscid_div_cand_desc(lx, pin, desc, sizeof(desc)));
    log_message(neko_log_buf);
    log_end_section();
    return pin;
  }

  for (int c = 0; c < INVISCID_DIV_CANDIDATES; c++) {
    time[c] = HUGE_VAL;
  }

  /* Verification on a synthetic state in (0.5, 1.5): the first call can be
     on a zero field, which every kernel reproduces. The buffers are shared
     with the host, so the inputs are filled and the outputs compared
     there. */
  if (neko_tune_verify_on()) {
    id<MTLDevice> device = neko_metal_device();
    id<MTLBuffer> synth[9];
    void *vbufs[27];
    float *ref = (float *) malloc(5 * (size_t) n * sizeof(float));
    const double tol = 1.0e-4;

    memcpy(vbufs, bufs, sizeof(vbufs));
    for (int i = 0; i < 9; i++) {
      synth[i] = [device newBufferWithLength:(NSUInteger) n * sizeof(float)
                                     options:MTLResourceStorageModeShared];
      float *x = (float *) synth[i].contents;
      for (int j = 0; j < n; j++) {
        unsigned int h = (unsigned int) j * 2654435761u + (i + 1) * 40503u;
        h ^= h >> 15;
        h *= 2246822519u;
        h ^= h >> 13;
        x[j] = 0.5f + (float) (h & 0xffffffu) / (float) 0x1000000u;
      }
      vbufs[5 + i] = (__bridge void *) synth[i];
    }

    sprintf(neko_log_buf, "Verify       : rel tol %.1e, synthetic input", tol);
    log_message(neko_log_buf);

    inviscid_div_run(lx, 0, vbufs, nel, 1);
    for (int q = 0; q < 5; q++) {
      memcpy(ref + (size_t) q * n, ((__bridge id<MTLBuffer>) bufs[q]).contents,
             (size_t) n * sizeof(float));
    }

    for (int c = 1; c < INVISCID_DIV_CANDIDATES; c++) {
      if (inviscid_div_cand_dup(lx, c)) {
        continue;
      }
      double maxdiff = 0.0, maxref = 0.0;
      if (inviscid_div_run(lx, c, vbufs, nel, 1) >= HUGE_VAL) {
        maxdiff = HUGE_VAL;
      }
      for (int q = 0; q < 5; q++) {
        const float *out = ((__bridge id<MTLBuffer>) bufs[q]).contents;
        for (int j = 0; j < n; j++) {
          const double y = ref[(size_t) q * n + j];
          const double d = fabs((double) out[j] - y);
          if (!(d < HUGE_VAL)) {
            maxdiff = HUGE_VAL;
          } else if (d > maxdiff) {
            maxdiff = d;
          }
          if (fabs(y) > maxref) {
            maxref = fabs(y);
          }
        }
      }
      rej[c] = !(maxdiff <= tol * maxref);
      if (rej[c]) {
        sprintf(neko_log_buf, "cand %2d REJECTED: |diff| %.1e, |ref| %.1e",
                c, maxdiff, maxref);
        log_message(neko_log_buf);
      }
    }
    free(ref);
  } else {
    sprintf(neko_log_buf, "Verify       : off (NEKO_TUNE_VERIFY=0)");
    log_message(neko_log_buf);
  }

  /* Warm every candidate before timing anything. Candidates that name the
     kernel of a lower one are skipped throughout */
  for (int c = 0; c < INVISCID_DIV_CANDIDATES; c++) {
    if (!inviscid_div_cand_dup(lx, c)) {
      inviscid_div_run(lx, c, bufs, nel, INVISCID_DIV_TUNE_WARMUP);
    }
  }

  /* Interleaved rounds, best time per candidate */
  for (int r = 0; r < rounds; r++) {
    for (int c = 0; c < INVISCID_DIV_CANDIDATES; c++) {
      if (!inviscid_div_cand_dup(lx, c)) {
        const double t = inviscid_div_run(lx, c, bufs, nel, iters);
        if (t < time[c]) {
          time[c] = t;
        }
      }
    }
  }

  for (int c = 0; c < INVISCID_DIV_CANDIDATES; c++) {
    if (time[c] >= HUGE_VAL) {
      continue;
    }
    sprintf(neko_log_buf, "cand %2d (%s): %9.2f us/call", c,
            inviscid_div_cand_desc(lx, c, desc, sizeof(desc)),
            time[c] * 1.0e6);
    log_message(neko_log_buf);
    if (rej[c]) {
      time[c] = HUGE_VAL;
    }
  }

  for (int c = 1; c < INVISCID_DIV_CANDIDATES; c++) {
    if (time[c] < time[best]) {
      best = c;
    }
  }

  inviscid_div_run(lx, best, bufs, nel, 1);

  sprintf(neko_log_buf, "Chose        : %d (%s)", best,
          inviscid_div_cand_desc(lx, best, desc, sizeof(desc)));
  log_message(neko_log_buf);
  log_end_section();
  return best;
}

/**
 * Fused inviscid divergence of the compressible residual on the Metal GPU:
 * one threadgroup per element, see compressible_res_inviscid_div_kernel.
 * The geometry is autotuned on the first call at each order.
 */
void compressible_res_inviscid_div_metal(void *rhs_rho, void *rhs_m_x,
                                         void *rhs_m_y, void *rhs_m_z,
                                         void *rhs_E, void *rho, void *m_x,
                                         void *m_y, void *m_z, void *E,
                                         void *p, void *u, void *v, void *w,
                                         void *drdx, void *dsdx, void *dtdx,
                                         void *drdy, void *dsdy, void *dtdy,
                                         void *drdz, void *dsdz, void *dtdz,
                                         void *dx, void *dy, void *dz,
                                         void *jacinv, int *nel, int *lx) {
  void *bufs[27] = { rhs_rho, rhs_m_x, rhs_m_y, rhs_m_z, rhs_E,
                     rho, m_x, m_y, m_z, E, p, u, v, w,
                     drdx, dsdx, dtdx, drdy, dsdy, dtdy, drdz, dsdz, dtdz,
                     dx, dy, dz, jacinv };

  if (*nel <= 0)
    return;

  if (*lx < 2 || *lx > 16) {
    fprintf(stderr, "%s: size not supported: %d\n", __FILE__, *lx);
    exit(1);
  }

  if (autotune_inviscid_div[*lx] == 0) {
    autotune_inviscid_div[*lx] = 1 + tune_inviscid_div(*lx, bufs, *nel);
  } else {
    inviscid_div_run(*lx, autotune_inviscid_div[*lx] - 1, bufs, *nel, 1);
  }
}

/**
 * Runge-Kutta stage summation on the Metal GPU.
 */
void compressible_res_part_rk_sum_metal(void *rho, void *m_x, void *m_y, void *m_z,
                                 void *E, void *k_rho_i, void *k_m_x_i,
                                 void *k_m_y_i, void *k_m_z_i, void *k_E_i,
                                 real *dt, real *c, int *n) {

  if (*n <= 0)
    return;

  if (pso_rk_sum == nil)
    pso_rk_sum =
      get_compressible_res_pipeline("compressible_res_part_rk_sum_kernel",
                                    NEKO_METAL_MAX_THREADS);

  id<MTLCommandBuffer> cmdBuf;
  id<MTLComputeCommandEncoder> enc = compressible_res_encoder(&cmdBuf);

  [enc setComputePipelineState:pso_rk_sum];

  [enc setBuffer:(__bridge id<MTLBuffer>)rho     offset:0 atIndex:0];
  [enc setBuffer:(__bridge id<MTLBuffer>)m_x     offset:0 atIndex:1];
  [enc setBuffer:(__bridge id<MTLBuffer>)m_y     offset:0 atIndex:2];
  [enc setBuffer:(__bridge id<MTLBuffer>)m_z     offset:0 atIndex:3];
  [enc setBuffer:(__bridge id<MTLBuffer>)E       offset:0 atIndex:4];
  [enc setBuffer:(__bridge id<MTLBuffer>)k_rho_i offset:0 atIndex:5];
  [enc setBuffer:(__bridge id<MTLBuffer>)k_m_x_i offset:0 atIndex:6];
  [enc setBuffer:(__bridge id<MTLBuffer>)k_m_y_i offset:0 atIndex:7];
  [enc setBuffer:(__bridge id<MTLBuffer>)k_m_z_i offset:0 atIndex:8];
  [enc setBuffer:(__bridge id<MTLBuffer>)k_E_i   offset:0 atIndex:9];
  [enc setBytes:dt length:sizeof(real) atIndex:10];
  [enc setBytes:c length:sizeof(real) atIndex:11];
  [enc setBytes:n length:sizeof(int) atIndex:12];

  compressible_res_dispatch(cmdBuf, enc, *n);
}
