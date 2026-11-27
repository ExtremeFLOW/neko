/*
 Copyright (c) 2025, The Neko Authors
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

#include <stdio.h>
#include <stdlib.h>
#include <device/device_config.h>
#include <device/cuda/check.h>
#include <math/bcknd/device/cuda/elem_block.h>
#include <math/bcknd/device/cuda/elem_block_tune.h>

#include "compressible_res_kernel.h"

extern "C" {
  #include <common/neko_log.h>
}

/* Candidate pinned by NEKO_INVISCID_DIV, or -1 to sweep, see neko_eb_pin() in
   elem_block_tune.h */
static int neko_inviscid_div_pin()
{
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

template <const int >
int tune_inviscid_div(void *rhs_rho, void *rhs_m_x, void *rhs_m_y,
                      void *rhs_m_z, void *rhs_E, void *rho, void *m_x,
                      void *m_y, void *m_z, void *E, void *p, void *u,
                      void *v, void *w, void *drdx, void *dsdx, void *dtdx,
                      void *drdy, void *dsdy, void *dtdy, void *drdz,
                      void *dsdz, void *dtdz, void *dx, void *dy, void *dz,
                      void *jacinv, int *nel, int *lx);

extern "C" {

  void compressible_res_part_visc_cuda(void *rhs_rho, void *rhs_m_x,
                                void *rhs_m_y, void *rhs_m_z, void *rhs_E,
                                void *visc_rho, void *visc_m_x, void *visc_m_y,
                                void *visc_m_z, void *visc_E,
                                void *mult, void *Binv, void *h1, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

    compressible_res_part_visc_kernel<real>
      <<<nblcks, nthrds, 0, stream>>>((real *) rhs_rho, (real *) rhs_m_x,
                                      (real *) rhs_m_y, (real *) rhs_m_z,
                                      (real *) rhs_E, (real *) visc_rho,
                                      (real *) visc_m_x, (real *) visc_m_y,
                                      (real *) visc_m_z, (real *) visc_E,
                                      (real *) mult, (real *) Binv,
                                      (real *) h1, *n);
    CUDA_CHECK(cudaGetLastError());
  }

  /**
   * Fused inviscid divergence of the compressible residual: one block per
   * element, see compressible_res_inviscid_div_kernel. The block geometry
   * is autotuned on the first call at each order.
   */
  void compressible_res_inviscid_div_cuda(void *rhs_rho, void *rhs_m_x,
                                          void *rhs_m_y, void *rhs_m_z,
                                          void *rhs_E, void *rho, void *m_x,
                                          void *m_y, void *m_z, void *E,
                                          void *p, void *u, void *v, void *w,
                                          void *drdx, void *dsdx, void *dtdx,
                                          void *drdy, void *dsdy, void *dtdy,
                                          void *drdz, void *dsdz, void *dtdz,
                                          void *dx, void *dy, void *dz,
                                          void *jacinv, int *nel, int *lx) {

    /* Candidate chosen by the tuner plus one, zero until it has run */
    static int autotune[17] = { 0 };

    const dim3 nblcks((*nel), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define CASE_CAND(LX, C)                                                      \
    compressible_res_inviscid_div_kernel<real, LX,                            \
                                         inviscid_div_chunks(LX, C),          \
                                         inviscid_div_stage(LX, C,            \
                                                            sizeof(real)),    \
                                         inviscid_div_regs(LX, C),            \
                                         inviscid_div_regm(LX, C)>            \
      <<<nblcks, inviscid_div_chunks(LX, C), 0, stream>>>                     \
        ((real *) rhs_rho, (real *) rhs_m_x, (real *) rhs_m_y,                \
         (real *) rhs_m_z, (real *) rhs_E, (real *) rho, (real *) m_x,        \
         (real *) m_y, (real *) m_z, (real *) E, (real *) p, (real *) u,      \
         (real *) v, (real *) w, (real *) drdx, (real *) dsdx,                \
         (real *) dtdx, (real *) drdy, (real *) dsdy, (real *) dtdy,          \
         (real *) drdz, (real *) dsdz, (real *) dtdz, (real *) dx,            \
         (real *) dy, (real *) dz, (real *) jacinv);                          \
    CUDA_CHECK(cudaGetLastError());

/* Runtime dispatch onto the tuned candidate */
#define CASE_CAND_SEL(LX, SEL)                                                \
    switch (SEL) {                                                            \
    case 0:  CASE_CAND(LX, 0); break;                                         \
    case 1:  CASE_CAND(LX, 1); break;                                         \
    case 2:  CASE_CAND(LX, 2); break;                                         \
    case 3:  CASE_CAND(LX, 3); break;                                         \
    case 4:  CASE_CAND(LX, 4); break;                                         \
    case 5:  CASE_CAND(LX, 5); break;                                         \
    case 6:  CASE_CAND(LX, 6); break;                                         \
    case 7:  CASE_CAND(LX, 7); break;                                         \
    case 8:  CASE_CAND(LX, 8); break;                                         \
    case 9:  CASE_CAND(LX, 9); break;                                         \
    case 10: CASE_CAND(LX, 10); break;                                        \
    case 11: CASE_CAND(LX, 11); break;                                        \
    case 12: CASE_CAND(LX, 12); break;                                        \
    case 13: CASE_CAND(LX, 13); break;                                        \
    case 14: CASE_CAND(LX, 14); break;                                        \
    default: CASE_CAND(LX, 15); break;                                        \
    }

/* The switch above has to name every candidate: one reaching the default
   arm would launch another geometry's kernel */
static_assert(INVISCID_DIV_CANDIDATES == 16,
              "CASE_CAND_SEL does not cover the candidate space");

#define CASE(LX)                                                              \
    case LX:                                                                  \
      if (autotune[LX] == 0) {                                                \
        autotune[LX] = 1 + tune_inviscid_div<LX>(rhs_rho, rhs_m_x, rhs_m_y,   \
                                                 rhs_m_z, rhs_E, rho, m_x,    \
                                                 m_y, m_z, E, p, u, v, w,     \
                                                 drdx, dsdx, dtdx, drdy,      \
                                                 dsdy, dtdy, drdz, dsdz,      \
                                                 dtdz, dx, dy, dz, jacinv,    \
                                                 nel, lx);                    \
      } else {                                                                \
        CASE_CAND_SEL(LX, autotune[LX] - 1);                                  \
      }                                                                       \
      break

    switch (*lx) {
      CASE(2);
      CASE(3);
      CASE(4);
      CASE(5);
      CASE(6);
      CASE(7);
      CASE(8);
      CASE(9);
      CASE(10);
      CASE(11);
      CASE(12);
      CASE(13);
      CASE(14);
      CASE(15);
      CASE(16);
    default:
      {
        fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
        exit(1);
      }
    }
#undef CASE
  }

  void compressible_res_part_rk_sum_cuda(void *rho, void *m_x, void *m_y,
                                  void *m_z, void *E,
                                  void *k_rho_i, void *k_m_x_i, void *k_m_y_i,
                                  void *k_m_z_i, void *k_E_i,
                                  real *dt, real *c, int *n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

    compressible_res_part_rk_sum_kernel<real>
      <<<nblcks, nthrds, 0, stream>>>((real *) rho, (real *) m_x, (real *) m_y,
                      (real *) m_z, (real *) E,
                      (real *) k_rho_i, (real *) k_m_x_i, (real *) k_m_y_i,
                      (real *) k_m_z_i, (real *) k_E_i,
                      *dt, *c, *n);
    CUDA_CHECK(cudaGetLastError());
  }
}

/*
 * Time every distinct candidate of the fused divergence at this order and
 * return the fastest, see the candidate description in
 * compressible_res_kernel.h. The sampling follows elem_block_tune.h:
 * interleaved rounds, min reduced, after a warm-up of every candidate.
 * NEKO_INVISCID_DIV pins a candidate, NEKO_TUNE_ROUNDS and NEKO_TUNE_ITERS
 * set the sampling.
 *
 * Every candidate computes the same sums in the same order, so the tuner
 * leaves the output of the chosen one in place merely by convention: it
 * stands in for the evaluation it was called for.
 */
template < const int LX >
int tune_inviscid_div(void *rhs_rho, void *rhs_m_x, void *rhs_m_y,
                      void *rhs_m_z, void *rhs_E, void *rho, void *m_x,
                      void *m_y, void *m_z, void *E, void *p, void *u,
                      void *v, void *w, void *drdx, void *dsdx, void *dtdx,
                      void *drdy, void *dsdy, void *dtdy, void *drdz,
                      void *dsdz, void *dtdz, void *dx, void *dy, void *dz,
                      void *jacinv, int *nel, int *lx) {
  cudaEvent_t start, stop;
  float time[INVISCID_DIV_CANDIDATES];
  const int rounds = neko_tune_rounds();
  const int iters = neko_tune_iters();
  const int pin = neko_inviscid_div_pin();
  const int tsize = (int) sizeof(real);
  int best = 0;

  const dim3 nblcks((*nel), 1, 1);
  const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

  char neko_log_buf[80];
  char desc[48];

  sprintf(neko_log_buf, "Autotune inviscid div (lx: %d)", *lx);
  log_section(neko_log_buf);

  if (pin >= 0) {
    CASE_CAND_SEL(LX, pin);
    sprintf(neko_log_buf, "Set by env   : %d (%s)", pin,
            inviscid_div_cand_desc(LX, pin, tsize, desc, sizeof(desc)));
    log_message(neko_log_buf);
    log_end_section();
    return pin;
  }

  for (int c = 0; c < INVISCID_DIV_CANDIDATES; c++) {
    time[c] = NEKO_TUNE_INIT;
  }

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  /* Warm every candidate before timing anything, see elem_block_tune.h.
     Candidates that instantiate the kernel of a lower one are skipped
     throughout, their time stays at NEKO_TUNE_INIT */
  for (int i = 0; i < NEKO_TUNE_WARMUP; i++) {
    for (int c = 0; c < INVISCID_DIV_CANDIDATES; c++) {
      if (!inviscid_div_cand_dup(LX, c, tsize)) {
        CASE_CAND_SEL(LX, c);
      }
    }
  }

  /* Interleaved rounds, best time per candidate */
  for (int r = 0; r < rounds; r++) {
    for (int c = 0; c < INVISCID_DIV_CANDIDATES; c++) {
      if (!inviscid_div_cand_dup(LX, c, tsize)) {
        NEKO_TUNE_TIME(time, CASE_CAND_SEL, LX, c, iters);
      }
    }
  }

  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  for (int c = 0; c < INVISCID_DIV_CANDIDATES; c++) {
    if (time[c] >= NEKO_TUNE_INIT) {
      continue;
    }
    sprintf(neko_log_buf, "cand %2d (%s): %9.2f us/call", c,
            inviscid_div_cand_desc(LX, c, tsize, desc, sizeof(desc)),
            NEKO_TUNE_US(time[c], iters));
    log_message(neko_log_buf);
  }

  NEKO_TUNE_BEST(time, best, INVISCID_DIV_CANDIDATES);

  CASE_CAND_SEL(LX, best);

  sprintf(neko_log_buf, "Chose        : %d (%s)", best,
          inviscid_div_cand_desc(LX, best, tsize, desc, sizeof(desc)));
  log_message(neko_log_buf);
  log_end_section();
  return best;
}
