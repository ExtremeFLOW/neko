/*
 Copyright (c) 2021-2026, The Neko Authors
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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "dudxyz_kernel.h"
#include "elem_block_tune.h"
#include <device/device_config.h>
#include <device/cuda/check.h>

extern "C" {
  #include <common/neko_log.h>
}

template <const int >
int tune_dudxyz(void *du, void *u,
                void *dr, void *ds, void *dt,
                void *dx, void *dy, void *dz,
                void *jacinv, int *nel, int *lx, int *eb_sel, int *ch_sel,
                int *nw_sel, int *tw_sel);

extern "C" {

  /**
   * Fortran wrapper for device cuda derivative kernels
   */
  void cuda_dudxyz(void *du, void *u,
                  void *dr, void *ds, void *dt,
                  void *dx, void *dy, void *dz,
                  void *jacinv, int *nel, int *lx) {

    static int autotune[16] = { 0 };
    /* elements per block candidate chosen by the tuner */
    static int autotune_eb[16] = { 0 };
    /* chunk candidate chosen for the 1d variant */
    static int autotune_ch[16] = { 0 };
    /* warps per block candidate chosen for the dmma variant */
    static int autotune_nw[16] = { 0 };
    /* warps per block candidate chosen for the tma staged dmma variant */
    static int autotune_tw[16] = { 0 };

    const dim3 nthrds_1d(1024, 1, 1);
    const dim3 nthrds_kstep((*lx), (*lx), 1);
    const dim3 nblcks((*nel), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define CASE_1D(LX, C)                                                          \
    dudxyz_kernel_1d<real, LX, NEKO_CHUNKS(LX, C)>                              \
      <<<nblcks, NEKO_CHUNKS_NTHRDS(LX, C), 0, stream>>>                        \
                              ((real *) du, (real *) u,                         \
                              (real *) dr, (real *) ds, (real *) dt,            \
                              (real *) dx, (real *) dy, (real *) dz,            \
                              (real *) jacinv);                                 \
    CUDA_CHECK(cudaGetLastError());

/* Runtime dispatch onto the tuned chunk candidate */
#define CASE_1D_SEL(LX, SEL)                                                    \
    switch (SEL) {                                                              \
    case 0:  CASE_1D(LX, 0); break;                                             \
    case 1:  CASE_1D(LX, 1); break;                                             \
    case 2:  CASE_1D(LX, 2); break;                                             \
    default: CASE_1D(LX, 3); break;                                             \
    }

#define CASE_KSTEP(LX, C)                                                       \
    dudxyz_kernel_kstep<real, LX, NEKO_EB(LX, C)>                               \
      <<<NEKO_EB_NBLCKS(*nel, LX, C), NEKO_EB_NTHRDS(LX, C), 0, stream>>>       \
                                ((real *) du, (real *) u,                       \
                                (real *) dr, (real *) ds, (real *) dt,          \
                                (real *) dx, (real *) dy, (real *) dz,          \
                                (real *) jacinv, *nel);                         \
      CUDA_CHECK(cudaGetLastError());

/* Runtime dispatch onto the tuned candidate */
#define CASE_KSTEP_SEL(LX, SEL)                                                 \
    switch (SEL) {                                                              \
    case 0:  CASE_KSTEP(LX, 0); break;                                          \
    case 1:  CASE_KSTEP(LX, 1); break;                                          \
    default: CASE_KSTEP(LX, 2); break;                                          \
    }

#define CASE_DMMA(LX, C)                                                        \
    dudxyz_kernel_dmma<real, LX, NEKO_DMMA_NW(C)>                               \
      <<<NEKO_DMMA_NBLCKS(*nel, LX), NEKO_DMMA_NTHRDS(C), 0, stream>>>          \
                                ((real *) du, (real *) u,                       \
                                 (real *) dr, (real *) ds, (real *) dt,         \
                                 (real *) dx, (real *) dy, (real *) dz,         \
                                 (real *) jacinv, *nel);                        \
      CUDA_CHECK(cudaGetLastError());

/* Runtime dispatch onto the tuned warps per block candidate */
#define CASE_DMMA_SEL(LX, SEL)                                                  \
    switch (SEL) {                                                              \
    case 0:  CASE_DMMA(LX, 0); break;                                           \
    case 1:  CASE_DMMA(LX, 1); break;                                           \
    default: CASE_DMMA(LX, 2); break;                                           \
    }

/*
 * The TMA staged dmma variant. Same grid and the same warps per block
 * candidates as CASE_DMMA -- at the one lx it supports NEKO_DMMA_PACK is 1,
 * so the grid is nel blocks -- but it takes no nel: there is no packed tail
 * to clamp when a block is exactly an element. See dmma_tma_kernel.h.
 */
#define CASE_DMMA_TMA(LX, C)                                                    \
    dudxyz_kernel_dmma_tma<real, LX, NEKO_DMMA_NW(C)>                           \
      <<<NEKO_DMMA_NBLCKS(*nel, LX), NEKO_DMMA_NTHRDS(C), 0, stream>>>          \
                                ((real *) du, (real *) u,                       \
                                 (real *) dr, (real *) ds, (real *) dt,         \
                                 (real *) dx, (real *) dy, (real *) dz,         \
                                 (real *) jacinv);                              \
      CUDA_CHECK(cudaGetLastError());

/* Runtime dispatch onto the tuned warps per block candidate */
#define CASE_DMMA_TMA_SEL(LX, SEL)                                              \
    switch (SEL) {                                                              \
    case 0:  CASE_DMMA_TMA(LX, 0); break;                                       \
    case 1:  CASE_DMMA_TMA(LX, 1); break;                                       \
    default: CASE_DMMA_TMA(LX, 2); break;                                       \
    }

 #define CASE(LX)                                                               \
    case LX:                                                                    \
      if(autotune[LX] == 0 ) {                                                  \
        autotune[LX]=tune_dudxyz<LX>(du, u,                                     \
                                     dr, ds, dt,                                \
                                     dx, dy, dz,                                \
                                     jacinv, nel, lx, &autotune_eb[LX],         \
                                     &autotune_ch[LX], &autotune_nw[LX],        \
                                     &autotune_tw[LX]);                         \
      } else if (autotune[LX] == 1 ) {                                          \
        CASE_1D_SEL(LX, autotune_ch[LX]);                                       \
      } else if (autotune[LX] == 2 ) {                                          \
        CASE_KSTEP_SEL(LX, autotune_eb[LX]);                                    \
      } else if (autotune[LX] == 3 ) {                                          \
        CASE_DMMA_SEL(LX, autotune_nw[LX]);                                     \
      } else if (autotune[LX] == 4 ) {                                          \
        CASE_DMMA_TMA_SEL(LX, autotune_tw[LX]);                                 \
      }                                                                         \
      break

#define CASE_LARGE(LX)                                                          \
    case LX:                                                                    \
      CASE_KSTEP(LX, 0);                                                        \
      break


    if ((*lx) < 11) {
      switch(*lx) {
        CASE(2);
        CASE(3);
        CASE(4);
        CASE(5);
        CASE(6);
        CASE(7);
        CASE(8);
        CASE(9);
        CASE(10);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
    }
    else {
      switch(*lx) {
        CASE_LARGE(11);
        CASE_LARGE(12);
        CASE_LARGE(13);
        CASE_LARGE(14);
        CASE_LARGE(15);
        CASE_LARGE(16);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
    }
  }
}

template < const int LX >
int tune_dudxyz(void *du, void *u,
                void *dr, void *ds, void *dt,
                void *dx, void *dy, void *dz,
                void *jacinv, int *nel, int *lx, int *eb_sel, int *ch_sel,
                int *nw_sel, int *tw_sel) {
  cudaEvent_t start,stop;
  float time1[NEKO_CHUNKS_CANDIDATES];
  int best1 = 0;
  float time2[NEKO_EB_CANDIDATES];
  float time3[NEKO_DMMA_CANDIDATES];
  int best3 = 0;
  float time4[NEKO_DMMA_CANDIDATES];
  int best4 = 0;
  const int rounds = neko_tune_rounds();
  const int iters = neko_tune_iters();
  const int sweep = neko_eb_sweep();
  const bool dmma = dmma_lx_supported<LX>() && cuda_have_dmma();
  /* Whether the pointers are bulk copy aligned is a property of this call
     rather than of the kernel, so it is checked rather than assumed, see
     dmma_tma_dudxyz_aligned() in dmma_tma_kernel.h */
  const bool tma = dmma_tma_dudxyz_lx_supported<LX>() && cuda_have_tma() &&
                   dmma_tma_dudxyz_aligned(du, u, dr, ds, dt, jacinv);
  int best = 0;
  int retval;

  for (int c = 0; c < NEKO_EB_CANDIDATES; c++) {
    time2[c] = NEKO_TUNE_INIT;
  }
  for (int c = 0; c < NEKO_CHUNKS_CANDIDATES; c++) {
    time1[c] = NEKO_TUNE_INIT;
  }
  for (int c = 0; c < NEKO_DMMA_CANDIDATES; c++) {
    time3[c] = NEKO_TUNE_INIT;
    time4[c] = NEKO_TUNE_INIT;
  }

  const dim3 nthrds_1d(1024, 1, 1);
  const dim3 nthrds_kstep((*lx), (*lx), 1);
  const dim3 nblcks((*nel), 1, 1);
  const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

  char *env_value = NULL;
  char neko_log_buf[80];

  env_value=getenv("NEKO_AUTOTUNE");

  sprintf(neko_log_buf, "Autotune dudxyz (lx: %d)", *lx);
  log_section(neko_log_buf);

  *eb_sel = 0;
  *ch_sel = 0;
  *nw_sel = 0;
  *tw_sel = 0;

  if(env_value) {
    if( !strcmp(env_value,"1D") ) {
      *ch_sel = neko_chunks_env();
      CASE_1D_SEL(LX, *ch_sel);
      sprintf(neko_log_buf,"Set by env   : 1 (1D, %d chunk)",
              NEKO_CHUNKS_SEL(LX, *ch_sel));
      log_message(neko_log_buf);
      log_end_section();
      return 1;
    } else if( !strcmp(env_value,"KSTEP") ) {
      *eb_sel = neko_eb_env();
      CASE_KSTEP_SEL(LX, *eb_sel);
      sprintf(neko_log_buf,"Set by env   : 2 (KSTEP, %d elem/block)",
              NEKO_EB_SEL(LX, *eb_sel));
      log_message(neko_log_buf);
      log_end_section();
      return 2;
    } else if( !strcmp(env_value,"DMMA") ) {
      if (dmma) {
        const int c = neko_dmma_env();
        *nw_sel = c;
        CASE_DMMA_SEL(LX, c);
        sprintf(neko_log_buf,"Set by env   : 3 (DMMA, %d warps)",
                NEKO_DMMA_NW(c));
        log_message(neko_log_buf);
        log_end_section();
        return 3;
      } else {
        sprintf(neko_log_buf, "DMMA strategy not available for this config");
        log_error(neko_log_buf);
      }
    } else if( !strcmp(env_value,"DMMA_TMA") ) {
      if (tma) {
        const int c = neko_dmma_tma_env();
        *tw_sel = c;
        CASE_DMMA_TMA_SEL(LX, c);
        sprintf(neko_log_buf,"Set by env   : 4 (DMMA_TMA, %d warps)",
                NEKO_DMMA_NW(c));
        log_message(neko_log_buf);
        log_end_section();
        return 4;
      } else {
        sprintf(neko_log_buf,
                "DMMA_TMA strategy not available for this config");
        log_error(neko_log_buf);
      }
    } else {
       sprintf(neko_log_buf, "Invalid value set for NEKO_AUTOTUNE");
       log_error(neko_log_buf);
    }
  }

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  /* Warm every variant before timing anything: each specialisation has to be
     resident and the clocks at steady state, or whichever is timed first is
     measured on a colder part */
  for (int i = 0; i < NEKO_TUNE_WARMUP; i++) {
    CASE_1D(LX, 0);
    CASE_1D(LX, 1);
    CASE_1D(LX, 2);
    CASE_1D(LX, 3);
    CASE_KSTEP(LX, 0);
    if (sweep) {
      CASE_KSTEP(LX, 1);
      CASE_KSTEP(LX, 2);
    }
    if (dmma) {
      CASE_DMMA(LX, 0);
      CASE_DMMA(LX, 1);
      CASE_DMMA(LX, 2);
    }
    if (tma) {
      CASE_DMMA_TMA(LX, 0);
      CASE_DMMA_TMA(LX, 1);
      CASE_DMMA_TMA(LX, 2);
    }
  }

  /* Interleaved rounds, best time per variant */
  for (int r = 0; r < rounds; r++) {
    NEKO_TUNE_TIME(time1, CASE_1D, LX, 0, iters);
    NEKO_TUNE_TIME(time1, CASE_1D, LX, 1, iters);
    NEKO_TUNE_TIME(time1, CASE_1D, LX, 2, iters);
    NEKO_TUNE_TIME(time1, CASE_1D, LX, 3, iters);
    NEKO_TUNE_TIME(time2, CASE_KSTEP, LX, 0, iters);
    if (sweep) {
      NEKO_TUNE_TIME(time2, CASE_KSTEP, LX, 1, iters);
      NEKO_TUNE_TIME(time2, CASE_KSTEP, LX, 2, iters);
    }
    if (dmma) {
      NEKO_TUNE_TIME(time3, CASE_DMMA, LX, 0, iters);
      NEKO_TUNE_TIME(time3, CASE_DMMA, LX, 1, iters);
      NEKO_TUNE_TIME(time3, CASE_DMMA, LX, 2, iters);
    }
    if (tma) {
      NEKO_TUNE_TIME(time4, CASE_DMMA_TMA, LX, 0, iters);
      NEKO_TUNE_TIME(time4, CASE_DMMA_TMA, LX, 1, iters);
      NEKO_TUNE_TIME(time4, CASE_DMMA_TMA, LX, 2, iters);
    }
  }

  NEKO_TUNE_LOG(LX, time1, time2);
  NEKO_TUNE_LOG_DMMA(LX, time3);
  NEKO_TUNE_LOG_DMMA_TMA(LX, time4);

  NEKO_TUNE_BEST(time1, best1, NEKO_CHUNKS_CANDIDATES);
  NEKO_TUNE_BEST(time2, best, NEKO_EB_CANDIDATES);
  NEKO_TUNE_BEST(time3, best3, NEKO_DMMA_CANDIDATES);
  NEKO_TUNE_BEST(time4, best4, NEKO_DMMA_CANDIDATES);
  *eb_sel = best;
  *ch_sel = best1;
  *nw_sel = best3;
  *tw_sel = best4;

  if (time1[best1] < time2[best]) {
    retval = 1;
  } else {
    retval = 2;
  }

  /* The dmma variants join the comparison only where they exist, their
     candidates are left at NEKO_TUNE_INIT otherwise */
  float tbest = (retval == 1) ? time1[best1] : time2[best];

  if (time3[best3] < tbest) {
    retval = 3;
    tbest = time3[best3];
  }
  if (time4[best4] < tbest) {
    retval = 4;
  }

  /* Leave the chosen kernel's output in place: the tuner stands in for a real
     evaluation and the variants do not sum in the same order */
  if (retval == 1) {
    CASE_1D_SEL(LX, best1);
  } else if (retval == 2) {
    CASE_KSTEP_SEL(LX, best);
  } else if (retval == 3) {
    CASE_DMMA_SEL(LX, best3);
  } else {
    CASE_DMMA_TMA_SEL(LX, best4);
  }

  if (retval == 1) {
    sprintf(neko_log_buf, "Chose        : 1 (1D, %d chunk)",
            NEKO_CHUNKS_SEL(LX, best1));
  } else if (retval == 2) {
    sprintf(neko_log_buf, "Chose        : 2 (KSTEP, %d elem/block)",
            NEKO_EB_SEL(LX, best));
  } else if (retval == 3) {
    sprintf(neko_log_buf, "Chose        : 3 (DMMA, %d warps, %d elem/blk)",
            NEKO_DMMA_NW(best3), NEKO_DMMA_PACK(LX));
  } else {
    sprintf(neko_log_buf, "Chose        : 4 (DMMA_TMA, %d warps)",
            NEKO_DMMA_NW(best4));
  }
  log_message(neko_log_buf);
  log_end_section();
  return retval;
}