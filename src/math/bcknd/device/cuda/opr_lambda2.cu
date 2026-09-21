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
#include "lambda2_kernel.h"
#include "elem_block_tune.h"
#include <device/device_config.h>
#include <device/cuda/check.h>

extern "C" {
  #include <common/neko_log.h>
}

template < const int >
int tune_lambda2(void *ux, void *uy, void *uz, void *u,
                void *dx, void *dy, void *dz,
                void *drdx, void *dsdx, void *dtdx,
                void *drdy, void *dsdy, void *dtdy,
                void *drdz, void *dsdz, void *dtdz,
                void *jacinv, int *nel, int *lx, int *eb_sel, int *ch_sel);

extern "C" {

  /** 
   * Fortran wrapper for device cuda convective terms
   */
  void cuda_lambda2(void *lambda2, void *u, void *v, void *w,
                   void *dx, void *dy, void *dz,
                   void *drdx, void *dsdx, void *dtdx,
                   void *drdy, void *dsdy, void *dtdy,
                   void *drdz, void *dsdz, void *dtdz,
                   void *jacinv, int *nel, int *lx) {

    static int autotune[17] = { 0 };
    /* elements per block candidate chosen by the tuner */
    static int autotune_eb[17] = { 0 };
    /* chunk candidate chosen for the 1d variant */
    static int autotune_ch[17] = { 0 };
    
    const dim3 nthrds_1d(1024, 1, 1);
    const dim3 nthrds_kstep((*lx), (*lx), 1);
    const dim3 nblcks((*nel), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;      

#define CASE_1D(LX, C)                                                          \
    lambda2_kernel_1d<real, LX, NEKO_CHUNKS(LX, C)>                             \
      <<<nblcks, NEKO_CHUNKS_NTHRDS(LX, C), 0, stream>>>                        \
      ((real *) lambda2, (real *) u, (real *) v, (real *) w,                    \
       (real *) dx, (real *) dy, (real *) dz,                                   \
       (real *) drdx, (real *) dsdx, (real *) dtdx,                             \
       (real *) drdy, (real *) dsdy, (real *) dtdy,                             \
       (real *) drdz, (real *) dsdz, (real *) dtdz,                             \
       (real *) jacinv);                                                        \
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
    lambda2_kernel_kstep<real, LX, NEKO_EB(LX, C)>                              \
    <<<NEKO_EB_NBLCKS(*nel, LX, C), NEKO_EB_NTHRDS(LX, C), 0, stream>>>         \
      ((real *) lambda2, (real *) u, (real *) v, (real *) w,                    \
       (real *) dx, (real *) dy, (real *) dz,                                   \
       (real *) drdx, (real *) dsdx, (real *) dtdx,                             \
       (real *) drdy, (real *) dsdy, (real *) dtdy,                             \
       (real *) drdz, (real *) dsdz, (real *) dtdz,                             \
       (real *) jacinv, *nel);                                                  \
    CUDA_CHECK(cudaGetLastError());                                             

/* Runtime dispatch onto the tuned candidate */
#define CASE_KSTEP_SEL(LX, SEL)                                                 \
    switch (SEL) {                                                              \
    case 0:  CASE_KSTEP(LX, 0); break;                                          \
    case 1:  CASE_KSTEP(LX, 1); break;                                          \
    default: CASE_KSTEP(LX, 2); break;                                          \
    }

#define CASE(LX)                                                                \
    case LX:                                                                    \
      if(autotune[LX] == 0 ) {                                                  \
        autotune[LX]=tune_lambda2<LX>(lambda2, u, v, w,                         \
                                     dx, dy, dz,                                \
                                     drdx, dsdx, dtdx,                          \
                                     drdy, dsdy, dtdy,                          \
                                     drdz, dsdz, dtdz,                          \
                                     jacinv, nel, lx, &autotune_eb[LX],         \
                                     &autotune_ch[LX]);                         \
      } else if (autotune[LX] == 1 ) {                                          \
        CASE_1D_SEL(LX, autotune_ch[LX]);                                       \
      } else if (autotune[LX] == 2 ) {                                          \
        CASE_KSTEP_SEL(LX, autotune_eb[LX]);                                    \
      }                                                                         \
      break
    
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
      CASE(11);
      CASE(12);
    default:
      {
        fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
        exit(1);
      }
    }
  } 
}

template < const int LX >
int tune_lambda2(void *lambda2, void *u, void *v, void *w,
                void *dx, void *dy, void *dz,
                void *drdx, void *dsdx, void *dtdx,
                void *drdy, void *dsdy, void *dtdy,
                void *drdz, void *dsdz, void *dtdz,
                void *jacinv, int *nel, int *lx, int *eb_sel, int *ch_sel) {
  cudaEvent_t start,stop;
  float time1[NEKO_CHUNKS_CANDIDATES];
  int best1 = 0;
  float time2[NEKO_EB_CANDIDATES];
  const int rounds = neko_tune_rounds();
  const int iters = neko_tune_iters();
  /* Candidates of the kstep sweep, one -- the unblocked shape -- when the
     elements per block sweep is off */
  const int eb_cand = neko_eb_sweep() ? NEKO_EB_CANDIDATES : 1;
  /* Geometry pinned by each formulation's own variable, -1 to sweep it */
  const int ch_pin = neko_chunks_pin();
  const int eb_pin = neko_eb_pin();
  /* Formulation pinned by NEKO_AUTOTUNE, as the identifier this returns */
  int strat = 0;
  int best = 0;
  int retval;

  for (int c = 0; c < NEKO_EB_CANDIDATES; c++) {
    time2[c] = NEKO_TUNE_INIT;
  }
  for (int c = 0; c < NEKO_CHUNKS_CANDIDATES; c++) {
    time1[c] = NEKO_TUNE_INIT;
  }

  const dim3 nthrds_1d(1024, 1, 1);
  const dim3 nthrds_kstep((*lx), (*lx), 1);
  const dim3 nblcks((*nel), 1, 1);
  const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;
  
  char *env_value = NULL;
  char neko_log_buf[80];
  
  env_value=getenv("NEKO_AUTOTUNE");

  sprintf(neko_log_buf, "Autotune lambda2 (lx: %d)", *lx);
  log_section(neko_log_buf);
  
  /*
   * NEKO_AUTOTUNE names a formulation, and that is all it does: the sweep
   * below is narrowed to that one kernel family, but its geometry -- the
   * chunk size, the elements per block -- is still measured candidate
   * against candidate, and reported.
   */
  if(env_value) {
    if( !strcmp(env_value,"1D") ) {
      strat = 1;
    } else if( !strcmp(env_value,"KSTEP") ) {
      strat = 2;
    } else {       
       sprintf(neko_log_buf, "Invalid value set for NEKO_AUTOTUNE");
       log_error(neko_log_buf);
    }
  }

  /* Geometry of the pinned formulation, if its own variable fixes that too.
     Both pinned leaves nothing to measure, so the kernel is launched once and
     reported, which is what pinning has always done */
  const int pin = (strat == 1) ? ch_pin : (strat == 2) ? eb_pin : -1;

  if (pin >= 0) {
    switch (strat) {
    case 1:
      *ch_sel = pin;
      CASE_1D_SEL(LX, pin);
      sprintf(neko_log_buf, "Set by env   : 1 (1D, %d chunk)",
              NEKO_CHUNKS_SEL(LX, pin));
      break;
    default:
      *eb_sel = pin;
      CASE_KSTEP_SEL(LX, pin);
      sprintf(neko_log_buf, "Set by env   : 2 (KSTEP, %d elem/block)",
              NEKO_EB_SEL(LX, pin));
      break;
    }
    log_message(neko_log_buf);
    log_end_section();
    return strat;
  }

  if (strat) {
    sprintf(neko_log_buf, "Set by env   : %d (%s)", strat, env_value);
    log_message(neko_log_buf);
  }

  /* Formulations the sweep considers, see NEKO_TUNE_FOR() */
  const bool try_1d = (strat == 0 || strat == 1);
  const bool try_kstep = (strat == 0 || strat == 2);

  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  /* Warm every variant before timing anything: each specialisation has to be
     resident and the clocks at steady state, or whichever is timed first is
     measured on a colder part */
  for (int i = 0; i < NEKO_TUNE_WARMUP; i++) {
    NEKO_TUNE_FOR(c, try_1d, ch_pin, NEKO_CHUNKS_CANDIDATES) {
      CASE_1D_SEL(LX, c);
    }
    NEKO_TUNE_FOR(c, try_kstep, eb_pin, eb_cand) {
      CASE_KSTEP_SEL(LX, c);
    }
  }

  /* Interleaved rounds, best time per variant */
  for (int r = 0; r < rounds; r++) {
    NEKO_TUNE_FOR(c, try_1d, ch_pin, NEKO_CHUNKS_CANDIDATES) {
      NEKO_TUNE_TIME(time1, CASE_1D_SEL, LX, c, iters);
    }
    NEKO_TUNE_FOR(c, try_kstep, eb_pin, eb_cand) {
      NEKO_TUNE_TIME(time2, CASE_KSTEP_SEL, LX, c, iters);
    }
  }

  NEKO_TUNE_LOG(LX, time1, time2);
  NEKO_TUNE_BEST(time1, best1, NEKO_CHUNKS_CANDIDATES);
  NEKO_TUNE_BEST(time2, best, NEKO_EB_CANDIDATES);
  *eb_sel = best;
  *ch_sel = best1;

  if (time1[best1] < time2[best]) {
    retval = 1;
  } else {
    retval = 2;
  }

  /* Leave the chosen kernel's output in place: the tuner stands in for a real
     evaluation and the variants do not sum in the same order */
  if (retval == 1) {
    CASE_1D_SEL(LX, best1);
  } else {
    CASE_KSTEP_SEL(LX, best);
  }

  if (retval == 1) {
    sprintf(neko_log_buf, "Chose        : 1 (1D, %d chunk)",
            NEKO_CHUNKS_SEL(LX, best1));
  } else {
    sprintf(neko_log_buf, "Chose        : 2 (KSTEP, %d elem/block)",
            NEKO_EB_SEL(LX, best));
  }
  log_message(neko_log_buf);
  log_end_section();
  return retval;
}

