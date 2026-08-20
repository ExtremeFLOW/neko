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
#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "opgrad_kernel.h"
#include "elem_block_tune.h"

extern "C" {
  #include <common/neko_log.h>
}

template < const int >
int tune_opgrad(void *ux, void *uy, void *uz, void *u,
                void *dx, void *dy, void *dz,
                void *drdx, void *dsdx, void *dtdx,
                void *drdy, void *dsdy, void *dtdy,
                void *drdz, void *dsdz, void *dtdz,
                void *w3, int *nel, int *lx, int *eb_sel, int *ch_sel);

extern "C" {

  /**
   * Fortran wrapper for device hip convective terms
   */
  void hip_opgrad(void *ux, void *uy, void *uz, void *u,
                   void *dx, void *dy, void *dz,
                   void *drdx, void *dsdx, void *dtdx,
                   void *drdy, void *dsdy, void *dtdy,
                   void *drdz, void *dsdz, void *dtdz,
                   void *w3, int *nel, int *lx) {

    static int autotune[19] = { 0 };
    /* elements per block candidate chosen by the tuner */
    static int autotune_eb[19] = { 0 };
    /* chunk candidate chosen for the 1d variant */
    static int autotune_ch[19] = { 0 };

    const dim3 nthrds_1d(1024, 1, 1);
    const dim3 nthrds_kstep((*lx), (*lx), 1);
    const dim3 nblcks((*nel), 1, 1);

#define CASE_1D(LX, C)                                                          \
    hipLaunchKernelGGL( HIP_KERNEL_NAME(                                        \
                        opgrad_kernel_1d<real, LX, NEKO_CHUNKS(LX, C)> ),       \
                        nblcks, NEKO_CHUNKS_NTHRDS(LX, C), 0,                   \
                        (hipStream_t) glb_cmd_queue,                            \
                        (real *) ux, (real *) uy, (real *) uz, (real *) u,      \
                        (real *) dx, (real *) dy, (real *) dz,                  \
                        (real *) drdx, (real *) dsdx, (real *) dtdx,            \
                        (real *) drdy, (real *) dsdy, (real *) dtdy,            \
                        (real *) drdz, (real *) dsdz, (real *) dtdz,            \
                        (real *) w3);                                           \
    HIP_CHECK(hipGetLastError());

/* Runtime dispatch onto the tuned chunk candidate */
#define CASE_1D_SEL(LX, SEL)                                                    \
    switch (SEL) {                                                              \
    case 0:  CASE_1D(LX, 0); break;                                             \
    case 1:  CASE_1D(LX, 1); break;                                             \
    case 2:  CASE_1D(LX, 2); break;                                             \
    default: CASE_1D(LX, 3); break;                                             \
    }


#define CASE_KSTEP(LX, C)                                                       \
    hipLaunchKernelGGL( HIP_KERNEL_NAME(                                        \
                        opgrad_kernel_kstep<real, LX, NEKO_EB(LX, C)> ),        \
                        NEKO_EB_NBLCKS(*nel, LX, C), NEKO_EB_NTHRDS(LX, C), 0,  \
                        (hipStream_t) glb_cmd_queue,                            \
                        (real *) ux, (real *) uy, (real *) uz, (real *) u,      \
                        (real *) dx, (real *) dy, (real *) dz,                  \
                        (real *) drdx, (real *) dsdx, (real *) dtdx,            \
                        (real *) drdy, (real *) dsdy, (real *) dtdy,            \
                        (real *) drdz, (real *) dsdz, (real *) dtdz,            \
                        (real *) w3, *nel);                                     \
    HIP_CHECK(hipGetLastError());

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
        autotune[LX]=tune_opgrad<LX>(ux, uy, uz, u,                             \
                                     dx, dy, dz,                                \
                                     drdx, dsdx, dtdx,                          \
                                     drdy, dsdy, dtdy,                          \
                                     drdz, dsdz, dtdz,                          \
                                     w3, nel, lx, &autotune_eb[LX],             \
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
      CASE(13);
      CASE(14);
      CASE(15);
      CASE(16);
      CASE(17);
      CASE(18);
    default:
      {
        fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
        exit(1);
      }
    }
  }
}

template < const int LX >
int tune_opgrad(void *ux, void *uy, void *uz, void *u,
                void *dx, void *dy, void *dz,
                void *drdx, void *dsdx, void *dtdx,
                void *drdy, void *dsdy, void *dtdy,
                void *drdz, void *dsdz, void *dtdz,
                void *w3, int *nel, int *lx, int *eb_sel, int *ch_sel) {
  hipEvent_t start,stop;
  float time1[NEKO_CHUNKS_CANDIDATES];
  int best1 = 0;
  float time2[NEKO_EB_CANDIDATES];
  const int rounds = neko_tune_rounds();
  const int iters = neko_tune_iters();
  const int sweep = neko_eb_sweep();
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
  const hipStream_t stream = (hipStream_t) glb_cmd_queue;

  char *env_value = NULL;
  char neko_log_buf[80];

  env_value=getenv("NEKO_AUTOTUNE");

  sprintf(neko_log_buf, "Autotune opgrad (lx: %d)", *lx);
  log_section(neko_log_buf);

  *eb_sel = 0;
  *ch_sel = 0;

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
      {
        const int c = neko_eb_env();
        *eb_sel = c;
        CASE_KSTEP_SEL(LX, c);
        sprintf(neko_log_buf,"Set by env   : 2 (KSTEP, %d elem/block)",
                NEKO_EB_SEL(LX, c));
      }
      log_message(neko_log_buf);
      log_end_section();
      return 2;
    } else {
       sprintf(neko_log_buf, "Invalid value set for NEKO_AUTOTUNE");
       log_error(neko_log_buf);
    }
  }

  HIP_CHECK(hipEventCreate(&start));
  HIP_CHECK(hipEventCreate(&stop));

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
