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
#include "ax_helm_kernel.h"
#include "elem_block_tune.h"

extern "C" {
  #include <common/neko_log.h>
}

template < const int>
int tune(void *w, void *u, void *dx, void *dy, void *dz,
         void *dxt, void *dyt, void *dzt, void *h1,
         void *g11, void *g22, void *g33, void *g12,
         void *g13, void *g23, int *nelv, int *lx, int *eb_sel, int *ch_sel);

template < const int>
int tune_padded(void *w, void *u, void *dx, void *dy, void *dz,
                void *dxt, void *dyt, void *dzt, void *h1,
                void *g11, void *g22, void *g33, void *g12,
                void *g13, void *g23, int *nelv, int *lx, int *eb_sel,
                int *ch_sel);

template < const int>
int tune_vector(void *au, void *av, void *aw, void *u, void *v, void *w,
                void *dx, void *dy, void *dz, void *h1,
                void *g11, void *g22, void *g33, void *g12,
                void *g13, void *g23, int *nelv, int *lx);

template < const int>
int tune_vector_padded(void *au, void *av, void *aw,
                       void *u, void *v, void *w,
                       void *dx, void *dy, void *dz, void *h1,
                       void *g11, void *g22, void *g33, void *g12,
                       void *g13, void *g23, int *nelv, int *lx);

extern "C" {

  /**
   * Fortran wrapper for device HIP Ax
   */
  void hip_ax_helm(void *w, void *u, void *dx, void *dy, void *dz,
                   void *dxt, void *dyt, void *dzt, void *h1,
                   void *g11, void *g22, void *g33, void *g12,
                   void *g13, void *g23, int *nelv, int *lx) {

    static int autotune[13] = { 0 };
    /* Elements per block candidate chosen by the autotuner, see
       elem_block<> in elem_block.h */
    static int autotune_eb[13] = { 0 };
    /* chunk candidate chosen for the 1d variant */
    static int autotune_ch[13] = { 0 };

    const dim3 nthrds_1d(1024, 1, 1);
    const dim3 nblcks_1d((*nelv), 1, 1);

#define CASE_1D(LX, C)                                                          \
    hipLaunchKernelGGL(HIP_KERNEL_NAME(                                         \
                       ax_helm_kernel_1d<real, LX, NEKO_CHUNKS(LX, C)> ),       \
                       nblcks_1d, NEKO_CHUNKS_NTHRDS(LX, C), 0,                 \
                       (hipStream_t) glb_cmd_queue,                             \
                       (real *) w, (real *) u,                                  \
                       (real *) dx, (real *) dy, (real *) dz,                   \
                       (real *) dxt, (real *) dyt, (real *) dzt, (real *) h1,   \
                       (real *) g11, (real *) g22, (real *) g33,                \
                       (real *) g12, (real *) g13, (real *) g23);               \
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
                          ax_helm_kernel_kstep<real, LX, NEKO_EB(LX, C)> ),     \
                        NEKO_EB_NBLCKS(*nelv, LX, C), NEKO_EB_NTHRDS(LX, C), 0, \
                        (hipStream_t) glb_cmd_queue,                            \
                        (real *) w, (real *) u,                                 \
                        (real *) dx, (real *) dy, (real *) dz,  (real *) h1,    \
                        (real *) g11, (real *) g22, (real *) g33,               \
                        (real *) g12, (real *) g13, (real *) g23, *nelv);       \
    HIP_CHECK(hipGetLastError());

#define CASE_KSTEP_PADDED(LX, C)                                                \
    hipLaunchKernelGGL( HIP_KERNEL_NAME(                                        \
                        ax_helm_kernel_kstep_padded<real, LX, NEKO_EB(LX, C)> ),\
                        NEKO_EB_NBLCKS(*nelv, LX, C), NEKO_EB_NTHRDS(LX, C), 0, \
                        (hipStream_t) glb_cmd_queue,                            \
                        (real *) w, (real *) u,                                 \
                        (real *) dx, (real *) dy, (real *) dz,  (real *) h1,    \
                        (real *) g11, (real *) g22, (real *) g33,               \
                        (real *) g12, (real *) g13, (real *) g23, *nelv);       \
    HIP_CHECK(hipGetLastError());

/* Runtime dispatch onto the tuned elements per block candidate */
#define CASE_KSTEP_SEL(LX, SEL)                                                 \
    switch (SEL) {                                                              \
    case 0:  CASE_KSTEP(LX, 0); break;                                          \
    case 1:  CASE_KSTEP(LX, 1); break;                                          \
    default: CASE_KSTEP(LX, 2); break;                                          \
    }

#define CASE_KSTEP_PADDED_SEL(LX, SEL)                                          \
    switch (SEL) {                                                              \
    case 0:  CASE_KSTEP_PADDED(LX, 0); break;                                   \
    case 1:  CASE_KSTEP_PADDED(LX, 1); break;                                   \
    default: CASE_KSTEP_PADDED(LX, 2); break;                                   \
    }

#define CASE(LX)                                                                \
    case LX:                                                                    \
      if(autotune[LX] == 0 ) {                                                  \
        autotune[LX]=tune<LX>( w,  u,                                           \
                               dx,  dy, dz,                                     \
                               dxt, dyt, dzt,h1,                                \
                               g11, g22, g33,                                   \
                               g12, g13, g23, nelv, lx,                         \
                               &autotune_eb[LX], &autotune_ch[LX]);             \
      } else if (autotune[LX] == 1 ) {                                          \
        CASE_1D_SEL(LX, autotune_ch[LX]);                                       \
      } else if (autotune[LX] == 2 ) {                                          \
        CASE_KSTEP_SEL(LX, autotune_eb[LX]);                                    \
      }                                                                         \
      break

#define CASE_PADDED(LX)                                                         \
    case LX:                                                                    \
      if(autotune[LX] == 0 ) {                                                  \
        autotune[LX]=tune_padded<LX>(w,  u,                                     \
                                     dx,  dy, dz,                               \
                                     dxt, dyt, dzt,h1,                          \
                                     g11, g22, g33,                             \
                                     g12, g13, g23,nelv,lx,                     \
                                     &autotune_eb[LX], &autotune_ch[LX]);       \
      } else if (autotune[LX] == 1 ) {                                          \
        CASE_1D_SEL(LX, autotune_ch[LX]);                                       \
      } else if (autotune[LX] == 2 ) {                                          \
        CASE_KSTEP_PADDED_SEL(LX, autotune_eb[LX]);                             \
      }                                                                         \
      break

/*
 * High order cases have no 1d variant to compare against (its LDS footprint
 * grows as LX^3), so they are not tuned and keep one element per block,
 * i.e. candidate 0
 */
#define CASE_LARGE(LX)                                                          \
    case LX:                                                                    \
      CASE_KSTEP(LX, 0);                                                        \
      break

#define CASE_LARGE_PADDED(LX)                                                   \
    case LX:                                                                    \
      CASE_KSTEP_PADDED(LX, 0);                                                 \
      break

    if ((*lx) < 13) {
      switch(*lx) {
        CASE(2);
        CASE(3);
        CASE(4);
        CASE(5);
        CASE(6);
        CASE(7);
        CASE_PADDED(8);
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
    else {
      switch(*lx) {
        CASE_LARGE(13);
        CASE_LARGE(14);
        CASE_LARGE(15);
        CASE_LARGE_PADDED(16);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
    }
  }
  /**
   * Fortran wrapper for device HIP Ax vector version
   */
  void hip_ax_helm_vector(void *au, void *av, void *aw,
                           void *u, void *v, void *w,
                           void *dx, void *dy, void *dz,
                           void *dxt, void *dyt, void *dzt,
                           void *h1, void *g11, void *g22,
                           void *g33, void *g12, void *g13,
                           void *g23, int *nelv, int *lx) {

    /* 0 = untuned, otherwise the chosen candidate plus one. Unlike the
       scalar Ax there is no 1d variant to fall back on, so the tuner only
       sweeps the elements per block candidates */
    static int autotune_v[17] = { 0 };

#define CASE_VECTOR_KSTEP(LX, C)                                                \
    hipLaunchKernelGGL( HIP_KERNEL_NAME(                                        \
                          ax_helm_kernel_vector_kstep<real, LX,                 \
                                                      NEKO_EB(LX, C)> ),        \
                        NEKO_EB_NBLCKS(*nelv, LX, C), NEKO_EB_NTHRDS(LX, C), 0, \
                        (hipStream_t) glb_cmd_queue,                            \
                        (real *) au, (real *) av, (real *) aw,                  \
                        (real *) u, (real *) v, (real *) w,                     \
                        (real *) dx, (real *) dy, (real *) dz,  (real *) h1,    \
                        (real *) g11, (real *) g22, (real *) g33,               \
                        (real *) g12, (real *) g13, (real *) g23, *nelv);       \
    HIP_CHECK(hipGetLastError());

#define CASE_VECTOR_KSTEP_PADDED(LX, C)                                         \
    hipLaunchKernelGGL( HIP_KERNEL_NAME(                                        \
                          ax_helm_kernel_vector_kstep_padded<real, LX,          \
                                                             NEKO_EB(LX, C)> ), \
                        NEKO_EB_NBLCKS(*nelv, LX, C), NEKO_EB_NTHRDS(LX, C), 0, \
                        (hipStream_t) glb_cmd_queue,                            \
                        (real *) au, (real *) av, (real *) aw,                  \
                        (real *) u, (real *) v, (real *) w,                     \
                        (real *) dx, (real *) dy, (real *) dz,  (real *) h1,    \
                        (real *) g11, (real *) g22, (real *) g33,               \
                        (real *) g12, (real *) g13, (real *) g23, *nelv);       \
    HIP_CHECK(hipGetLastError());

#define CASE_VECTOR_SEL(LX, SEL)                                                \
    switch (SEL) {                                                              \
    case 0:  CASE_VECTOR_KSTEP(LX, 0); break;                                   \
    case 1:  CASE_VECTOR_KSTEP(LX, 1); break;                                   \
    default: CASE_VECTOR_KSTEP(LX, 2); break;                                   \
    }

#define CASE_VECTOR_PADDED_SEL(LX, SEL)                                         \
    switch (SEL) {                                                              \
    case 0:  CASE_VECTOR_KSTEP_PADDED(LX, 0); break;                            \
    case 1:  CASE_VECTOR_KSTEP_PADDED(LX, 1); break;                            \
    default: CASE_VECTOR_KSTEP_PADDED(LX, 2); break;                            \
    }

#define CASE_VECTOR(LX)                                                         \
    case LX:                                                                    \
      if (autotune_v[LX] == 0) {                                                \
        autotune_v[LX] = 1 + tune_vector<LX>(au, av, aw, u, v, w,               \
                                             dx, dy, dz, h1,                    \
                                             g11, g22, g33,                     \
                                             g12, g13, g23, nelv, lx);          \
      } else {                                                                  \
        CASE_VECTOR_SEL(LX, autotune_v[LX] - 1);                                \
      }                                                                         \
      break

#define CASE_VECTOR_PADDED(LX)                                                  \
    case LX:                                                                    \
      if (autotune_v[LX] == 0) {                                                \
        autotune_v[LX] = 1 + tune_vector_padded<LX>(au, av, aw, u, v, w,        \
                                                    dx, dy, dz, h1,             \
                                                    g11, g22, g33,              \
                                                    g12, g13, g23, nelv, lx);   \
      } else {                                                                  \
        CASE_VECTOR_PADDED_SEL(LX, autotune_v[LX] - 1);                         \
      }                                                                         \
      break

    switch(*lx) {
      CASE_VECTOR(2);
      CASE_VECTOR(3);
      CASE_VECTOR_PADDED(4);
      CASE_VECTOR(5);
      CASE_VECTOR(6);
      CASE_VECTOR(7);
      CASE_VECTOR_PADDED(8);
      CASE_VECTOR(9);
      CASE_VECTOR(10);
      CASE_VECTOR(11);
      CASE_VECTOR(12);
      CASE_VECTOR(13);
      CASE_VECTOR(14);
      CASE_VECTOR(15);
      CASE_VECTOR_PADDED(16);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
  }

  /**
   * Fortran wrapper for device HIP Ax vector version part2
   */
  void hip_ax_helm_vector_part2(void *au, void *av, void *aw,
                                 void *u, void *v, void *w,
                                 void *h2, void *B, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;

    hipLaunchKernelGGL( HIP_KERNEL_NAME( ax_helm_kernel_vector_part2<real> ),
                        nblcks, nthrds, 0, stream,
                        (real *) au, (real *) av, (real *) aw,
                        (real *) u, (real *) v, (real *) w,
                        (real *) h2, (real *) B, *n);
  }

}

template < const int LX >
int tune(void *w, void *u, void *dx, void *dy, void *dz,
         void *dxt, void *dyt, void *dzt, void *h1,
         void *g11, void *g22, void *g33, void *g12,
         void *g13, void *g23, int *nelv, int *lx, int *eb_sel, int *ch_sel) {
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
  const dim3 nblcks_1d((*nelv), 1, 1);
  const hipStream_t stream = (hipStream_t) glb_cmd_queue;

  char *env_value = NULL;
  char neko_log_buf[80];

  env_value=getenv("NEKO_AUTOTUNE");

  sprintf(neko_log_buf, "Autotune Ax helm (lx: %d)", *lx);
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
      const int c = neko_eb_env();
      *eb_sel = c;
      CASE_KSTEP_SEL(LX, c);
      sprintf(neko_log_buf,"Set by env   : 2 (KSTEP, %d elem/block)",
              NEKO_EB_SEL(LX, c));
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

  /* Warm every variant before timing anything: each specialisation has to
     be resident and the clocks at steady state, or whichever variant is
     timed first is measured on a colder part */
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

  /* Interleaved rounds, best time per variant: timing them one after another
     in a fixed order lets clock drift bias the comparison by position */
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

  HIP_CHECK(hipEventDestroy(start));
  HIP_CHECK(hipEventDestroy(stop));

  /* The tuner stands in for a real Ax evaluation, and the variants do not
     sum in the same order, so leave the output of the chosen kernel in w */
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

template < const int LX >
int tune_padded(void *w, void *u, void *dx, void *dy, void *dz,
                void *dxt, void *dyt, void *dzt, void *h1,
                void *g11, void *g22, void *g33, void *g12,
                void *g13, void *g23, int *nelv, int *lx, int *eb_sel, int *ch_sel) {
  hipEvent_t start, stop;
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
  const dim3 nblcks_1d((*nelv), 1, 1);
  const hipStream_t stream = (hipStream_t) glb_cmd_queue;

  char *env_value = NULL;
  char neko_log_buf[80];

  env_value=getenv("NEKO_AUTOTUNE");

  sprintf(neko_log_buf, "Autotune Ax helm (lx: %d)", *lx);
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
      const int c = neko_eb_env();
      *eb_sel = c;
      CASE_KSTEP_PADDED_SEL(LX, c);
      sprintf(neko_log_buf,"Set by env   : 2 (KSTEP, %d elem/block)",
              NEKO_EB_SEL(LX, c));
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

  /* Warm every variant before timing anything: each specialisation has to
     be resident and the clocks at steady state, or whichever variant is
     timed first is measured on a colder part */
  for (int i = 0; i < NEKO_TUNE_WARMUP; i++) {
    CASE_1D(LX, 0);
    CASE_1D(LX, 1);
    CASE_1D(LX, 2);
    CASE_1D(LX, 3);
    CASE_KSTEP_PADDED(LX, 0);
    if (sweep) {
      CASE_KSTEP_PADDED(LX, 1);
      CASE_KSTEP_PADDED(LX, 2);
    }
  }

  /* Interleaved rounds, best time per variant: timing them one after another
     in a fixed order lets clock drift bias the comparison by position */
  for (int r = 0; r < rounds; r++) {
    NEKO_TUNE_TIME(time1, CASE_1D, LX, 0, iters);
    NEKO_TUNE_TIME(time1, CASE_1D, LX, 1, iters);
    NEKO_TUNE_TIME(time1, CASE_1D, LX, 2, iters);
    NEKO_TUNE_TIME(time1, CASE_1D, LX, 3, iters);
    NEKO_TUNE_TIME(time2, CASE_KSTEP_PADDED, LX, 0, iters);
    if (sweep) {
      NEKO_TUNE_TIME(time2, CASE_KSTEP_PADDED, LX, 1, iters);
      NEKO_TUNE_TIME(time2, CASE_KSTEP_PADDED, LX, 2, iters);
    }
  }

  NEKO_TUNE_LOG(LX, time1, time2);

  NEKO_TUNE_BEST(time1, best1, NEKO_CHUNKS_CANDIDATES);
  NEKO_TUNE_BEST(time2, best, NEKO_EB_CANDIDATES);
  *eb_sel = best;
  *ch_sel = best1;

  if (time1[best1] < time2[best]) {
    retval=1;
  } else {
    retval=2;
  }

  HIP_CHECK(hipEventDestroy(start));
  HIP_CHECK(hipEventDestroy(stop));

  /* Leave the output of the chosen kernel in w, see tune() */
  if (retval == 1) {
    CASE_1D_SEL(LX, best1);
  } else {
    CASE_KSTEP_PADDED_SEL(LX, best);
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

/*
 * The vector Ax has no 1d variant, so these only sweep the elements per
 * block candidates. Setting NEKO_AUTOTUNE to anything pins candidate 0,
 * i.e. the one element per block shape these kernels have always had,
 * which is the A/B baseline for the blocking.
 */
template < const int LX >
int tune_vector(void *au, void *av, void *aw, void *u, void *v, void *w,
                void *dx, void *dy, void *dz, void *h1,
                void *g11, void *g22, void *g33, void *g12,
                void *g13, void *g23, int *nelv, int *lx) {
  hipEvent_t start, stop;
  float time[NEKO_EB_CANDIDATES];
  const int rounds = neko_tune_rounds();
  const int iters = neko_tune_iters();
  const int sweep = neko_eb_sweep();
  int best = 0;

  for (int c = 0; c < NEKO_EB_CANDIDATES; c++) {
    time[c] = NEKO_TUNE_INIT;
  }

  const hipStream_t stream = (hipStream_t) glb_cmd_queue;
  char neko_log_buf[80];

  sprintf(neko_log_buf, "Autotune Ax helm vector (lx: %d)", *lx);
  log_section(neko_log_buf);

  /*
   * Default to one element per block, as on the CUDA backend. These kernels
   * carry six T[LX] register arrays against the scalar kernels' two, so they
   * reach the register file limit -- and start spilling to scratch -- well
   * before a wider block can buy anything. Sweeping them by default was a
   * regression. Set NEKO_EB_TUNE=1 to opt back into the sweep
   */
  if (getenv("NEKO_AUTOTUNE") ||
      getenv("NEKO_EB") || !neko_eb_sweep()) {
    const int c = neko_eb_env();
    CASE_VECTOR_SEL(LX, c);
    sprintf(neko_log_buf, "Vector Ax    : %d elem/block (untuned)",
            NEKO_EB_SEL(LX, c));
    log_message(neko_log_buf);
    log_end_section();
    return c;
  }

  HIP_CHECK(hipEventCreate(&start));
  HIP_CHECK(hipEventCreate(&stop));

  for (int i = 0; i < NEKO_TUNE_WARMUP; i++) {
    CASE_VECTOR_KSTEP(LX, 0);
    CASE_VECTOR_KSTEP(LX, 1);
    CASE_VECTOR_KSTEP(LX, 2);
  }

  for (int r = 0; r < rounds; r++) {
    NEKO_TUNE_TIME(time, CASE_VECTOR_KSTEP, LX, 0, iters);
    NEKO_TUNE_TIME(time, CASE_VECTOR_KSTEP, LX, 1, iters);
    NEKO_TUNE_TIME(time, CASE_VECTOR_KSTEP, LX, 2, iters);
  }

  for (int c = 0; c < NEKO_EB_CANDIDATES; c++) {
    sprintf(neko_log_buf, "eb=%-4d      : %9.2f us/call",
            NEKO_EB_SEL(LX, c), time[c] * 10.0);
    log_message(neko_log_buf);
  }

  NEKO_TUNE_BEST(time, best, NEKO_EB_CANDIDATES);

  HIP_CHECK(hipEventDestroy(start));
  HIP_CHECK(hipEventDestroy(stop));

  /* Leave the output of the chosen kernel in au/av/aw */
  CASE_VECTOR_SEL(LX, best);

  sprintf(neko_log_buf, "Chose        : %d elem/block", NEKO_EB_SEL(LX, best));
  log_message(neko_log_buf);
  log_end_section();
  return best;
}

template < const int LX >
int tune_vector_padded(void *au, void *av, void *aw,
                       void *u, void *v, void *w,
                       void *dx, void *dy, void *dz, void *h1,
                       void *g11, void *g22, void *g33, void *g12,
                       void *g13, void *g23, int *nelv, int *lx) {
  hipEvent_t start, stop;
  float time[NEKO_EB_CANDIDATES];
  const int rounds = neko_tune_rounds();
  const int iters = neko_tune_iters();
  const int sweep = neko_eb_sweep();
  int best = 0;

  for (int c = 0; c < NEKO_EB_CANDIDATES; c++) {
    time[c] = NEKO_TUNE_INIT;
  }

  const hipStream_t stream = (hipStream_t) glb_cmd_queue;
  char neko_log_buf[80];

  sprintf(neko_log_buf, "Autotune Ax helm vector (lx: %d)", *lx);
  log_section(neko_log_buf);

  /*
   * Default to one element per block, as on the CUDA backend. These kernels
   * carry six T[LX] register arrays against the scalar kernels' two, so they
   * reach the register file limit -- and start spilling to scratch -- well
   * before a wider block can buy anything. Sweeping them by default was a
   * regression. Set NEKO_EB_TUNE=1 to opt back into the sweep
   */
  if (getenv("NEKO_AUTOTUNE") ||
      getenv("NEKO_EB") || !neko_eb_sweep()) {
    const int c = neko_eb_env();
    CASE_VECTOR_PADDED_SEL(LX, c);
    sprintf(neko_log_buf, "Vector Ax    : %d elem/block (untuned)",
            NEKO_EB_SEL(LX, c));
    log_message(neko_log_buf);
    log_end_section();
    return c;
  }

  HIP_CHECK(hipEventCreate(&start));
  HIP_CHECK(hipEventCreate(&stop));

  for (int i = 0; i < NEKO_TUNE_WARMUP; i++) {
    CASE_VECTOR_KSTEP_PADDED(LX, 0);
    CASE_VECTOR_KSTEP_PADDED(LX, 1);
    CASE_VECTOR_KSTEP_PADDED(LX, 2);
  }

  for (int r = 0; r < rounds; r++) {
    NEKO_TUNE_TIME(time, CASE_VECTOR_KSTEP_PADDED, LX, 0, iters);
    NEKO_TUNE_TIME(time, CASE_VECTOR_KSTEP_PADDED, LX, 1, iters);
    NEKO_TUNE_TIME(time, CASE_VECTOR_KSTEP_PADDED, LX, 2, iters);
  }

  for (int c = 0; c < NEKO_EB_CANDIDATES; c++) {
    sprintf(neko_log_buf, "eb=%-4d      : %9.2f us/call",
            NEKO_EB_SEL(LX, c), time[c] * 10.0);
    log_message(neko_log_buf);
  }

  NEKO_TUNE_BEST(time, best, NEKO_EB_CANDIDATES);

  HIP_CHECK(hipEventDestroy(start));
  HIP_CHECK(hipEventDestroy(stop));

  /* Leave the output of the chosen kernel in au/av/aw */
  CASE_VECTOR_PADDED_SEL(LX, best);

  sprintf(neko_log_buf, "Chose        : %d elem/block", NEKO_EB_SEL(LX, best));
  log_message(neko_log_buf);
  log_end_section();
  return best;
}
