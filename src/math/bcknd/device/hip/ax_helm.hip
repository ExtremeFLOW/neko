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
         void *g13, void *g23, int *nelv, int *lx, int *eb_sel, int *ch_sel,
         int *nwf_sel);

template < const int>
int tune_padded(void *w, void *u, void *dx, void *dy, void *dz,
                void *dxt, void *dyt, void *dzt, void *h1,
                void *g11, void *g22, void *g33, void *g12,
                void *g13, void *g23, int *nelv, int *lx, int *eb_sel,
                int *ch_sel, int *nwf_sel);

template < const int>
int tune_vector(void *au, void *av, void *aw, void *u, void *v, void *w,
                void *dx, void *dy, void *dz, void *h1,
                void *g11, void *g22, void *g33, void *g12,
                void *g13, void *g23, int *nelv, int *lx, int *eb_sel,
                int *nwf_sel);

template < const int>
int tune_vector_padded(void *au, void *av, void *aw,
                       void *u, void *v, void *w,
                       void *dx, void *dy, void *dz, void *h1,
                       void *g11, void *g22, void *g33, void *g12,
                       void *g13, void *g23, int *nelv, int *lx, int *eb_sel,
                       int *nwf_sel);

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
    /* wavefronts per block candidate chosen for the mfma variant */
    static int autotune_nwf[13] = { 0 };

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

#define CASE_MFMA(LX, C)                                                        \
    hipLaunchKernelGGL( HIP_KERNEL_NAME(                                        \
                          ax_helm_kernel_mfma<real, LX, NEKO_MFMA_NWF(C),       \
                                              NEKO_MFMA_TILE(C)> ),             \
                        NEKO_MFMA_NBLCKS(*nelv, LX, C), NEKO_MFMA_NTHRDS(C), 0, \
                        (hipStream_t) glb_cmd_queue,                            \
                        (real *) w, (real *) u,                                 \
                        (real *) dx, (real *) dy, (real *) dz,  (real *) h1,    \
                        (real *) g11, (real *) g22, (real *) g33,               \
                        (real *) g12, (real *) g13, (real *) g23, *nelv);       \
    HIP_CHECK(hipGetLastError());

/* Runtime dispatch onto the tuned candidate, which names both the wavefronts
   per block and the matrix core tile, see NEKO_MFMA_TILE in mfma_kernel.h */
#define CASE_MFMA_SEL(LX, SEL)                                                  \
    switch (SEL) {                                                              \
    case 0:  CASE_MFMA(LX, 0); break;                                           \
    case 1:  CASE_MFMA(LX, 1); break;                                           \
    case 2:  CASE_MFMA(LX, 2); break;                                           \
    case 3:  CASE_MFMA(LX, 3); break;                                           \
    case 4:  CASE_MFMA(LX, 4); break;                                           \
    case 5:  CASE_MFMA(LX, 5); break;                                           \
    case 6:  CASE_MFMA(LX, 6); break;                                           \
    default: CASE_MFMA(LX, 7); break;                                           \
    }

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
                               &autotune_eb[LX], &autotune_ch[LX],              \
                               &autotune_nwf[LX]);                              \
      } else if (autotune[LX] == 1 ) {                                          \
        CASE_1D_SEL(LX, autotune_ch[LX]);                                       \
      } else if (autotune[LX] == 2 ) {                                          \
        CASE_KSTEP_SEL(LX, autotune_eb[LX]);                                    \
      } else if (autotune[LX] == 3 ) {                                          \
        CASE_MFMA_SEL(LX, autotune_nwf[LX]);                                    \
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
                                     &autotune_eb[LX], &autotune_ch[LX],        \
                                     &autotune_nwf[LX]);                        \
      } else if (autotune[LX] == 1 ) {                                          \
        CASE_1D_SEL(LX, autotune_ch[LX]);                                       \
      } else if (autotune[LX] == 2 ) {                                          \
        CASE_KSTEP_PADDED_SEL(LX, autotune_eb[LX]);                             \
      } else if (autotune[LX] == 3 ) {                                          \
        CASE_MFMA_SEL(LX, autotune_nwf[LX]);                                    \
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

    /* Strategy chosen by the autotuner: 0 untuned, 2 kstep, 3 mfma. There
       is no 1d variant of the vector operator, so the kstep family is the
       baseline and the matrix core variant the only alternative */
    static int autotune_v[17] = { 0 };
    /* elements per block candidate chosen for the kstep variant */
    static int autotune_v_eb[17] = { 0 };
    /* wavefronts per block candidate chosen for the mfma variant */
    static int autotune_v_nwf[17] = { 0 };

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

/* The matrix core variant: same grid and wavefronts per block candidates as
   the scalar CASE_MFMA, see mfma_kernel.h */
#define CASE_VECTOR_MFMA(LX, C)                                                 \
    hipLaunchKernelGGL( HIP_KERNEL_NAME(                                        \
                          ax_helm_kernel_mfma_vector<real, LX,                  \
                                                     NEKO_MFMA_NWF(C),          \
                                                     NEKO_MFMA_TILE(C)> ),      \
                        NEKO_MFMA_NBLCKS(*nelv, LX, C), NEKO_MFMA_NTHRDS(C), 0, \
                        (hipStream_t) glb_cmd_queue,                            \
                        (real *) au, (real *) av, (real *) aw,                  \
                        (real *) u, (real *) v, (real *) w,                     \
                        (real *) dx, (real *) dy, (real *) dz,  (real *) h1,    \
                        (real *) g11, (real *) g22, (real *) g33,               \
                        (real *) g12, (real *) g13, (real *) g23, *nelv);       \
    HIP_CHECK(hipGetLastError());

#define CASE_VECTOR_KSTEP_SEL(LX, SEL)                                          \
    switch (SEL) {                                                              \
    case 0:  CASE_VECTOR_KSTEP(LX, 0); break;                                   \
    case 1:  CASE_VECTOR_KSTEP(LX, 1); break;                                   \
    default: CASE_VECTOR_KSTEP(LX, 2); break;                                   \
    }

#define CASE_VECTOR_KSTEP_PADDED_SEL(LX, SEL)                                   \
    switch (SEL) {                                                              \
    case 0:  CASE_VECTOR_KSTEP_PADDED(LX, 0); break;                            \
    case 1:  CASE_VECTOR_KSTEP_PADDED(LX, 1); break;                            \
    default: CASE_VECTOR_KSTEP_PADDED(LX, 2); break;                            \
    }

/* Runtime dispatch onto the tuned candidate, see CASE_MFMA_SEL */
#define CASE_VECTOR_MFMA_SEL(LX, SEL)                                           \
    switch (SEL) {                                                              \
    case 0:  CASE_VECTOR_MFMA(LX, 0); break;                                    \
    case 1:  CASE_VECTOR_MFMA(LX, 1); break;                                    \
    case 2:  CASE_VECTOR_MFMA(LX, 2); break;                                    \
    case 3:  CASE_VECTOR_MFMA(LX, 3); break;                                    \
    case 4:  CASE_VECTOR_MFMA(LX, 4); break;                                    \
    case 5:  CASE_VECTOR_MFMA(LX, 5); break;                                    \
    case 6:  CASE_VECTOR_MFMA(LX, 6); break;                                    \
    default: CASE_VECTOR_MFMA(LX, 7); break;                                    \
    }

#define CASE_VECTOR(LX)                                                         \
    case LX:                                                                    \
      if (autotune_v[LX] == 0) {                                                \
        autotune_v[LX] = tune_vector<LX>(au, av, aw, u, v, w,                   \
                                         dx, dy, dz, h1,                        \
                                         g11, g22, g33,                         \
                                         g12, g13, g23, nelv, lx,               \
                                         &autotune_v_eb[LX],                    \
                                         &autotune_v_nwf[LX]);                  \
      } else if (autotune_v[LX] == 2) {                                         \
        CASE_VECTOR_KSTEP_SEL(LX, autotune_v_eb[LX]);                           \
      } else if (autotune_v[LX] == 3) {                                         \
        CASE_VECTOR_MFMA_SEL(LX, autotune_v_nwf[LX]);                           \
      }                                                                         \
      break

#define CASE_VECTOR_PADDED(LX)                                                  \
    case LX:                                                                    \
      if (autotune_v[LX] == 0) {                                                \
        autotune_v[LX] = tune_vector_padded<LX>(au, av, aw, u, v, w,            \
                                                dx, dy, dz, h1,                 \
                                                g11, g22, g33,                  \
                                                g12, g13, g23, nelv, lx,        \
                                                &autotune_v_eb[LX],             \
                                                &autotune_v_nwf[LX]);           \
      } else if (autotune_v[LX] == 2) {                                         \
        CASE_VECTOR_KSTEP_PADDED_SEL(LX, autotune_v_eb[LX]);                    \
      } else if (autotune_v[LX] == 3) {                                         \
        CASE_VECTOR_MFMA_SEL(LX, autotune_v_nwf[LX]);                           \
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
         void *g13, void *g23, int *nelv, int *lx, int *eb_sel, int *ch_sel,
         int *nwf_sel) {
  hipEvent_t start,stop;
  float time1[NEKO_CHUNKS_CANDIDATES];
  int best1 = 0;
  float time2[NEKO_EB_CANDIDATES];
  float time3[NEKO_MFMA_CANDIDATES];
  int best3 = 0;
  const int rounds = neko_tune_rounds();
  const int iters = neko_tune_iters();
  /* Candidates of the kstep sweep, one -- the unblocked shape -- when the
     elements per block sweep is off */
  const int eb_cand = neko_eb_sweep() ? NEKO_EB_CANDIDATES : 1;
  /* Candidates of the matrix core sweep, see neko_mfma_candidates() */
  const int mfma_cand = neko_mfma_candidates();
  /* Geometry pinned by each formulation's own variable, -1 to sweep it */
  const int ch_pin = neko_chunks_pin();
  const int eb_pin = neko_eb_pin();
  const int nwf_pin = neko_mfma_pin();
  /* Formulation pinned by NEKO_AUTOTUNE, as the identifier tune() returns */
  int strat = 0;
  /* Hardware capability, which the explicit NEKO_AUTOTUNE=MFMA pin needs */
  const bool mfma = mfma_lx_supported<LX>() && hip_have_mfma();
  int best = 0;
  int retval;

  for (int c = 0; c < NEKO_EB_CANDIDATES; c++) {
    time2[c] = NEKO_TUNE_INIT;
  }
  for (int c = 0; c < NEKO_CHUNKS_CANDIDATES; c++) {
    time1[c] = NEKO_TUNE_INIT;
  }
  for (int c = 0; c < NEKO_MFMA_CANDIDATES; c++) {
    time3[c] = NEKO_TUNE_INIT;
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
  *nwf_sel = 0;

  /*
   * NEKO_AUTOTUNE names a formulation, and that is all it does: the sweep
   * below is narrowed to that one kernel family, but its geometry -- the
   * chunk size, the elements per block, the wavefronts per block -- is still
   * measured candidate against candidate. A formulation this build or this
   * device does not have is reported and ignored, leaving the full sweep.
   */
  if(env_value) {
    if( !strcmp(env_value,"1D") ) {
      strat = 1;
    } else if( !strcmp(env_value,"KSTEP") ) {
      strat = 2;
    } else if( !strcmp(env_value,"MFMA") ) {
      if (mfma) {
        strat = 3;
      } else {
        sprintf(neko_log_buf, "MFMA strategy not available for this config");
        log_error(neko_log_buf);
      }
    } else {
       sprintf(neko_log_buf, "Invalid value set for NEKO_AUTOTUNE");
       log_error(neko_log_buf);
    }
  }

  /* Geometry of the pinned formulation, if its own variable fixes that too.
     Both pinned leaves nothing to measure, so the kernel is launched once and
     reported, which is what pinning has always done */
  const int pin = (strat == 1) ? ch_pin : (strat == 2) ? eb_pin :
                  (strat == 3) ? nwf_pin : -1;

  if (pin >= 0) {
    switch (strat) {
    case 1:
      *ch_sel = pin;
      CASE_1D_SEL(LX, pin);
      sprintf(neko_log_buf, "Set by env   : 1 (1D, %d chunk)",
              NEKO_CHUNKS_SEL(LX, pin));
      break;
    case 2:
      *eb_sel = pin;
      CASE_KSTEP_SEL(LX, pin);
      sprintf(neko_log_buf, "Set by env   : 2 (KSTEP, %d elem/block)",
              NEKO_EB_SEL(LX, pin));
      break;
    default:
      *nwf_sel = pin;
      CASE_MFMA_SEL(LX, pin);
      sprintf(neko_log_buf, "Set by env   : 3 (MFMA %s, %d wf, %d elem/blk)",
              mfma_tile_name(pin), NEKO_MFMA_NWF(pin), NEKO_MFMA_EB(LX, pin));
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

  /* Formulations the sweep considers, see NEKO_TUNE_FOR(). NEKO_MFMA_TUNE
     keeps the matrix core variant out of an unpinned sweep, but naming it
     explicitly still measures it */
  const bool try_1d = (strat == 0 || strat == 1);
  const bool try_kstep = (strat == 0 || strat == 2);
  const bool try_mfma = mfma && ((strat == 3) ||
                                 (strat == 0 && neko_mfma_sweep()));

  HIP_CHECK(hipEventCreate(&start));
  HIP_CHECK(hipEventCreate(&stop));

  /* Warm every variant before timing anything: each specialisation has to
     be resident and the clocks at steady state, or whichever variant is
     timed first is measured on a colder part */
  for (int i = 0; i < NEKO_TUNE_WARMUP; i++) {
    NEKO_TUNE_FOR(c, try_1d, ch_pin, NEKO_CHUNKS_CANDIDATES) {
      CASE_1D_SEL(LX, c);
    }
    NEKO_TUNE_FOR(c, try_kstep, eb_pin, eb_cand) {
      CASE_KSTEP_SEL(LX, c);
    }
    NEKO_TUNE_FOR(c, try_mfma, nwf_pin, mfma_cand) {
      CASE_MFMA_SEL(LX, c);
    }
  }

  /* Interleaved rounds, best time per variant: timing them one after another
     in a fixed order lets clock drift bias the comparison by position */
  for (int r = 0; r < rounds; r++) {
    NEKO_TUNE_FOR(c, try_1d, ch_pin, NEKO_CHUNKS_CANDIDATES) {
      NEKO_TUNE_TIME(time1, CASE_1D_SEL, LX, c, iters);
    }
    NEKO_TUNE_FOR(c, try_kstep, eb_pin, eb_cand) {
      NEKO_TUNE_TIME(time2, CASE_KSTEP_SEL, LX, c, iters);
    }
    NEKO_TUNE_FOR(c, try_mfma, nwf_pin, mfma_cand) {
      NEKO_TUNE_TIME(time3, CASE_MFMA_SEL, LX, c, iters);
    }
  }

  NEKO_TUNE_LOG(LX, time1, time2);
  NEKO_TUNE_LOG_MFMA(LX, time3);

  NEKO_TUNE_BEST(time1, best1, NEKO_CHUNKS_CANDIDATES);
  NEKO_TUNE_BEST(time2, best, NEKO_EB_CANDIDATES);
  NEKO_TUNE_BEST(time3, best3, NEKO_MFMA_CANDIDATES);
  *eb_sel = best;
  *ch_sel = best1;
  *nwf_sel = best3;

  if (time1[best1] < time2[best]) {
    retval = 1;
  } else {
    retval = 2;
  }

  /* The mfma variant joins the comparison only where it exists, its
     candidates are left at NEKO_TUNE_INIT otherwise */
  if (time3[best3] < ((retval == 1) ? time1[best1] : time2[best])) {
    retval = 3;
  }

  HIP_CHECK(hipEventDestroy(start));
  HIP_CHECK(hipEventDestroy(stop));

  /* The tuner stands in for a real Ax evaluation, and the variants do not
     sum in the same order, so leave the output of the chosen kernel in w */
  if (retval == 1) {
    CASE_1D_SEL(LX, best1);
  } else if (retval == 2) {
    CASE_KSTEP_SEL(LX, best);
  } else {
    CASE_MFMA_SEL(LX, best3);
  }

  if (retval == 1) {
    sprintf(neko_log_buf, "Chose        : 1 (1D, %d chunk)",
            NEKO_CHUNKS_SEL(LX, best1));
  } else if (retval == 2) {
    sprintf(neko_log_buf, "Chose        : 2 (KSTEP, %d elem/block)",
            NEKO_EB_SEL(LX, best));
  } else {
    sprintf(neko_log_buf, "Chose        : 3 (MFMA %s, %d wf, %d elem/blk)",
            mfma_tile_name(best3), NEKO_MFMA_NWF(best3),
            NEKO_MFMA_EB(LX, best3));
  }
  log_message(neko_log_buf);
  log_end_section();
  return retval;
}

template < const int LX >
int tune_padded(void *w, void *u, void *dx, void *dy, void *dz,
                void *dxt, void *dyt, void *dzt, void *h1,
                void *g11, void *g22, void *g33, void *g12,
                void *g13, void *g23, int *nelv, int *lx,
                int *eb_sel, int *ch_sel, int *nwf_sel) {
  hipEvent_t start, stop;
  float time1[NEKO_CHUNKS_CANDIDATES];
  int best1 = 0;
  float time2[NEKO_EB_CANDIDATES];
  float time3[NEKO_MFMA_CANDIDATES];
  int best3 = 0;
  const int rounds = neko_tune_rounds();
  const int iters = neko_tune_iters();
  /* Candidates of the kstep sweep, see tune() */
  const int eb_cand = neko_eb_sweep() ? NEKO_EB_CANDIDATES : 1;
  /* Candidates of the matrix core sweep, see neko_mfma_candidates() */
  const int mfma_cand = neko_mfma_candidates();
  /* Geometry pinned by each formulation's own variable, -1 to sweep it */
  const int ch_pin = neko_chunks_pin();
  const int eb_pin = neko_eb_pin();
  const int nwf_pin = neko_mfma_pin();
  /* Formulation pinned by NEKO_AUTOTUNE, as the identifier this returns */
  int strat = 0;
  /* Hardware capability, which the explicit NEKO_AUTOTUNE=MFMA pin needs */
  const bool mfma = mfma_lx_supported<LX>() && hip_have_mfma();
  int best = 0;
  int retval;

  for (int c = 0; c < NEKO_EB_CANDIDATES; c++) {
    time2[c] = NEKO_TUNE_INIT;
  }
  for (int c = 0; c < NEKO_CHUNKS_CANDIDATES; c++) {
    time1[c] = NEKO_TUNE_INIT;
  }
  for (int c = 0; c < NEKO_MFMA_CANDIDATES; c++) {
    time3[c] = NEKO_TUNE_INIT;
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
  *nwf_sel = 0;

  /* NEKO_AUTOTUNE names a formulation and nothing more, see tune() */
  if(env_value) {
    if( !strcmp(env_value,"1D") ) {
      strat = 1;
    } else if( !strcmp(env_value,"KSTEP") ) {
      strat = 2;
    } else if( !strcmp(env_value,"MFMA") ) {
      if (mfma) {
        strat = 3;
      } else {
        sprintf(neko_log_buf, "MFMA strategy not available for this config");
        log_error(neko_log_buf);
      }
    } else {
      sprintf(neko_log_buf, "Invalid value set for NEKO_AUTOTUNE");
      log_error(neko_log_buf);
    }
  }

  /* Formulation and geometry both pinned leaves nothing to measure */
  const int pin = (strat == 1) ? ch_pin : (strat == 2) ? eb_pin :
                  (strat == 3) ? nwf_pin : -1;

  if (pin >= 0) {
    switch (strat) {
    case 1:
      *ch_sel = pin;
      CASE_1D_SEL(LX, pin);
      sprintf(neko_log_buf, "Set by env   : 1 (1D, %d chunk)",
              NEKO_CHUNKS_SEL(LX, pin));
      break;
    case 2:
      *eb_sel = pin;
      CASE_KSTEP_PADDED_SEL(LX, pin);
      sprintf(neko_log_buf, "Set by env   : 2 (KSTEP, %d elem/block)",
              NEKO_EB_SEL(LX, pin));
      break;
    default:
      *nwf_sel = pin;
      CASE_MFMA_SEL(LX, pin);
      sprintf(neko_log_buf, "Set by env   : 3 (MFMA %s, %d wf, %d elem/blk)",
              mfma_tile_name(pin), NEKO_MFMA_NWF(pin), NEKO_MFMA_EB(LX, pin));
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
  const bool try_mfma = mfma && ((strat == 3) ||
                                 (strat == 0 && neko_mfma_sweep()));

  HIP_CHECK(hipEventCreate(&start));
  HIP_CHECK(hipEventCreate(&stop));

  /* Warm every variant before timing anything: each specialisation has to
     be resident and the clocks at steady state, or whichever variant is
     timed first is measured on a colder part */
  for (int i = 0; i < NEKO_TUNE_WARMUP; i++) {
    NEKO_TUNE_FOR(c, try_1d, ch_pin, NEKO_CHUNKS_CANDIDATES) {
      CASE_1D_SEL(LX, c);
    }
    NEKO_TUNE_FOR(c, try_kstep, eb_pin, eb_cand) {
      CASE_KSTEP_PADDED_SEL(LX, c);
    }
    NEKO_TUNE_FOR(c, try_mfma, nwf_pin, mfma_cand) {
      CASE_MFMA_SEL(LX, c);
    }
  }

  /* Interleaved rounds, best time per variant: timing them one after another
     in a fixed order lets clock drift bias the comparison by position */
  for (int r = 0; r < rounds; r++) {
    NEKO_TUNE_FOR(c, try_1d, ch_pin, NEKO_CHUNKS_CANDIDATES) {
      NEKO_TUNE_TIME(time1, CASE_1D_SEL, LX, c, iters);
    }
    NEKO_TUNE_FOR(c, try_kstep, eb_pin, eb_cand) {
      NEKO_TUNE_TIME(time2, CASE_KSTEP_PADDED_SEL, LX, c, iters);
    }
    NEKO_TUNE_FOR(c, try_mfma, nwf_pin, mfma_cand) {
      NEKO_TUNE_TIME(time3, CASE_MFMA_SEL, LX, c, iters);
    }
  }

  NEKO_TUNE_LOG(LX, time1, time2);
  NEKO_TUNE_LOG_MFMA(LX, time3);

  NEKO_TUNE_BEST(time1, best1, NEKO_CHUNKS_CANDIDATES);
  NEKO_TUNE_BEST(time2, best, NEKO_EB_CANDIDATES);
  NEKO_TUNE_BEST(time3, best3, NEKO_MFMA_CANDIDATES);
  *eb_sel = best;
  *ch_sel = best1;
  *nwf_sel = best3;

  if (time1[best1] < time2[best]) {
    retval = 1;
  } else {
    retval = 2;
  }

  /* The mfma variant joins the comparison only where it exists, its
     candidates are left at NEKO_TUNE_INIT otherwise */
  if (time3[best3] < ((retval == 1) ? time1[best1] : time2[best])) {
    retval = 3;
  }

  HIP_CHECK(hipEventDestroy(start));
  HIP_CHECK(hipEventDestroy(stop));

  /* Leave the output of the chosen kernel in w, see tune() */
  if (retval == 1) {
    CASE_1D_SEL(LX, best1);
  } else if (retval == 2) {
    CASE_KSTEP_PADDED_SEL(LX, best);
  } else {
    CASE_MFMA_SEL(LX, best3);
  }

  if (retval == 1) {
    sprintf(neko_log_buf, "Chose        : 1 (1D, %d chunk)",
            NEKO_CHUNKS_SEL(LX, best1));
  } else if (retval == 2) {
    sprintf(neko_log_buf, "Chose        : 2 (KSTEP, %d elem/block)",
            NEKO_EB_SEL(LX, best));
  } else {
    sprintf(neko_log_buf, "Chose        : 3 (MFMA %s, %d wf, %d elem/blk)",
            mfma_tile_name(best3), NEKO_MFMA_NWF(best3),
            NEKO_MFMA_EB(LX, best3));
  }
  log_message(neko_log_buf);
  log_end_section();
  return retval;
}

/*
 * The vector Ax has no 1d variant, so the tuner picks between the kstep
 * family -- one element per block, plus the blocked candidates when the
 * elements per block sweep is on -- and the matrix core variant, where the
 * device and the polynomial order allow it. NEKO_AUTOTUNE=KSTEP (or 1D,
 * there being nothing else to pin it to) narrows the sweep to the kstep
 * variant, NEKO_AUTOTUNE=MFMA to the matrix core one, and in both cases the
 * geometry candidates of the named variant are still measured -- it takes
 * NEKO_EB or NEKO_MFMA_NWF to fix that too, as for the scalar operator.
 *
 * @note The section title has to fit the 30 character header field in
 * log_section(), which truncates rather than wraps -- hence "Ax vector"
 * rather than the "Ax helm vector" it reads as. It matches the CUDA copy.
 */
template < const int LX >
int tune_vector(void *au, void *av, void *aw, void *u, void *v, void *w,
                void *dx, void *dy, void *dz, void *h1,
                void *g11, void *g22, void *g33, void *g12,
                void *g13, void *g23, int *nelv, int *lx, int *eb_sel,
                int *nwf_sel) {
  hipEvent_t start, stop;
  float time2[NEKO_EB_CANDIDATES];
  int best2 = 0;
  float time3[NEKO_MFMA_CANDIDATES];
  int best3 = 0;
  const int rounds = neko_tune_rounds();
  const int iters = neko_tune_iters();
  const int sweep = neko_eb_sweep();
  /* Candidates of the kstep sweep, one -- the unblocked shape -- when the
     elements per block sweep is off */
  const int eb_cand = sweep ? NEKO_EB_CANDIDATES : 1;
  /* Candidates of the matrix core sweep, see neko_mfma_candidates() */
  const int mfma_cand = neko_mfma_candidates();
  /* Geometry pinned by each formulation's own variable, -1 to sweep it */
  const int eb_pin = neko_eb_pin();
  const int nwf_pin = neko_mfma_pin();
  /* Formulation pinned by NEKO_AUTOTUNE, as the identifier this returns */
  int strat = 0;
  /* Hardware capability, which the explicit NEKO_AUTOTUNE=MFMA pin needs */
  const bool mfma = mfma_lx_supported<LX>() && hip_have_mfma();
  int retval;

  for (int c = 0; c < NEKO_EB_CANDIDATES; c++) {
    time2[c] = NEKO_TUNE_INIT;
  }
  for (int c = 0; c < NEKO_MFMA_CANDIDATES; c++) {
    time3[c] = NEKO_TUNE_INIT;
  }

  const hipStream_t stream = (hipStream_t) glb_cmd_queue;

  char *env_value = NULL;
  char neko_log_buf[80];

  env_value = getenv("NEKO_AUTOTUNE");

  sprintf(neko_log_buf, "Autotune Ax vector (lx: %d)", *lx);
  log_section(neko_log_buf);

  *eb_sel = 0;
  *nwf_sel = 0;

  /* NEKO_AUTOTUNE names a formulation and nothing more, see tune() */
  if (env_value) {
    if ( !strcmp(env_value, "KSTEP") || !strcmp(env_value, "1D") ) {
      /* No 1d variant here, so either pin lands on kstep */
      strat = 2;
    } else if ( !strcmp(env_value, "MFMA") ) {
      if (mfma) {
        strat = 3;
      } else {
        sprintf(neko_log_buf, "MFMA strategy not available for this config");
        log_error(neko_log_buf);
      }
    } else {
      sprintf(neko_log_buf, "Invalid value set for NEKO_AUTOTUNE");
      log_error(neko_log_buf);
    }
  }

  /* Formulation and geometry both pinned leaves nothing to measure */
  const int pin = (strat == 2) ? eb_pin : (strat == 3) ? nwf_pin : -1;

  if (pin >= 0) {
    if (strat == 2) {
      *eb_sel = pin;
      CASE_VECTOR_KSTEP_SEL(LX, pin);
      sprintf(neko_log_buf, "Set by env   : 2 (KSTEP, %d elem/block)",
              NEKO_EB_SEL(LX, pin));
    } else {
      *nwf_sel = pin;
      CASE_VECTOR_MFMA_SEL(LX, pin);
      sprintf(neko_log_buf, "Set by env   : 3 (MFMA %s, %d wf, %d elem/blk)",
              mfma_tile_name(pin), NEKO_MFMA_NWF(pin), NEKO_MFMA_EB(LX, pin));
    }
    log_message(neko_log_buf);
    log_end_section();
    return strat;
  }

  if (strat) {
    sprintf(neko_log_buf, "Set by env   : %d (%s)", strat, env_value);
    log_message(neko_log_buf);
  }

  /* Formulations the sweep considers, see NEKO_TUNE_FOR(). NEKO_MFMA_TUNE
     keeps the matrix core variant out of an unpinned sweep, but naming it
     explicitly still measures it */
  const bool try_kstep = (strat == 0 || strat == 2);
  const bool try_mfma = mfma && ((strat == 3) ||
                                 (strat == 0 && neko_mfma_sweep()));

  /* With neither the blocking sweep nor a matrix core variant there is one
     candidate, which is what these kernels have always run */
  if (!sweep && !try_mfma) {
    CASE_VECTOR_KSTEP(LX, 0);
    sprintf(neko_log_buf, "Vector Ax    : %d elem/block (untuned)",
            NEKO_EB_SEL(LX, 0));
    log_message(neko_log_buf);
    log_end_section();
    return 2;
  }

  HIP_CHECK(hipEventCreate(&start));
  HIP_CHECK(hipEventCreate(&stop));

  /* Warm every variant before timing anything, see tune() */
  for (int i = 0; i < NEKO_TUNE_WARMUP; i++) {
    NEKO_TUNE_FOR(c, try_kstep, eb_pin, eb_cand) {
      CASE_VECTOR_KSTEP_SEL(LX, c);
    }
    NEKO_TUNE_FOR(c, try_mfma, nwf_pin, mfma_cand) {
      CASE_VECTOR_MFMA_SEL(LX, c);
    }
  }

  /* Interleaved rounds, best time per variant, see tune() */
  for (int r = 0; r < rounds; r++) {
    NEKO_TUNE_FOR(c, try_kstep, eb_pin, eb_cand) {
      NEKO_TUNE_TIME(time2, CASE_VECTOR_KSTEP_SEL, LX, c, iters);
    }
    NEKO_TUNE_FOR(c, try_mfma, nwf_pin, mfma_cand) {
      NEKO_TUNE_TIME(time3, CASE_VECTOR_MFMA_SEL, LX, c, iters);
    }
  }

  NEKO_TUNE_LOG_KSTEP(LX, time2);
  NEKO_TUNE_LOG_MFMA_VEC(LX, time3);

  NEKO_TUNE_BEST(time2, best2, NEKO_EB_CANDIDATES);
  NEKO_TUNE_BEST(time3, best3, NEKO_MFMA_CANDIDATES);
  *eb_sel = best2;
  *nwf_sel = best3;

  /* The mfma variant joins the comparison only where it exists, its
     candidates are left at NEKO_TUNE_INIT otherwise */
  retval = (time3[best3] < time2[best2]) ? 3 : 2;

  HIP_CHECK(hipEventDestroy(start));
  HIP_CHECK(hipEventDestroy(stop));

  /* The tuner stands in for a real Ax evaluation, and the variants do not
     sum in the same order, so leave the output of the chosen kernel in
     au, av and aw */
  if (retval == 2) {
    CASE_VECTOR_KSTEP_SEL(LX, best2);
    sprintf(neko_log_buf, "Chose        : 2 (KSTEP, %d elem/block)",
            NEKO_EB_SEL(LX, best2));
  } else {
    CASE_VECTOR_MFMA_SEL(LX, best3);
    sprintf(neko_log_buf, "Chose        : 3 (MFMA %s, %d wf, %d elem/blk)",
            mfma_tile_name(best3), NEKO_MFMA_NWF(best3),
            NEKO_MFMA_EB(LX, best3));
  }
  log_message(neko_log_buf);
  log_end_section();
  return retval;
}

/* Padded variant of tune_vector(), for the orders whose kstep kernel stages
   with a padded stride. The matrix core variant is the same either way */
template < const int LX >
int tune_vector_padded(void *au, void *av, void *aw,
                       void *u, void *v, void *w,
                       void *dx, void *dy, void *dz, void *h1,
                       void *g11, void *g22, void *g33, void *g12,
                       void *g13, void *g23, int *nelv, int *lx, int *eb_sel,
                       int *nwf_sel) {
  hipEvent_t start, stop;
  float time2[NEKO_EB_CANDIDATES];
  int best2 = 0;
  float time3[NEKO_MFMA_CANDIDATES];
  int best3 = 0;
  const int rounds = neko_tune_rounds();
  const int iters = neko_tune_iters();
  const int sweep = neko_eb_sweep();
  /* Candidates of the kstep sweep, one -- the unblocked shape -- when the
     elements per block sweep is off */
  const int eb_cand = sweep ? NEKO_EB_CANDIDATES : 1;
  /* Candidates of the matrix core sweep, see neko_mfma_candidates() */
  const int mfma_cand = neko_mfma_candidates();
  /* Geometry pinned by each formulation's own variable, -1 to sweep it */
  const int eb_pin = neko_eb_pin();
  const int nwf_pin = neko_mfma_pin();
  /* Formulation pinned by NEKO_AUTOTUNE, as the identifier this returns */
  int strat = 0;
  /* Hardware capability, which the explicit NEKO_AUTOTUNE=MFMA pin needs */
  const bool mfma = mfma_lx_supported<LX>() && hip_have_mfma();
  int retval;

  for (int c = 0; c < NEKO_EB_CANDIDATES; c++) {
    time2[c] = NEKO_TUNE_INIT;
  }
  for (int c = 0; c < NEKO_MFMA_CANDIDATES; c++) {
    time3[c] = NEKO_TUNE_INIT;
  }

  const hipStream_t stream = (hipStream_t) glb_cmd_queue;

  char *env_value = NULL;
  char neko_log_buf[80];

  env_value = getenv("NEKO_AUTOTUNE");

  sprintf(neko_log_buf, "Autotune Ax vector (lx: %d)", *lx);
  log_section(neko_log_buf);

  *eb_sel = 0;
  *nwf_sel = 0;

  /* NEKO_AUTOTUNE names a formulation and nothing more, see tune() */
  if (env_value) {
    if ( !strcmp(env_value, "KSTEP") || !strcmp(env_value, "1D") ) {
      /* No 1d variant here, so either pin lands on kstep */
      strat = 2;
    } else if ( !strcmp(env_value, "MFMA") ) {
      if (mfma) {
        strat = 3;
      } else {
        sprintf(neko_log_buf, "MFMA strategy not available for this config");
        log_error(neko_log_buf);
      }
    } else {
      sprintf(neko_log_buf, "Invalid value set for NEKO_AUTOTUNE");
      log_error(neko_log_buf);
    }
  }

  /* Formulation and geometry both pinned leaves nothing to measure */
  const int pin = (strat == 2) ? eb_pin : (strat == 3) ? nwf_pin : -1;

  if (pin >= 0) {
    if (strat == 2) {
      *eb_sel = pin;
      CASE_VECTOR_KSTEP_PADDED_SEL(LX, pin);
      sprintf(neko_log_buf, "Set by env   : 2 (KSTEP, %d elem/block)",
              NEKO_EB_SEL(LX, pin));
    } else {
      *nwf_sel = pin;
      CASE_VECTOR_MFMA_SEL(LX, pin);
      sprintf(neko_log_buf, "Set by env   : 3 (MFMA %s, %d wf, %d elem/blk)",
              mfma_tile_name(pin), NEKO_MFMA_NWF(pin), NEKO_MFMA_EB(LX, pin));
    }
    log_message(neko_log_buf);
    log_end_section();
    return strat;
  }

  if (strat) {
    sprintf(neko_log_buf, "Set by env   : %d (%s)", strat, env_value);
    log_message(neko_log_buf);
  }

  /* Formulations the sweep considers, see NEKO_TUNE_FOR(). NEKO_MFMA_TUNE
     keeps the matrix core variant out of an unpinned sweep, but naming it
     explicitly still measures it */
  const bool try_kstep = (strat == 0 || strat == 2);
  const bool try_mfma = mfma && ((strat == 3) ||
                                 (strat == 0 && neko_mfma_sweep()));

  /* One candidate, see tune_vector() */
  if (!sweep && !try_mfma) {
    CASE_VECTOR_KSTEP_PADDED(LX, 0);
    sprintf(neko_log_buf, "Vector Ax    : %d elem/block (untuned)",
            NEKO_EB_SEL(LX, 0));
    log_message(neko_log_buf);
    log_end_section();
    return 2;
  }

  HIP_CHECK(hipEventCreate(&start));
  HIP_CHECK(hipEventCreate(&stop));

  /* Warm every variant before timing anything, see tune() */
  for (int i = 0; i < NEKO_TUNE_WARMUP; i++) {
    NEKO_TUNE_FOR(c, try_kstep, eb_pin, eb_cand) {
      CASE_VECTOR_KSTEP_PADDED_SEL(LX, c);
    }
    NEKO_TUNE_FOR(c, try_mfma, nwf_pin, mfma_cand) {
      CASE_VECTOR_MFMA_SEL(LX, c);
    }
  }

  /* Interleaved rounds, best time per variant, see tune() */
  for (int r = 0; r < rounds; r++) {
    NEKO_TUNE_FOR(c, try_kstep, eb_pin, eb_cand) {
      NEKO_TUNE_TIME(time2, CASE_VECTOR_KSTEP_PADDED_SEL, LX, c, iters);
    }
    NEKO_TUNE_FOR(c, try_mfma, nwf_pin, mfma_cand) {
      NEKO_TUNE_TIME(time3, CASE_VECTOR_MFMA_SEL, LX, c, iters);
    }
  }

  NEKO_TUNE_LOG_KSTEP(LX, time2);
  NEKO_TUNE_LOG_MFMA_VEC(LX, time3);

  NEKO_TUNE_BEST(time2, best2, NEKO_EB_CANDIDATES);
  NEKO_TUNE_BEST(time3, best3, NEKO_MFMA_CANDIDATES);
  *eb_sel = best2;
  *nwf_sel = best3;

  /* The mfma variant joins the comparison only where it exists, its
     candidates are left at NEKO_TUNE_INIT otherwise */
  retval = (time3[best3] < time2[best2]) ? 3 : 2;

  HIP_CHECK(hipEventDestroy(start));
  HIP_CHECK(hipEventDestroy(stop));

  /* Leave the output of the chosen kernel in au, av and aw, see
     tune_vector() */
  if (retval == 2) {
    CASE_VECTOR_KSTEP_PADDED_SEL(LX, best2);
    sprintf(neko_log_buf, "Chose        : 2 (KSTEP, %d elem/block)",
            NEKO_EB_SEL(LX, best2));
  } else {
    CASE_VECTOR_MFMA_SEL(LX, best3);
    sprintf(neko_log_buf, "Chose        : 3 (MFMA %s, %d wf, %d elem/blk)",
            mfma_tile_name(best3), NEKO_MFMA_NWF(best3),
            NEKO_MFMA_EB(LX, best3));
  }
  log_message(neko_log_buf);
  log_end_section();
  return retval;
}
