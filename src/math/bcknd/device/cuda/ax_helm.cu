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
#include "ax_helm_kernel.h"
#include "elem_block_tune.h"
#include <device/device_config.h>
#include <device/cuda/check.h>

extern "C" {
  #include <common/neko_log.h>
}

template < const int>
int tune(void *w, void *u, void *dx, void *dy, void *dz,
         void *dxt, void *dyt, void *dzt, void *h1,
         void *g11, void *g22, void *g33, void *g12,
         void *g13, void *g23, int *nelv, int *lx,
         int *eb_sel, int *ch_sel, int *nw_sel, int *tw_sel);

template < const int>
int tune_padded(void *w, void *u, void *dx, void *dy, void *dz,
                void *dxt, void *dyt, void *dzt, void *h1,
                void *g11, void *g22, void *g33, void *g12,
                void *g13, void *g23, int *nelv, int *lx,
                int *eb_sel, int *ch_sel, int *nw_sel, int *tw_sel);

template < const int>
int tune_vector(void *au, void *av, void *aw, void *u, void *v, void *w,
                void *dx, void *dy, void *dz, void *h1,
                void *g11, void *g22, void *g33, void *g12,
                void *g13, void *g23, int *nelv, int *lx,
                int *eb_sel, int *nw_sel, int *tw_sel,
                int *bw_sel);

template < const int>
int tune_vector_padded(void *au, void *av, void *aw, void *u, void *v, void *w,
                       void *dx, void *dy, void *dz, void *h1,
                       void *g11, void *g22, void *g33, void *g12,
                       void *g13, void *g23, int *nelv, int *lx,
                int *eb_sel, int *nw_sel, int *tw_sel,
                int *bw_sel);

extern "C" {

  /**
   * Fortran wrapper for device CUDA Ax
   */
  void cuda_ax_helm(void *w, void *u, void *dx, void *dy, void *dz,
                    void *dxt, void *dyt, void *dzt, void *h1,
                    void *g11, void *g22, void *g33, void *g12,
                    void *g13, void *g23, int *nelv, int *lx) {

    static int autotune[17] = { 0 };
    /* Elements per block candidate chosen by the autotuner, see
       elem_block<> in elem_block.h */
    static int autotune_eb[17] = { 0 };
    /* chunk candidate chosen for the 1d variant */
    static int autotune_ch[17] = { 0 };
    /* warps per block candidate chosen for the dmma variant */
    static int autotune_nw[17] = { 0 };
    /* warps per block candidate chosen for the tma staged dmma variant */
    static int autotune_tw[17] = { 0 };

    const dim3 nblcks_1d((*nelv), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define CASE_1D(LX, C)                                                          \
    ax_helm_kernel_1d<real, LX, NEKO_CHUNKS(LX, C)>                             \
      <<<nblcks_1d, NEKO_CHUNKS_NTHRDS(LX, C), 0, stream>>>                     \
                         ((real *) w, (real *) u,                               \
                          (real *) dx, (real *) dy, (real *) dz,                \
                          (real *) dxt, (real *) dyt, (real *) dzt, (real *) h1,\
                          (real *) g11, (real *) g22, (real *) g33,             \
                          (real *) g12, (real *) g13, (real *) g23);            \
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
    ax_helm_kernel_kstep<real, LX, NEKO_EB(LX, C)>                              \
      <<<NEKO_EB_NBLCKS(*nelv, LX, C), NEKO_EB_NTHRDS(LX, C), 0, stream>>>      \
                          ((real *) w, (real *) u,                              \
                           (real *) dx, (real *) dy, (real *) dz, (real *) h1,  \
                           (real *) g11, (real *) g22, (real *) g33,            \
                           (real *) g12, (real *) g13, (real *) g23, *nelv);    \
      CUDA_CHECK(cudaGetLastError());

#define CASE_KSTEP_PADDED(LX, C)                                                \
    ax_helm_kernel_kstep_padded<real, LX, NEKO_EB(LX, C)>                       \
      <<<NEKO_EB_NBLCKS(*nelv, LX, C), NEKO_EB_NTHRDS(LX, C), 0, stream>>>      \
                          ((real *) w, (real *) u,                              \
                           (real *) dx, (real *) dy, (real *) dz, (real *) h1,  \
                           (real *) g11, (real *) g22, (real *) g33,            \
                           (real *) g12, (real *) g13, (real *) g23, *nelv);    \
      CUDA_CHECK(cudaGetLastError());

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

#define CASE_DMMA(LX, C)                                                        \
    ax_helm_kernel_dmma<real, LX, NEKO_DMMA_NW(C)>                              \
      <<<NEKO_DMMA_NBLCKS(*nelv, LX), NEKO_DMMA_NTHRDS(C), 0, stream>>>         \
                          ((real *) w, (real *) u,                              \
                           (real *) dx, (real *) dy, (real *) dz, (real *) h1,  \
                           (real *) g11, (real *) g22, (real *) g33,            \
                           (real *) g12, (real *) g13, (real *) g23, *nelv);    \
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
 * so the grid is nelv blocks -- but it takes no nelv: there is no packed tail
 * to clamp when a block is exactly an element. See dmma_tma_kernel.h.
 */
#define CASE_DMMA_TMA(LX, C)                                                    \
    ax_helm_kernel_dmma_tma<real, LX, NEKO_DMMA_NW(C)>                          \
      <<<NEKO_DMMA_NBLCKS(*nelv, LX), NEKO_DMMA_NTHRDS(C), 0, stream>>>         \
                          ((real *) w, (real *) u,                              \
                           (real *) dx, (real *) dy, (real *) dz, (real *) h1,  \
                           (real *) g11, (real *) g22, (real *) g33,            \
                           (real *) g12, (real *) g13, (real *) g23);           \
      CUDA_CHECK(cudaGetLastError());

/* Runtime dispatch onto the tuned warps per block candidate */
#define CASE_DMMA_TMA_SEL(LX, SEL)                                              \
    switch (SEL) {                                                              \
    case 0:  CASE_DMMA_TMA(LX, 0); break;                                       \
    case 1:  CASE_DMMA_TMA(LX, 1); break;                                       \
    default: CASE_DMMA_TMA(LX, 2); break;                                       \
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
                               &autotune_nw[LX], &autotune_tw[LX]);             \
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

#define CASE_PADDED(LX)                                                         \
    case LX:                                                                    \
      if(autotune[LX] == 0 ) {                                                  \
        autotune[LX]=tune_padded<LX>(w,  u,                                     \
                                     dx,  dy, dz,                               \
                                     dxt, dyt, dzt,h1,                          \
                                     g11, g22, g33,                             \
                                     g12, g13, g23,nelv,lx,                     \
                                     &autotune_eb[LX], &autotune_ch[LX],        \
                                     &autotune_nw[LX], &autotune_tw[LX]);       \
      } else if (autotune[LX] == 1 ) {                                          \
        CASE_1D_SEL(LX, autotune_ch[LX]);                                       \
      } else if (autotune[LX] == 2 ) {                                          \
        CASE_KSTEP_PADDED_SEL(LX, autotune_eb[LX]);                             \
      } else if (autotune[LX] == 3 ) {                                          \
        CASE_DMMA_SEL(LX, autotune_nw[LX]);                                     \
      } else if (autotune[LX] == 4 ) {                                          \
        CASE_DMMA_TMA_SEL(LX, autotune_tw[LX]);                                 \
      }                                                                         \
      break

/*
 * High order cases have no 1d variant to compare against (its shared memory
 * footprint grows as LX^3), so they are not tuned and keep one element per
 * block, i.e. candidate 0
 */
#define CASE_LARGE(LX)                                                          \
    case LX:                                                                    \
      CASE_KSTEP(LX, 0);                                                        \
      break

#define CASE_LARGE_PADDED(LX)                                                   \
    case LX:                                                                    \
      CASE_KSTEP_PADDED(LX, 0);                                                 \
      break

    if ((*lx) < 12) {
      switch(*lx) {
        CASE(2);
        CASE(3);
        CASE_PADDED(4);
        CASE(5);
        CASE(6);
        CASE(7);
        CASE_PADDED(8);
        CASE(9);
        CASE(10);
        CASE(11);
      default:
        {
          fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
          exit(1);
        }
      }
    }
    else {
      switch(*lx) {
        CASE_LARGE(12);
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
   * Fortran wrapper for device CUDA Ax vector version
   */
  void cuda_ax_helm_vector(void *au, void *av, void *aw,
                           void *u, void *v, void *w,
                           void *dx, void *dy, void *dz,
                           void *dxt, void *dyt, void *dzt,
                           void *h1, void *g11, void *g22,
                           void *g33, void *g12, void *g13,
                           void *g23, int *nelv, int *lx) {

    /* Strategy chosen by the autotuner: 2 kstep, 3 dmma */
    static int autotune_vec[17] = { 0 };
    /* elements per block candidate chosen for the kstep variant */
    static int autotune_vec_eb[17] = { 0 };
    /* warps per block candidate chosen for the dmma variant */
    static int autotune_vec_nw[17] = { 0 };
    /* warps per block candidate chosen for the tma staged dmma variant */
    static int autotune_vec_tw[17] = { 0 };
    /* warps per block candidate chosen for the batched tma variant */
    static int autotune_vec_bw[17] = { 0 };

    const dim3 nblcks_dmma((*nelv), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

/*
 * The vector kstep variant sweeps the same elements per block candidates as
 * the scalar one, and the tuner picks between it and the dmma variant. The
 * vector kernels hold three times the register blocked state -- six T[LX]
 * arrays rather than two -- and measure at 254-255 registers on sm_90 for
 * every lx from 8 up, spilling at 10, 12, 14 and 16, so blocking is not
 * expected to buy much here: it does not change registers per thread, only
 * threads per block, so it saves the derivative matrix loads and nothing
 * else. It is swept rather than assumed because that is a measurement, and
 * because elem_block<>'s thread clamp already keeps every candidate inside
 * the shared memory budget at every lx.
 */
#define CASE_VECTOR_KSTEP(LX, C)                                               \
    ax_helm_kernel_vector_kstep<real, LX, NEKO_EB(LX, C)>                      \
    <<<NEKO_EB_NBLCKS(*nelv, LX, C), NEKO_EB_NTHRDS(LX, C), 0, stream>>>       \
                                    ((real *) au, (real *) av, (real *) aw,    \
                                     (real *) u, (real *) v, (real *) w,       \
                                     (real *) dx, (real *) dy, (real *) dz,    \
                                     (real *) h1, (real *) g11, (real *) g22,  \
                                     (real *) g33, (real *) g12, (real *) g13, \
                                     (real *) g23, *nelv);                     \
    CUDA_CHECK(cudaGetLastError());

#define CASE_VECTOR_KSTEP_PADDED(LX, C)                                        \
    ax_helm_kernel_vector_kstep_padded<real, LX, NEKO_EB(LX, C)>               \
    <<<NEKO_EB_NBLCKS(*nelv, LX, C), NEKO_EB_NTHRDS(LX, C), 0, stream>>>       \
                                    ((real *) au, (real *) av, (real *) aw,    \
                                     (real *) u, (real *) v, (real *) w,       \
                                     (real *) dx, (real *) dy, (real *) dz,    \
                                     (real *) h1, (real *) g11, (real *) g22,  \
                                     (real *) g33, (real *) g12, (real *) g13, \
                                     (real *) g23, *nelv);                     \
    CUDA_CHECK(cudaGetLastError());

/* Runtime dispatch onto the tuned elements per block candidate */
#define CASE_VECTOR_KSTEP_SEL(LX, SEL)                                         \
    switch (SEL) {                                                             \
    case 0:  CASE_VECTOR_KSTEP(LX, 0); break;                                  \
    case 1:  CASE_VECTOR_KSTEP(LX, 1); break;                                  \
    default: CASE_VECTOR_KSTEP(LX, 2); break;                                  \
    }

#define CASE_VECTOR_KSTEP_PADDED_SEL(LX, SEL)                                  \
    switch (SEL) {                                                             \
    case 0:  CASE_VECTOR_KSTEP_PADDED(LX, 0); break;                           \
    case 1:  CASE_VECTOR_KSTEP_PADDED(LX, 1); break;                           \
    default: CASE_VECTOR_KSTEP_PADDED(LX, 2); break;                           \
    }

#define CASE_VECTOR_DMMA(LX, C)                                                \
    ax_helm_kernel_dmma_vector<real, LX, NEKO_DMMA_NW(C)>                      \
    <<<nblcks_dmma, NEKO_DMMA_NTHRDS(C), 0, stream>>>                          \
                                    ((real *) au, (real *) av, (real *) aw,    \
                                     (real *) u, (real *) v, (real *) w,       \
                                     (real *) dx, (real *) dy, (real *) dz,    \
                                     (real *) h1, (real *) g11, (real *) g22,  \
                                     (real *) g33, (real *) g12, (real *) g13, \
                                     (real *) g23);                            \
    CUDA_CHECK(cudaGetLastError());

/* Runtime dispatch onto the tuned warps per block candidate */
#define CASE_VECTOR_DMMA_SEL(LX, SEL)                                          \
    switch (SEL) {                                                             \
    case 0:  CASE_VECTOR_DMMA(LX, 0); break;                                   \
    case 1:  CASE_VECTOR_DMMA(LX, 1); break;                                   \
    default: CASE_VECTOR_DMMA(LX, 2); break;                                   \
    }

/* The TMA staged vector variant, same grid and candidates as
   CASE_VECTOR_DMMA; see dmma_tma_kernel.h */
#define CASE_VECTOR_DMMA_TMA(LX, C)                                            \
    ax_helm_kernel_dmma_tma_vector<real, LX, NEKO_DMMA_NW(C)>                  \
    <<<nblcks_dmma, NEKO_DMMA_NTHRDS(C), 0, stream>>>                          \
                                    ((real *) au, (real *) av, (real *) aw,    \
                                     (real *) u, (real *) v, (real *) w,       \
                                     (real *) dx, (real *) dy, (real *) dz,    \
                                     (real *) h1, (real *) g11, (real *) g22,  \
                                     (real *) g33, (real *) g12, (real *) g13, \
                                     (real *) g23);                            \
    CUDA_CHECK(cudaGetLastError());

/* Runtime dispatch onto the tuned warps per block candidate */
#define CASE_VECTOR_DMMA_TMA_SEL(LX, SEL)                                      \
    switch (SEL) {                                                             \
    case 0:  CASE_VECTOR_DMMA_TMA(LX, 0); break;                               \
    case 1:  CASE_VECTOR_DMMA_TMA(LX, 1); break;                               \
    default: CASE_VECTOR_DMMA_TMA(LX, 2); break;                               \
    }

/*
 * The batched TMA variant. Unlike every other launch here it carries a dynamic
 * shared memory allocation -- 54800 B, past the 48 kB a block gets for free --
 * which the kernel has to opt into once before its first launch; see
 * ax_helm_dmma_tma_batch_optin(), which is a predictable branch after that.
 */
#define CASE_VECTOR_DMMA_TMA_BATCH(LX, C)                                      \
    (void) ax_helm_dmma_tma_batch_optin<real, LX, NEKO_DMMA_NW(C)>();          \
    ax_helm_kernel_dmma_tma_batch<real, LX, NEKO_DMMA_NW(C)>                   \
    <<<nblcks_dmma, NEKO_DMMA_NTHRDS(C),                                       \
       NEKO_DMMA_TMA_BATCH_SMEM, stream>>>                                     \
                                    ((real *) au, (real *) av, (real *) aw,    \
                                     (real *) u, (real *) v, (real *) w,       \
                                     (real *) dx, (real *) dy, (real *) dz,    \
                                     (real *) h1, (real *) g11, (real *) g22,  \
                                     (real *) g33, (real *) g12, (real *) g13, \
                                     (real *) g23);                            \
    CUDA_CHECK(cudaGetLastError());

/* Runtime dispatch onto the tuned warps per block candidate */
#define CASE_VECTOR_DMMA_TMA_BATCH_SEL(LX, SEL)                                \
    switch (SEL) {                                                             \
    case 0:  CASE_VECTOR_DMMA_TMA_BATCH(LX, 0); break;                         \
    case 1:  CASE_VECTOR_DMMA_TMA_BATCH(LX, 1); break;                         \
    default: CASE_VECTOR_DMMA_TMA_BATCH(LX, 2); break;                         \
    }

#define CASE_VECTOR(LX)                                                        \
    case LX:                                                                   \
      if (autotune_vec[LX] == 0 ) {                                            \
        autotune_vec[LX] = tune_vector<LX>(au, av, aw, u, v, w,                \
                                           dx, dy, dz, h1,                     \
                                           g11, g22, g33, g12, g13, g23,       \
                                           nelv, lx, &autotune_vec_eb[LX],     \
                                           &autotune_vec_nw[LX],               \
                                           &autotune_vec_tw[LX],               \
                                           &autotune_vec_bw[LX]);              \
      } else if (autotune_vec[LX] == 2 ) {                                     \
        CASE_VECTOR_KSTEP_SEL(LX, autotune_vec_eb[LX]);                        \
      } else if (autotune_vec[LX] == 3 ) {                                     \
        CASE_VECTOR_DMMA_SEL(LX, autotune_vec_nw[LX]);                         \
      } else if (autotune_vec[LX] == 4 ) {                                     \
        CASE_VECTOR_DMMA_TMA_SEL(LX, autotune_vec_tw[LX]);                     \
      } else if (autotune_vec[LX] == 5 ) {                                     \
        CASE_VECTOR_DMMA_TMA_BATCH_SEL(LX, autotune_vec_bw[LX]);               \
      }                                                                        \
       break

#define CASE_VECTOR_PADDED(LX)                                                 \
    case LX:                                                                   \
      if (autotune_vec[LX] == 0 ) {                                            \
        autotune_vec[LX] = tune_vector_padded<LX>(au, av, aw, u, v, w,         \
                                                  dx, dy, dz, h1,              \
                                                  g11, g22, g33,               \
                                                  g12, g13, g23,               \
                                                  nelv, lx,                    \
                                                  &autotune_vec_eb[LX],        \
                                                  &autotune_vec_nw[LX],        \
                                                  &autotune_vec_tw[LX],        \
                                                  &autotune_vec_bw[LX]);       \
      } else if (autotune_vec[LX] == 2 ) {                                     \
        CASE_VECTOR_KSTEP_PADDED_SEL(LX, autotune_vec_eb[LX]);                 \
      } else if (autotune_vec[LX] == 3 ) {                                     \
        CASE_VECTOR_DMMA_SEL(LX, autotune_vec_nw[LX]);                         \
      } else if (autotune_vec[LX] == 4 ) {                                     \
        CASE_VECTOR_DMMA_TMA_SEL(LX, autotune_vec_tw[LX]);                     \
      } else if (autotune_vec[LX] == 5 ) {                                     \
        CASE_VECTOR_DMMA_TMA_BATCH_SEL(LX, autotune_vec_bw[LX]);               \
      }                                                                        \
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
   * Fortran wrapper for device CUDA Ax vector version part2
   */
  void cuda_ax_helm_vector_part2(void *au, void *av, void *aw,
                                 void *u, void *v, void *w,
                                 void *h2, void *B, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

    ax_helm_kernel_vector_part2<real>
      <<<nblcks, nthrds, 0, stream>>> ((real *) au, (real *) av, (real *) aw,
                                       (real *) u, (real *) v, (real *) w,
                                       (real *) h2, (real *) B, *n);
  }
}

template < const int LX >
int tune(void *w, void *u, void *dx, void *dy, void *dz,
         void *dxt, void *dyt, void *dzt, void *h1,
         void *g11, void *g22, void *g33, void *g12,
         void *g13, void *g23, int *nelv, int *lx,
         int *eb_sel, int *ch_sel, int *nw_sel, int *tw_sel) {
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
     dmma_tma_aligned() in dmma_tma_kernel.h */
  const bool tma = dmma_tma_lx_supported<LX>() && cuda_have_tma() &&
                   dmma_tma_aligned(w, u, h1, g11, g22, g33, g12, g13, g23);
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

  const dim3 nblcks_1d((*nelv), 1, 1);
  const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

  char *env_value = NULL;
  char neko_log_buf[80];

  env_value=getenv("NEKO_AUTOTUNE");

  sprintf(neko_log_buf, "Autotune Ax helm (lx: %d)", *lx);
  log_section(neko_log_buf);

  *eb_sel = 0;
  *ch_sel = 0;
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
      const int c = neko_eb_env();
      *eb_sel = c;
      CASE_KSTEP_SEL(LX, c);
      sprintf(neko_log_buf,"Set by env   : 2 (KSTEP, %d elem/block)",
              NEKO_EB_SEL(LX, c));
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

  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  /* The tuner stands in for a real Ax evaluation, and the variants do not
     sum in the same order, so leave the output of the chosen kernel in w */
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

template < const int LX >
int tune_padded(void *w, void *u, void *dx, void *dy, void *dz,
                void *dxt, void *dyt, void *dzt, void *h1,
                void *g11, void *g22, void *g33, void *g12,
                void *g13, void *g23, int *nelv, int *lx,
                int *eb_sel, int *ch_sel, int *nw_sel, int *tw_sel) {
  cudaEvent_t start, stop;
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
     dmma_tma_aligned() in dmma_tma_kernel.h */
  const bool tma = dmma_tma_lx_supported<LX>() && cuda_have_tma() &&
                   dmma_tma_aligned(w, u, h1, g11, g22, g33, g12, g13, g23);
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

  const dim3 nblcks_1d((*nelv), 1, 1);
  const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

  char *env_value = NULL;
  char neko_log_buf[80];

  env_value=getenv("NEKO_AUTOTUNE");

  sprintf(neko_log_buf, "Autotune Ax helm (lx: %d)", *lx);
  log_section(neko_log_buf);

  *eb_sel = 0;
  *ch_sel = 0;
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
      const int c = neko_eb_env();
      *eb_sel = c;
      CASE_KSTEP_PADDED_SEL(LX, c);
      sprintf(neko_log_buf,"Set by env   : 2 (KSTEP, %d elem/block)",
              NEKO_EB_SEL(LX, c));
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

  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  /* Leave the output of the chosen kernel in w, see tune() */
  if (retval == 1) {
    CASE_1D_SEL(LX, best1);
  } else if (retval == 2) {
    CASE_KSTEP_PADDED_SEL(LX, best);
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


/*
 * Autotuner for the vector Ax
 *
 * The vector operator has no 1d variant, so the sweep is the elements per
 * block candidates of the kstep variant against the warps per block
 * candidates of the dmma variant -- the same two dimensions the scalar tuner
 * sweeps, and the same shape as the HIP vector tuner. See the note at the top
 * of ax_helm_kernel.h for why blocking is not expected to pay here, and why
 * it is measured rather than assumed. Sampling follows the scalar tuner:
 * interleaved rounds, best round per candidate.
 *
 * @note The section title has to fit the 30 character header field in
 * log_section(), which truncates rather than wraps -- hence "Ax vector"
 * rather than the "Ax helm vector" it reads as.
 */
template < const int LX >
int tune_vector(void *au, void *av, void *aw, void *u, void *v, void *w,
                void *dx, void *dy, void *dz, void *h1,
                void *g11, void *g22, void *g33, void *g12,
                void *g13, void *g23, int *nelv, int *lx,
                int *eb_sel, int *nw_sel, int *tw_sel,
                int *bw_sel) {
  cudaEvent_t start, stop;
  float time2[NEKO_EB_CANDIDATES];
  int best2 = 0;
  float time3[NEKO_DMMA_CANDIDATES];
  int best3 = 0;
  float time4[NEKO_DMMA_CANDIDATES];
  int best4 = 0;
  float time5[NEKO_DMMA_CANDIDATES];
  int best5 = 0;
  const int rounds = neko_tune_rounds();
  const int iters = neko_tune_iters();
  const int sweep = neko_eb_sweep();
  const bool dmma = dmma_vector_lx_supported<LX>() && cuda_have_dmma();
  /* Whether the pointers are bulk copy aligned is a property of this call
     rather than of the kernel, see dmma_tma_vector_aligned() */
  const bool tma = dmma_tma_vector_lx_supported<LX>() && cuda_have_tma() &&
                   dmma_tma_vector_aligned(au, av, aw, u, v, w, h1,
                                           g11, g22, g33, g12, g13, g23);
  /* Same pointers, and additionally the device has to grant the block the
     dynamic allocation, see cuda_have_tma_batch() */
  const bool batch = dmma_tma_batch_lx_supported<LX>() && cuda_have_tma_batch()
                     && dmma_tma_vector_aligned(au, av, aw, u, v, w, h1,
                                                g11, g22, g33,
                                                g12, g13, g23);
  int retval;

  for (int c = 0; c < NEKO_EB_CANDIDATES; c++) {
    time2[c] = NEKO_TUNE_INIT;
  }
  for (int c = 0; c < NEKO_DMMA_CANDIDATES; c++) {
    time3[c] = NEKO_TUNE_INIT;
    time4[c] = NEKO_TUNE_INIT;
    time5[c] = NEKO_TUNE_INIT;
  }

  const dim3 nblcks_dmma((*nelv), 1, 1);
  const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

  char *env_value = NULL;
  char neko_log_buf[80];

  env_value=getenv("NEKO_AUTOTUNE");

  sprintf(neko_log_buf, "Autotune Ax vector (lx: %d)", *lx);
  log_section(neko_log_buf);

  *nw_sel = 0;
  *tw_sel = 0;
  *bw_sel = 0;

  if(env_value) {
    if( !strcmp(env_value,"KSTEP") || !strcmp(env_value,"1D") ) {
      /* No 1d variant here, so either pin lands on kstep */
      const int c = neko_eb_env();
      *eb_sel = c;
      CASE_VECTOR_KSTEP_SEL(LX, c);
      sprintf(neko_log_buf,"Set by env   : 2 (KSTEP, %d elem/block)",
              NEKO_EB_SEL(LX, c));
      log_message(neko_log_buf);
      log_end_section();
      return 2;
    } else if( !strcmp(env_value,"DMMA") ) {
      if (dmma) {
        const int c = neko_dmma_env();
        *nw_sel = c;
        CASE_VECTOR_DMMA_SEL(LX, c);
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
        CASE_VECTOR_DMMA_TMA_SEL(LX, c);
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
    } else if( !strcmp(env_value,"DMMA_TMA_BATCH") ) {
      if (batch) {
        const int c = neko_dmma_tma_env();
        *bw_sel = c;
        CASE_VECTOR_DMMA_TMA_BATCH_SEL(LX, c);
        sprintf(neko_log_buf,"Set by env   : 5 (DMMA_TMA_BATCH, %d warps)",
                NEKO_DMMA_NW(c));
        log_message(neko_log_buf);
        log_end_section();
        return 5;
      } else {
        sprintf(neko_log_buf,
                "DMMA_TMA_BATCH strategy not available for this config");
        log_error(neko_log_buf);
      }
    } else {
      sprintf(neko_log_buf, "Invalid value set for NEKO_AUTOTUNE");
      log_error(neko_log_buf);
    }
  }

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  /* Warm every variant before timing anything, see tune() */
  for (int i = 0; i < NEKO_TUNE_WARMUP; i++) {
    CASE_VECTOR_KSTEP(LX, 0);
    if (sweep) {
      CASE_VECTOR_KSTEP(LX, 1);
      CASE_VECTOR_KSTEP(LX, 2);
    }
    if (dmma) {
      CASE_VECTOR_DMMA(LX, 0);
      CASE_VECTOR_DMMA(LX, 1);
      CASE_VECTOR_DMMA(LX, 2);
    }
    if (tma) {
      CASE_VECTOR_DMMA_TMA(LX, 0);
      CASE_VECTOR_DMMA_TMA(LX, 1);
      CASE_VECTOR_DMMA_TMA(LX, 2);
    }
    if (batch) {
      CASE_VECTOR_DMMA_TMA_BATCH(LX, 0);
      CASE_VECTOR_DMMA_TMA_BATCH(LX, 1);
      CASE_VECTOR_DMMA_TMA_BATCH(LX, 2);
    }
  }

  for (int r = 0; r < rounds; r++) {
    NEKO_TUNE_TIME(time2, CASE_VECTOR_KSTEP, LX, 0, iters);
    if (sweep) {
      NEKO_TUNE_TIME(time2, CASE_VECTOR_KSTEP, LX, 1, iters);
      NEKO_TUNE_TIME(time2, CASE_VECTOR_KSTEP, LX, 2, iters);
    }
    if (dmma) {
      NEKO_TUNE_TIME(time3, CASE_VECTOR_DMMA, LX, 0, iters);
      NEKO_TUNE_TIME(time3, CASE_VECTOR_DMMA, LX, 1, iters);
      NEKO_TUNE_TIME(time3, CASE_VECTOR_DMMA, LX, 2, iters);
    }
    if (tma) {
      NEKO_TUNE_TIME(time4, CASE_VECTOR_DMMA_TMA, LX, 0, iters);
      NEKO_TUNE_TIME(time4, CASE_VECTOR_DMMA_TMA, LX, 1, iters);
      NEKO_TUNE_TIME(time4, CASE_VECTOR_DMMA_TMA, LX, 2, iters);
    }
    if (batch) {
      NEKO_TUNE_TIME(time5, CASE_VECTOR_DMMA_TMA_BATCH, LX, 0, iters);
      NEKO_TUNE_TIME(time5, CASE_VECTOR_DMMA_TMA_BATCH, LX, 1, iters);
      NEKO_TUNE_TIME(time5, CASE_VECTOR_DMMA_TMA_BATCH, LX, 2, iters);
    }
  }

  for (int c = 0; c < NEKO_EB_CANDIDATES; c++) {
    if (time2[c] >= NEKO_TUNE_INIT) { continue; }
    sprintf(neko_log_buf, "KSTEP eb=%-4d: %9.2f us/call",
            NEKO_EB_SEL(LX, c), NEKO_TUNE_US(time2[c], iters));
    log_message(neko_log_buf);
  }
  NEKO_TUNE_LOG_DMMA(LX, time3);
  NEKO_TUNE_LOG_DMMA_TMA(LX, time4);
  NEKO_TUNE_LOG_DMMA_TMA_BATCH(LX, time5);

  NEKO_TUNE_BEST(time2, best2, NEKO_EB_CANDIDATES);
  NEKO_TUNE_BEST(time3, best3, NEKO_DMMA_CANDIDATES);
  NEKO_TUNE_BEST(time4, best4, NEKO_DMMA_CANDIDATES);
  NEKO_TUNE_BEST(time5, best5, NEKO_DMMA_CANDIDATES);
  *eb_sel = best2;
  *nw_sel = best3;
  *tw_sel = best4;
  *bw_sel = best5;

  /* The dmma variants join the comparison only where they exist, their
     candidates are left at NEKO_TUNE_INIT otherwise */
  float tbest = time2[best2];

  retval = 2;
  if (time3[best3] < tbest) { retval = 3; tbest = time3[best3]; }
  if (time4[best4] < tbest) { retval = 4; tbest = time4[best4]; }
  if (time5[best5] < tbest) { retval = 5; }

  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  /* The tuner stands in for a real Ax evaluation, and the variants do not
     sum in the same order, so leave the output of the chosen kernel in
     au, av and aw */
  if (retval == 2) {
    CASE_VECTOR_KSTEP_SEL(LX, best2);
    sprintf(neko_log_buf, "Chose        : 2 (KSTEP, %d elem/block)",
            NEKO_EB_SEL(LX, best2));
  } else if (retval == 3) {
    CASE_VECTOR_DMMA_SEL(LX, best3);
    sprintf(neko_log_buf, "Chose        : 3 (DMMA, %d warps, %d elem/blk)",
            NEKO_DMMA_NW(best3), NEKO_DMMA_PACK(LX));
  } else if (retval == 4) {
    CASE_VECTOR_DMMA_TMA_SEL(LX, best4);
    sprintf(neko_log_buf, "Chose        : 4 (DMMA_TMA, %d warps)",
            NEKO_DMMA_NW(best4));
  } else {
    CASE_VECTOR_DMMA_TMA_BATCH_SEL(LX, best5);
    sprintf(neko_log_buf, "Chose        : 5 (DMMA_TMA_BATCH, %d warps)",
            NEKO_DMMA_NW(best5));
  }
  log_message(neko_log_buf);
  log_end_section();
  return retval;
}

/* Padded variant of tune_vector(), see the note on bank conflicts above */
template < const int LX >
int tune_vector_padded(void *au, void *av, void *aw, void *u, void *v, void *w,
                       void *dx, void *dy, void *dz, void *h1,
                       void *g11, void *g22, void *g33, void *g12,
                       void *g13, void *g23, int *nelv, int *lx,
                int *eb_sel, int *nw_sel, int *tw_sel,
                int *bw_sel) {
  cudaEvent_t start, stop;
  float time2[NEKO_EB_CANDIDATES];
  int best2 = 0;
  float time3[NEKO_DMMA_CANDIDATES];
  int best3 = 0;
  float time4[NEKO_DMMA_CANDIDATES];
  int best4 = 0;
  float time5[NEKO_DMMA_CANDIDATES];
  int best5 = 0;
  const int rounds = neko_tune_rounds();
  const int iters = neko_tune_iters();
  const int sweep = neko_eb_sweep();
  const bool dmma = dmma_vector_lx_supported<LX>() && cuda_have_dmma();
  /* Whether the pointers are bulk copy aligned is a property of this call
     rather than of the kernel, see dmma_tma_vector_aligned() */
  const bool tma = dmma_tma_vector_lx_supported<LX>() && cuda_have_tma() &&
                   dmma_tma_vector_aligned(au, av, aw, u, v, w, h1,
                                           g11, g22, g33, g12, g13, g23);
  /* Same pointers, and additionally the device has to grant the block the
     dynamic allocation, see cuda_have_tma_batch() */
  const bool batch = dmma_tma_batch_lx_supported<LX>() && cuda_have_tma_batch()
                     && dmma_tma_vector_aligned(au, av, aw, u, v, w, h1,
                                                g11, g22, g33,
                                                g12, g13, g23);
  int retval;

  for (int c = 0; c < NEKO_EB_CANDIDATES; c++) {
    time2[c] = NEKO_TUNE_INIT;
  }
  for (int c = 0; c < NEKO_DMMA_CANDIDATES; c++) {
    time3[c] = NEKO_TUNE_INIT;
    time4[c] = NEKO_TUNE_INIT;
    time5[c] = NEKO_TUNE_INIT;
  }

  const dim3 nblcks_dmma((*nelv), 1, 1);
  const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

  char *env_value = NULL;
  char neko_log_buf[80];

  env_value=getenv("NEKO_AUTOTUNE");

  sprintf(neko_log_buf, "Autotune Ax vector (lx: %d)", *lx);
  log_section(neko_log_buf);

  *nw_sel = 0;
  *tw_sel = 0;
  *bw_sel = 0;

  if(env_value) {
    if( !strcmp(env_value,"KSTEP") || !strcmp(env_value,"1D") ) {
      /* No 1d variant here, so either pin lands on kstep */
      const int c = neko_eb_env();
      *eb_sel = c;
      CASE_VECTOR_KSTEP_PADDED_SEL(LX, c);
      sprintf(neko_log_buf,"Set by env   : 2 (KSTEP, %d elem/block)",
              NEKO_EB_SEL(LX, c));
      log_message(neko_log_buf);
      log_end_section();
      return 2;
    } else if( !strcmp(env_value,"DMMA") ) {
      if (dmma) {
        const int c = neko_dmma_env();
        *nw_sel = c;
        CASE_VECTOR_DMMA_SEL(LX, c);
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
        CASE_VECTOR_DMMA_TMA_SEL(LX, c);
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
    } else if( !strcmp(env_value,"DMMA_TMA_BATCH") ) {
      if (batch) {
        const int c = neko_dmma_tma_env();
        *bw_sel = c;
        CASE_VECTOR_DMMA_TMA_BATCH_SEL(LX, c);
        sprintf(neko_log_buf,"Set by env   : 5 (DMMA_TMA_BATCH, %d warps)",
                NEKO_DMMA_NW(c));
        log_message(neko_log_buf);
        log_end_section();
        return 5;
      } else {
        sprintf(neko_log_buf,
                "DMMA_TMA_BATCH strategy not available for this config");
        log_error(neko_log_buf);
      }
    } else {
      sprintf(neko_log_buf, "Invalid value set for NEKO_AUTOTUNE");
      log_error(neko_log_buf);
    }
  }

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  /* Warm every variant before timing anything, see tune() */
  for (int i = 0; i < NEKO_TUNE_WARMUP; i++) {
    CASE_VECTOR_KSTEP_PADDED(LX, 0);
    if (sweep) {
      CASE_VECTOR_KSTEP_PADDED(LX, 1);
      CASE_VECTOR_KSTEP_PADDED(LX, 2);
    }
    if (dmma) {
      CASE_VECTOR_DMMA(LX, 0);
      CASE_VECTOR_DMMA(LX, 1);
      CASE_VECTOR_DMMA(LX, 2);
    }
    if (tma) {
      CASE_VECTOR_DMMA_TMA(LX, 0);
      CASE_VECTOR_DMMA_TMA(LX, 1);
      CASE_VECTOR_DMMA_TMA(LX, 2);
    }
    if (batch) {
      CASE_VECTOR_DMMA_TMA_BATCH(LX, 0);
      CASE_VECTOR_DMMA_TMA_BATCH(LX, 1);
      CASE_VECTOR_DMMA_TMA_BATCH(LX, 2);
    }
  }

  for (int r = 0; r < rounds; r++) {
    NEKO_TUNE_TIME(time2, CASE_VECTOR_KSTEP_PADDED, LX, 0, iters);
    if (sweep) {
      NEKO_TUNE_TIME(time2, CASE_VECTOR_KSTEP_PADDED, LX, 1, iters);
      NEKO_TUNE_TIME(time2, CASE_VECTOR_KSTEP_PADDED, LX, 2, iters);
    }
    if (dmma) {
      NEKO_TUNE_TIME(time3, CASE_VECTOR_DMMA, LX, 0, iters);
      NEKO_TUNE_TIME(time3, CASE_VECTOR_DMMA, LX, 1, iters);
      NEKO_TUNE_TIME(time3, CASE_VECTOR_DMMA, LX, 2, iters);
    }
    if (tma) {
      NEKO_TUNE_TIME(time4, CASE_VECTOR_DMMA_TMA, LX, 0, iters);
      NEKO_TUNE_TIME(time4, CASE_VECTOR_DMMA_TMA, LX, 1, iters);
      NEKO_TUNE_TIME(time4, CASE_VECTOR_DMMA_TMA, LX, 2, iters);
    }
    if (batch) {
      NEKO_TUNE_TIME(time5, CASE_VECTOR_DMMA_TMA_BATCH, LX, 0, iters);
      NEKO_TUNE_TIME(time5, CASE_VECTOR_DMMA_TMA_BATCH, LX, 1, iters);
      NEKO_TUNE_TIME(time5, CASE_VECTOR_DMMA_TMA_BATCH, LX, 2, iters);
    }
  }

  for (int c = 0; c < NEKO_EB_CANDIDATES; c++) {
    if (time2[c] >= NEKO_TUNE_INIT) { continue; }
    sprintf(neko_log_buf, "KSTEP eb=%-4d: %9.2f us/call",
            NEKO_EB_SEL(LX, c), NEKO_TUNE_US(time2[c], iters));
    log_message(neko_log_buf);
  }
  NEKO_TUNE_LOG_DMMA(LX, time3);
  NEKO_TUNE_LOG_DMMA_TMA(LX, time4);
  NEKO_TUNE_LOG_DMMA_TMA_BATCH(LX, time5);

  NEKO_TUNE_BEST(time2, best2, NEKO_EB_CANDIDATES);
  NEKO_TUNE_BEST(time3, best3, NEKO_DMMA_CANDIDATES);
  NEKO_TUNE_BEST(time4, best4, NEKO_DMMA_CANDIDATES);
  NEKO_TUNE_BEST(time5, best5, NEKO_DMMA_CANDIDATES);
  *eb_sel = best2;
  *nw_sel = best3;
  *tw_sel = best4;
  *bw_sel = best5;

  /* The dmma variants join the comparison only where they exist, their
     candidates are left at NEKO_TUNE_INIT otherwise */
  float tbest = time2[best2];

  retval = 2;
  if (time3[best3] < tbest) { retval = 3; tbest = time3[best3]; }
  if (time4[best4] < tbest) { retval = 4; tbest = time4[best4]; }
  if (time5[best5] < tbest) { retval = 5; }

  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  /* Leave the output of the chosen kernel in au, av and aw, see tune() */
  if (retval == 2) {
    CASE_VECTOR_KSTEP_PADDED_SEL(LX, best2);
    sprintf(neko_log_buf, "Chose        : 2 (KSTEP, %d elem/block)",
            NEKO_EB_SEL(LX, best2));
  } else if (retval == 3) {
    CASE_VECTOR_DMMA_SEL(LX, best3);
    sprintf(neko_log_buf, "Chose        : 3 (DMMA, %d warps, %d elem/blk)",
            NEKO_DMMA_NW(best3), NEKO_DMMA_PACK(LX));
  } else if (retval == 4) {
    CASE_VECTOR_DMMA_TMA_SEL(LX, best4);
    sprintf(neko_log_buf, "Chose        : 4 (DMMA_TMA, %d warps)",
            NEKO_DMMA_NW(best4));
  } else {
    CASE_VECTOR_DMMA_TMA_BATCH_SEL(LX, best5);
    sprintf(neko_log_buf, "Chose        : 5 (DMMA_TMA_BATCH, %d warps)",
            NEKO_DMMA_NW(best5));
  }
  log_message(neko_log_buf);
  log_end_section();
  return retval;
}
