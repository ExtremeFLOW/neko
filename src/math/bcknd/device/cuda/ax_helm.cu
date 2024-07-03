/*
 Copyright (c) 2021-2024, The Neko Authors
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
#include <device/device_config.h>
#include <device/cuda/check.h>

extern "C" {
  #include <common/neko_log.h>
}

template < const int>
int tune(void *w, void *u, void *dx, void *dy, void *dz,
         void *dxt, void *dyt, void *dzt, void *h1,
         void *g11, void *g22, void *g33, void *g12,
         void *g13, void *g23, int *nelv, int *lx);

template < const int>
int tune_padded(void *w, void *u, void *dx, void *dy, void *dz,
                void *dxt, void *dyt, void *dzt, void *h1,
                void *g11, void *g22, void *g33, void *g12,
                void *g13, void *g23, int *nelv, int *lx);

extern "C" {

  /**
   * Fortran wrapper for device CUDA Ax
   */
  void cuda_ax_helm(void *w, void *u, void *dx, void *dy, void *dz,
                    void *dxt, void *dyt, void *dzt, void *h1,
                    void *g11, void *g22, void *g33, void *g12,
                    void *g13, void *g23, int *nelv, int *lx) {

    static int autotune[17] = { 0 };

    const dim3 nthrds_1d(1024, 1, 1);
    const dim3 nblcks_1d((*nelv), 1, 1);
    const dim3 nthrds_kstep((*lx), (*lx), 1);
    const dim3 nblcks_kstep((*nelv), 1, 1);

    const dim3 nthrds((*lx), (*lx), 1);
    const dim3 nblcks((*nelv), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define CASE_1D(LX)                                                             \
    ax_helm_kernel_1d<real, LX, 1024>                                           \
      <<<nblcks_1d, nthrds_1d, 0, stream>>>((real *) w, (real *) u,             \
                          (real *) dx, (real *) dy, (real *) dz,                \
                          (real *) dxt, (real *) dyt, (real *) dzt, (real *) h1,\
                          (real *) g11, (real *) g22, (real *) g33,             \
                          (real *) g12, (real *) g13, (real *) g23);            \
      CUDA_CHECK(cudaGetLastError());

#define CASE_KSTEP(LX)                                                          \
    ax_helm_kernel_kstep<real, LX>                                              \
      <<<nblcks_kstep, nthrds_kstep, 0, stream>>>((real *) w, (real *) u,       \
                           (real *) dx, (real *) dy, (real *) dz, (real *) h1,  \
                           (real *) g11, (real *) g22, (real *) g33,            \
                           (real *) g12, (real *) g13, (real *) g23);           \
      CUDA_CHECK(cudaGetLastError());

#define CASE_KSTEP_PADDED(LX)                                                   \
    ax_helm_kernel_kstep_padded<real, LX>                                       \
    <<<nblcks_kstep, nthrds_kstep, 0, stream>>>((real *) w, (real *) u,         \
                           (real *) dx, (real *) dy, (real *) dz, (real *) h1,  \
                           (real *) g11, (real *) g22, (real *) g33,            \
                           (real *) g12, (real *) g13, (real *) g23);           \
      CUDA_CHECK(cudaGetLastError());

#define CASE(LX)                                                                \
    case LX:                                                                    \
      if(autotune[LX] == 0 ) {                                                  \
        autotune[LX]=tune<LX>( w,  u,                                           \
                               dx,  dy, dz,                                     \
                               dxt, dyt, dzt,h1,                                \
                               g11, g22, g33,                                   \
                               g12, g13, g23, nelv, lx);                        \
      } else if (autotune[LX] == 1 ) {                                          \
        CASE_1D(LX);                                                            \
      } else if (autotune[LX] == 2 ) {                                          \
        CASE_KSTEP(LX);                                                         \
      }                                                                         \
      break


#define CASE_PADDED(LX)                                                         \
    case LX:                                                                    \
      if(autotune[LX] == 0 ) {                                                  \
        autotune[LX]=tune_padded<LX>(w,  u,                                     \
                                     dx,  dy, dz,                               \
                                     dxt, dyt, dzt,h1,                          \
                                     g11, g22, g33,                             \
                                     g12, g13, g23,nelv,lx);                    \
      } else if (autotune[LX] == 1 ) {                                          \
        CASE_1D(LX);                                                            \
      } else if (autotune[LX] == 2 ) {                                          \
        CASE_KSTEP_PADDED(LX);                                                  \
      }                                                                         \
      break

#define CASE_LARGE(LX)                                                          \
    case LX:                                                                    \
      CASE_KSTEP(LX);                                                           \
      break

#define CASE_LARGE_PADDED(LX)                                                   \
    case LX:                                                                    \
      CASE_KSTEP_PADDED(LX);                                                    \
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

    const dim3 nthrds((*lx), (*lx), 1);
    const dim3 nblcks((*nelv), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define CASE_VECTOR_KSTEP(LX)                                                  \
    ax_helm_kernel_vector_kstep<real, LX>                                      \
    <<<nblcks, nthrds, 0, stream>>> ((real *) au, (real *) av, (real *) aw,    \
                                     (real *) u, (real *) v, (real *) w,       \
                                     (real *) dx, (real *) dy, (real *) dz,    \
                                     (real *) h1, (real *) g11, (real *) g22,  \
                                     (real *) g33, (real *) g12, (real *) g13, \
                                     (real *) g23);                            \
    CUDA_CHECK(cudaGetLastError());

#define CASE_VECTOR_KSTEP_PADDED(LX)                                           \
    ax_helm_kernel_vector_kstep_padded<real, LX>                               \
    <<<nblcks, nthrds, 0, stream>>> ((real *) au, (real *) av, (real *) aw,    \
                                     (real *) u, (real *) v, (real *) w,       \
                                     (real *) dx, (real *) dy, (real *) dz,    \
                                     (real *) h1, (real *) g11, (real *) g22,  \
                                     (real *) g33, (real *) g12, (real *) g13, \
                                     (real *) g23);                            \
    CUDA_CHECK(cudaGetLastError());


#define CASE_VECTOR(LX)                                                        \
    case LX:                                                                   \
      CASE_VECTOR_KSTEP(LX);                                                   \
       break

#define CASE_VECTOR_PADDED(LX)                                                 \
    case LX:                                                                   \
      CASE_VECTOR_KSTEP_PADDED(LX);                                            \
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
         void *g13, void *g23, int *nelv, int *lx) {
  cudaEvent_t start,stop;
  float time1,time2;
  int retval;

  const dim3 nthrds_1d(1024, 1, 1);
  const dim3 nblcks_1d((*nelv), 1, 1);
  const dim3 nthrds_kstep((*lx), (*lx), 1);
  const dim3 nblcks_kstep((*nelv), 1, 1);
 const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

  char *env_value = NULL;
  char neko_log_buf[80];

  env_value=getenv("NEKO_AUTOTUNE");

  sprintf(neko_log_buf, "Autotune Ax helm (lx: %d)", *lx);
  log_section(neko_log_buf);

  if(env_value) {
    if( !strcmp(env_value,"1D") ) {
      CASE_1D(LX);
      sprintf(neko_log_buf,"Set by env : 1 (1D)");
      log_message(neko_log_buf);
      log_end_section();
      return 1;
    } else if( !strcmp(env_value,"KSTEP") ) {
      CASE_KSTEP(LX);
      sprintf(neko_log_buf,"Set by env : 2 (KSTEP)");
      log_message(neko_log_buf);
      log_end_section();
      return 2;
    } else {
       sprintf(neko_log_buf, "Invalid value set for NEKO_AUTOTUNE");
       log_error(neko_log_buf);
    }
  }

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start,0);

  for(int i = 0; i < 100; i++) {
    CASE_1D(LX);
  }

  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time1, start, stop);

  cudaEventRecord(start,0);

  for(int i = 0; i < 100; i++) {
     CASE_KSTEP(LX);
   }

  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time2, start, stop);

  if(time1 < time2) {
     retval = 1;
  } else {
    retval = 2;
  }

  sprintf(neko_log_buf, "Chose      : %d (%s)", retval,
          (retval > 1 ? "KSTEP" : "1D"));
  log_message(neko_log_buf);
  log_end_section();
  return retval;
}

template < const int LX >
int tune_padded(void *w, void *u, void *dx, void *dy, void *dz,
                void *dxt, void *dyt, void *dzt, void *h1,
                void *g11, void *g22, void *g33, void *g12,
                void *g13, void *g23, int *nelv, int *lx) {
  cudaEvent_t start, stop;
  float time1, time2;
  int retval;

  const dim3 nthrds_1d(1024, 1, 1);
  const dim3 nblcks_1d((*nelv), 1, 1);
  const dim3 nthrds_kstep((*lx), (*lx), 1);
  const dim3 nblcks_kstep((*nelv), 1, 1);
 const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

  char *env_value = NULL;
  char neko_log_buf[80];

  env_value=getenv("NEKO_AUTOTUNE");

  sprintf(neko_log_buf, "Autotune Ax helm (lx: %d)", *lx);
  log_section(neko_log_buf);

  if(env_value) {
    if( !strcmp(env_value,"1D") ) {
      CASE_1D(LX);
      sprintf(neko_log_buf,"Set by env : 1 (1D)");
      log_message(neko_log_buf);
      log_end_section();
      return 1;
     } else if( !strcmp(env_value,"KSTEP") ) {
      CASE_KSTEP(LX);
      sprintf(neko_log_buf,"Set by env : 2 (KSTEP)");
      log_message(neko_log_buf);
      log_end_section();
      return 2;
    } else {
      sprintf(neko_log_buf, "Invalid value set for NEKO_AUTOTUNE");
      log_error(neko_log_buf);
    }
  }

  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start,0);

  for(int i = 0; i < 100; i++) {
    CASE_1D(LX);
  }

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time1, start, stop);

  cudaEventRecord(start, 0);

  for(int i = 0; i < 100; i++) {
    CASE_KSTEP_PADDED(LX);
  }

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time2, start, stop);

  if(time1 < time2) {
    retval=1;
  } else {
    retval=2;
  }

  sprintf(neko_log_buf, "Chose      : %d (%s)", retval,
          (retval > 1 ? "KSTEP" : "1D"));
  log_message(neko_log_buf);
  log_end_section();
  return retval;
}
