/*
 Copyright (c) 2021-2023, The Neko Authors
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
#include "cdtp_kernel.h"
#include <device/device_config.h>
#include <device/cuda/check.h>

extern "C" {
  #include <common/neko_log.h>
}

template < const int >
int tune_cdtp(void *dtx, void *x,
              void *dr, void *ds, void *dt,
              void *dxt, void *dyt, void *dzt,
              void *B, void *jac, int *nel, int *lx);

extern "C" {

  /** 
   * Fortran wrapper for device cuda \f$ D^T X \f$
   */
  void cuda_cdtp(void *dtx, void *x,
                 void *dr, void *ds, void *dt,
                 void *dxt, void *dyt, void *dzt,
                 void *B, void *jac, int *nel, int *lx) {
    
    static int autotune[17] = { 0 };
    
    const dim3 nthrds_1d(1024, 1, 1);
    const dim3 nthrds_kstep((*lx), (*lx), 1);
    const dim3 nblcks((*nel), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define CASE_1D(LX)                                                             \
    cdtp_kernel_1d<real, LX, 1024>                                              \
      <<<nblcks, nthrds_1d, 0, stream>>>((real *) dtx, (real *) x,              \
                              (real *) dr, (real *) ds, (real *) dt,            \
                              (real *) dxt, (real *) dyt, (real *) dzt,         \
                              (real *) B, (real *) jac);                        \
    CUDA_CHECK(cudaGetLastError());

#define CASE_KSTEP(LX)                                                          \
    cdtp_kernel_kstep<real, LX>                                                 \
      <<<nblcks, nthrds_kstep, 0, stream>>>((real *) dtx, (real *) x,           \
                                 (real *) dr, (real *) ds, (real *) dt,         \
                                 (real *) dxt, (real *) dyt, (real *) dzt,      \
                                 (real *) B, (real *) jac);                     \
    CUDA_CHECK(cudaGetLastError());

#define CASE(LX)                                                                \
    case LX:                                                                    \
      if(autotune[LX] == 0 ) {                                                  \
        autotune[LX]=tune_cdtp<LX>(dtx, x,                                      \
                                   dr, ds, dt,                                  \
                                   dxt, dyt, dzt,                               \
                                   B, jac, nel, lx);                            \
      } else if (autotune[LX] == 1 ) {                                          \
        CASE_1D(LX);                                                            \
      } else if (autotune[LX] == 2 ) {                                          \
        CASE_KSTEP(LX);                                                         \
      }                                                                         \
      break

#define CASE_LARGE(LX)                                                          \
    case LX:                                                                    \
      CASE_KSTEP(LX);                                                           \
      break


    if ((*lx) < 13) {      
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
    else {
      switch(*lx) {
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
int tune_cdtp(void *dtx, void *x,
              void *dr, void *ds, void *dt,
              void *dxt, void *dyt, void *dzt,
              void *B, void *jac, int *nel, int *lx) {
  cudaEvent_t start,stop;
  float time1,time2;
  int retval;

  const dim3 nthrds_1d(1024, 1, 1);
  const dim3 nthrds_kstep((*lx), (*lx), 1);
  const dim3 nblcks((*nel), 1, 1);
  const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;
  
  char *env_value = NULL;
  char neko_log_buf[80];
  
  env_value=getenv("NEKO_AUTOTUNE");

  sprintf(neko_log_buf, "Autotune cdtp (lx: %d)", *lx);
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
