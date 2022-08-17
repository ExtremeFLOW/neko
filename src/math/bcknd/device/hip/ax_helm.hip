/*
 Copyright (c) 2021-2022, The Neko Authors
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
   * Fortran wrapper for device HIP Ax
   */
  void hip_ax_helm(void *w, void *u, void *dx, void *dy, void *dz,
                   void *dxt, void *dyt, void *dzt, void *h1,
                   void *g11, void *g22, void *g33, void *g12,
                   void *g13, void *g23, int *nelv, int *lx) {

    static int autotune[13] = { 0 };

    const dim3 nthrds_1d(1024, 1, 1);
    const dim3 nblcks_1d((*nelv), 1, 1);
    const dim3 nthrds_kstep((*lx), (*lx), 1);
    const dim3 nblcks_kstep((*nelv), 1, 1);

#define CASE_1D(LX)                                                             \
      hipLaunchKernelGGL(HIP_KERNEL_NAME( ax_helm_kernel_1d<real, LX, 1024> ),  \
                         nblcks_1d, nthrds_1d, 0, 0,                            \
                         (real *) w, (real *) u,                                \
                         (real *) dx, (real *) dy, (real *) dz,                 \
                         (real *) dxt, (real *) dyt, (real *) dzt, (real *) h1, \
                         (real *) g11, (real *) g22, (real *) g33,              \
                         (real *) g12, (real *) g13, (real *) g23);             \
      HIP_CHECK(hipGetLastError());

#define CASE_KSTEP(LX)                                                               \
      hipLaunchKernelGGL( HIP_KERNEL_NAME( ax_helm_kernel_kstep<real, LX> ),         \
                          nblcks_kstep, nthrds_kstep, 0, 0,                          \
                          (real *) w, (real *) u,                                    \
                          (real *) dx, (real *) dy, (real *) dz,  (real *) h1,       \
                          (real *) g11, (real *) g22, (real *) g33,                  \
                          (real *) g12, (real *) g13, (real *) g23);                 \
      HIP_CHECK(hipGetLastError());


 
#define CASE_KSTEP_PADDED(LX)                                                        \
      hipLaunchKernelGGL( HIP_KERNEL_NAME( ax_helm_kernel_kstep_padded<real, LX> ),  \
                          nblcks_kstep, nthrds_kstep, 0, 0,                          \
                          (real *) w, (real *) u,                                    \
                          (real *) dx, (real *) dy, (real *) dz,  (real *) h1,       \
                          (real *) g11, (real *) g22, (real *) g33,                  \
                          (real *) g12, (real *) g13, (real *) g23);                 \
      HIP_CHECK(hipGetLastError());

#define CASE(LX)                                                                     \
     case LX:                                                                        \
       if(autotune[LX] == 0 ) {                                                      \
         autotune[LX]=tune<LX>( w,  u,                                               \
                            dx,  dy, dz,                                             \
                            dxt, dyt, dzt,h1,                                        \
                            g11, g22, g33,                                           \
                            g12, g13, g23, nelv, lx);                                \
       } else if (autotune[LX] == 1 ) {                                              \
          CASE_1D(LX);                                                               \
       } else if (autotune[LX] == 2 ) {                                              \
          CASE_KSTEP(LX);                                                            \
       }                                                                             \
     break


#define CASE_PADDED(LX)                                                              \
     case LX:                                                                        \
       if(autotune[LX] == 0 ) {                                                      \
         autotune[LX]=tune_padded<LX>(w,  u,                                         \
                            dx,  dy, dz,                                             \
                            dxt, dyt, dzt,h1,                                        \
                            g11, g22, g33,                                           \
                            g12, g13, g23,nelv,lx);                                  \
       } else if (autotune[LX] == 1 ) {                                              \
          CASE_1D(LX);                                                               \
       } else if (autotune[LX] == 2 ) {                                              \
          CASE_KSTEP_PADDED(LX);                                                     \
       }                                                                             \
     break


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

}

template < const int LX > 
int tune(void *w, void *u, void *dx, void *dy, void *dz,
         void *dxt, void *dyt, void *dzt, void *h1,
         void *g11, void *g22, void *g33, void *g12,
         void *g13, void *g23, int *nelv, int *lx) {
   hipEvent_t start,stop;
   float time1,time2;
   int retval;

   const dim3 nthrds_1d(1024, 1, 1);
   const dim3 nblcks_1d((*nelv), 1, 1);
   const dim3 nthrds_kstep((*lx), (*lx), 1);
   const dim3 nblcks_kstep((*nelv), 1, 1);

   char *env_value=NULL;

   env_value=getenv("NEKO_AUTOTUNE");

   if(env_value) {
     if( !strcmp(env_value,"1D") ) {
       CASE_1D(LX);
       printf("Autotune %d, set by env to 1\n",*lx);
       return 1;
     } else if( !strcmp(env_value,"KSTEP") ) {
       CASE_KSTEP(LX);
       printf("Autotune %d, set by env to 2\n",*lx);
       return 2;
     } else {
       printf("ERROR : Invalid value set for NEKO_AUTOTUNE\n");
     }
   }
   
   hipEventCreate(&start);
   hipEventCreate(&stop);

   hipEventRecord(start,0);
   
   for(int i=0;i<5;i++) {
     CASE_1D(LX);
   }

   hipEventRecord(stop,0); 
   hipEventSynchronize(stop);
   hipEventElapsedTime(&time1, start, stop);

   hipEventRecord(start,0);

   for(int i=0;i<5;i++) {
     CASE_KSTEP(LX);
   }

   hipEventRecord(stop,0); 
   hipEventSynchronize(stop);
   hipEventElapsedTime(&time2, start, stop);

   if(time1<time2) {
     retval=1;
   } else {
     retval=2;
   }

   printf("Autotune %d, chose %d\n",*lx,retval);
   return retval;
}

template < const int LX >
int tune_padded(void *w, void *u, void *dx, void *dy, void *dz,
                void *dxt, void *dyt, void *dzt, void *h1,
                void *g11, void *g22, void *g33, void *g12,
                void *g13, void *g23, int *nelv, int *lx) {
   hipEvent_t start,stop;
   float time1,time2;
   int retval;

   const dim3 nthrds_1d(1024, 1, 1);
   const dim3 nblcks_1d((*nelv), 1, 1);
   const dim3 nthrds_kstep((*lx), (*lx), 1);
   const dim3 nblcks_kstep((*nelv), 1, 1);

   char *env_value=NULL;

   env_value=getenv("NEKO_AUTOTUNE");

   if(env_value) {
     if( !strcmp(env_value,"1D") ) {
       CASE_1D(LX);
       printf("Autotune %d, set by env to 1\n",*lx);
       return 1;
     } else if( !strcmp(env_value,"KSTEP") ) {
       CASE_KSTEP(LX);
       printf("Autotune %d, set by env to 2\n",*lx);
       return 2;
     } else {
       printf("ERROR : Invalid value set for NEKO_AUTOTUNE\n");
     }
   }


   hipEventCreate(&start);
   hipEventCreate(&stop);

   hipEventRecord(start,0);

   for(int i=0;i<5;i++) {
    CASE_1D(LX);
   }

   hipEventRecord(stop,0);
   hipEventSynchronize(stop);
   hipEventElapsedTime(&time1, start, stop);

   hipEventRecord(start,0);

   for(int i=0;i<5;i++) {
    CASE_KSTEP_PADDED(LX);
   }

   hipEventRecord(stop,0);
   hipEventSynchronize(stop);
   hipEventElapsedTime(&time2, start, stop);

   if(time1<time2) {
     retval=1;
   } else {
     retval=2;
   }

   printf("Autotune %d, chose %d\n",*lx,retval);
   return retval;
}

