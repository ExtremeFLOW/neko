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

#include <stdio.h>
#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "ax_helm_kernel.h"

extern "C" {

  /** 
   * Fortran wrapper for device HIP Ax
   */
  void hip_ax_helm(void *w, void *u, void *dx, void *dy, void *dz,
                   void *dxt, void *dyt, void *dzt, void *h1,
                   void *g11, void *g22, void *g33, void *g12,
                   void *g13, void *g23, int *nelv, int *lx) {

    const dim3 nthrds_1d(1024, 1, 1);
    const dim3 nblcks_1d((*nelv), 1, 1);
    const dim3 nthrds_kstep((*lx), (*lx), 1);
    const dim3 nblcks_kstep((*nelv), 1, 1);

#define CASE_1D(LX)                                                             \
    case LX:                                                                    \
      hipLaunchKernelGGL(HIP_KERNEL_NAME( ax_helm_kernel_1d<real, LX, 1024> ),  \
                         nblcks_1d, nthrds_1d, 0, 0,                            \
                         (real *) w, (real *) u,                                \
                         (real *) dx, (real *) dy, (real *) dz,                 \
                         (real *) dxt, (real *) dyt, (real *) dzt, (real *) h1, \
                         (real *) g11, (real *) g22, (real *) g33,              \
                         (real *) g12, (real *) g13, (real *) g23);             \
      HIP_CHECK(hipGetLastError());                                             \
      break

#define CASE_KSTEP(LX)                                                               \
    case LX:                                                                         \
      hipLaunchKernelGGL( HIP_KERNEL_NAME( ax_helm_kernel_kstep<real, LX> ),         \
                          nblcks_kstep, nthrds_kstep, 0, 0,                          \
                          (real *) w, (real *) u,                                    \
                          (real *) dx, (real *) dy, (real *) dz,  (real *) h1,       \
                          (real *) g11, (real *) g22, (real *) g33,                  \
                          (real *) g12, (real *) g13, (real *) g23);                 \
      HIP_CHECK(hipGetLastError());                                                  \
      break;


 
#define CASE_KSTEP_PADDED(LX)                                                        \
    case LX:                                                                         \
      hipLaunchKernelGGL( HIP_KERNEL_NAME( ax_helm_kernel_kstep_padded<real, LX> ),  \
                          nblcks_kstep, nthrds_kstep, 0, 0,                          \
                          (real *) w, (real *) u,                                    \
                          (real *) dx, (real *) dy, (real *) dz,  (real *) h1,       \
                          (real *) g11, (real *) g22, (real *) g33,                  \
                          (real *) g12, (real *) g13, (real *) g23);                 \
      HIP_CHECK(hipGetLastError());                                                  \
      break;

    switch(*lx) {
      CASE_KSTEP(2);
      CASE_KSTEP(3);
      CASE_KSTEP(4);
      CASE_KSTEP(5);
      CASE_KSTEP(6);
      CASE_KSTEP(7);
      CASE_KSTEP_PADDED(8);
      CASE_KSTEP(9);
      CASE_KSTEP(10);
      CASE_1D(11);
      CASE_1D(12);
    default:
      {
        fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
        exit(1);
      }
    }
  }
}
