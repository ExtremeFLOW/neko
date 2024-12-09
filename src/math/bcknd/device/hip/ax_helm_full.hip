/*
 Copyright (c) 2024, The Neko Authors
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
#include "ax_helm_full_kernel.h"

extern "C" {
  #include <common/neko_log.h>
}

extern "C" {

  /**
   * Fortran wrapper for device HIP Ax vector version
   */
  void hip_ax_helm_stress_vector(void *au, void *av, void *aw,
                           void *u, void *v, void *w,
                           void *dx, void *dy, void *dz,
                           void *dxt, void *dyt, void *dzt,
                           void *h1,
                           void *drdx, void *drdy, void *drdz,
                           void *dsdx, void *dsdy, void *dsdz,
                           void *dtdx, void *dtdy, void *dtdz,
                           void *jacinv, void *w3, int *nelv, int *lx) {

    const dim3 nthrds((*lx), (*lx), 1);
    const dim3 nblcks((*nelv), 1, 1);

#define CASE_VECTOR_KSTEP(LX)                                                   \
    hipLaunchKernelGGL( HIP_KERNEL_NAME( ax_helm_stress_kernel_vector_kstep<real, LX> ), \
                        nblcks, nthrds, 0,                                      \
                        (hipStream_t) glb_cmd_queue,                            \
                        (real *) au, (real *) av, (real *) aw,                  \
                        (real *) u, (real *) v, (real *) w,                     \
                        (real *) dx, (real *) dy, (real *) dz,  (real *) h1,    \
                        (real *) drdx, (real *) drdy, (real *) drdz,            \
                        (real *) dsdx, (real *) dsdy, (real *) dsdz,            \
                        (real *) dtdx, (real *) dtdy, (real *) dtdz,            \
                        (real *) jacinv, (real *) w3);              \
    HIP_CHECK(hipGetLastError());

#define CASE_VECTOR_KSTEP_PADDED(LX)                                            \
    hipLaunchKernelGGL( HIP_KERNEL_NAME( ax_helm_stress_kernel_vector_kstep_padded<real, LX> ), \
                        nblcks, nthrds, 0,                                      \
                        (hipStream_t) glb_cmd_queue,                            \
                        (real *) au, (real *) av, (real *) aw,                  \
                        (real *) u, (real *) v, (real *) w,                     \
                        (real *) dx, (real *) dy, (real *) dz,  (real *) h1,    \
                        (real *) drdx, (real *) drdy, (real *) drdz,            \
                        (real *) dsdx, (real *) dsdy, (real *) dsdz,            \
                        (real *) dtdx, (real *) dtdy, (real *) dtdz,            \
                        (real *) jacinv, (real *) w3);              \
    HIP_CHECK(hipGetLastError());

#define CASE_VECTOR(LX)                                                         \
    case LX:                                                                    \
      CASE_VECTOR_KSTEP(LX);                                                    \
       break

#define CASE_VECTOR_PADDED(LX)                                                  \
    case LX:                                                                    \
      CASE_VECTOR_KSTEP_PADDED(LX);                                             \
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
  void hip_ax_helm_stress_vector_part2(void *au, void *av, void *aw,
                                 void *u, void *v, void *w,
                                 void *h2, void *B, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;

    hipLaunchKernelGGL( HIP_KERNEL_NAME( ax_helm_stress_kernel_vector_part2<real> ),
                        nblcks, nthrds, 0, stream,
                        (real *) au, (real *) av, (real *) aw,
                        (real *) u, (real *) v, (real *) w,
                        (real *) h2, (real *) B, *n);
  }

}
