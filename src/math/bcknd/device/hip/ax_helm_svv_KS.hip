/*
 Copyright (c) 2025, The Neko Authors
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
#include <stdlib.h>
#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "ax_helm_svv_KS_kernel.h"

extern "C" {

/**
 * Fortran wrapper for the fused HIP SVV Helmholtz operator.
 */
void hip_ax_helm_svv_KS(
    void *w, void *u, void *dx, void *dy, void *dz, void *h1,
    void *drdx, void *drdy, void *drdz,
    void *dsdx, void *dsdy, void *dsdz,
    void *dtdx, void *dtdy, void *dtdz,
    void *jacinv, void *w3, void *h1_svv,
    void *filter_r, void *filter_s, void *filter_t,
    int *nelv, int *lx) {

  const dim3 threads(*lx, *lx, 1);
  const dim3 blocks(*nelv, 1, 1);
  const size_t shared_size =
      2 * (*lx) * (*lx) * (*lx) * sizeof(real);

#define CASE(LX)                                                               \
  case LX:                                                                     \
    hipLaunchKernelGGL(HIP_KERNEL_NAME(ax_helm_svv_KS_kernel<real, LX>),          \
        blocks, threads, shared_size, (hipStream_t) glb_cmd_queue,             \
        (real *) w, (real *) u,                                                \
        (real *) dx, (real *) dy, (real *) dz, (real *) h1,                   \
        (real *) drdx, (real *) drdy, (real *) drdz,                          \
        (real *) dsdx, (real *) dsdy, (real *) dsdz,                          \
        (real *) dtdx, (real *) dtdy, (real *) dtdz,                          \
        (real *) jacinv, (real *) w3, (real *) h1_svv,                        \
        (real *) filter_r, (real *) filter_s, (real *) filter_t);             \
    HIP_CHECK(hipGetLastError());                                              \
    break

  switch (*lx) {
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
    default:
      fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
      exit(1);
  }
}

}
