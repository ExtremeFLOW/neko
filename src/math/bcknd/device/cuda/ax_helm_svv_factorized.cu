/*
 Copyright (c) 2026, The Neko Authors
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
#include <device/device_config.h>
#include <device/cuda/check.h>
#include "ax_helm_kernel.h"

extern "C" {

/** Fortran wrapper for the fused factorized CUDA SVV Helmholtz operator. */
void cuda_ax_helm_svv_factorized(
    void *w, void *u, void *dx, void *dy, void *dz, void *h1,
    void *Br, void *Bs, void *Bt, void *h1_svv,
    void *g11, void *g22, void *g33,
    void *g12, void *g13, void *g23,
    int *nelv, int *lx) {

  const dim3 threads(*lx, *lx, 1);
  const dim3 blocks(*nelv, 1, 1);
  const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

#define LAUNCH(LX) \
    ax_helm_kernel_kstep<real, LX, 1, false> \
        <<<blocks, threads, 0, stream>>>( \
        (real *) w, (real *) u, \
        (real *) dx, (real *) dy, (real *) dz, (real *) h1, \
        (real *) g11, (real *) g22, (real *) g33, \
        (real *) g12, (real *) g13, (real *) g23, *nelv); \
    CUDA_CHECK(cudaGetLastError()); \
    ax_helm_kernel_kstep<real, LX, 1, true> \
        <<<blocks, threads, 0, stream>>>( \
        (real *) w, (real *) u, \
        (real *) Br, (real *) Bs, (real *) Bt, (real *) h1_svv, \
        (real *) g11, (real *) g22, (real *) g33, \
        (real *) g12, (real *) g13, (real *) g23, *nelv); \
    CUDA_CHECK(cudaGetLastError())

#define LAUNCH_PADDED(LX) \
    ax_helm_kernel_kstep_padded<real, LX, 1, false> \
        <<<blocks, threads, 0, stream>>>( \
        (real *) w, (real *) u, \
        (real *) dx, (real *) dy, (real *) dz, (real *) h1, \
        (real *) g11, (real *) g22, (real *) g33, \
        (real *) g12, (real *) g13, (real *) g23, *nelv); \
    CUDA_CHECK(cudaGetLastError()); \
    ax_helm_kernel_kstep_padded<real, LX, 1, true> \
        <<<blocks, threads, 0, stream>>>( \
        (real *) w, (real *) u, \
        (real *) Br, (real *) Bs, (real *) Bt, (real *) h1_svv, \
        (real *) g11, (real *) g22, (real *) g33, \
        (real *) g12, (real *) g13, (real *) g23, *nelv); \
    CUDA_CHECK(cudaGetLastError())

#define CASE(LX) \
  case LX: \
    LAUNCH(LX); \
    break

#define CASE_PADDED(LX) \
  case LX: \
    LAUNCH_PADDED(LX); \
    break

  switch (*lx) {
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
    CASE(12);
    CASE(13);
    CASE(14);
    CASE(15);
    CASE_PADDED(16);
    default:
      fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
      exit(1);
  }
}

}
