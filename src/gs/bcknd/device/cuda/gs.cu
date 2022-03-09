/*
 Copyright (c) 2021, The Neko Authors
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

#include <climits>
#include <cstdio>
#include <device/device_config.h>
#include <device/cuda/check.h>
#include "gs_kernels.h"

#define GS_OP_ADD  1
#define GS_OP_MUL  2
#define GS_OP_MIN  3
#define GS_OP_MAX  4

extern "C" {

  /** 
   * Fortran wrapper for device gather kernels
   */
  void cuda_gather_kernel(void *v, int *m, int *o, void *dg,
			  void *u, int *n, void *gd, int *nb,
			  void *b, void *bo, int *op) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+ 1024 - 1)/ 1024, 1, 1);

    switch (*op) {
    case GS_OP_ADD:
      gather_kernel_add<real>
	<<<nblcks, nthrds>>>((real *) v, *m, *o, (int *) dg,
			     (real *) u, *n, (int *) gd,
			     *nb, (int *) b, (int *) bo);
      CUDA_CHECK(cudaGetLastError());
      break;
    case GS_OP_MUL:
      gather_kernel_mul<real>
	<<<nblcks, nthrds>>>((real *) v, *m, *o, (int *) dg,
			     (real *) u, *n, (int *) gd,
			     *nb, (int *) b, (int *) bo);
      CUDA_CHECK(cudaGetLastError());
      break;
    case GS_OP_MIN:
      gather_kernel_min<real>
	<<<nblcks, nthrds>>>((real *) v, *m, *o, (int *) dg,
			     (real *) u, *n, (int *) gd,
			     *nb, (int *) b, (int *) bo);
      CUDA_CHECK(cudaGetLastError());
      break;
    case GS_OP_MAX:
      gather_kernel_max<real>
	<<<nblcks, nthrds>>>((real *) v, *m, *o, (int *) dg,
			     (real *) u, *n, (int *) gd,
			     *nb, (int *) b, (int *) bo);
      CUDA_CHECK(cudaGetLastError());
      break;
    }
  }

  /**
   * Fortran wrapper for device scatter kernel
   */
  void cuda_scatter_kernel(void *v, int *m, void *dg,
			   void *u, int *n, void *gd,
			   int *nb, void *b, void *bo) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    scatter_kernel<real>
      <<<nblcks, nthrds>>>((real *) v, *m, (int *) dg,
			   (real *) u, *n, (int *) gd,
			   *nb, (int *) b, (int *) bo);
    CUDA_CHECK(cudaGetLastError());
  }
}
