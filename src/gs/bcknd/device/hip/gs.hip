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

#include <hip/hip_runtime.h>
#include <climits>
#include <cstdio>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "gs_kernels.h"

#define GS_OP_ADD  1
#define GS_OP_MUL  2
#define GS_OP_MIN  3
#define GS_OP_MAX  4

extern "C" {

  /** 
   * Fortran wrapper for device gather kernels
   */
  void hip_gather_kernel(void *v, int *m, int *o, void *dg,
			 void *u, int *n, void *gd, int *nb,
			 void *b, void *bo, int *op) {

    if ((*m) == 0) return;
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+ 1024 - 1)/ 1024, 1, 1);

    switch (*op) {
    case GS_OP_ADD:
      hipLaunchKernelGGL(HIP_KERNEL_NAME(gather_kernel_add<real>),
			 nblcks, nthrds, 0, 0,
			 (real *) v, *m, *o, (int *) dg,
			 (real *) u, *n, (int *) gd,
			 *nb, (int *) b, (int *) bo);
      HIP_CHECK(hipGetLastError());
      break;
    case GS_OP_MUL:
      hipLaunchKernelGGL(HIP_KERNEL_NAME(gather_kernel_mul<real>),
			 nblcks, nthrds, 0, 0,
			 (real *) v, *m, *o, (int *) dg,
			 (real *) u, *n, (int *) gd,
			 *nb, (int *) b, (int *) bo);
      HIP_CHECK(hipGetLastError());
      break;
    case GS_OP_MIN:
      hipLaunchKernelGGL(HIP_KERNEL_NAME(gather_kernel_min<real>),
			 nblcks, nthrds, 0, 0,
			 (real *) v, *m, *o, (int *) dg,
			 (real *) u, *n, (int *) gd,
			 *nb, (int *) b, (int *) bo);
      HIP_CHECK(hipGetLastError());
      break;
    case GS_OP_MAX:
      hipLaunchKernelGGL(HIP_KERNEL_NAME(gather_kernel_max<real>),
			 nblcks, nthrds, 0, 0,
			 (real *) v, *m, *o, (int *) dg,
			 (real *) u, *n, (int *) gd,
			 *nb, (int *) b, (int *) bo);
      HIP_CHECK(hipGetLastError());
      break;
    }
  }

  /**
   * Fortran wrapper for device scatter kernel
   */
  void hip_scatter_kernel(void *v, int *m, void *dg,
			  void *u, int *n, void *gd,
			  int *nb, void *b, void *bo) {

    if ((*m) == 0) return;
	
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(scatter_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (real *) v, *m, (int *) dg,
		       (real *) u, *n, (int *) gd,
		       *nb, (int *) b, (int *) bo);
    HIP_CHECK(hipGetLastError());
  }

  void hip_gs_pack(void *u_d, void *buf_d, void *dof_d, int n) {

    const int nthrds = 1024;
    const int nblcks = (n + nthrds - 1) / nthrds;

    hipLaunchKernelGGL(HIP_KERNEL_NAME(gs_pack_kernel<real>),
		       nblcks, nthrds, 0, 0,
		       (real *) u_d, (real *) buf_d,
		       (int *) dof_d, n);
    HIP_CHECK(hipGetLastError());
  }

  void hip_gs_unpack(void *u_d, int op, void *buf_d, int *dof_d, int n) {

    const int nthrds = 1024;
    const int nblcks = (n + nthrds - 1) / nthrds;

    switch (op) {
    case GS_OP_ADD:
      hipLaunchKernelGGL(HIP_KERNEL_NAME(gs_unpack_add_kernel<real>),
			 nblcks, nthrds, 0, 0,
			 (real *) u_d, (real *) buf_d,
			 (int *) dof_d, n);
      break;
    default:
      printf("%s: unknown gs op %d\n", __FILE__, op);
      abort();
    }

    HIP_CHECK(hipGetLastError());
  }
}
