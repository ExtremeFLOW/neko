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

#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>

#include "constrain_mixed_bc_kernel.h"

extern "C" {

  void hip_constrain_mixed_bc_zero(void *mixed_msk, void *x, void *y,
                                   void *z, int *constraint_n,
                                   int *constraint_t1, int *constraint_t2,
                                   void *n, void *t1, void *t2, int *m,
                                   hipStream_t strm) {

    if (*m <= 0)
      return;

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(constrain_mixed_bc_zero_kernel<real>),
                       nblcks, nthrds, 0, strm,
                       (const int *) mixed_msk,
                       (real *) x,
                       (real *) y,
                       (real *) z,
                       *constraint_n,
                       *constraint_t1,
                       *constraint_t2,
                       (const real *) n,
                       (const real *) t1,
                       (const real *) t2,
                       *m);
    HIP_CHECK(hipGetLastError());
  }

  void hip_constrain_mixed_bc_set(void *mixed_msk, void *x, void *y,
                                  void *z, int *constraint_n,
                                  int *constraint_t1, int *constraint_t2,
                                  void *n, void *t1, void *t2,
                                  void *values_n, void *values_t1,
                                  void *values_t2, int *m,
                                  hipStream_t strm) {

    if (*m <= 0)
      return;

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(constrain_mixed_bc_set_kernel<real>),
                       nblcks, nthrds, 0, strm,
                       (const int *) mixed_msk,
                       (real *) x,
                       (real *) y,
                       (real *) z,
                       *constraint_n,
                       *constraint_t1,
                       *constraint_t2,
                       (const real *) n,
                       (const real *) t1,
                       (const real *) t2,
                       (const real *) values_n,
                       (const real *) values_t1,
                       (const real *) values_t2,
                       *m);
    HIP_CHECK(hipGetLastError());
  }

  void hip_constrain_mixed_bc_set_const(void *mixed_msk, void *x, void *y,
                                        void *z, int *constraint_n,
                                        int *constraint_t1,
                                        int *constraint_t2, void *n,
                                        void *t1, void *t2, real *value_n,
                                        real *value_t1, real *value_t2,
                                        int *m, hipStream_t strm) {

    if (*m <= 0)
      return;

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(constrain_mixed_bc_set_const_kernel<real>),
                       nblcks, nthrds, 0, strm,
                       (const int *) mixed_msk,
                       (real *) x,
                       (real *) y,
                       (real *) z,
                       *constraint_n,
                       *constraint_t1,
                       *constraint_t2,
                       (const real *) n,
                       (const real *) t1,
                       (const real *) t2,
                       *value_n,
                       *value_t1,
                       *value_t2,
                       *m);
    HIP_CHECK(hipGetLastError());
  }

}
