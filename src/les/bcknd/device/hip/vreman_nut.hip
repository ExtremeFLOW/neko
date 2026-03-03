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

#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "vreman_nut_kernel.h"

extern "C" {
  void hip_vreman_nut_compute(void *a11, void *a12, void *a13,
                             void *a21, void *a22, void *a23,
                             void *a31, void *a32, void *a33,
                             void *delta, void *nut, void *mult,
                             real *c, real *eps, int * n){

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    hipLaunchKernelGGL(HIP_KERNEL_NAME(vreman_nut_compute<real>),
                                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                                       (real *) a11, (real *) a12, (real *) a13,
                                       (real *) a21, (real *) a22, (real *) a23,
                                       (real *) a31, (real *) a32, (real *) a33,
                                       (real *) delta, (real *) nut, (real *) mult,
                                       *c, * eps, * n);
    HIP_CHECK(hipGetLastError());
  }
}


extern "C" {
  void hip_vreman_nut_compute_buoy(void *a11, void *a12, void *a13,
                                   void *a21, void *a22, void *a23,
                                   void *a31, void *a32, void *a33,
                                   void *delta, void *nut, void *mult,
                                   real *c, real *eps, int *n,
                                   void *dTdx, void *dTdy, void *dTdz,
                                   real *g, real *ri_c, real *ref_temp)
  {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    // Normalize g on host (as CPU and CUDA do)
    real g1 = g[0], g2 = g[1], g3 = g[2];
    real gmag = sqrt(g1*g1 + g2*g2 + g3*g3);

    real n1 = g1 / gmag;
    real n2 = g2 / gmag;
    real n3 = g3 / gmag;

    hipLaunchKernelGGL(
        HIP_KERNEL_NAME(vreman_nut_compute_buoy<real>),
        nblcks, nthrds, 0, (hipStream_t)glb_cmd_queue,
        (real*)a11, (real*)a12, (real*)a13,
        (real*)a21, (real*)a22, (real*)a23,
        (real*)a31, (real*)a32, (real*)a33,
        (real*)delta, (real*)nut, (real*)mult,
        *c, *eps, *n,
        (real*)dTdx, (real*)dTdy, (real*)dTdz,
        n1, n2, n3,      // normalized gravity direction
        g1, g2, g3,      // original gravity vector
        *ri_c, *ref_temp
    );

    HIP_CHECK(hipGetLastError());
  }
}
