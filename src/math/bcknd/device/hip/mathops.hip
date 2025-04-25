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
#include <device/device_config.h>
#include <device/hip/check.h>
#include "mathops_kernel.h"

extern "C" {

  /** Fortran wrapper for opchsign \f$ a = -a \f$ */
  void hip_opchsign(void *a1, void *a2, void *a3, int *gdim, int *n) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(opchsign_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a1, (real *) a2, (real *) a3,
                       *gdim, *n);
    HIP_CHECK(hipGetLastError());
  }

  /** Fortran wrapper for opcolv \f$ a = a * c \f$ */
  void hip_opcolv(void *a1, void *a2, void *a3, void *c, int *gdim, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(opcolv_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a1, (real *) a2, (real *) a3, (real *) c,
                       *gdim, *n);
    HIP_CHECK(hipGetLastError());
    
  }

  /** Fortran wrapper for opcolv3c \f$ a(i) = b(i) * c(i) * d \f$ */
  void hip_opcolv3c(void *a1, void *a2, void *a3, void *b1, void *b2, void *b3,
                    void *c, real *d, int *gdim, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(opcolv3c_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a1, (real *) a2, (real *) a3,
                       (real *) b1, (real *) b2, (real *) b3,
                       (real *) c, *d, *gdim, *n);
    HIP_CHECK(hipGetLastError());
    
  }

  /** Fortran wrapper for opadd2cm \f$ a(i) = a + b(i) * c \f$ */
  void hip_opadd2cm(void *a1, void *a2, void *a3, 
                    void *b1, void *b2, void *b3, real *c, int *gdim, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(opadd2cm_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a1, (real *) a2, (real *) a3,
                       (real *) b1, (real *) b2, (real *) b3,
                       *c, *gdim, *n);
    HIP_CHECK(hipGetLastError());
    
  }

  /** Fortran wrapper for opadd2col \f$ a(i) = a + b(i) * c(i) \f$ */
  void hip_opadd2col(void *a1, void *a2, void *a3, 
                     void *b1, void *b2, void *b3, void *c, int *gdim, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(opadd2col_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) a1, (real *) a2, (real *) a3,
                       (real *) b1, (real *) b2, (real *) b3,
                       (real *) c, *gdim, *n);
    HIP_CHECK(hipGetLastError());
    
  }

}
