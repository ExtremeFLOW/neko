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
#include "pipecg_kernel.h"

/**
 * @todo Make sure that this gets deleted at some point...
 */
real *buf = NULL;
real *buf_d1 = NULL;
real *buf_d2 = NULL;
real *buf_d3 = NULL;
int buf_len = 0;

extern "C" {
  
  void hip_cg_update_xp(void *x, void *p, void *u, void *alpha, void *beta,
                        int *p_cur, int *p_space, int *n) {
        
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(cg_update_xp_kernel<real>),
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,
                       (real *) x, (real *) p, (real **) u, (real *) alpha,
                       (real *) beta, *p_cur, *p_space, *n);
    HIP_CHECK(hipGetLastError());
  }

  void hip_pipecg_vecops(void *p, void *q, void *r, void *s, void *u1,
                            void *u2, void *w, void *z, 
                            void *ni, void *mi, real *alpha, 
                            real *beta, void *mult, 
                            real *reduction, int *n) {
        
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;

    if (buf != NULL && buf_len < nb) {
      HIP_CHECK(hipHostFree(buf));
      HIP_CHECK(hipFree(buf_d1));
      HIP_CHECK(hipFree(buf_d2));
      HIP_CHECK(hipFree(buf_d3));
      buf = NULL;
    }

    if (buf == NULL){
      HIP_CHECK(hipHostMalloc(&buf, nb*sizeof(real)));
      HIP_CHECK(hipMalloc(&buf_d1, nb*sizeof(real)));
      HIP_CHECK(hipMalloc(&buf_d2, nb*sizeof(real)));
      HIP_CHECK(hipMalloc(&buf_d3, nb*sizeof(real)));
      buf_len = nb;
    }
     
    hipLaunchKernelGGL(HIP_KERNEL_NAME(pipecg_vecops_kernel<real>),
                       nblcks, nthrds, 0, stream,
                       (real *) p, (real *) q, (real *) r, (real *) s,
                       (real *) u1, (real *) u2, (real *) w, (real *) z,
                       (real *) ni, (real *) mi, *alpha, *beta,
                       (real *) mult, buf_d1, buf_d2, buf_d3, *n);
    HIP_CHECK(hipGetLastError());

    hipLaunchKernelGGL(HIP_KERNEL_NAME(reduce_kernel<real>),
                       1, 1024, 0 , stream, buf_d1, nb);
    HIP_CHECK(hipGetLastError());
    hipLaunchKernelGGL(HIP_KERNEL_NAME(reduce_kernel<real>),
                       1, 1024, 0 , stream, buf_d2, nb);
    HIP_CHECK(hipGetLastError());
        hipLaunchKernelGGL(HIP_KERNEL_NAME(reduce_kernel<real>),
                       1, 1024, 0 , stream, buf_d3, nb);
    HIP_CHECK(hipGetLastError());

    HIP_CHECK(hipMemcpyAsync(&buf[0], buf_d1, sizeof(real),
                             hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipGetLastError());
    HIP_CHECK(hipMemcpyAsync(&buf[1], buf_d2, sizeof(real),
                             hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipGetLastError());
    HIP_CHECK(hipMemcpyAsync(&buf[2], buf_d3, sizeof(real),
                             hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipGetLastError());

    hipStreamSynchronize(stream);
        
    reduction[0] = buf[0];
    reduction[1] = buf[1];
    reduction[2] = buf[2];

  }
}
