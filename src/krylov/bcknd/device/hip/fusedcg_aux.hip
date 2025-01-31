/*
 Copyright (c) 2021-2024, The Neko Authors
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
#include "fusedcg_kernel.h"
#include <device/device_config.h>
#include <device/hip/check.h>

/**
 * @todo Make sure that this gets deleted at some point...
 */
real *fusedcg_buf = NULL;
real *fusedcg_buf_d = NULL;
int fusedcg_buf_len = 0;

extern "C" {

#include <math/bcknd/device/device_mpi_reduce.h>
#include <math/bcknd/device/device_mpi_op.h>


  
  void hip_fusedcg_update_p(void *p, void *z, void *po, real *beta, int *n) {
        
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(fusedcg_update_p_kernel<real>),
                       nblcks, nthrds, 0, stream,
                       (real *) p, (real *) z, (real *) po, *beta, *n);
    HIP_CHECK(hipGetLastError());

  }

  void hip_fusedcg_update_x(void *x, void *p, void *alpha, int *p_cur, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(fusedcg_update_x_kernel<real>),
                       nblcks, nthrds, 0, stream,
                       (real *) x, (const real **) p,
                       (const real *) alpha, *p_cur, *n);
    HIP_CHECK(hipGetLastError());
  }


  real hip_fusedcg_part2(void *a, void *b, void *c,
                          void *alpha_d , real *alpha, int *p_cur, int * n) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;

    if (fusedcg_buf != NULL && fusedcg_buf_len < nb) {
      HIP_CHECK(hipHostFree(fusedcg_buf));
      HIP_CHECK(hipFree(fusedcg_buf_d));
      fusedcg_buf = NULL;
    }
    
    if (fusedcg_buf == NULL){
      HIP_CHECK(hipHostMalloc(&fusedcg_buf, 2*sizeof(real)));
      HIP_CHECK(hipMalloc(&fusedcg_buf_d, nb*sizeof(real)));
      fusedcg_buf_len = nb;
    }

    /* Store alpha(p_cur) in pinned memory */
    fusedcg_buf[1] = (*alpha);

    /* Update alpha_d(p_cur) = alpha(p_cur) */
    real *alpha_d_p_cur = ((real *) alpha_d) + ((*p_cur - 1));   
    HIP_CHECK(hipMemcpyAsync(alpha_d_p_cur, &fusedcg_buf[1], sizeof(real),
                               hipMemcpyHostToDevice, stream));
   
    hipLaunchKernelGGL(HIP_KERNEL_NAME(fusedcg_part2_kernel<real>),
                       nblcks, nthrds, 0, stream,
                       (real *) a, (real *) b, (real *) c,
                       *alpha, fusedcg_buf_d, * n);
    HIP_CHECK(hipGetLastError());

    hipLaunchKernelGGL(HIP_KERNEL_NAME(reduce_kernel<real>),
                       1, 1024, 0, stream, fusedcg_buf_d, nb);
    HIP_CHECK(hipGetLastError());

#ifdef HAVE_DEVICE_MPI
    HIP_CHECK(hipStreamSynchronize(stream));
    device_mpi_allreduce(fusedcg_buf_d, fusedcg_buf, 1,
                         sizeof(real), DEVICE_MPI_SUM);
#else
    HIP_CHECK(hipMemcpyAsync(fusedcg_buf, fusedcg_buf_d, sizeof(real),
                               hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipStreamSynchronize(stream));
#endif

    return fusedcg_buf[0];
  }
}

