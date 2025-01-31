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

#include "fusedcg_cpld_kernel.h"
#include <device/device_config.h>
#include <device/cuda/check.h>

/**
 * @todo Make sure that this gets deleted at some point...
 */
real *fusedcg_cpld_buf = NULL;
real *fusedcg_cpld_buf_d = NULL;
int fusedcg_cpld_buf_len = 0;

extern "C" {

#include <math/bcknd/device/device_mpi_reduce.h>
#include <math/bcknd/device/device_mpi_op.h>

  void cuda_fusedcg_cpld_part1(void *a1, void *a2, void *a3,
                               void *b1, void *b2, void *b3,
                               void *tmp, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;
    
    fusedcg_cpld_part1_kernel<real>
      <<<nblcks, nthrds, 0, stream>>>((real *) a1, (real *) a2, (real *) a3,
                                      (real *) b1, (real *) b2, (real *) b3,
                                      (real *) tmp, *n);
    CUDA_CHECK(cudaGetLastError());
  }
  
  void cuda_fusedcg_cpld_update_p(void *p1, void *p2, void *p3,
                                  void *z1, void *z2, void *z3,
                                  void *po1, void *po2, void *po3,
                                  real *beta, int *n) {
        
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;
    
    fusedcg_cpld_update_p_kernel<real>
      <<<nblcks, nthrds, 0, stream>>>((real *) p1, (real *) p2, (real *) p3,
                                      (real *) z1, (real *) z2, (real *) z3, 
                                      (real *) po1, (real *) po2, (real *) po3,
                                      *beta, *n);
    CUDA_CHECK(cudaGetLastError());

  }

  void cuda_fusedcg_cpld_update_x(void *x1, void *x2, void *x3,
                                  void *p1, void *p2, void *p3,
                                  void *alpha, int *p_cur, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;
    
    fusedcg_cpld_update_x_kernel<real>
      <<<nblcks, nthrds, 0, stream>>>((real *) x1, (real *) x2, (real *) x3,
                                      (const real **) p1, (const real **) p2,
                                      (const real **) p3, (const real *) alpha,
                                      *p_cur, *n);
    CUDA_CHECK(cudaGetLastError());
  }

  real cuda_fusedcg_cpld_part2(void *a1, void *a2, void *a3, void *b,
                               void *c1, void *c2, void *c3, void *alpha_d ,
                               real *alpha, int *p_cur, int * n) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

    if (fusedcg_cpld_buf != NULL && fusedcg_cpld_buf_len < nb) {
      CUDA_CHECK(cudaFreeHost(fusedcg_cpld_buf));
      CUDA_CHECK(cudaFree(fusedcg_cpld_buf_d));
      fusedcg_cpld_buf = NULL;
    }
    
    if (fusedcg_cpld_buf == NULL){
      CUDA_CHECK(cudaMallocHost(&fusedcg_cpld_buf, 2*sizeof(real)));
      CUDA_CHECK(cudaMalloc(&fusedcg_cpld_buf_d, nb*sizeof(real)));
      fusedcg_cpld_buf_len = nb;
    }

    /* Store alpha(p_cur) in pinned memory */
    fusedcg_cpld_buf[1] = (*alpha);

    /* Update alpha_d(p_cur) = alpha(p_cur) */
    real *alpha_d_p_cur = ((real *) alpha_d) + ((*p_cur - 1));   
    CUDA_CHECK(cudaMemcpyAsync(alpha_d_p_cur, &fusedcg_cpld_buf[1],
                               sizeof(real), cudaMemcpyHostToDevice,
                               stream));

   
    fusedcg_cpld_part2_kernel<real>
      <<<nblcks, nthrds, 0, stream>>>((real *) a1, (real *) a2, (real *) a3,
                                      (real *) b, (real *) c1, (real *) c2,
                                      (real *) c3, *alpha, fusedcg_cpld_buf_d,
                                      *n);
    CUDA_CHECK(cudaGetLastError());

    reduce_kernel<real><<<1, 1024, 0, stream>>>(fusedcg_cpld_buf_d, nb);
    CUDA_CHECK(cudaGetLastError());

#ifdef HAVE_DEVICE_MPI
    cudaStreamSynchronize(stream);
    device_mpi_allreduce(fusedcg_cpld_buf_d, fusedcg_cpld_buf, 1,
                         sizeof(real), DEVICE_MPI_SUM);
#else
    CUDA_CHECK(cudaMemcpyAsync(fusedcg_cpld_buf, fusedcg_cpld_buf_d, sizeof(real),
                               cudaMemcpyDeviceToHost, stream));
    cudaStreamSynchronize(stream);
#endif

    return fusedcg_cpld_buf[0];
  }
}

