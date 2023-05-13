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

#include "pipecg_kernel.h"
#include <device/device_config.h>
#include <device/cuda/check.h>

/**
 * @todo Make sure that this gets deleted at some point...
 */
real *buf1 = NULL;
real *buf2 = NULL;
real *buf3 = NULL;
real *buf_d1 = NULL;
real *buf_d2 = NULL;
real *buf_d3 = NULL;
int buf_len = 0;

extern "C" {
  
  void cuda_cg_update_xp(void *x, void *p, void *u, void *alpha, void *beta,
                         int *p_cur, int *p_space, int *n) {
        
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    cg_update_xp_kernel<real>
      <<<nblcks, nthrds>>>((real *) x, (real *) p,(real **) u, (real *) alpha,
                           (real *) beta, *p_cur, *p_space, *n);
    CUDA_CHECK(cudaGetLastError());

  }

  void cuda_pipecg_vecops(void *p, void *q, void *r, void *s, void *u1,
                            void *u2, void *w, void *z, 
                            void *ni, void *mi, real *alpha, 
                            real *beta, void *mult, 
                            real *reduction, int *n) {
        
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;

    if (buf1 != NULL && buf_len < nb) {
      CUDA_CHECK(cudaFreeHost(buf1));
      CUDA_CHECK(cudaFreeHost(buf2));
      CUDA_CHECK(cudaFreeHost(buf3));
      CUDA_CHECK(cudaFree(buf_d1));
      CUDA_CHECK(cudaFree(buf_d2));
      CUDA_CHECK(cudaFree(buf_d3));
      buf1 = NULL;
    }

    if (buf1 == NULL){
      CUDA_CHECK(cudaMallocHost(&buf1, nb*sizeof(real)));
      CUDA_CHECK(cudaMallocHost(&buf2, nb*sizeof(real)));
      CUDA_CHECK(cudaMallocHost(&buf3, nb*sizeof(real)));
      CUDA_CHECK(cudaMalloc(&buf_d1, nb*sizeof(real)));
      CUDA_CHECK(cudaMalloc(&buf_d2, nb*sizeof(real)));
      CUDA_CHECK(cudaMalloc(&buf_d3, nb*sizeof(real)));
      buf_len = nb;
    }
     
    pipecg_vecops_kernel<real>
      <<<nblcks, nthrds>>>((real *) p, (real *) q,
                           (real *) r, (real *) s,
                           (real *) u1, (real *) u2,
                           (real *) w, (real *) z,
                           (real *) ni, (real *) mi, 
                           *alpha, *beta, (real *)mult, 
                           buf_d1, buf_d2, buf_d3, *n);
    
    CUDA_CHECK(cudaGetLastError());
    reduce_kernel<real><<<1, 1024>>>(buf_d1, nb);
    CUDA_CHECK(cudaGetLastError());
    reduce_kernel<real><<<1, 1024>>>(buf_d2, nb);
    CUDA_CHECK(cudaGetLastError());
    reduce_kernel<real><<<1, 1024>>>(buf_d3, nb);
    CUDA_CHECK(cudaGetLastError());
    
    CUDA_CHECK(cudaMemcpy(&reduction[0], buf_d1, sizeof(real),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpy(&reduction[1], buf_d2, sizeof(real),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpy(&reduction[2], buf_d3, sizeof(real),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaGetLastError());

  }
}

