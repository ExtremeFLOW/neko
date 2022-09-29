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

#include "cpr_kernel.h"
#include <device/device_config.h>
#include <device/cuda/check.h>

extern "C" {
  
  /**
   * Fortran wrapper that  
   * Weighted inner product \f$ a^T b c \f$
   */
  void cuda_glsc3_elem(void *res, void *a, void *b, void *c,int *lx, int *nelv) {
    const dim3 nthrds((*lx)*(*lx)*(*lx), 1, 1);
    const dim3 nblcks((*nelv), 1, 1);
    const int nb = (*nelv);
    //const int nxyz = (*lx)*(*lx)*(*lx);
    
    real * buf = (real *) malloc(nb * sizeof(real));
    real * buf_d;
    CUDA_CHECK(cudaMalloc(&buf_d, nb*sizeof(real)));
     
    glsc3_elem_kernel<real, 512><<<nblcks, nthrds>>>((real *) a, (real *) b,
                                             (real *) c, buf_d);
    CUDA_CHECK(cudaGetLastError());

    //CUDA_CHECK(cudaMemcpy(buf, buf_d, nb * sizeof(real),
    //                      cudaMemcpyDeviceToHost));

    CUDA_CHECK(cudaMemcpyAsync(res, buf_d, nb * sizeof(real),cudaMemcpyDeviceToDevice));    

    free(buf);
    CUDA_CHECK(cudaFree(buf_d));

  }
  



}
