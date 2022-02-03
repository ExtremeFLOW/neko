/*
 Copyright (c) 2022, The Neko Authors
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

#include <stdio.h>
#include "gmres_kernel.h"
#include <device/device_config.h>
#include <device/cuda/check.h>

extern "C" {

  real * gmres_bf1;
  real * gmres_bfd1;
  
  real cuda_gmres_part2(void *w, void *v, void *h, void * mult, int *j, int *n) {
        
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    
    if (!gmres_bf1){
      gmres_bf1 = (real *) malloc(nb * sizeof(real));
      CUDA_CHECK(cudaMalloc(&gmres_bfd1, nb*sizeof(real)));
    }
     
    gmres_part2_kernel<real><<<nblcks, nthrds>>>((real *) w, (real **) v,
                                                 (real *) mult,(real *) h,
                                                 gmres_bfd1, *j, *n);
    CUDA_CHECK(cudaGetLastError());

    CUDA_CHECK(cudaMemcpy(gmres_bf1, gmres_bfd1, nb * sizeof(real),
                          cudaMemcpyDeviceToHost));

    real res1 = 0.0;
    for (int i = 0; i < nb; i++) {
      res1 += gmres_bf1[i];
    }

    return res1;
  }
}
