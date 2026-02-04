/*
 Copyright (c) 2025, The Neko Authors
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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <device/device_config.h>
#include <device/cuda/check.h>
#include "deardorff_nut_kernel.h"

extern "C" {
  #include <common/neko_log.h>
}

extern "C" {
  void cuda_deardorff_nut_compute(void *TKE, 
                             void *dTdx, void *dTdy, void *dTdz,
                             void *a11, void *a12, void *a13,
                             void *a21, void *a22, void *a23,
                             void *a31, void *a32, void *a33, 
                             void *delta, void *nut, void *temperature_alphat,
                             void *TKE_alphat, void *TKE_source,
                             real *c_k, real *T0, 
                             real *g1, real *g2, real *g3, real *eps, int * n){

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

    deardorff_nut_compute<real>
    <<<nblcks, nthrds, 0, stream>>>((real *) TKE, 
                                    (real *) dTdx, (real *) dTdy, (real *) dTdz,
                                    (real *) a11, (real *) a12, (real *) a13,
                                    (real *) a21, (real *) a22, (real *) a23,
                                    (real *) a31, (real *) a32, (real *) a33,
                                    (real *) delta, (real *) nut, 
                                    (real *) temperature_alphat,
                                    (real *) TKE_alphat, (real *) TKE_source,
                                    *c_k, *T0, *g1, *g2, *g3, *eps, *n);
    CUDA_CHECK(cudaGetLastError());
  }
}