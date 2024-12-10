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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <device/device_config.h>
#include <device/cuda/check.h>
#include "smagorinsky_nut_kernel.h"

extern "C" {
  #include <common/neko_log.h>
}

extern "C" {
  void cuda_smagorinsky_nut_compute(void *s11, void *s22, void *s33,
                             void *s12, void *s13, void *s23,
                             void *delta, void *nut, void *mult, 
                             real *c_s, int * n){
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

    smagorinsky_nut_compute<real>
    <<<nblcks, nthrds, 0, stream>>>((real *) s11, (real *) s22, (real *) s33,
                                    (real *) s12, (real *) s13, (real *) s23,
                                    (real *) delta, (real *) nut, (real *) mult, 
                                    *c_s, * n);
    CUDA_CHECK(cudaGetLastError());
  }
}