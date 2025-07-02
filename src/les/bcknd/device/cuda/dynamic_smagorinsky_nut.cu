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
#include "dynamic_smagorinsky_nut_kernel.h"

extern "C" {
  #include <common/neko_log.h>
}

extern "C" {
  void cuda_s_abs_compute(void *s_abs, void *s11, void *s22, void *s33,
                         void *s12, void *s13, void *s23,
                         int * n){
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

    s_abs_compute<real>
    <<<nblcks, nthrds, 0, stream>>>((real *) s_abs,
                                    (real *) s11, (real *) s22, (real *) s33,
                                    (real *) s12, (real *) s13, (real *) s23,
                                    * n);
    CUDA_CHECK(cudaGetLastError());
  }
  
  void cuda_lij_compute_part1(void *l11, void *l22, void *l33,
                             void *l12, void *l13, void *l23,
                             void *u, void *v, void *w,
                             void *fu, void *fv, void *fw,
                             void *fuu, void *fvv, void *fww,
                             void *fuv, void *fuw, void *fvw,
                             int * n){
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

    lij_compute_part1<real>
    <<<nblcks, nthrds, 0, stream>>>((real *) l11, (real *) l22, (real *) l33,
                                    (real *) l12, (real *) l13, (real *) l23,
                                    (real *) u, (real *) v, (real *) w,
                                    (real *) fu, (real *) fv, (real *) fw,
                                    (real *) fuu, (real *) fvv, (real *) fww,
                                    (real *) fuv, (real *) fuw, (real *) fvw,
                                    * n);
    CUDA_CHECK(cudaGetLastError());
  }

  void cuda_lij_compute_part2(void *l11, void *l22, void *l33,
                             void *l12, void *l13, void *l23,
                             void *fuu, void *fvv, void *fww,
                             void *fuv, void *fuw, void *fvw,
                             int * n){
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;
    
    lij_compute_part2<real>
    <<<nblcks, nthrds, 0, stream>>>((real *) l11, (real *) l22, (real *) l33,
                                    (real *) l12, (real *) l13, (real *) l23,
                                    (real *) fuu, (real *) fvv, (real *) fww,
                                    (real *) fuv, (real *) fuw, (real *) fvw,
                                    * n);
    CUDA_CHECK(cudaGetLastError());
  }

  void cuda_mij_compute_part1(void *m11, void *m22, void *m33,
                              void *m12, void *m13, void *m23,
                              void *s_abs, void *s11, void *s22, void *s33,
                              void *s12, void *s13, void *s23,
                              void *fs_abs, void *fs11, void *fs22, void *fs33,
                              void *fs12, void *fs13, void *fs23,
                              void *fsabss11, void *fsabss22, void *fsabss33,
                              void *fsabss12, void *fsabss13, void *fsabss23,
                              real * delta_ratio2, int * n){
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

    mij_compute_part1<real>
    <<<nblcks, nthrds, 0, stream>>>((real *) m11, (real *) m22, (real *) m33,
                                    (real *) m12, (real *) m13, (real *) m23,
                                    (real *) s_abs, (real *) s11, (real *) s22, 
                                    (real *) s33, (real *) s12, (real *) s13, (real *) s23,
                                    (real *) fs_abs, (real *) fs11, (real *) fs22, 
                                    (real *) fs33, (real *) fs12, (real *) fs13, (real *) fs23,
                                    (real *) fsabss11, (real *) fsabss22, (real *) fsabss33,
                                    (real *) fsabss12, (real *) fsabss13, (real *) fsabss23,
                                    * delta_ratio2, * n);
    CUDA_CHECK(cudaGetLastError());
  }

  void cuda_mij_nut_compute_part2(void *m11, void *m22, void *m33,
                             void *m12, void *m13, void *m23,
                             void *l11, void *l22, void *l33,
                             void *l12, void *l13, void *l23,
                             void *fsabss11, void *fsabss22, void *fsabss33,
                             void *fsabss12, void *fsabss13, void *fsabss23,
                             void *num, void *den, void *c_dyn, void *delta,
                             void *s_abs, void *nut,
                             real * alpha, int * n){
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

    mij_nut_compute_part2<real>
    <<<nblcks, nthrds, 0, stream>>>((real *) m11, (real *) m22, (real *) m33,
                                    (real *) m12, (real *) m13, (real *) m23,
                                    (real *) l11, (real *) l22, (real *) l33,
                                    (real *) l12, (real *) l13, (real *) l23,
                                    (real *) fsabss11, (real *) fsabss22, (real *) fsabss33,
                                    (real *) fsabss12, (real *) fsabss13, (real *) fsabss23,
                                    (real *) num, (real *) den, (real *) c_dyn,
                                    (real *) delta, (real *) s_abs, (real *) nut,
                                    * alpha, * n);
    CUDA_CHECK(cudaGetLastError());
  }
}