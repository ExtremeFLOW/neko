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

#include <device/device_config.h>
#include <device/cuda/check.h>

#include "sumab_kernel.h"
#include "makeabf_kernel.h"
#include "makebdf_kernel.h"

extern "C" {
  void fluid_sumab_cuda(void *u, void *v, void *w,
                        void *uu, void *vv, void *ww,
                        void *ulag1, void *ulag2, void *vlag1,
                        void *vlag2, void *wlag1, void *wlag2,
                        real *ab1, real *ab2, real *ab3, int *nab, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    sumab_kernel<real>
      <<<nblcks, nthrds>>>((real *) u, (real *) v, (real *) w,
                           (real *) uu, (real *) vv, (real *) ww,
                           (real *) ulag1, (real *) ulag2, (real *) vlag1,
                           (real *) vlag2, (real *) wlag1, (real *) wlag2,
                           *ab1, *ab2, *ab3, *nab, *n);
    CUDA_CHECK(cudaGetLastError());
  }

  void fluid_makeabf_cuda(void *ta1, void *ta2, void *ta3,
                          void *abx1, void *aby1, void *abz1, 
                          void *abx2, void *aby2, void *abz2,
                          void *bfx, void *bfy, void *bfz,
                          real *rho, real *ab1, real *ab2, real *ab3, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    
    makeabf_kernel<real>
      <<<nblcks, nthrds>>>((real *) ta1, (real *) ta2, (real *) ta3,
                           (real *) abx1, (real *) aby1, (real *) abz1, 
                           (real *) abx2, (real *) aby2, (real *) abz2,
                           (real *) bfx, (real *) bfy, (real *) bfz,
                           *rho, *ab1, *ab2, *ab3, *n);      
      CUDA_CHECK(cudaGetLastError());
  }
  
  void fluid_makebdf_cuda(void *ta1, void *ta2, void *ta3,
                          void *tb1, void *tb2, void *tb3,
                          void *ulag1, void *ulag2, void *vlag1,
                          void *vlag2, void *wlag1, void *wlag2, 
                          void *bfx, void *bfy, void *bfz,
                          void *u, void *v, void *w, void *B, 
                          real *rho, real *dt, real *bd2,
                          real *bd3, real *bd4, int *nbd, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    makebdf_kernel<real>
      <<<nblcks, nthrds>>>((real *) ta1, (real *) ta2, (real *) ta3,
                           (real *) tb1, (real *) tb2, (real *) tb3,
                           (real *) ulag1, (real *) ulag2, (real *) vlag1,
                           (real *) vlag2, (real *) wlag1, (real *) wlag2, 
                           (real *) bfx, (real *) bfy, (real *) bfz,
                           (real *) u, (real *) v, (real *) w, (real *) B, 
                           *rho, *dt, *bd2, *bd3, *bd4, *nbd,  *n);
    CUDA_CHECK(cudaGetLastError());
  }
}

