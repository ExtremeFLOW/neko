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

#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>

#include "sumab_kernel.h"
#include "makeext_kernel.h"
#include "makebdf_kernel.h"

extern "C" {
  void rhs_maker_sumab_hip(void *u, void *v, void *w,
                       void *uu, void *vv, void *ww,
                       void *ulag1, void *ulag2, void *vlag1,
                       void *vlag2, void *wlag1, void *wlag2,
                       real *ab1, real *ab2, real *ab3, int *nab, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);


    hipLaunchKernelGGL(HIP_KERNEL_NAME( sumab_kernel<real> ),
                       nblcks, nthrds, 0, 0,
                       (real *) u, (real *) v, (real *) w,
                       (real *) uu, (real *) vv, (real *) ww,
                       (real *) ulag1, (real *) ulag2, (real *) vlag1,
                       (real *) vlag2, (real *) wlag1, (real *) wlag2,
                       *ab1, *ab2, *ab3, *nab, *n);
    HIP_CHECK(hipGetLastError());
  }

  void rhs_maker_ext_hip(void *abx1, void *aby1, void *abz1, 
                         void *abx2, void *aby2, void *abz2,
                         void *bfx, void *bfy, void *bfz,
                         real *rho, real *ab1, real *ab2, real *ab3, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME( makeext_kernel<real> ),
                       nblcks, nthrds, 0, 0,
                       (real *) abx1, (real *) aby1, (real *) abz1, 
                       (real *) abx2, (real *) aby2, (real *) abz2,
                       (real *) bfx, (real *) bfy, (real *) bfz,
                       *rho, *ab1, *ab2, *ab3, *n);      
    HIP_CHECK(hipGetLastError());
  }

  void scalar_rhs_maker_ext_hip(void *fs_lag, void *fs_laglag, void *fs, 
                                real *rho, real *ext1, real *ext2, real *ext3,
                                int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME( scalar_makeext_kernel<real> ),
                       nblcks, nthrds, 0, 0,
                       (real *) fs_lag, (real *) fs_laglag, (real *) fs, 
                       (real *) fs, *rho, *ext1, *ext2, *ext3, *n);      
    HIP_CHECK(hipGetLastError());
  }
  
  void rhs_maker_bdf_hip(void *ulag1, void *ulag2, void *vlag1,
                         void *vlag2, void *wlag1, void *wlag2, 
                         void *bfx, void *bfy, void *bfz,
                         void *u, void *v, void *w, void *B, 
                         real *rho, real *dt, real *bd2,
                         real *bd3, real *bd4, int *nbd, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME( makebdf_kernel<real> ),
                       nblcks, nthrds, 0, 0,
                       (real *) ulag1, (real *) ulag2, (real *) vlag1,
                       (real *) vlag2, (real *) wlag1, (real *) wlag2, 
                       (real *) bfx, (real *) bfy, (real *) bfz,
                       (real *) u, (real *) v, (real *) w, (real *) B, 
                       *rho, *dt, *bd2, *bd3, *bd4, *nbd,  *n);
    HIP_CHECK(hipGetLastError());
  }

  void scalar_rhs_maker_bdf_hip(void *s_lag, void *s_laglag, void *fs,
                                void *s, void *B, real *rho, real *dt,
                                real *bd2, real *bd3, real *bd4,
                                int *nbd, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME( scalar_makebdf_kernel<real> ),
                       nblcks, nthrds, 0, 0,
                       (real *) s_lag, (real *) s_laglag, (real *) fs,
                       (real *) s, (real *) B, *rho, *dt, *bd2, *bd3,
                       *bd4, *nbd,  *n);
    HIP_CHECK(hipGetLastError());
  }
}

