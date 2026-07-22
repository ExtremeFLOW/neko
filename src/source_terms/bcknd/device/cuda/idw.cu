/*
 Copyright (c) 2024-2026, The Neko Authors
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
#include "idw_kernel.h"
#include <stdio.h>
#include <stdlib.h>

extern "C" {

  /**
   * Fortran wrapper for the atomic-free IDW gather kernel.
   */
  void cuda_idw_gather_one_sided(void *fu, void *fv, void *fw,
                                 void *fu_ib, void *fv_ib, void *fw_ib,
                                 void *fum_ib, void *fvm_ib, void *fwm_ib,
                                 void *x, void *y, void *z, void *ds,
                                 void *pmsk, void *w, void *wm,
                                 void *lpx, void *lpy, void *lpz,
                                 void *active_el, void *el_off, void *el_lag,
                                 int *n_active, int *lx3, real *dt,
                                 real *rmax, real *pwr, real *eps, real *wtol) {

    const int total = (*n_active) * (*lx3);
    if (total < 1)
      return;

    const dim3 nthrds(256, 1, 1);
    const dim3 nblcks((total + 256 - 1) / 256, 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

    idw_gather_one_sided_kernel<real>
      <<<nblcks, nthrds, 0, stream>>>((real *) fu, (real *) fv, (real *) fw,
                                      (real *) fu_ib, (real *) fv_ib,
                                      (real *) fw_ib, (real *) fum_ib,
                                      (real *) fvm_ib, (real *) fwm_ib,
                                      (real *) x, (real *) y, (real *) z,
                                      (real *) ds, (real *) pmsk, (real *) w,
                                      (real *) wm, (real *) lpx, (real *) lpy,
                                      (real *) lpz, (int *) active_el,
                                      (int *) el_off, (int *) el_lag,
                                      *n_active, *lx3,
                                      *dt, *rmax, *pwr, *eps, *wtol);
    CUDA_CHECK(cudaGetLastError());
  }

}
