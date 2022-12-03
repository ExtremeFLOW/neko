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
#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "cfl_kernel.h"

extern "C" {

  /**
   * @todo Make sure that this gets deleted at some point...
   */
  real *cfl_d = NULL;
  
  /** 
   * Fortran wrapper for device hip CFL
   */
  real hip_cfl(real *dt, void *u, void *v, void *w,
	       void *drdx, void *dsdx, void *dtdx,
	       void *drdy, void *dsdy, void *dtdy,
	       void *drdz, void *dsdz, void *dtdz,
	       void *dr_inv, void *ds_inv, void *dt_inv,
	       void *jacinv, int *nel, int *lx) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);
  
    if (cfl_d == NULL) {
      HIP_CHECK(hipMalloc(&cfl_d, (*nel) * sizeof(real)));
    }

#define CASE(LX)                                                                \
    case LX:                                                                    \
        hipLaunchKernelGGL(HIP_KERNEL_NAME( cfl_kernel<real, LX, 1024> ),       \
                           nblcks, nthrds, 0, 0,                                \
                           *dt, (real *) u, (real *) v, (real *) w,             \
                           (real *) drdx, (real *) dsdx, (real *) dtdx,         \
                           (real *) drdy, (real *) dsdy, (real *) dtdy,         \
                           (real *) drdz, (real *) dsdz, (real *) dtdz,         \
                           (real *) dr_inv, (real *) ds_inv, (real *) dt_inv,   \
                           (real *) jacinv, (real *) cfl_d);                    \
      HIP_CHECK(hipGetLastError());                                             \
      break
      
    switch(*lx) {
      CASE(2);
      CASE(3);
      CASE(4);
      CASE(5);
      CASE(6);
      CASE(7);
      CASE(8);
      CASE(9);
      CASE(10);
    default:
      {
        fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
        exit(1);
      }
    }


    hipLaunchKernelGGL(HIP_KERNEL_NAME(cfl_reduce_kernel<real>),
                       1, 1024, 0, 0,
                       cfl_d, (*nel));
    HIP_CHECK(hipGetLastError());

    real cfl = 0.0;
    HIP_CHECK(hipMemcpy(&cfl, cfl_d, sizeof(real),
			hipMemcpyDeviceToHost));
    


    return cfl;
  } 
}
