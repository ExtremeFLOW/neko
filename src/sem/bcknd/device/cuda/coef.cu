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
#include "coef_kernel.h"
#include <device/device_config.h>
#include <device/cuda/check.h>

extern "C" {
  
  /** 
   * Fortran wrapper for generating geometric factors
   */
  void cuda_coef_generate_geo(void *G11, void *G12, void *G13, 
                              void *G22, void *G23, void *G33, 
                              void *drdx, void *drdy, void *drdz,
                              void *dsdx, void *dsdy, void *dsdz, 
                              void *dtdx, void *dtdy, void *dtdz, 
                              void *jacinv, void *w3, int *nel, 
                              int *lx, int *gdim) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);

#define GEO_CASE(LX)                                                            \
    case LX:                                                                    \
      coef_generate_geo_kernel<real, LX, 1024>                                  \
        <<<nblcks, nthrds>>>((real *) G11, (real *) G12, (real *) G13,          \
                             (real *) G22, (real *) G23, (real *) G33,          \
                             (real *) drdx, (real *) drdy, (real *) drdz,       \
                             (real *) dsdx, (real *) dsdy, (real *) dsdz,       \
                             (real *) dtdx, (real *) dtdy, (real *) dtdz,       \
                             (real *) jacinv, (real *) w3, *gdim);              \
      CUDA_CHECK(cudaGetLastError());                                           \
      break
    
    switch(*lx) {
      GEO_CASE(2);
      GEO_CASE(3);
      GEO_CASE(4);
      GEO_CASE(5);
      GEO_CASE(6);
      GEO_CASE(7);
      GEO_CASE(8);
      GEO_CASE(9);
      GEO_CASE(10);
      GEO_CASE(11);
      GEO_CASE(12);
      GEO_CASE(13);
      GEO_CASE(14);
      GEO_CASE(15);
      GEO_CASE(16);
    default:
      {
        fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
        exit(1);
      }
    }
  }

  /** 
   * Fortran wrapper for generating geometric factors
   */
  void cuda_coef_generate_dxyzdrst(void *drdx, void *drdy, void *drdz, 
				   void *dsdx, void *dsdy, void *dsdz, 
				   void *dtdx, void *dtdy, void *dtdz, 
				   void *dxdr, void *dydr, void *dzdr, 
				   void *dxds, void *dyds, void *dzds, 
				   void *dxdt, void *dydt, void *dzdt,
				   void *dx, void *dy, void *dz, 
				   void *x, void *y, void *z,
				   void *jacinv, void *jac,
				   int *lx, int *nel)  {

    const int n = (*nel) * (*lx) * (*lx) * (*lx);
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks_dxyz((*nel), 1, 1);
    const dim3 nblcks_drst((n + 1024 - 1)/ 1024, 1, 1);

#define DXYZDRST_CASE(LX)					               \
    case LX:								       \
      coef_generate_dxyz_kernel<real, LX, 1024>                                \
	<<<nblcks_dxyz, nthrds>>>((real *) dxdr, (real *) dydr, (real *) dzdr, \
				  (real *) dxds, (real *) dyds, (real *) dzds, \
				  (real *) dxdt, (real *) dydt, (real *) dzdt, \
				  (real *) dx, (real *) dy, (real *) dz,       \
				  (real *) x, (real *) y, (real *) z);	       \
      CUDA_CHECK(cudaGetLastError());					       \
      break

    switch(*lx) {
      DXYZDRST_CASE(2);
      DXYZDRST_CASE(3);
      DXYZDRST_CASE(4);
      DXYZDRST_CASE(5);
      DXYZDRST_CASE(6);
      DXYZDRST_CASE(7);
      DXYZDRST_CASE(8);
      DXYZDRST_CASE(9);
      DXYZDRST_CASE(10);
      DXYZDRST_CASE(11);
      DXYZDRST_CASE(12);
      DXYZDRST_CASE(13);
      DXYZDRST_CASE(14);
      DXYZDRST_CASE(15);
      DXYZDRST_CASE(16);
    default:
      {
        fprintf(stderr, __FILE__ ": size not supported: %d\n", *lx);
        exit(1);
      }
    }

    coef_generate_drst_kernel<real>
      <<<nblcks_drst, nthrds>>>((real *) jac, (real *) jacinv, 
                                (real *) drdx, (real *) drdy, (real *) drdz,
                                (real *) dsdx, (real *) dsdy, (real *) dsdz,
                                (real *) dtdx, (real *) dtdy, (real *) dtdz,
                                (real *) dxdr, (real *) dydr, (real *) dzdr, 
                                (real *) dxds, (real *) dyds, (real *) dzds, 
                                (real *) dxdt, (real *) dydt, (real *) dzdt, n);
    CUDA_CHECK(cudaGetLastError());

  }
}




