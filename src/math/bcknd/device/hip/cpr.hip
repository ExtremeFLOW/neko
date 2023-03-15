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

#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "cpr_kernel.h"

extern "C" {

#include <math/bcknd/device/device_mpi_reduce.h>
  
  /*
   * Reduction buffer
   */
  int red_a = 0;
  real * buf = NULL;
  real * buf_d = NULL;
   
  /**
   * Fortran wrapper lcsc3
   * Weighted inner product \f$ a^T b c \f$
   */
  void hip_lcsc3(void *res, void *a, void *b, void *c,int *lx, int *nelv) {
    
    const dim3 nthrds((*lx)*(*lx)*(*lx), 1, 1);
    const dim3 nblcks((*nelv), 1, 1);
    const int nb = (*nelv);
         
    if ( nb > red_a){
      if (buf != NULL) {
        HIP_CHECK(hipHostFree(buf));
        HIP_CHECK(hipFree(buf_d));        
      }
      HIP_CHECK(hipHostMalloc(&buf,nb*sizeof(real),hipHostMallocDefault));
      HIP_CHECK(hipMalloc(&buf_d, nb*sizeof(real)));
    }
     
    hipLaunchKernelGGL(HIP_KERNEL_NAME(lcsc3_kernel<real,1024>),
                       nblcks, nthrds, 0, 0,
                       buf_d, (real *) a, (real *) b,
                       (real *) c);
    HIP_CHECK(hipGetLastError());

    HIP_CHECK(hipMemcpy(res, buf_d, nb * sizeof(real),
                        hipMemcpyDeviceToDevice));
  }


  /**
   * Fortran wrapper
   * local sum
   */
  void hip_lcsum(void *res, void *a, int *lx, int *nelv) {
    
    const dim3 nthrds((*lx)*(*lx)*(*lx), 1, 1);
    const dim3 nblcks((*nelv), 1, 1);
    const int nb = (*nelv);
         
    if ( nb > red_a){
      if (buf != NULL) {
        HIP_CHECK(hipHostFree(buf));
        HIP_CHECK(hipFree(buf_d));        
      }
      HIP_CHECK(hipHostMalloc(&buf,nb*sizeof(real),hipHostMallocDefault));
      HIP_CHECK(hipMalloc(&buf_d, nb*sizeof(real)));
    }
     
    hipLaunchKernelGGL(HIP_KERNEL_NAME(lcsum_kernel<real,1024>),
                       nblcks, nthrds, 0, 0,
                       buf_d, (real *) a);
    HIP_CHECK(hipGetLastError());

    HIP_CHECK(hipMemcpy(res, buf_d, nb * sizeof(real),
                        hipMemcpyDeviceToDevice));
  }


  /**
   * Fortran wrapper
   * local min
   */
  void hip_lcmin(void *res, void *a, int *lx, int *nelv) {
    
    const dim3 nthrds((*lx)*(*lx)*(*lx), 1, 1);
    const dim3 nblcks((*nelv), 1, 1);
    const int nb = (*nelv);
         
    if ( nb > red_a){
      if (buf != NULL) {
        HIP_CHECK(hipHostFree(buf));
        HIP_CHECK(hipFree(buf_d));        
      }
      HIP_CHECK(hipHostMalloc(&buf,nb*sizeof(real),hipHostMallocDefault));
      HIP_CHECK(hipMalloc(&buf_d, nb*sizeof(real)));
    }
     
    hipLaunchKernelGGL(HIP_KERNEL_NAME(lcmin_kernel<real,1024>),
                       nblcks, nthrds, 0, 0,
                       buf_d, (real *) a);
    HIP_CHECK(hipGetLastError());

    HIP_CHECK(hipMemcpy(res, buf_d, nb * sizeof(real),
                        hipMemcpyDeviceToDevice));
  }


  /**
   * Fortran wrapper
   * local max
   */
  void hip_lcmax(void *res, void *a, int *lx, int *nelv) {
    
    const dim3 nthrds((*lx)*(*lx)*(*lx), 1, 1);
    const dim3 nblcks((*nelv), 1, 1);
    const int nb = (*nelv);
         
    if ( nb > red_a){
      if (buf != NULL) {
        HIP_CHECK(hipHostFree(buf));
        HIP_CHECK(hipFree(buf_d));        
      }
      HIP_CHECK(hipHostMalloc(&buf,nb*sizeof(real),hipHostMallocDefault));
      HIP_CHECK(hipMalloc(&buf_d, nb*sizeof(real)));
    }
     
    hipLaunchKernelGGL(HIP_KERNEL_NAME(lcmax_kernel<real,1024>),
                       nblcks, nthrds, 0, 0,
                       buf_d, (real *) a);
    HIP_CHECK(hipGetLastError());

    HIP_CHECK(hipMemcpy(res, buf_d, nb * sizeof(real),
                        hipMemcpyDeviceToDevice));
  }
  

  /**
   * Fortran wrapper 
   * local sort with absolute values
   */
  void hip_lcsort_abs(void *asort, void *keysort, void *a, void *key,int *lx, int *nelv) {
    
    const dim3 nthrds((*lx)*(*lx)*(*lx), 1, 1);
    const dim3 nblcks((*nelv), 1, 1);
          
    hipLaunchKernelGGL(HIP_KERNEL_NAME(lcsort_abs_kernel<real,1024>),
                       nblcks, nthrds, 0, 0,
                       (real *) asort, (int *) keysort,
                       (real *) a, (int *) key);
    HIP_CHECK(hipGetLastError());

  }


  /**
   * Fortran wrapper 
   * local sort by key
   */
  void hip_lcsort_bykey(void *asort, void *keysort, void *a, void *key,int *lx, int *nelv) {
    
    const dim3 nthrds((*lx)*(*lx)*(*lx), 1, 1);
    const dim3 nblcks((*nelv), 1, 1);
          
    hipLaunchKernelGGL(HIP_KERNEL_NAME(lcsort_bykey_kernel<real,1024>),
                       nblcks, nthrds, 0, 0,
                       (real *) asort, (int *) keysort,
                       (real *) a, (int *) key);
    HIP_CHECK(hipGetLastError());

  }


  /**
   * Fortran wrapper 
   * lctnsr3d
   */ 
  void hip_lctnsr3d(void *v, int *nv, void *u, int *nu,
		   void *A, void *Bt, void *Ct,
		   void *bp, void*bp_key,  int *nel) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);
          
    hipLaunchKernelGGL(HIP_KERNEL_NAME(lctnsr3d_kernel<real>),
                       nblcks, nthrds, 0, 0,
                       (real *) v, *nv, (real *) u,
		       *nu, (real *) A, (real *) Bt, (real *) Ct,
		       (real *) bp, (real *) bp_key);
    HIP_CHECK(hipGetLastError());

  }














}



