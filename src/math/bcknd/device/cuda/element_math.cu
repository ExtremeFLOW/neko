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

#include "element_math_kernel.h"
#include <device/device_config.h>
#include <device/cuda/check.h>

extern "C" {
  
  /**
   * Fortran wrapper that  
   * Weighted inner product \f$ a^T b c \f$
   */
  void cuda_lcsc3(void *res, void *a, void *b, void *c,int *lx, int *nelv) {
    const dim3 nthrds((*lx)*(*lx)*(*lx), 1, 1);
    const dim3 nblcks((*nelv), 1, 1);
    const int nb = (*nelv);
    //const int nxyz = (*lx)*(*lx)*(*lx);
    
    real * buf = (real *) malloc(nb * sizeof(real));
    real * buf_d;
    CUDA_CHECK(cudaMalloc(&buf_d, nb*sizeof(real)));
     
    lcsc3_kernel<real, 1024><<<nblcks, nthrds>>>(buf_d, (real *) a, 
		                               (real *) b,
                                               (real *) c);
    CUDA_CHECK(cudaGetLastError());

    CUDA_CHECK(cudaMemcpyAsync(res, buf_d, nb * sizeof(real),cudaMemcpyDeviceToDevice));    

    free(buf);
    CUDA_CHECK(cudaFree(buf_d));

  }
  

  /**
   * Fortran wrapper that  
   * local sum$
   */
  void cuda_lcsum(void *res, void *a,int *lx, int *nelv) {
    const dim3 nthrds((*lx)*(*lx)*(*lx), 1, 1);
    const dim3 nblcks((*nelv), 1, 1);
    const int nb = (*nelv);
    //const int nxyz = (*lx)*(*lx)*(*lx);
    
    real * buf = (real *) malloc(nb * sizeof(real));
    real * buf_d;
    CUDA_CHECK(cudaMalloc(&buf_d, nb*sizeof(real)));
     
    lcsum_kernel<real, 1024><<<nblcks, nthrds>>>(buf_d, (real *) a);
    CUDA_CHECK(cudaGetLastError());

    CUDA_CHECK(cudaMemcpyAsync(res, buf_d, nb * sizeof(real),cudaMemcpyDeviceToDevice));    

    free(buf);
    CUDA_CHECK(cudaFree(buf_d));

  }


  /**
   * Fortran wrapper that  
   * local min$
   */
  void cuda_lcmin(void *res, void *a,int *lx, int *nelv) {
    const dim3 nthrds((*lx)*(*lx)*(*lx), 1, 1);
    const dim3 nblcks((*nelv), 1, 1);
    const int nb = (*nelv);
    //const int nxyz = (*lx)*(*lx)*(*lx);
    
    real * buf = (real *) malloc(nb * sizeof(real));
    real * buf_d;
    CUDA_CHECK(cudaMalloc(&buf_d, nb*sizeof(real)));
     
    lcmin_kernel<real, 1024><<<nblcks, nthrds>>>(buf_d, (real *) a);
    CUDA_CHECK(cudaGetLastError());

    CUDA_CHECK(cudaMemcpyAsync(res, buf_d, nb * sizeof(real),cudaMemcpyDeviceToDevice));    

    free(buf);
    CUDA_CHECK(cudaFree(buf_d));

  }


  /**
   * Fortran wrapper that  
   * local max$
   */
  void cuda_lcmax(void *res, void *a,int *lx, int *nelv) {
    const dim3 nthrds((*lx)*(*lx)*(*lx), 1, 1);
    const dim3 nblcks((*nelv), 1, 1);
    const int nb = (*nelv);
    //const int nxyz = (*lx)*(*lx)*(*lx);
    
    real * buf = (real *) malloc(nb * sizeof(real));
    real * buf_d;
    CUDA_CHECK(cudaMalloc(&buf_d, nb*sizeof(real)));
     
    lcmax_kernel<real, 1024><<<nblcks, nthrds>>>(buf_d, (real *) a);
    CUDA_CHECK(cudaGetLastError());

    CUDA_CHECK(cudaMemcpyAsync(res, buf_d, nb * sizeof(real),cudaMemcpyDeviceToDevice));    

    free(buf);
    CUDA_CHECK(cudaFree(buf_d));

  }


  /**
   * Fortran wrapper that  
   * local sort decending order$
   */

  void cuda_lcsort_abs(void *asort,void *keysort, void *a, void *key,int *lx, int *nelv) {
  /*  const dim3 nthrds((*lx)*(*lx)*(*lx), 1, 1);
    const dim3 nblcks((*nelv), 1, 1);
    //const int nb = (*nelv);
    //const int nxyz = (*lx)*(*lx)*(*lx);
     
    lcsort_abs_kernel<real, 1024><<<nblcks, nthrds>>>((real *) asort, (int *) keysort, (real *) a, (int *) key);
  */
    const dim3 nthrds((*lx)*(*lx)*(*lx)/2, 1, 1);
    const dim3 nblcks(4320, 1, 1);
    const int n = *nelv;
    const int d = (*lx)*(*lx)*(*lx);
    //const int nb = (*nelv);
    //const int nxyz = (*lx)*(*lx)*(*lx);
     
    lcsort_abs_kernel<real><<<nblcks, nthrds,d*(sizeof(real)+sizeof(int)),0>>>((real *) asort, (int *) keysort, (real *) a, (int *) key, n, d);
    CUDA_CHECK(cudaGetLastError());

  }

  /**
   * Fortran wrapper that  
   * local sort by key$
   */
  void cuda_lcsort_bykey(void *asort,void *keysort, void *a, void *key,int *lx, int *nelv) {
    const dim3 nthrds((*lx)*(*lx)*(*lx), 1, 1);
    const dim3 nblcks((*nelv), 1, 1);
    //const int nb = (*nelv);
    //const int nxyz = (*lx)*(*lx)*(*lx);
     
    lcsort_bykey_kernel<real, 1024><<<nblcks, nthrds>>>((real *) asort, (int *) keysort, (real *) a, (int *) key);
    CUDA_CHECK(cudaGetLastError());

  }


  /** Fortran wrapper for lctnsr3d **/
  void cuda_lctnsr3d(void *v, int *nv, void *u, int *nu,
		   void *A, void *Bt, void *Ct,
		   void *bp, void*bp_key,  int *nel) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(*nel, 1, 1);
    
    lctnsr3d_kernel<real>
      <<<nblcks, nthrds>>>((real *) v, *nv, (real *) u, 
			   *nu, (real *) A, (real *) Bt, (real *) Ct,
			   (real *) bp, (real *) bp_key);
    CUDA_CHECK(cudaGetLastError());
    
  }



}
