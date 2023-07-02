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

#include "projection_kernel.h"
#include <math/bcknd/device/cuda/math_kernel.h>


/*
 * Reduction buffer
 */
int proj_red_s = 0;
real * proj_bufred_d = NULL;

extern "C" {

#include <math/bcknd/device/device_mpi_reduce.h>
#include <math/bcknd/device/device_mpi_op.h>

  void cuda_project_on(void *alpha, void * b, void *xx, void *bb, void *mult,
                       void *xbar, int *j, int *n){ 
    
    int pow2 = 1;
    while(pow2 < (*j)){
      pow2 = 2*pow2;
    }
    const int nt = 1024/pow2;   
    const dim3 glsc3_nthrds(pow2, nt, 1);
    const dim3 glsc3_nblcks(((*n)+nt - 1)/nt, 1, 1);
    const int glsc3_nb = ((*n) + nt - 1)/nt;
    if((*j)*glsc3_nb>proj_red_s){
      proj_red_s = (*j)*glsc3_nb;
      if (proj_bufred_d != NULL) {
	CUDA_CHECK(cudaFree(proj_bufred_d));
      }
      CUDA_CHECK(cudaMalloc(&proj_bufred_d, (*j)*glsc3_nb*sizeof(real)));
    }

    /* First glsc3_many call */
    glsc3_many_kernel<real>
      <<<glsc3_nblcks, glsc3_nthrds>>>((const real *) b,
                                       (const real **) xx,
                                       (const real *) mult,
                                       proj_bufred_d, *j, *n);
    CUDA_CHECK(cudaGetLastError());
    glsc3_reduce_kernel<real><<<(*j),1024>>> (proj_bufred_d, glsc3_nb, *j);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpy(alpha, proj_bufred_d, (*j) * sizeof(real),
                          cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemsetAsync(xbar, 0, (*n) * sizeof(real)));

    cudaDeviceSynchronize();
    device_mpi_allreduce_inplace(alpha, (*j), sizeof(real), DEVICE_MPI_SUM);

    const dim3 vec_nthrds(1024, 1, 1);
    const dim3 vec_nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    /* First vector operation block */
    project_on_vec_kernel<real><<<vec_nblcks, vec_nthrds>>>((real *) xbar,
                                                      (const real **) xx,
                                                      (real *) b,
                                                      (const real **) bb,
                                                      (const real *) alpha,
                                                      *j, *n);
    /* Second glsc3_many call */
    glsc3_many_kernel<real>
      <<<glsc3_nblcks, glsc3_nthrds>>>((const real *) b,
                                       (const real **) xx,
                                       (const real *) mult,
                                       proj_bufred_d, *j, *n);
    CUDA_CHECK(cudaGetLastError());
    glsc3_reduce_kernel<real><<<(*j),1024>>> (proj_bufred_d, glsc3_nb, *j);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpy(alpha, proj_bufred_d, (*j) * sizeof(real),
                          cudaMemcpyDeviceToDevice));

    cudaDeviceSynchronize();
    device_mpi_allreduce_inplace(alpha, (*j), sizeof(real), DEVICE_MPI_SUM);

    /* Second vector operation block */
    project_on_vec_kernel<real><<<vec_nblcks, vec_nthrds>>>((real *) xbar,
                                                            (const real **) xx,
                                                            (real *) b,
                                                            (const real **) bb,
                                                            (const real *) alpha,
                                                            *j, *n);
  }

  void cuda_project_ortho(void *alpha, void * b, void *xx, void *bb,
                          void *w, void *xm, int *j, int *n, real *nrm){

    int pow2 = 1;
    while(pow2 < (*j)){
      pow2 = 2*pow2;
    }
    const int nt = 1024/pow2;   
    const dim3 glsc3_nthrds(pow2, nt, 1);
    const dim3 glsc3_nblcks(((*n)+nt - 1)/nt, 1, 1);
    const int glsc3_nb = ((*n) + nt - 1)/nt;
    if((*j)*glsc3_nb>proj_red_s){
      proj_red_s = (*j)*glsc3_nb;
      if (proj_bufred_d != NULL) {
	CUDA_CHECK(cudaFree(proj_bufred_d));
      }
      CUDA_CHECK(cudaMalloc(&proj_bufred_d, (*j)*glsc3_nb*sizeof(real)));
    }

    /* First glsc3_many call */
    glsc3_many_kernel<real>
      <<<glsc3_nblcks, glsc3_nthrds>>>((const real *) b,
                                       (const real **) xx,
                                       (const real *) w,
                                       proj_bufred_d, *j, *n);
    CUDA_CHECK(cudaGetLastError());
    glsc3_reduce_kernel<real><<<(*j),1024>>> (proj_bufred_d, glsc3_nb, *j);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpy(alpha, proj_bufred_d, (*j) * sizeof(real),
                          cudaMemcpyDeviceToDevice));

    cudaDeviceSynchronize();
    device_mpi_allreduce_inplace(alpha, (*j), sizeof(real), DEVICE_MPI_SUM);

    CUDA_CHECK(cudaMemcpy(nrm, (real *) alpha + (*j - 1),
                          sizeof(real), cudaMemcpyDeviceToHost));
    (*nrm) = sqrt(*nrm);


    const dim3 vec_nthrds(1024, 1, 1);
    const dim3 vec_nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    /* First vector operation block */
    project_ortho_vec_kernel<real>
      <<<vec_nblcks, vec_nthrds>>>((real *) xm, (const real **) xx,
                                   (real *) b, (const real **) bb,
                                   (const real *) alpha, *j, *n);

    /* Second glsc3_many call */
    glsc3_many_kernel<real>
      <<<glsc3_nblcks, glsc3_nthrds>>>((const real *) b,
                                       (const real **) xx,
                                       (const real *) w,
                                       proj_bufred_d, *j, *n);
    CUDA_CHECK(cudaGetLastError());
    glsc3_reduce_kernel<real><<<(*j),1024>>> (proj_bufred_d, glsc3_nb, *j);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpy(alpha, proj_bufred_d, (*j) * sizeof(real),
                          cudaMemcpyDeviceToDevice));

    cudaDeviceSynchronize();
    device_mpi_allreduce_inplace(alpha, (*j), sizeof(real), DEVICE_MPI_SUM);

    /* Second vector operation block */
    project_ortho_vec_kernel<real>
      <<<vec_nblcks, vec_nthrds>>>((real *) xm, (const real **) xx,
                                   (real *) b, (const real **) bb,
                                   (const real *) alpha, *j, *n);
    
  }
      
}

