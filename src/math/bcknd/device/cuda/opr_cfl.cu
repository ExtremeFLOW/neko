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
#include "cfl_kernel.h"
#include <device/device_config.h>
#include <device/cuda/check.h>


extern "C" {

#include <math/bcknd/device/device_mpi_reduce.h>
#include <math/bcknd/device/device_mpi_op.h>

  /**
   * @todo Make sure that this gets deleted at some point...
   */
  real *cfl_d = NULL;
  
  /** 
   * Fortran wrapper for device cuda CFL
   */
  real cuda_cfl(real *dt, void *u, void *v, void *w,
                void *drdx, void *dsdx, void *dtdx,
                void *drdy, void *dsdy, void *dtdy,
                void *drdz, void *dsdz, void *dtdz,
                void *dr_inv, void *ds_inv, void *dt_inv,
                void *jacinv, int *nel, int *lx) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;      

    if (cfl_d == NULL) {
      CUDA_CHECK(cudaMalloc(&cfl_d, (*nel) * sizeof(real)));
    }



#define CASE(LX)                                                                \
    case LX:                                                                    \
      cfl_kernel<real, LX, 1024>                                                \
        <<<nblcks, nthrds, 0, stream>>>                                         \
        (*dt, (real *) u, (real *) v, (real *) w,                               \
         (real *) drdx, (real *) dsdx, (real *) dtdx,                           \
         (real *) drdy, (real *) dsdy, (real *) dtdy,                           \
         (real *) drdz, (real *) dsdz, (real *) dtdz,                           \
         (real *) dr_inv, (real *) ds_inv, (real *) dt_inv,                     \
         (real *) jacinv, (real *) cfl_d);                                      \
      CUDA_CHECK(cudaGetLastError());                                           \
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

    cfl_reduce_kernel<real><<<1, 1024>>> (cfl_d, (*nel));
    CUDA_CHECK(cudaGetLastError());

    real cfl;
#ifdef HAVE_DEVICE_MPI
    cudaStreamSynchronize(stream);
    device_mpi_allreduce(cfl_d, &cfl, 1, sizeof(real), DEVICE_MPI_MAX);
#else
    CUDA_CHECK(cudaMemcpyAsync(&cfl, cfl_d, sizeof(real),
                               cudaMemcpyDeviceToHost, stream));
#endif
    
    return cfl;
  } 
}
