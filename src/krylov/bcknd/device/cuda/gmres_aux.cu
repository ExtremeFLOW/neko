/*
 Copyright (c) 2022-2025, The Neko Authors
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
#include "gmres_kernel.h"
#include <device/device_config.h>
#include <device/cuda/check.h>

#ifdef HAVE_NVSHMEM
#include <nvshmem.h>
#include <nvshmemx.h>
#endif

extern "C" {

#include <math/bcknd/device/device_mpi_reduce.h>
#include <math/bcknd/device/device_mpi_op.h>

#ifdef HAVE_NCCL
#include <math/bcknd/device/device_nccl_reduce.h>
#include <math/bcknd/device/device_nccl_op.h>
#endif

  /**
   * @todo Make sure that this gets deleted at some point...
   */
  real * gmres_bf1 = NULL;
  real * gmres_bfd1 = NULL;
  int gmres_bf_len = 0;

  real cuda_gmres_part2(void *w, void *v, void *h, void * mult, int *j, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;

    if (gmres_bf1 != NULL && gmres_bf_len < nb) {
      CUDA_CHECK(cudaFreeHost(gmres_bf1));
#ifdef HAVE_NVSHMEM
      nvshmem_free(gmres_bfd1);
#else
      CUDA_CHECK(cudaFree(gmres_bfd1));
#endif
      gmres_bf1 = NULL;
    }

    if (gmres_bf1 == NULL){
      CUDA_CHECK(cudaMallocHost(&gmres_bf1, nb*sizeof(real)));
#ifdef HAVE_NVSHMEM
      gmres_bfd1 = (real *) nvshmem_malloc(nb*sizeof(real));
#else
      CUDA_CHECK(cudaMalloc(&gmres_bfd1, nb*sizeof(real)));
#endif
      gmres_bf_len = nb;
    }

    gmres_part2_kernel<real>
      <<<nblcks, nthrds, 0, stream>>>((real *) w, (real **) v,
                                      (real *) mult,(real *) h,
                                      (real *) gmres_bfd1, *j, *n);
    CUDA_CHECK(cudaGetLastError());
    reduce_kernel<real><<<1, 1024, 0, stream>>>((real *) gmres_bfd1, nb);
    CUDA_CHECK(cudaGetLastError());

#ifdef HAVE_NCCL
    device_nccl_allreduce(gmres_bfd1, gmres_bfd1, 1, sizeof(real),
                          DEVICE_NCCL_SUM, stream);
    CUDA_CHECK(cudaMemcpyAsync(gmres_bf1, gmres_bfd1, sizeof(real),
                               cudaMemcpyDeviceToHost, stream));
    cudaStreamSynchronize(stream);
#elif HAVE_NVSHMEM
    if (sizeof(real) == sizeof(float)) {
      nvshmemx_float_sum_reduce_on_stream(NVSHMEM_TEAM_WORLD,
                                           (float *) gmres_bfd1,
                                           (float *) gmres_bfd1, 1, stream);
    }
    else if (sizeof(real) == sizeof(double)) {
      nvshmemx_double_sum_reduce_on_stream(NVSHMEM_TEAM_WORLD,
                                           (double *) gmres_bfd1,
                                           (double *) gmres_bfd1, 1, stream);
    }
    CUDA_CHECK(cudaMemcpyAsync(gmres_bf1, gmres_bfd1, sizeof(real),
                               cudaMemcpyDeviceToHost, stream));
    cudaStreamSynchronize(stream);
#elif HAVE_DEVICE_MPI
    cudaStreamSynchronize(stream);
    device_mpi_allreduce(gmres_bfd1, gmres_bf1, 1, sizeof(real), DEVICE_MPI_SUM);
#else

    CUDA_CHECK(cudaMemcpyAsync(gmres_bf1, gmres_bfd1, sizeof(real),
                               cudaMemcpyDeviceToHost, stream));
    cudaStreamSynchronize(stream);
#endif
    return gmres_bf1[0];
  }
}
