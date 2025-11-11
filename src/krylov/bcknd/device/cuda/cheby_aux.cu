/*
 Copyright (c) 2025, The Neko Authors
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


template< typename T >
__global__ void cheby_part1(T * __restrict__ d,
                            T * __restrict__ x,
                            const T inv_tha,
                            const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    const real dt = d[i] * inv_tha;
    d[i] = dt;
    x[i] = x[i] + dt;    
  }  
}

template< typename T >
__global__ void cheby_part2(T * __restrict__ d,
                            T * __restrict__ w,
                            T * __restrict__ x,
                            const T tmp1,
                            const T tmp2,
                            const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    const real dt = tmp1 * d[i] + tmp2 * w[i];
    d[i] = dt;
    x[i] = x[i] + dt; 
  }  
}

extern "C" {

  /**
   * Fortran wrapper for part1 of cheby
   */
  void cuda_cheby_part1(void *d, void *x,
                       real *inv_tha, int *n, cudaStream_t strm) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    cheby_part1<real><<<nblcks, nthrds, 0, strm>>>((real *) d, (real *) x, *inv_tha, *n);
    CUDA_CHECK(cudaGetLastError());
    
  }
  
  /**
   * Fortran wrapper for part2 of cheby
   */
  void cuda_cheby_part2(void *d, void *w, void *x,
                       real *tmp1, real *tmp2, int *n, cudaStream_t strm) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    
    cheby_part2<real>
      <<<nblcks, nthrds, 0, strm>>>((real *) d, (real *) w, (real *) x, *tmp1, *tmp2, *n);
    CUDA_CHECK(cudaGetLastError());
    
  }
}
