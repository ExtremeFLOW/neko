/*
 Copyright (c) 2024-2025, The Neko Authors
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

#ifndef __GS_NVSHMEM_KERNELS__
#define __GS_NVHSMEM_KERNELS__

#include <nvshmemx.h>


template< typename T >
__global__ void pack_pushShmemKernel(const T * __restrict__ u,
                                     T * dest,
                                     T * __restrict__ src,
                                     const int * __restrict__ dof,
                                     const int destRank,
                                     const int srcRank,
                                     const int n,
                                     uint64_t counter,
                                     uint64_t * notifyDone,
                                     uint64_t * notifyReady);

template<>
__global__ void pack_pushShmemKernel(const float * __restrict__ u,
                                     float * dest,
                                     float * __restrict__ src,
                                     const int * __restrict__ dof,
                                     const int destRank,
                                     const int srcRank,
                                     const int n,
                                     uint64_t counter,
                                     uint64_t * notifyDone,
                                     uint64_t * notifyReady)
{


  const int j = threadIdx.x + blockDim.x * blockIdx.x;

  if (j <  n) {
    src[j] = u[dof[j]-1];
  }
  __syncthreads();
  
  //TO DO: 1 block transfers seem best from initial investigations, check this more thoroughly
  size_t numBlocksForTransfer = 1; 
  if(blockIdx.x < numBlocksForTransfer)
  {
    size_t n_per_block = n/numBlocksForTransfer;
    size_t block_offset = n_per_block*blockIdx.x;
    size_t dataSize = blockIdx.x != (numBlocksForTransfer - 1) ?
      n_per_block : max(n - block_offset, n_per_block);
    
    // Notify ready to sending rank, and wait until recieving rank is ready
    if (threadIdx.x == 0) {
      nvshmemx_signal_op(notifyReady, counter, NVSHMEM_SIGNAL_SET, srcRank);
      nvshmem_signal_wait_until(notifyReady, NVSHMEM_CMP_EQ, counter);
    }
    __syncthreads();
    
    // Push data
    nvshmemx_float_put_signal_nbi_block(dest + block_offset, src +
                                        block_offset, dataSize,
                                        notifyDone, counter,
                                        NVSHMEM_SIGNAL_SET, destRank);
  }
}

template<>
__global__ void pack_pushShmemKernel(const double * __restrict__ u,
                                     double * dest,
                                     double * __restrict__ src,
                                     const int * __restrict__ dof,
                                     const int destRank,
                                     const int srcRank,
                                     const int n,
                                     uint64_t counter,
                                     uint64_t * notifyDone,
                                     uint64_t * notifyReady)
{


  const int j = threadIdx.x + blockDim.x * blockIdx.x;

  if (j <  n) {
    src[j] = u[dof[j]-1];
  }
  __syncthreads();
  
  //TO DO: 1 block transfers seem best from initial investigations, check this more thoroughly
  size_t numBlocksForTransfer = 1; 
  if(blockIdx.x < numBlocksForTransfer)
  {
    size_t n_per_block = n/numBlocksForTransfer;
    size_t block_offset = n_per_block*blockIdx.x;
    size_t dataSize = blockIdx.x != (numBlocksForTransfer - 1) ?
      n_per_block : max(n - block_offset, n_per_block);
    
    // Notify ready to sending rank, and wait until recieving rank is ready
    if (threadIdx.x == 0) {
      nvshmemx_signal_op(notifyReady, counter, NVSHMEM_SIGNAL_SET, srcRank);
      nvshmem_signal_wait_until(notifyReady, NVSHMEM_CMP_EQ, counter);
    }
    __syncthreads();
    
    // Push data
    nvshmemx_double_put_signal_nbi_block(dest + block_offset, src +
                                         block_offset, dataSize,
                                         notifyDone, counter,
                                         NVSHMEM_SIGNAL_SET, destRank);
  }
}

__global__ void pushShmemKernelWait(uint64_t counter,
                                    uint64_t *notifyDone) 
{
  // Notify done to receiving rank, and wait for data from sending rank
  if (blockIdx.x==0 && threadIdx.x == 0) {
    nvshmem_signal_wait_until(notifyDone, NVSHMEM_CMP_EQ, counter);
  }
}




#endif
