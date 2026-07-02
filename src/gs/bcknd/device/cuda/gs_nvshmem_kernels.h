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
#define __GS_NVSHMEM_KERNELS__

#include <nvshmemx.h>

/*
 * Fused pack-and-push kernels with rank-indexed signaling.
 *
 * Signaling uses two symmetric arrays of pe_size slots, indexed by the
 * REMOTE PE's rank (so peer lists need not be uniform in length or aligned
 * in order across ranks, and the symmetric allocations are collective-safe):
 *  - doneSig  = &done_sig[my_rank] on the destination: set to iter by our
 *    put_signal when our slab has landed there.
 *  - readySlot = &ready_sig[destRank] locally: set to iter by the
 *    destination once it has consumed our round-iter slab.
 * The round counter iter advances once per gs op (lockstep across ranks),
 * and all waits use CMP_GE, so no cross-rank counter matching is needed.
 *
 * These are SINGLE-BLOCK kernels (launched with one block): the pack is a
 * block-stride loop, so the block-local __syncthreads() suffices to order
 * the pack before the put. A multi-block pack would race with the push --
 * there is no grid-wide barrier between pack and push. Single-block
 * transfers were also found to perform best.
 *
 * The ready wait comes BEFORE the pack: the destination posts ready(iter-1)
 * only after consuming our round iter-1 slab, which implies our previous
 * (non-blocking) put has fully drained src -- so the pack below can safely
 * overwrite it. Packing before this wait would race with a still-in-flight
 * nbi put on proxy-based transports.
 */

template< typename T >
__global__ void pack_pushShmemKernel(const T * __restrict__ u,
                                     T * dest,
                                     T * __restrict__ src,
                                     const int * __restrict__ dof,
                                     const int destRank,
                                     const int n,
                                     uint64_t iter,
                                     uint64_t * doneSig,
                                     uint64_t * readySlot);

template<>
__global__ void pack_pushShmemKernel(const float * __restrict__ u,
                                     float * dest,
                                     float * __restrict__ src,
                                     const int * __restrict__ dof,
                                     const int destRank,
                                     const int n,
                                     uint64_t iter,
                                     uint64_t * doneSig,
                                     uint64_t * readySlot)
{

  /* Wait until destRank has consumed our previous round (see note above) */
  if (threadIdx.x == 0) {
    nvshmem_signal_wait_until(readySlot, NVSHMEM_CMP_GE, iter - 1);
  }
  __syncthreads();

  /* Pack with a block-stride loop (single-block kernel) */
  for (int j = threadIdx.x; j < n; j += blockDim.x) {
    src[j] = u[dof[j]-1];
  }
  __syncthreads();

  /* Push data and set done_sig[my_rank] = iter on the destination */
  nvshmemx_float_put_signal_nbi_block(dest, src, n,
                                      doneSig, iter,
                                      NVSHMEM_SIGNAL_SET, destRank);
}

template<>
__global__ void pack_pushShmemKernel(const double * __restrict__ u,
                                     double * dest,
                                     double * __restrict__ src,
                                     const int * __restrict__ dof,
                                     const int destRank,
                                     const int n,
                                     uint64_t iter,
                                     uint64_t * doneSig,
                                     uint64_t * readySlot)
{

  /* Wait until destRank has consumed our previous round (see note above) */
  if (threadIdx.x == 0) {
    nvshmem_signal_wait_until(readySlot, NVSHMEM_CMP_GE, iter - 1);
  }
  __syncthreads();

  /* Pack with a block-stride loop (single-block kernel) */
  for (int j = threadIdx.x; j < n; j += blockDim.x) {
    src[j] = u[dof[j]-1];
  }
  __syncthreads();

  /* Push data and set done_sig[my_rank] = iter on the destination */
  nvshmemx_double_put_signal_nbi_block(dest, src, n,
                                       doneSig, iter,
                                       NVSHMEM_SIGNAL_SET, destRank);
}

/* Wait until the slab from a recv peer has landed (doneSlot is our local
   done_sig[src] slot, set by the peer's put_signal) */
__global__ void pushShmemKernelWait(uint64_t iter,
                                    uint64_t *doneSlot)
{
  if (blockIdx.x==0 && threadIdx.x == 0) {
    nvshmem_signal_wait_until(doneSlot, NVSHMEM_CMP_GE, iter);
  }
}

/* Post our ready signal to the peer we receive from: sets
   ready_sig[my_rank] = iter on srcRank, allowing it to overwrite its send
   slab for the next round. Launched after the unpack on the same stream. */
__global__ void postReadyShmemKernel(uint64_t *readySlot,
                                     uint64_t iter,
                                     const int srcRank)
{
  if (blockIdx.x==0 && threadIdx.x == 0) {
    nvshmemx_signal_op(readySlot, iter, NVSHMEM_SIGNAL_SET, srcRank);
  }
}

/*
 * Fused nc-component pack-and-push (single-block, rank-indexed signaling,
 * see the scalar kernels above). u is the compact shared buffer,
 * component-outer with per-component stride ns; src/dest are interleaved
 * (nc per position). The pushed payload is nc*n elements.
 */
template< typename T >
__global__ void pack_pushShmemKernel_vec(const T * __restrict__ u,
                                         T * dest,
                                         T * __restrict__ src,
                                         const int * __restrict__ dof,
                                         const int destRank,
                                         const int n,
                                         const int nc,
                                         const int ns,
                                         uint64_t iter,
                                         uint64_t * doneSig,
                                         uint64_t * readySlot);

template<>
__global__ void pack_pushShmemKernel_vec(const float * __restrict__ u,
                                         float * dest,
                                         float * __restrict__ src,
                                         const int * __restrict__ dof,
                                         const int destRank,
                                         const int n,
                                         const int nc,
                                         const int ns,
                                         uint64_t iter,
                                         uint64_t * doneSig,
                                         uint64_t * readySlot)
{

  /* Wait until destRank has consumed our previous round */
  if (threadIdx.x == 0) {
    nvshmem_signal_wait_until(readySlot, NVSHMEM_CMP_GE, iter - 1);
  }
  __syncthreads();

  /* Pack with a block-stride loop (single-block kernel) */
  for (int j = threadIdx.x; j < n; j += blockDim.x) {
    const int idx = dof[j] - 1;
    for (int c = 0; c < nc; c++)
      src[nc*j + c] = u[c*ns + idx];
  }
  __syncthreads();

  /* Push data and set done_sig[my_rank] = iter on the destination */
  nvshmemx_float_put_signal_nbi_block(dest, src, (size_t) nc * n,
                                      doneSig, iter,
                                      NVSHMEM_SIGNAL_SET, destRank);
}

template<>
__global__ void pack_pushShmemKernel_vec(const double * __restrict__ u,
                                         double * dest,
                                         double * __restrict__ src,
                                         const int * __restrict__ dof,
                                         const int destRank,
                                         const int n,
                                         const int nc,
                                         const int ns,
                                         uint64_t iter,
                                         uint64_t * doneSig,
                                         uint64_t * readySlot)
{

  /* Wait until destRank has consumed our previous round */
  if (threadIdx.x == 0) {
    nvshmem_signal_wait_until(readySlot, NVSHMEM_CMP_GE, iter - 1);
  }
  __syncthreads();

  /* Pack with a block-stride loop (single-block kernel) */
  for (int j = threadIdx.x; j < n; j += blockDim.x) {
    const int idx = dof[j] - 1;
    for (int c = 0; c < nc; c++)
      src[nc*j + c] = u[c*ns + idx];
  }
  __syncthreads();

  /* Push data and set done_sig[my_rank] = iter on the destination */
  nvshmemx_double_put_signal_nbi_block(dest, src, (size_t) nc * n,
                                       doneSig, iter,
                                       NVSHMEM_SIGNAL_SET, destRank);
}

#endif
