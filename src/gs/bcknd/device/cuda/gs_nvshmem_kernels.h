/*
 Copyright (c) 2024-2026, The Neko Authors
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
 * Push kernels with rank-indexed signaling.
 *
 * Signaling uses two symmetric arrays of pe_size slots, indexed by the
 * REMOTE PE's rank (so peer lists need not be uniform in length across
 * ranks, and the symmetric allocations are collective-safe):
 *  - doneSig  = &done_sig[my_rank] on the destination: set to iter by our
 *    put_signal when our slab has landed there.
 *  - readySlot = &ready_sig[destRank] locally: set to iter by the
 *    destination once it has consumed our round-iter slab.
 * The round counter iter advances once per gs op (lockstep across ranks),
 * and all waits use CMP_GE, so no cross-rank counter matching is needed.
 *
 * Packing is NOT done here: all peer slabs are packed by one bulk kernel
 * on the main stream in nbsend, into the round's parity half of the
 * double-buffered send buffer (see gs_device_shmem.F90). That orders every
 * read of the shared buffer u before any unpack writes u. Packing inside
 * the push kernel, gated on the remote ready signal, let an unpack from a
 * fast peer modify u before the pack for a slow peer had read it -- for a
 * dof shared with both peers the slow peer then received an already
 * partially reduced value (observed as divergence at large rank counts,
 * where multi-peer dofs and round-level skew are common).
 *
 * The ready wait below therefore gates only the put: the destination posts
 * ready(iter-1) once it has consumed our round iter-1 slab, so its recv
 * slab may be overwritten. The pack needs no remote gate; see the nbsend
 * comment in gs_device_shmem.F90 for why the parity slab has always
 * drained by the time it is repacked.
 *
 * These are SINGLE-BLOCK kernels (launched with one block); single-block
 * transfers were found to perform best at gs slab sizes.
 */

template< typename T >
__global__ void pushShmemKernel(T * dest,
                                const T * __restrict__ src,
                                const size_t n,
                                const int destRank,
                                uint64_t iter,
                                uint64_t * doneSig,
                                uint64_t * readySlot);

template<>
__global__ void pushShmemKernel(float * dest,
                                const float * __restrict__ src,
                                const size_t n,
                                const int destRank,
                                uint64_t iter,
                                uint64_t * doneSig,
                                uint64_t * readySlot)
{

  /* Wait until destRank has consumed our previous round (see note above) */
  if (threadIdx.x == 0) {
    nvshmem_signal_wait_until(readySlot, NVSHMEM_CMP_GE, iter - 1);
  }
  __syncthreads();

  /* Push data and set done_sig[my_rank] = iter on the destination */
  nvshmemx_float_put_signal_nbi_block(dest, src, n,
                                      doneSig, iter,
                                      NVSHMEM_SIGNAL_SET, destRank);
}

template<>
__global__ void pushShmemKernel(double * dest,
                                const double * __restrict__ src,
                                const size_t n,
                                const int destRank,
                                uint64_t iter,
                                uint64_t * doneSig,
                                uint64_t * readySlot)
{

  /* Wait until destRank has consumed our previous round (see note above) */
  if (threadIdx.x == 0) {
    nvshmem_signal_wait_until(readySlot, NVSHMEM_CMP_GE, iter - 1);
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
   ready_sig[my_rank] = iter on srcRank, allowing it to put its next round
   into our recv slab. Launched after the unpack on the same stream. */
__global__ void postReadyShmemKernel(uint64_t *readySlot,
                                     uint64_t iter,
                                     const int srcRank)
{
  if (blockIdx.x==0 && threadIdx.x == 0) {
    nvshmemx_signal_op(readySlot, iter, NVSHMEM_SIGNAL_SET, srcRank);
  }
}

#endif
