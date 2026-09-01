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

#include <device/device_config.h>
#include <device/cuda/check.h>
#ifdef HAVE_NVSHMEM
#include <nvshmem.h>
#include "gs_kernels.h"
#include "gs_nvshmem_kernels.h"
#endif

extern "C" {

#ifdef HAVE_NVSHMEM
  void cudamalloc_nvshmem(void** ptr, size_t size)
  {
    *ptr = nvshmem_malloc(size);
    CUDA_CHECK(cudaGetLastError());
    cudaMemset(*ptr, 0, size);
    CUDA_CHECK(cudaGetLastError());
  }

  void cudafree_nvshmem(void** ptr, size_t size)
  {
    nvshmem_free(*ptr);
    CUDA_CHECK(cudaGetLastError());
  }

  /**
   * Bulk-pack every peer's slab into one parity half of the double-buffered
   * send buffer. boffset is the parity slab offset (in elements); the dof
   * map always starts at 0 and covers all peers. Launched on the main
   * stream in nbsend, before any per-peer work, so every read of u
   * precedes any unpack write to u (see gs_nvshmem_kernels.h).
   */
  void cuda_gs_nvshmem_pack(void *u_d, void *buf_d, void *dof_d,
                            int boffset, int n, cudaStream_t stream)
  {
    if (n == 0) return;

    const int nthrds = 1024;
    const int nblcks = (n + nthrds - 1) / nthrds;

    gs_pack_kernel<real>
      <<<nblcks, nthrds, 0, stream>>>((real *) u_d, (real *) buf_d + boffset,
                                      (int *) dof_d, n);
    CUDA_CHECK(cudaGetLastError());
  }

  /**
   * Fused nc-component bulk pack (see cuda_gs_nvshmem_pack). u_d is the
   * compact shared buffer (component-outer, per-component stride ns);
   * buf_d is interleaved (nc per position).
   */
  void cuda_gs_nvshmem_pack_vec(void *u_d, void *buf_d, void *dof_d,
                                int boffset, int n, int nc, int ns,
                                cudaStream_t stream)
  {
    if (n == 0) return;

    const int nthrds = 1024;
    const int nblcks = (n + nthrds - 1) / nthrds;

    gs_pack_vec_kernel<real>
      <<<nblcks, nthrds, 0, stream>>>((real *) u_d, (real *) buf_d + boffset,
                                      (int *) dof_d, n, nc, ns);
    CUDA_CHECK(cudaGetLastError());
  }

  /**
   * Push one packed peer slab with rank-indexed signaling. soffset/roffset/n
   * are in elements (the fused vector path passes nc-scaled values); roffset
   * is the offset in the destination's recv buffer where our slab lands;
   * done_d and ready_d are the symmetric rank-indexed signal arrays; mype is
   * our rank. Single-block launch (see gs_nvshmem_kernels.h).
   */
  void cuda_gs_push(void *sbuf_d, int soffset, int n, cudaStream_t stream,
                    int destRank, void *rbuf_d, int roffset, int iter,
                    void *done_d, void *ready_d, int mype)
  {
    const int nthrds = 1024;

    pushShmemKernel<real>
      <<<1,nthrds,0,stream>>>((real *) rbuf_d + roffset,
                              (real *) sbuf_d + soffset,
                              (size_t) n, destRank, (uint64_t) iter,
                              (uint64_t *) done_d + mype,
                              (uint64_t *) ready_d + destRank);
    CUDA_CHECK(cudaGetLastError());
  }

  /**
   * Wait until the slab from srcRank has landed (our local done_sig[srcRank]
   * reaches iter).
   */
  void cuda_gs_push_wait(cudaStream_t stream, int iter,
                         void *done_d, int srcRank)
  {
    pushShmemKernelWait<<<1,1,0,stream>>>((uint64_t) iter,
                                          (uint64_t *) done_d + srcRank);
    CUDA_CHECK(cudaGetLastError());
  }

  /**
   * Post our ready signal to srcRank (sets ready_sig[mype] = iter there),
   * allowing it to overwrite its send slab for the next round. Launch after
   * the unpack on the same stream.
   */
  void cuda_gs_post_ready(cudaStream_t stream, int iter, void *ready_d,
                          int mype, int srcRank)
  {
    postReadyShmemKernel<<<1,1,0,stream>>>((uint64_t *) ready_d + mype,
                                           (uint64_t) iter, srcRank);
    CUDA_CHECK(cudaGetLastError());
  }
#endif
}
