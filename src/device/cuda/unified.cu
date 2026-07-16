/*
 Copyright (c) 2026, The Neko Authors
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

/**
 * Zero-copy mapping of host arrays for unified memory architectures.
 *
 * On coherent CPU-GPU platforms (e.g. Grace Hopper), system
 * allocations are directly and coherently accessible from device
 * kernels via NVLink-C2C and ATS.  Mapped host arrays then alias
 * their host allocation instead of being replicated in a separate
 * cudaMalloc buffer, and host-device copies degenerate to no-ops.
 *
 * Unlike a single-memory APU, Grace Hopper has two physical
 * memories (LPDDR on the CPU, HBM on the GPU) and host first-touch
 * places mapped arrays in LPDDR.  Aliased ranges are therefore
 * prefetched to device memory and advised to stay there, giving
 * device-resident data that remains host-visible over the coherent
 * interconnect.
 */

#include <stdlib.h>
#include <string.h>
#include <cuda_runtime.h>
#include <device/cuda/unified.h>

extern "C" {

/**
 * Whether mapped arrays may alias host allocations (zero-copy).
 *
 * Zero-copy is only possible on a cache-coherent platform where the
 * device accesses pageable system memory through the host's page
 * tables (ATS), e.g. Grace Hopper:
 *
 *  - Pageable memory access reflects whether system allocations are
 *    accessible from the device at all.
 *  - Pageable memory access using host page tables identifies a
 *    hardware-coherent platform.  Without it (e.g. x86 + HMM),
 *    pageable access means page fault-and-migrate, which would
 *    thrash host arrays back and forth -- such systems always use
 *    replicated buffers.
 *
 * On capable devices, zero-copy can be disabled with
 * NEKO_CUDA_ZEROCOPY=0.
 */
int cuda_zerocopy(void)
{
  static int zerocopy = -1;
  if (zerocopy < 0) {
    int dev = 0, pageable = 0, hostpt = 0;
    const char *env = getenv("NEKO_CUDA_ZEROCOPY");

    if (cudaGetDevice(&dev) == cudaSuccess &&
        cudaDeviceGetAttribute(&pageable,
                               cudaDevAttrPageableMemoryAccess,
                               dev) == cudaSuccess &&
        cudaDeviceGetAttribute(&hostpt,
                          cudaDevAttrPageableMemoryAccessUsesHostPageTables,
                               dev) == cudaSuccess)
      zerocopy = (pageable && hostpt);
    else
      zerocopy = 0;
    (void) cudaGetLastError();

    if (zerocopy && env != NULL && atoi(env) == 0)
      zerocopy = 0;
  }
  return zerocopy;
}

/**
 * Map \p s bytes of host memory at \p ptr_h to the device.
 *
 * On coherent platforms (see cuda_zerocopy) the device pointer
 * aliases the host allocation; otherwise a separate replicated
 * buffer is allocated with cudaMalloc.  Aliased ranges are advised
 * to prefer device memory and prefetched there, counteracting the
 * host-side first touch which would otherwise leave them in the
 * (slower) CPU memory.
 */
cudaError_t cuda_map(void **ptr_d, void *ptr_h, size_t s)
{
  if (cuda_zerocopy() && ptr_h != NULL) {
    int dev = 0;
    *ptr_d = ptr_h;
    /* Placement only; ignore failures */
    if (cudaGetDevice(&dev) == cudaSuccess) {
#if CUDART_VERSION >= 13000
      /* CUDA 13 dropped the int-device overloads */
      struct cudaMemLocation loc;
      loc.type = cudaMemLocationTypeDevice;
      loc.id = dev;
      (void) cudaMemAdvise(ptr_h, s, cudaMemAdviseSetPreferredLocation, loc);
      (void) cudaMemPrefetchAsync(ptr_h, s, loc, 0, 0);
#else
      (void) cudaMemAdvise(ptr_h, s, cudaMemAdviseSetPreferredLocation, dev);
      (void) cudaMemPrefetchAsync(ptr_h, s, dev, 0);
#endif
    }
    (void) cudaGetLastError();
    return cudaSuccess;
  }

  return cudaMalloc(ptr_d, s);
}

/**
 * Whether \p ptr_d is a device allocation (cudaMalloc); zero-copy
 * mappings are unregistered host pointers.
 */
static int cuda_is_device_alloc(void *ptr_d)
{
  struct cudaPointerAttributes attr;
  if (cudaPointerGetAttributes(&attr, ptr_d) != cudaSuccess) {
    (void) cudaGetLastError();
    return 0;
  }
  return (attr.type == cudaMemoryTypeDevice);
}

/**
 * Free a device pointer obtained from cuda_map (or cudaMalloc).
 *
 * Pointers aliasing host memory are left untouched; their
 * allocation is owned by the host side.
 */
cudaError_t cuda_map_free(void *ptr_d)
{
  if (cuda_zerocopy() && !cuda_is_device_alloc(ptr_d))
    return cudaSuccess;

  return cudaFree(ptr_d);
}

/**
 * Memset on a device pointer obtained from cuda_map (or cudaMalloc).
 *
 * cudaMemset only accepts device allocations, so pointers aliasing
 * host memory are set on the host instead, ordered after any
 * in-flight device work on \p stream.
 */
cudaError_t cuda_map_memset(void *ptr_d, int value, size_t s, void *stream)
{
  if (cuda_zerocopy() && !cuda_is_device_alloc(ptr_d)) {
    cudaError_t err = cudaStreamSynchronize((cudaStream_t) stream);
    if (err != cudaSuccess)
      return err;
    memset(ptr_d, value, s);
    return cudaSuccess;
  }

  return cudaMemsetAsync(ptr_d, value, s, (cudaStream_t) stream);
}

} /* extern "C" */
