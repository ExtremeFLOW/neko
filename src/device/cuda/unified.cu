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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <cuda_runtime.h>
#include <device/cuda/unified.h>

/**
 * Prefetch a range to device memory (see cuda_map, cuda_map_prefetch).
 * Placement only, but warn on the first failure: if migration is
 * unavailable, zero-copy arrays stay CPU-resident and kernels run
 * at interconnect bandwidth.
 */
static void cuda_unified_prefetch(void *ptr, size_t s, cudaStream_t stream)
{
  static int warned = 0;
  int dev = 0;
  cudaError_t err;

  if (cudaGetDevice(&dev) != cudaSuccess) {
    (void) cudaGetLastError();
    return;
  }

#if CUDART_VERSION >= 13000
  struct cudaMemLocation loc;
  loc.type = cudaMemLocationTypeDevice;
  loc.id = dev;
  err = cudaMemPrefetchAsync(ptr, s, loc, 0, stream);
#else
  err = cudaMemPrefetchAsync(ptr, s, dev, stream);
#endif
  if (err != cudaSuccess && !warned) {
    fprintf(stderr,
            "Neko: warning: prefetch of zero-copy mapping failed (%s); "
            "mapped arrays may stay CPU-resident (slow kernels)\n",
            cudaGetErrorString(err));
    warned = 1;
  }
  (void) cudaGetLastError();
}

/**
 * Grid-stride byte memset kernel for zero-copy mappings.
 */
__global__ void cuda_unified_memset_kernel(unsigned char *a,
                                           unsigned char v, size_t n)
{
  const size_t idx = blockIdx.x * (size_t) blockDim.x + threadIdx.x;
  const size_t str = (size_t) blockDim.x * gridDim.x;
  for (size_t i = idx; i < n; i += str)
    a[i] = v;
}

/**
 * Grid-stride copy kernels for zero-copy mappings (word and byte
 * variants).
 */
__global__ void cuda_unified_copy8_kernel(unsigned long long * __restrict__ a,
                                          const unsigned long long
                                          * __restrict__ b, size_t n)
{
  const size_t idx = blockIdx.x * (size_t) blockDim.x + threadIdx.x;
  const size_t str = (size_t) blockDim.x * gridDim.x;
  for (size_t i = idx; i < n; i += str)
    a[i] = b[i];
}

__global__ void cuda_unified_copy1_kernel(unsigned char * __restrict__ a,
                                          const unsigned char * __restrict__ b,
                                          size_t n)
{
  const size_t idx = blockIdx.x * (size_t) blockDim.x + threadIdx.x;
  const size_t str = (size_t) blockDim.x * gridDim.x;
  for (size_t i = idx; i < n; i += str)
    a[i] = b[i];
}

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
  /* Cached on first use; device selection must already have
     happened (device_init runs before any map/memset/copy) */
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
    /* Placement only; prefer device memory and migrate any pages the
       host has already touched (failures ignored, see helper) */
    if (cudaGetDevice(&dev) == cudaSuccess) {
#if CUDART_VERSION >= 13000
      /* CUDA 13 dropped the int-device overloads */
      struct cudaMemLocation loc;
      loc.type = cudaMemLocationTypeDevice;
      loc.id = dev;
      (void) cudaMemAdvise(ptr_h, s, cudaMemAdviseSetPreferredLocation, loc);
#else
      (void) cudaMemAdvise(ptr_h, s, cudaMemAdviseSetPreferredLocation, dev);
#endif
    }
    (void) cudaGetLastError();
    cuda_unified_prefetch(ptr_h, s, 0);
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
 * Pointers aliasing host memory are left untouched; their allocation
 * is owned by the host side.  cudaFree implicitly synchronizes the
 * device, and the host typically deallocates right after an unmap,
 * so preserve that barrier for aliased mappings such that the
 * allocation cannot be released under in-flight device work.
 */
cudaError_t cuda_map_free(void *ptr_d)
{
  if (cuda_zerocopy() && !cuda_is_device_alloc(ptr_d))
    return cudaDeviceSynchronize();

  return cudaFree(ptr_d);
}

/**
 * Memcpy between device pointers obtained from cuda_map (or
 * cudaMalloc) and/or host pointers.
 *
 * Under zero-copy either side may alias pageable host memory, which
 * cudaMemcpy treats as a staged pageable copy (slow, serializing);
 * copy with a kernel instead, which runs at full memory bandwidth
 * and stays stream-ordered (on a coherent platform any pointer is
 * dereferenceable from the device).
 */
cudaError_t cuda_map_memcpy(void *dst, void *src, size_t s,
                            int kind, void *stream)
{
  if (cuda_zerocopy()) {
    if (s == 0)
      return cudaSuccess;
    if ((((uintptr_t) dst | (uintptr_t) src | s) & 7) == 0) {
      const size_t n = s / 8;
      const size_t nb = (n + 1024 - 1) / 1024;
      const unsigned int nblcks = (unsigned int) (nb > 65535 ? 65535 : nb);
      cuda_unified_copy8_kernel<<<nblcks, 1024, 0, (cudaStream_t) stream>>>
        ((unsigned long long *) dst, (const unsigned long long *) src, n);
    }
    else {
      const size_t nb = (s + 1024 - 1) / 1024;
      const unsigned int nblcks = (unsigned int) (nb > 65535 ? 65535 : nb);
      cuda_unified_copy1_kernel<<<nblcks, 1024, 0, (cudaStream_t) stream>>>
        ((unsigned char *) dst, (const unsigned char *) src, s);
    }
    return cudaGetLastError();
  }

  return cudaMemcpyAsync(dst, src, s, (enum cudaMemcpyKind) kind,
                         (cudaStream_t) stream);
}

/**
 * Memset on a device pointer obtained from cuda_map (or cudaMalloc).
 *
 * cudaMemset only accepts device allocations, so pointers aliasing
 * host memory are set with a kernel instead.  This keeps the
 * operation stream-ordered, and on platforms with separate physical
 * memories (e.g. Grace Hopper) a device-side first touch places the
 * pages in device memory, where mapped arrays should reside.
 */
cudaError_t cuda_map_memset(void *ptr_d, int value, size_t s, void *stream)
{
  if (cuda_zerocopy() && !cuda_is_device_alloc(ptr_d)) {
    const size_t nb = (s + 1024 - 1) / 1024;
    const unsigned int nblcks = (unsigned int) (nb > 65535 ? 65535 :
                                                (nb == 0 ? 1 : nb));
    cuda_unified_memset_kernel<<<nblcks, 1024, 0, (cudaStream_t) stream>>>
      ((unsigned char *) ptr_d, (unsigned char) value, s);
    return cudaGetLastError();
  }

  return cudaMemsetAsync(ptr_d, value, s, (cudaStream_t) stream);
}

/**
 * Prefetch a zero-copy mapping to device memory.
 *
 * Stands in for the host-to-device copy of an aliased mapping: on
 * platforms with separate physical memories (e.g. Grace Hopper),
 * host-side writes populate pages in CPU memory, and this migrates
 * them to device memory at the point where the copy semantics say
 * the data moves.  Placement only; failures are ignored.
 */
cudaError_t cuda_map_prefetch(void *ptr_d, size_t s, void *stream)
{
  if (cuda_zerocopy())
    cuda_unified_prefetch(ptr_d, s, (cudaStream_t) stream);

  return cudaSuccess;
}

} /* extern "C" */
