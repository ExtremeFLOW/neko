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
 * On APUs with a single physical memory (e.g. MI300A) and XNACK
 * enabled (HSA_XNACK=1), system allocations are directly and
 * coherently accessible from device kernels.  Mapped host arrays
 * then alias their host allocation instead of being replicated in a
 * separate hipMalloc buffer, and host-device copies degenerate to
 * no-ops.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#ifdef __linux__
#include <unistd.h>
#include <sys/mman.h>
#endif
#include <hip/hip_runtime.h>
#include <device/hip/unified.h>

/**
 * Grid-stride byte memset kernel for zero-copy mappings.
 */
__global__ void hip_unified_memset_kernel(unsigned char *a,
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
__global__ void hip_unified_copy8_kernel(unsigned long long * __restrict__ a,
                                         const unsigned long long
                                         * __restrict__ b, size_t n)
{
  const size_t idx = blockIdx.x * (size_t) blockDim.x + threadIdx.x;
  const size_t str = (size_t) blockDim.x * gridDim.x;
  for (size_t i = idx; i < n; i += str)
    a[i] = b[i];
}

__global__ void hip_unified_copy1_kernel(unsigned char * __restrict__ a,
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
 * Zero-copy is OFF by default and must be opted into with
 * NEKO_HIP_ZEROCOPY=1: it is correct on a unified-memory APU but
 * currently slower than replicated buffers on MI300A (the gather-
 * scatter exchange in particular), so replication is the default
 * until that overhead is understood.  It also halves the memory
 * footprint of mapped data, which is its reason to exist as an
 * opt-in capacity mode.
 *
 * When opted in, zero-copy is still only possible on an APU with
 * unified physical memory and XNACK enabled, such that system
 * allocations are directly and coherently accessible from the
 * device:
 *
 *  - The integrated attribute identifies an APU (e.g. MI300A, as
 *    opposed to the discrete MI300X, which shares the same gfx942
 *    architecture).  Discrete GPUs (MI250X, MI300X) always use
 *    replicated buffers; with XNACK they also report pageable
 *    memory access, but there it means page migration, which
 *    would thrash host arrays back and forth.
 *  - Pageable memory access reflects the effective XNACK state of
 *    the process (HSA_XNACK=1, kernel support and xnack-capable
 *    code objects); on an APU with XNACK off it reports 0 and we
 *    fall back to replicated buffers.
 */
int hip_zerocopy(void)
{
  static int zerocopy = -1;
  if (zerocopy < 0) {
    const char *env = getenv("NEKO_HIP_ZEROCOPY");

    if (env != NULL && atoi(env) != 0) {
      int dev = 0, integrated = 0, pageable = 0;
      if (hipGetDevice(&dev) == hipSuccess &&
          hipDeviceGetAttribute(&integrated,
                                hipDeviceAttributeIntegrated,
                                dev) == hipSuccess &&
          hipDeviceGetAttribute(&pageable,
                                hipDeviceAttributePageableMemoryAccess,
                                dev) == hipSuccess)
        zerocopy = (integrated && pageable);
      else
        zerocopy = 0;
      (void) hipGetLastError();

      if (!zerocopy)
        fprintf(stderr,
                "Neko: warning: NEKO_HIP_ZEROCOPY=1 requested but this "
                "device does not support zero-copy mapping (needs a "
                "unified-memory APU with XNACK); using replicated "
                "buffers\n");
    } else {
      zerocopy = 0;
    }
  }
  return zerocopy;
}

/**
 * Map \p s bytes of host memory at \p ptr_h to the device.
 *
 * On unified memory (see hip_zerocopy) the device pointer aliases
 * the host allocation; otherwise a separate replicated buffer is
 * allocated with hipMalloc.  Aliased ranges are advised to
 * coarse-grained coherence, since mapped arrays are never accessed
 * concurrently from host and device and fine-grained system memory
 * can be slower in kernels.
 */
hipError_t hip_map(void **ptr_d, void *ptr_h, size_t s)
{
  if (hip_zerocopy() && ptr_h != NULL) {
    static int warned = 0;
    int dev = 0;
    *ptr_d = ptr_h;
    /* Advice only, but warn on the first failure: without coarse
       grain, kernels on fine-grained system memory can be slower */
    if (hipGetDevice(&dev) == hipSuccess) {
      if (hipMemAdvise(ptr_h, s, hipMemAdviseSetCoarseGrain, dev)
          != hipSuccess && !warned) {
        fprintf(stderr,
                "Neko: warning: coarse-grain advise of zero-copy mapping "
                "failed; mapped arrays stay fine-grained (slower kernels)\n");
        warned = 1;
      }
    }
    (void) hipGetLastError();
#ifdef __linux__
    /* Best effort: back the mapping with transparent huge pages,
       reducing GPU TLB pressure (system pages are much smaller than
       the large pages backing hipMalloc allocations) */
    {
      const uintptr_t pgsz = (uintptr_t) sysconf(_SC_PAGESIZE);
      const uintptr_t beg = ((uintptr_t) ptr_h + pgsz - 1) & ~(pgsz - 1);
      const uintptr_t end = ((uintptr_t) ptr_h + s) & ~(pgsz - 1);
      if (end > beg)
        (void) madvise((void *) beg, end - beg, MADV_HUGEPAGE);
    }
#endif
    return hipSuccess;
  }

  return hipMalloc(ptr_d, s);
}

/**
 * Whether \p ptr_d is a device allocation (hipMalloc); zero-copy
 * mappings are unregistered host pointers.
 */
static int hip_is_device_alloc(void *ptr_d)
{
  hipPointerAttribute_t attr;
  if (hipPointerGetAttributes(&attr, ptr_d) != hipSuccess) {
    (void) hipGetLastError();
    return 0;
  }
#if HIP_VERSION_MAJOR >= 6
  return (attr.type == hipMemoryTypeDevice);
#else
  return (attr.memoryType == hipMemoryTypeDevice);
#endif
}

/**
 * Free a device pointer obtained from hip_map (or hipMalloc).
 *
 * Pointers aliasing host memory are left untouched; their allocation
 * is owned by the host side.  hipFree implicitly synchronizes the
 * device, and the host typically deallocates right after an unmap,
 * so preserve that barrier for aliased mappings such that the
 * allocation cannot be released under in-flight device work.
 */
hipError_t hip_map_free(void *ptr_d)
{
  if (hip_zerocopy() && !hip_is_device_alloc(ptr_d))
    return hipDeviceSynchronize();

  return hipFree(ptr_d);
}

/**
 * Memset on a device pointer obtained from hip_map (or hipMalloc).
 *
 * hipMemset only accepts device allocations, so pointers aliasing
 * host memory are set with a kernel instead, which also keeps the
 * operation stream-ordered.
 */
hipError_t hip_map_memset(void *ptr_d, int value, size_t s, void *stream)
{
  if (hip_zerocopy() && !hip_is_device_alloc(ptr_d)) {
    const size_t nb = (s + 1024 - 1) / 1024;
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((unsigned int) (nb > 65535 ? 65535 :
                                      (nb == 0 ? 1 : nb)), 1, 1);
    hipLaunchKernelGGL(hip_unified_memset_kernel, nblcks, nthrds, 0,
                       (hipStream_t) stream,
                       (unsigned char *) ptr_d, (unsigned char) value, s);
    return hipGetLastError();
  }

  return hipMemsetAsync(ptr_d, value, s, (hipStream_t) stream);
}

/**
 * Memcpy between device pointers obtained from hip_map (or
 * hipMalloc) and/or host pointers.
 *
 * Under zero-copy either side may alias pageable host memory, which
 * hipMemcpy treats as a staged pageable copy (slow, serializing);
 * copy with a kernel instead, which runs at full memory bandwidth
 * and stays stream-ordered (on a unified memory APU any pointer is
 * dereferenceable from the device).
 */
hipError_t hip_map_memcpy(void *dst, void *src, size_t s,
                          int kind, void *stream)
{
  if (hip_zerocopy()) {
    if (s == 0)
      return hipSuccess;
    if ((((uintptr_t) dst | (uintptr_t) src | s) & 7) == 0) {
      const size_t n = s / 8;
      const size_t nb = (n + 1024 - 1) / 1024;
      const dim3 nthrds(1024, 1, 1);
      const dim3 nblcks((unsigned int) (nb > 65535 ? 65535 : nb), 1, 1);
      hipLaunchKernelGGL(hip_unified_copy8_kernel, nblcks, nthrds, 0,
                         (hipStream_t) stream,
                         (unsigned long long *) dst,
                         (const unsigned long long *) src, n);
    }
    else {
      const size_t nb = (s + 1024 - 1) / 1024;
      const dim3 nthrds(1024, 1, 1);
      const dim3 nblcks((unsigned int) (nb > 65535 ? 65535 : nb), 1, 1);
      hipLaunchKernelGGL(hip_unified_copy1_kernel, nblcks, nthrds, 0,
                         (hipStream_t) stream,
                         (unsigned char *) dst,
                         (const unsigned char *) src, s);
    }
    return hipGetLastError();
  }

  return hipMemcpyAsync(dst, src, s, (hipMemcpyKind) kind,
                        (hipStream_t) stream);
}

} /* extern "C" */
