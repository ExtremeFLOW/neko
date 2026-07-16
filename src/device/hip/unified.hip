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

#include <stdlib.h>
#include <string.h>
#include <hip/hip_runtime.h>
#include <device/hip/unified.h>

extern "C" {

/**
 * Whether mapped arrays may alias host allocations (zero-copy).
 *
 * Zero-copy is only possible on an APU with unified physical memory
 * and XNACK enabled, such that system allocations are directly and
 * coherently accessible from the device:
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
 *
 * On capable devices, zero-copy can be disabled with
 * NEKO_HIP_ZEROCOPY=0.
 */
int hip_zerocopy(void)
{
  static int zerocopy = -1;
  if (zerocopy < 0) {
    int dev = 0, integrated = 0, pageable = 0;
    const char *env = getenv("NEKO_HIP_ZEROCOPY");

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

    if (zerocopy && env != NULL && atoi(env) == 0)
      zerocopy = 0;
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
    int dev = 0;
    *ptr_d = ptr_h;
    /* Advice only; ignore failures */
    if (hipGetDevice(&dev) == hipSuccess)
      (void) hipMemAdvise(ptr_h, s, hipMemAdviseSetCoarseGrain, dev);
    (void) hipGetLastError();
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
 * Pointers aliasing host memory are left untouched; their
 * allocation is owned by the host side.
 */
hipError_t hip_map_free(void *ptr_d)
{
  if (hip_zerocopy() && !hip_is_device_alloc(ptr_d))
    return hipSuccess;

  return hipFree(ptr_d);
}

/**
 * Memset on a device pointer obtained from hip_map (or hipMalloc).
 *
 * hipMemset only accepts device allocations, so pointers aliasing
 * host memory are set on the host instead, ordered after any
 * in-flight device work on \p stream.
 */
hipError_t hip_map_memset(void *ptr_d, int value, size_t s, void *stream)
{
  if (hip_zerocopy() && !hip_is_device_alloc(ptr_d)) {
    hipError_t err = hipStreamSynchronize((hipStream_t) stream);
    if (err != hipSuccess)
      return err;
    memset(ptr_d, value, s);
    return hipSuccess;
  }

  return hipMemsetAsync(ptr_d, value, s, (hipStream_t) stream);
}

} /* extern "C" */
