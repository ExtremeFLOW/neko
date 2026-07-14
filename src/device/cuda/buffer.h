#ifndef __CUDA_BUFFER_H
#define __CUDA_BUFFER_H
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

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Device buffer: a lazily grown pair of a pinned host buffer and a
 * device buffer, owned by the device layer. A buffer is linked into
 * a registry on its first reservation, and cuda_buffer_free_all()
 * (called from cuda_finalize) releases all registered buffers and
 * resets them to an empty state, from which they grow again on
 * demand after a subsequent device init
 */
struct cuda_buffer {
  void *host;               /**< Pinned host buffer (NULL if device-only) */
  void *dev;                /**< Device buffer */
  size_t size;              /**< Capacity in bytes */
  int symm;                 /**< Device buffer on NVSHMEM symmetric heap */
  int dev_only;             /**< Device buffer only, no pinned host buffer */
  int registered;           /**< Buffer linked into the registry */
  struct cuda_buffer *next; /**< Registry link */
};
typedef struct cuda_buffer cuda_buffer_t;

#define CUDA_BUFFER_INIT {NULL, NULL, 0, 0, 0, 0, NULL}
#define CUDA_BUFFER_INIT_SYMM {NULL, NULL, 0, 1, 0, 0, NULL}
#define CUDA_BUFFER_INIT_DEV {NULL, NULL, 0, 0, 1, 0, NULL}

/**
 * Ensure that a buffer holds at least @a size bytes on both the
 * host and the device side (device side only for device-only
 * buffers), growing it if necessary
 */
void cuda_buffer_reserve(cuda_buffer_t *buf, size_t size);

/**
 * Free all registered buffers and reset them to an empty state
 */
void cuda_buffer_free_all(void);

#ifdef __cplusplus
}
#endif

#endif /* __CUDA_BUFFER_H */
