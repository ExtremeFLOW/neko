#ifndef __OPENCL_BUFFER_H
#define __OPENCL_BUFFER_H
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

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

/**
 * Device buffer: a lazily grown pair of a host buffer and a device
 * buffer, owned by the device layer. A buffer is linked into a
 * registry on its first reservation, and opencl_buffer_free_all()
 * (called from opencl_finalize before the context is released)
 * releases all registered buffers and resets them to an empty
 * state, from which they grow again on demand after a subsequent
 * device init
 */
struct opencl_buffer {
  void *host;                 /**< Host buffer */
  cl_mem dev;                 /**< Device buffer */
  size_t size;                /**< Capacity in bytes */
  int registered;             /**< Buffer linked into the registry */
  struct opencl_buffer *next; /**< Registry link */
};
typedef struct opencl_buffer opencl_buffer_t;

#define OPENCL_BUFFER_INIT {NULL, NULL, 0, 0, NULL}

/**
 * Ensure that a buffer holds at least @a size bytes on both the
 * host and the device side, growing it if necessary
 */
void opencl_buffer_reserve(opencl_buffer_t *buf, size_t size);

/**
 * Free all registered buffers and reset them to an empty state
 */
void opencl_buffer_free_all(void);

#endif /* __OPENCL_BUFFER_H */
