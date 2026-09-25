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

#include <stdlib.h>

#include <device/device_config.h>
#include <device/opencl/check.h>
#include <device/opencl/buffer.h>

/**
 * Registry of all reserved device buffers
 */
static opencl_buffer_t *buffer_registry = NULL;

void opencl_buffer_reserve(opencl_buffer_t *buf, size_t size) {
  cl_int err;

  if (size <= buf->size)
    return;

  if (buf->host != NULL) {
    free(buf->host);
    CL_CHECK(clReleaseMemObject(buf->dev));
  }

  buf->host = malloc(size);
  buf->dev = clCreateBuffer(glb_ctx, CL_MEM_READ_WRITE, size, NULL, &err);
  CL_CHECK(err);
  buf->size = size;

  if (!buf->registered) {
    buf->next = buffer_registry;
    buffer_registry = buf;
    buf->registered = 1;
  }
}

void opencl_buffer_free_all(void) {
  opencl_buffer_t *buf = buffer_registry;

  while (buf != NULL) {
    opencl_buffer_t *next = buf->next;
    if (buf->host != NULL) {
      free(buf->host);
      buf->host = NULL;
      CL_CHECK(clReleaseMemObject(buf->dev));
      buf->dev = NULL;
    }
    buf->size = 0;
    buf->registered = 0;
    buf->next = NULL;
    buf = next;
  }
  buffer_registry = NULL;
}
