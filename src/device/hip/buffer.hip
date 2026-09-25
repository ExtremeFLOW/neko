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

#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include <device/hip/buffer.h>

extern "C" {

  /**
   * Registry of all reserved device buffers
   */
  static hip_buffer_t *buffer_registry = NULL;

  void hip_buffer_reserve(hip_buffer_t *buf, size_t size) {

    if (size <= buf->size)
      return;

    if (buf->dev != NULL) {
      if (buf->host != NULL)
        HIP_CHECK(hipHostFree(buf->host));
      HIP_CHECK(hipFree(buf->dev));
    }

    if (!buf->dev_only)
      HIP_CHECK(hipHostMalloc(&buf->host, size, hipHostMallocDefault));
    HIP_CHECK(hipMalloc(&buf->dev, size));
    buf->size = size;

    if (!buf->registered) {
      buf->next = buffer_registry;
      buffer_registry = buf;
      buf->registered = 1;
    }
  }

  void hip_buffer_free_all(void) {
    hip_buffer_t *buf = buffer_registry;

    while (buf != NULL) {
      hip_buffer_t *next = buf->next;
      if (buf->host != NULL) {
        HIP_CHECK(hipHostFree(buf->host));
        buf->host = NULL;
      }
      if (buf->dev != NULL) {
        HIP_CHECK(hipFree(buf->dev));
        buf->dev = NULL;
      }
      buf->size = 0;
      buf->registered = 0;
      buf->next = NULL;
      buf = next;
    }
    buffer_registry = NULL;
  }

}
