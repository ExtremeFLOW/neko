/*
 Copyright (c) 2021-2026, The Neko Authors
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

#ifndef __DEVICE_METAL_CHECK_H__
#define __DEVICE_METAL_CHECK_H__

#include <stdio.h>
#include <stdlib.h>

/**
 * Check a Metal/NSError status and abort on failure
 */
#define METAL_CHECK(err)                                                   \
  do {                                                                      \
    if ((err) != nil) {                                                     \
      fprintf(stderr, "%s:%d Metal error: %s\n", __FILE__, __LINE__,       \
              [[err localizedDescription] UTF8String]);                     \
      exit(EXIT_FAILURE);                                                   \
    }                                                                       \
  } while (0)

#ifdef __OBJC__
#import <Metal/Metal.h>

/** Threadgroup size every Neko compute pipeline is compiled to support */
#define NEKO_METAL_MAX_THREADS 1024

/**
 * Create a compute pipeline state for @a fn. With @a max_threads > 0 the
 * pipeline is built from a descriptor whose maxTotalThreadsPerThreadgroup
 * is @a max_threads, so the compiler bounds register use to guarantee that
 * threadgroup size (Apple's compiler otherwise lowers the limit of
 * register-heavy kernels, e.g. to 768, and a larger dispatch is rejected
 * silently). With 0 the compiler picks the limit and the wrapper has to
 * respect pso.maxTotalThreadsPerThreadgroup itself.
 */
id<MTLComputePipelineState> neko_metal_new_pipeline(id<MTLFunction> fn,
                                                    NSUInteger max_threads,
                                                    NSError **error);

/** Abort if a completed command buffer failed on the GPU */
#define METAL_CHECK_CB(cb)                                                  \
  do {                                                                      \
    if ([(cb) status] == MTLCommandBufferStatusError) {                     \
      NSError *cb_err__ = [(cb) error];                                     \
      fprintf(stderr, "%s:%d Metal command buffer failed: %s\n",            \
              __FILE__, __LINE__,                                           \
              cb_err__ ? [[cb_err__ localizedDescription] UTF8String]       \
                       : "unknown error");                                  \
      exit(EXIT_FAILURE);                                                   \
    }                                                                       \
  } while (0)
#endif /* __OBJC__ */

#endif /* __DEVICE_METAL_CHECK_H__ */
