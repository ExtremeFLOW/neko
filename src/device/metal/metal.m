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
 * C-callable wrappers around Metal API for the Neko device layer.
 *
 * On Apple Silicon the CPU and GPU share physical memory.  We use
 * MTLResourceStorageModeShared so that a single MTLBuffer is
 * accessible from both sides.  The opaque c_ptr that Fortran stores
 * is the ARC-retained MTLBuffer pointer.
 */

#ifdef __APPLE__

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#include <string.h>

/*
 *  Module-local Metal device and default command queue
 */
static id<MTLDevice>       _mtl_device    = nil;
static id<MTLCommandQueue> _mtl_glb_queue = nil;
static id<MTLCommandQueue> _mtl_aux_queue = nil;

/*
 * Public accessor so that kernel dispatch code (ax_helm.m, etc.)
 * can reach the device and the default library.
 */
id<MTLDevice> neko_metal_device(void) { return _mtl_device; }

/*
 * Embedded Metal shader library
 */

/* Generated at build time by xxd -i from neko_kernels.metallib */
extern const unsigned char neko_kernels_metallib[];
extern const unsigned int  neko_kernels_metallib_len;

static id<MTLLibrary> _neko_mtl_library = nil;

/**
 * Return the shared MTLLibrary built from all embedded Metal kernels.
 * The library is compiled into the binary at build time and loaded on
 * first call via newLibraryWithData:.
 */
id<MTLLibrary> neko_metal_library(void)
{
    if (_neko_mtl_library == nil) {
        NSError *error = nil;
        dispatch_data_t data = dispatch_data_create(
            neko_kernels_metallib,
            neko_kernels_metallib_len,
            NULL,
            DISPATCH_DATA_DESTRUCTOR_DEFAULT);
        _neko_mtl_library = [_mtl_device newLibraryWithData:data error:&error];
        if (_neko_mtl_library == nil) {
            fprintf(stderr, "Metal: failed to load embedded shader library\n");
            if (error)
                fprintf(stderr, "  %s\n",
                        [[error localizedDescription] UTF8String]);
            exit(EXIT_FAILURE);
        }
    }
    return _neko_mtl_library;
}

/*
 * Init / finalize
 */

void metal_init(void **glb_queue, void **aux_queue)
{
    _mtl_device = MTLCreateSystemDefaultDevice();
    if (!_mtl_device) {
        fprintf(stderr, "Metal: failed to create system default device\n");
        exit(1);
    }
    _mtl_glb_queue = [_mtl_device newCommandQueue];
    _mtl_aux_queue = [_mtl_device newCommandQueue];

    /* Return the queues as opaque pointers. */
    *glb_queue = (__bridge_retained void *)_mtl_glb_queue;
    *aux_queue = (__bridge_retained void *)_mtl_aux_queue;
}

void metal_finalize(void *glb_queue, void *aux_queue)
{
    /* Transfer ownership back to ARC so the queues are released. */
    if (glb_queue)
        (void)(__bridge_transfer id<MTLCommandQueue>)glb_queue;
    if (aux_queue)
        (void)(__bridge_transfer id<MTLCommandQueue>)aux_queue;

    _mtl_glb_queue = nil;
    _mtl_aux_queue = nil;
    _mtl_device    = nil;
}

/*
 * Device queries
 */

void metal_device_name(char *name, int len)
{
    if (!_mtl_device) {
        strncpy(name, "No Metal device", len);
        return;
    }
    const char *n = [[_mtl_device name] UTF8String];
    strncpy(name, n, len);
    name[len - 1] = '\0';
}

int metal_device_count(void)
{
    /* On macOS with Apple Silicon there is exactly one GPU. */
    return (_mtl_device != nil) ? 1 : (MTLCreateSystemDefaultDevice() ? 1 : 0);
}

/*
 * Memory allocation / deallocation
 */

/**
 * Allocate a Metal buffer of \p s bytes.
 * Returns 0 on success, non-zero on failure.
 * The returned pointer is an ARC-retained MTLBuffer.
 */
int metal_alloc(void **ptr, size_t s)
{
    id<MTLBuffer> buf =
        [_mtl_device newBufferWithLength:s
                                options:MTLResourceStorageModeShared];
    if (!buf) {
        *ptr = NULL;
        return -1;
    }
    *ptr = (__bridge_retained void *)buf;
    return 0;
}

/**
 * Free a previously allocated Metal buffer.
 * Returns 0 on success.
 */
int metal_free(void *ptr)
{
    if (ptr) {
        /*
         * Transfer ownership back to ARC — ARC releases when buf
         * goes out of scope.
         */
        id<MTLBuffer> buf __attribute__((unused)) =
            (__bridge_transfer id<MTLBuffer>)ptr;
    }
    return 0;
}

/*
 * Memcpy
 */

/**
 * Copy \p s bytes from host memory \p src to Metal buffer \p dst_d.
 *
 * Because we use shared-mode buffers the copy is a plain memcpy.
 * A blit encoder could be used for managed-mode on discrete GPUs,
 * but Apple Silicon only has shared memory.
 */
int metal_memcpy_htod(void *dst_d, const void *src, size_t s)
{
    id<MTLBuffer> buf = (__bridge id<MTLBuffer>)dst_d;
    memcpy([buf contents], src, s);
    return 0;
}

/**
 * Copy \p s bytes from Metal buffer \p src_d to host memory \p dst.
 */
int metal_memcpy_dtoh(void *dst, void *src_d, size_t s)
{
    id<MTLBuffer> buf = (__bridge id<MTLBuffer>)src_d;
    memcpy(dst, [buf contents], s);
    return 0;
}

/**
 * Copy \p s bytes between two Metal buffers.
 */
int metal_memcpy_dtod(void *dst_d, void *src_d, size_t s)
{
    id<MTLBuffer> dst_buf = (__bridge id<MTLBuffer>)dst_d;
    id<MTLBuffer> src_buf = (__bridge id<MTLBuffer>)src_d;
    memcpy([dst_buf contents], [src_buf contents], s);
    return 0;
}

/*
 * Memset
 */

int metal_memset(void *buf_d, int value, size_t s)
{
    id<MTLBuffer> buf = (__bridge id<MTLBuffer>)buf_d;
    memset([buf contents], value, s);
    return 0;
}

/*
 * Synchronization
 */

/**
 * Wait until all commands in the default queue have completed.
 *
 * Metal doesn't expose a queue-level sync directly.  We commit
 * an empty command buffer and wait on it.
 */
int metal_device_sync(void)
{
    id<MTLCommandBuffer> cb = [_mtl_glb_queue commandBuffer];
    [cb commit];
    [cb waitUntilCompleted];
    return 0;
}

/**
 * Wait until all commands in the given queue have completed.
 */
int metal_stream_sync(void *queue)
{
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)queue;
    id<MTLCommandBuffer> cb = [q commandBuffer];
    [cb commit];
    [cb waitUntilCompleted];
    return 0;
}

/*
 * Streams (command queues)
 */

int metal_stream_create(void **stream)
{
    id<MTLCommandQueue> q = [_mtl_device newCommandQueue];
    if (!q) return -1;
    *stream = (__bridge_retained void *)q;
    return 0;
}

int metal_stream_destroy(void *stream)
{
    if (stream) {
        id<MTLCommandQueue> q __attribute__((unused)) =
            (__bridge_transfer id<MTLCommandQueue>)stream;
    }
    return 0;
}

/*
 * Events
 */

int metal_event_create(void **event)
{
    id<MTLSharedEvent> e = [_mtl_device newSharedEvent];
    if (!e) return -1;
    *event = (__bridge_retained void *)e;
    return 0;
}

int metal_event_destroy(void *event)
{
    if (event) {
        id<MTLSharedEvent> e __attribute__((unused)) =
            (__bridge_transfer id<MTLSharedEvent>)event;
    }
    return 0;
}

int metal_event_record(void *event, void *queue)
{
    id<MTLSharedEvent> e = (__bridge id<MTLSharedEvent>)event;
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)queue;
    uint64_t val = e.signaledValue + 1;
    id<MTLCommandBuffer> cb = [q commandBuffer];
    [cb encodeSignalEvent:e value:val];
    [cb commit];
    return 0;
}

int metal_event_sync(void *event)
{
    id<MTLSharedEvent> e = (__bridge id<MTLSharedEvent>)event;
    /* Spin-wait until the GPU signals the event. */
    uint64_t target = e.signaledValue;
    while (e.signaledValue < target) {
        /*
         * On Apple Silicon this is essentially a no-op since
         * event_record commits the command buffer.
         */
    }
    return 0;
}

int metal_stream_wait_event(void *queue, void *event)
{
    id<MTLSharedEvent> e = (__bridge id<MTLSharedEvent>)event;
    id<MTLCommandQueue> q = (__bridge id<MTLCommandQueue>)queue;
    uint64_t val = e.signaledValue;
    id<MTLCommandBuffer> cb = [q commandBuffer];
    [cb encodeWaitForEvent:e value:val];
    [cb commit];
    return 0;
}

#endif /* __APPLE__ */
