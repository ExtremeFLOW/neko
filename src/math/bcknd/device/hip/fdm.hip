/*
 Copyright (c) 2022, The Neko Authors
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
#include "fdm_kernel.h"

extern "C" {

  /** Fortran wrapper for tnsr3d **/
  void hip_fdm_do_fast(void *e, void *r, void *s, void *d, int *nl, int *nel, hipStream_t stream) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(*nel, 1, 1);

#define CASE(NL)                                                                 \
    case NL:                                                                     \
    hipLaunchKernelGGL(HIP_KERNEL_NAME(fdm_do_fast_kernel<real,NL>),             \
                       nblcks, nthrds, 0, stream,           \
                       (real *) e, (real *) r,                                   \
                       (real *) s,(real *) d);                                   \
    HIP_CHECK(hipGetLastError());                                                \
    break;

    switch(*nl) {
       CASE(2);
       CASE(3);
       CASE(4);
       CASE(5);
       CASE(6);
       CASE(7);
       CASE(8);
       CASE(9);
       CASE(10);
       CASE(11);
       CASE(12);
       CASE(13);
       CASE(14);
       CASE(15);
     default:
       {
         fprintf(stderr, __FILE__ ": size not supported: %d\n", *nl);
         exit(1);
       }
    }
  }
}
