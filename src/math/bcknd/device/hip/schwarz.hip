/*
 Copyright (c) 2021-2022, The Neko Authors
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
#include "schwarz_kernel.h"
#include <stdio.h>

extern "C" {
  /** 
   * Fortran wrapper for device extrude
   */
  void hip_schwarz_extrude(void *arr1, int * l1, real * f1,
                            void *arr2, int * l2, real * f2,
                            int * nx, int * nel, hipStream_t stream) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);
#define CASE(NX)                                                          \
    case NX:                                                              \
    hipLaunchKernelGGL(HIP_KERNEL_NAME(schwarz_extrude_kernel<real,NX>),  \
<<<<<<< HEAD
                       nblcks, nthrds, 0, stream,                              \
=======
                       nblcks, nthrds, 0, (hipStream_t) glb_cmd_queue,    \
>>>>>>> 24580ed7fce3071907a6d0ffea3b9170379e35e0
                       (real *) arr1,* l1, * f1,                          \
                       (real *) arr2, *l2, *f2);                          \
    HIP_CHECK(hipGetLastError());                                         \
    break;

    switch(*nx) {
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
      CASE(16);
    default:
      {
        fprintf(stderr, __FILE__ ": size not supported: %d\n", *nx);
        exit(1);
      }
   }

  } 

  void hip_schwarz_toext3d(void *a, void *b,int * nx, int * nel, hipStream_t stream){
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(schwarz_toext3d_kernel<real>),
                                       nblcks, nthrds, 0,
                                       stream,
                                       (real *) a,(real *) b, * nx);
    HIP_CHECK(hipGetLastError());
  } 

  void hip_schwarz_toreg3d(void *b, void *a,int * nx, int * nel, hipStream_t stream){
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);
    hipLaunchKernelGGL(HIP_KERNEL_NAME(schwarz_toreg3d_kernel<real>),
                                       nblcks, nthrds, 0, stream,
                                       (real *) b,(real *) a, * nx);
    HIP_CHECK(hipGetLastError());
  } 

}
