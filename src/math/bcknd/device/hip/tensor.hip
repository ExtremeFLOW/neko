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
#include "tensor_kernel.h"
#include <stdio.h>

#define max(a,b) \
  ({ __typeof__ (a) _a = (a);  \
     __typeof__ (b) _b = (b);   \
     _a > _b ? _a : _b; })


/*
 * Threads per block for tnsr3d.
 *
 * One thread per output point of the widest of the three phases, rounded up
 * to whole wavefronts and capped at the block maximum. The three phases
 * produce nu*nu*nv, nu*nv*nv and nv*nv*nv points, so the widest is nu*nu*nv
 * when restricting and nv*nv*nv when prolonging.
 *
 * This used to be a flat 1024 regardless, which on a p-transfer left most of
 * the block with nothing to do: an 8 -> 4 restriction has 256, 128 and 64
 * points in its phases, so three quarters to fifteen sixteenths of the
 * threads idled. p-multigrid transfers between adjacent levels constantly,
 * which is where that showed up.
 *
 * Any block size is correct here: the phases are grid stride loops and the
 * barriers between them sit outside the loop bodies, so a thread with no
 * iterations still reaches every barrier.
 */
static int hip_tnsr3d_nthrds(const int nu, const int nv)
{
  const int w1 = nu * nu * nv;
  const int w3 = nv * nv * nv;
  int work = (w1 > w3) ? w1 : w3;

  work = ((work + 64 - 1) / 64) * 64;
  if (work > 1024) work = 1024;
  if (work < 64) work = 64;
  return work;
}

extern "C" {

  /** Fortran wrapper for tnsr3d **/
  void hip_tnsr3d(void *v, int *nv, void *u, int *nu,
		  void *A, void *Bt, void *Ct, int *nel) {
    
    const dim3 nthrds(hip_tnsr3d_nthrds(*nu, *nv), 1, 1);
    const dim3 nblcks(*nel, 1, 1);

    int n = max(*nu,*nv);
#define CASE(N)                                                     \
    case N:                                                         \
    hipLaunchKernelGGL(HIP_KERNEL_NAME(tnsr3d_kernel<real, N>),     \
                       nblcks, nthrds, 0,                           \
                       (hipStream_t) glb_cmd_queue,                 \
                       (real *) v, *nv,                             \
                       (real *) u, *nu,                             \
                       (real *) A, (real *) Bt, (real *) Ct);       \
    HIP_CHECK(hipGetLastError());                                   \
    break

    switch(n) {
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
        fprintf(stderr, __FILE__ ": size not supported: %d\n", n);
        exit(1);
      }
    }
  }

 /** Fortran wrapper for tnsr3d **/
  void hip_tnsr3d_el_list(void *v, int *nv, void *u, int *nu,
		   void *A, void *Bt, void *Ct, int * elements, int* n_points) {
    const dim3 nthrds(hip_tnsr3d_nthrds(*nu, *nv), 1, 1);
    const dim3 nblcks(*n_points, 1, 1);

    int n = max(*nu,*nv);
#define CASE2(N)                                                    \
    case N:                                                         \
    hipLaunchKernelGGL(HIP_KERNEL_NAME(tnsr3d_el_kernel<real, N>),  \
                       nblcks, nthrds, 0,                           \
                       (hipStream_t) glb_cmd_queue,                 \
                       (real *) v, *nv,                             \
                       (real *) u, *nu,                             \
                       (real *) A, (real *) Bt, (real *) Ct,        \
                       (int *) elements, *n_points);                \
    HIP_CHECK(hipGetLastError());                                   \
    break

    switch(n) {
      CASE2(2);
      CASE2(3);
      CASE2(4);
      CASE2(5);
      CASE2(6);
      CASE2(7);
      CASE2(8);
      CASE2(9);
      CASE2(10);
      CASE2(11);
      CASE2(12);
      CASE2(13);
      CASE2(14);
    default:
      {
        fprintf(stderr, __FILE__ ": size not supported: %d\n", n);
        exit(1);
      }
    }
  }


}
