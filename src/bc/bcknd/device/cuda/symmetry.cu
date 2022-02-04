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

#include <algorithm>
#include <device/device_config.h>
#include <device/cuda/check.h>

#include "symmetry_kernel.h"


extern "C" {

  /** 
   * Fortran wrapper for device symmetry apply vector
   */
  void cuda_symmetry_apply_vector(void *xmsk, void *ymsk, void *zmsk,
                                 void *x, void *y, void *z,
                                 int *m, int *n, int *l) {

    const int max_len = std::max(std::max(*m, *n), *l);
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((max_len) + 1024 - 1)/ 1024, 1, 1);

    symmetry_apply_vector_kernel<real>
      <<<nblcks, nthrds>>>((int *) xmsk, (int *) ymsk, (int *) zmsk,
                           (real *) x, (real *) y, (real *) z, *m, *n, *l);
    CUDA_CHECK(cudaGetLastError());
  }
 
}
