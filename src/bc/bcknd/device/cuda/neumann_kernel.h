/*
 Copyright (c) 2025, The Neko Authors
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

#ifndef __BC_NEUMANN_KERNEL__
#define __BC_NEUMANN_KERNEL__

#include "bc_utils.h"

/**
 * Computes the linear index for area and normal arrays
 * @note Fortran indexing input, C indexing output
 */

#define coef_normal_area_idx(i, j, k, l, lx, nf) \
  (((i) + (lx) * (((j) - 1) + (lx) * (((k) - 1) + (nf) * (((l) - 1))))) - 1)


/**
 * Device kernel for neumann scalar boundary condition
 */
template< typename T >
__global__
void neumann_apply_scalar_kernel(const int * __restrict__ msk,
                                 const int * __restrict__ facet,
                                 T * __restrict__ x,
                                 const T * __restrict__ flux,
                                 const T * __restrict__ area,
                                 const int lx,
                                 const int m) {
  int index[4];
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = (idx + 1); i < m; i += str) {
    const int k = (msk[i] - 1);
    const int f = (facet[i]);
    nonlinear_index(msk[i], lx, index);

    switch(f) {
    case 1:
    case 2:
      {
        const int na_idx = coef_normal_area_idx(index[1], index[2],
                                                f, index[3], lx, 6);
        x[k] += flux[i-1] * area[na_idx];
        break;
      }
    case 3:
    case 4:
      {
        const int na_idx = coef_normal_area_idx(index[0], index[2],
                                                f, index[3], lx, 6);
        x[k] += flux[i-1] * area[na_idx];
        break;
      }
    case 5:
    case 6:
      {
        const int na_idx = coef_normal_area_idx(index[0], index[1],
                                                f, index[3], lx, 6);
        x[k] += flux[i-1] * area[na_idx];
        break;
      }
    }
  }
}

#endif // __BC_NEUMANN_KERNEL__
