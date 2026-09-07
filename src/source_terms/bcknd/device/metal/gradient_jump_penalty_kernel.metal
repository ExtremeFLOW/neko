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
 * Metal compute kernels for the gradient jump penalty source term.
 *
 * @note Apple GPUs do not support FP64; all arithmetic uses float.
 */

#include <metal_stdlib>
using namespace metal;

/**
 * Device kernel for pick_facet_value_hex
 * @note One threadgroup per element.
 */
kernel void pick_facet_value_hex_kernel(device float *b [[ buffer(0) ]],
                                        device const float *a [[ buffer(1) ]],
                                        constant int &nx [[ buffer(2) ]],
                                        uint e [[ threadgroup_position_in_grid ]],
                                        uint idx [[ thread_position_in_threadgroup ]],
                                        uint nthrds [[ threads_per_threadgroup ]]) {
  const int nx2 = nx + 2;
  const int el2 = (int) e * nx2 * nx2 * nx2;
  const int el = (int) e * nx * nx * nx;

  for (int ijk = (int) idx; ijk < nx * nx * nx; ijk += (int) nthrds) {
    const int jk = ijk / nx;
    const int i = ijk - jk * nx;
    const int k = jk / nx;
    const int j = jk - k * nx;
    if (i == 0) {
      b[0 + (j + 1) * nx2 + (k + 1) * nx2 * nx2 + el2] = a[ijk + el];
    }
    if (i == nx - 1) {
      b[nx2 - 1 + (j + 1) * nx2 + (k + 1) * nx2 * nx2 + el2] = a[ijk + el];
    }
    if (j == 0) {
      b[(i + 1) + 0 * nx2 + (k + 1) * nx2 * nx2 + el2] = a[ijk + el];
    }
    if (j == nx - 1) {
      b[(i + 1) + (nx2 - 1) * nx2 + (k + 1) * nx2 * nx2 + el2] = a[ijk + el];
    }
    if (k == 0) {
      b[(i + 1) + (j + 1) * nx2 + 0 * nx2 * nx2 + el2] = a[ijk + el];
    }
    if (k == nx - 1) {
      b[(i + 1) + (j + 1) * nx2 + (nx2 - 1) * nx2 * nx2 + el2] = a[ijk + el];
    }
  }
}

/**
 * Device kernel for gradient_jump_penalty_finalize
 * @note One threadgroup per element.
 */
kernel void gradient_jump_penalty_finalize_kernel(
    device float *penalty_d [[ buffer(0) ]],
    device const float *penalty_facet_d [[ buffer(1) ]],
    device const float *dphidxi_d [[ buffer(2) ]],
    constant int &nx [[ buffer(3) ]],
    uint e [[ threadgroup_position_in_grid ]],
    uint idx [[ thread_position_in_threadgroup ]],
    uint nthrds [[ threads_per_threadgroup ]]) {
  const int nx2 = nx + 2;
  const int el2 = (int) e * nx2 * nx2 * nx2;
  const int el = (int) e * nx * nx * nx;

  for (int ijk = (int) idx; ijk < nx * nx * nx; ijk += (int) nthrds) {
    const int jk = ijk / nx;
    const int i = ijk - jk * nx;
    const int k = jk / nx;
    const int j = jk - k * nx;
    penalty_d[ijk + el] =
      penalty_facet_d[0 + (j + 1) * nx2 + (k + 1) * nx2 * nx2 + el2]
      * dphidxi_d[0 + i * nx] +
      penalty_facet_d[(nx2 - 1) + (j + 1) * nx2 + (k + 1) * nx2 * nx2 + el2]
      * dphidxi_d[nx - 1 + i * nx] +
      penalty_facet_d[(i + 1) + 0 * nx2 + (k + 1) * nx2 * nx2 + el2]
      * dphidxi_d[0 + j * nx] +
      penalty_facet_d[(i + 1) + (nx2 - 1) * nx2 + (k + 1) * nx2 * nx2 + el2]
      * dphidxi_d[nx - 1 + j * nx] +
      penalty_facet_d[(i + 1) + (j + 1) * nx2 + 0 * nx2 * nx2 + el2]
      * dphidxi_d[0 + k * nx] +
      penalty_facet_d[(i + 1) + (j + 1) * nx2 + (nx2 - 1) * nx2 * nx2 + el2]
      * dphidxi_d[nx - 1 + k * nx];
  }
}
