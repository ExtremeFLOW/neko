#ifndef __SOURCE_TERMS_GRADIENT_JUMP_PENALTY_KERNEL_CL__
#define __SOURCE_TERMS_GRADIENT_JUMP_PENALTY_KERNEL_CL__
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
 * Device kernel for pick_facet_value_hex
 * @note One work-group per element.
 */
__kernel void pick_facet_value_hex_kernel(__global real * __restrict__ b,
                                          __global const real * __restrict__ a,
                                          const int nx) {

  const int idx = get_local_id(0);
  const int nx2 = nx + 2;
  const int el2 = get_group_id(0) * nx2 * nx2 * nx2;
  const int el = get_group_id(0) * nx * nx * nx;

  for (int ijk = idx; ijk < nx * nx * nx; ijk += get_local_size(0)) {
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
 * @note One work-group per element.
 */
__kernel void gradient_jump_penalty_finalize_kernel(
    __global real * __restrict__ penalty_d,
    __global const real * __restrict__ penalty_facet_d,
    __global const real * __restrict__ dphidxi_d,
    const int nx) {

  const int idx = get_local_id(0);
  const int nx2 = nx + 2;
  const int el2 = get_group_id(0) * nx2 * nx2 * nx2;
  const int el = get_group_id(0) * nx * nx * nx;

  for (int ijk = idx; ijk < nx * nx * nx; ijk += get_local_size(0)) {
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

#endif // __SOURCE_TERMS_GRADIENT_JUMP_PENALTY_KERNEL_CL__
