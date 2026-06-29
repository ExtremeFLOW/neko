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

#include <metal_stdlib>
#include "bc_shader_utils.h"
using namespace metal;

kernel void facet_normal_apply_surfvec_kernel(
        device const int *msk[[ buffer(0) ]],
        device const int *facet[[ buffer(1) ]],
        device float *x[[ buffer(2) ]],
        device float *y[[ buffer(3) ]],
        device float *z[[ buffer(4) ]],
        device const float *u[[ buffer(5) ]],
        device const float *v[[ buffer(6) ]],
        device const float *w[[ buffer(7) ]],
        device const float *nx[[ buffer(8) ]],
        device const float *ny[[ buffer(9) ]],
        device const float *nz[[ buffer(10) ]],
        device const float *area[[ buffer(11) ]],
        constant int &lx[[ buffer(12) ]],
        constant int &m[[ buffer(13) ]],
        uint tid [[ thread_position_in_grid ]])
{
    const int i = (int)tid + 1;
    if (i >= m) return;
    const int k = msk[i] - 1;
    const int f = facet[i];

    int index[4];
    nonlinear_index(msk[i], lx, index);

    int na_idx;
    switch (f) {
    case 1: case 2:
        na_idx = coef_normal_area_idx(index[1], index[2], f, index[3], lx, 6);
        break;
    case 3: case 4:
        na_idx = coef_normal_area_idx(index[0], index[2], f, index[3], lx, 6);
        break;
    case 5: case 6:
        na_idx = coef_normal_area_idx(index[0], index[1], f, index[3], lx, 6);
        break;
    default:
        return;
    }

    x[k] = u[k] * nx[na_idx] * area[na_idx];
    y[k] = v[k] * ny[na_idx] * area[na_idx];
    z[k] = w[k] * nz[na_idx] * area[na_idx];
}
