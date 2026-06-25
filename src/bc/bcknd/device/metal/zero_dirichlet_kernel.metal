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
using namespace metal;

/**
 * Device kernel for scalar apply for a zero-Dirichlet condition
 */
kernel void zero_dirichlet_apply_scalar_kernel(
        device const int *msk[[ buffer(0) ]],
        device float *x[[ buffer(1) ]],
        constant int &m[[ buffer(2) ]],
        uint tid [[ thread_position_in_grid ]])
{
    const int i = (int)tid + 1;
    if (i >= m) return;
    const int k = msk[i] - 1;
    x[k] = 0.0f;
}

/**
 * Device kernel for vector apply for a zero-Dirichlet condition
 */
kernel void zero_dirichlet_apply_vector_kernel(
        device const int *msk[[ buffer(0) ]],
        device float *x[[ buffer(1) ]],
        device float *y[[ buffer(2) ]],
        device float *z[[ buffer(3) ]],
        constant int &m[[ buffer(4) ]],
        uint tid [[ thread_position_in_grid ]])
{
    const int i = (int)tid + 1;
    if (i >= m) return;
    const int k = msk[i] - 1;
    x[k] = 0.0f;
    y[k] = 0.0f;
    z[k] = 0.0f;
}
