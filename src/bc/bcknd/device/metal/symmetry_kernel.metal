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

kernel void symmetry_aligned_apply_vector_kernel(
        device const int *xmsk[[ buffer(0) ]],
        device const int *ymsk[[ buffer(1) ]],
        device const int *zmsk[[ buffer(2) ]],
        device float *x[[ buffer(3) ]],
        device float *y[[ buffer(4) ]],
        device float *z[[ buffer(5) ]],
        constant int &m[[ buffer(6) ]],
        constant int &n[[ buffer(7) ]],
        constant int &l[[ buffer(8) ]],
        uint tid [[ thread_position_in_grid ]])
{
    const int i = (int)tid;
    if (i < m) {
        const int k = xmsk[i + 1] - 1;
        x[k] = 0.0f;
    }
    if (i < n) {
        const int k = ymsk[i + 1] - 1;
        y[k] = 0.0f;
    }
    if (i < l) {
        const int k = zmsk[i + 1] - 1;
        z[k] = 0.0f;
    }
}
