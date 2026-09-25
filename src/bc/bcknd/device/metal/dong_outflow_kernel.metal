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

kernel void dong_outflow_apply_scalar_kernel(
        device const int *msk[[ buffer(0) ]],
        device float *x[[ buffer(1) ]],
        device const float *normal_x[[ buffer(2) ]],
        device const float *normal_y[[ buffer(3) ]],
        device const float *normal_z[[ buffer(4) ]],
        device const float *u[[ buffer(5) ]],
        device const float *v[[ buffer(6) ]],
        device const float *w[[ buffer(7) ]],
        constant float &uinf[[ buffer(8) ]],
        constant float &delta[[ buffer(9) ]],
        constant int &m[[ buffer(10) ]],
        uint tid [[ thread_position_in_grid ]])
{
    const int i = (int)tid;
    if (i >= m) return;
    const int k = msk[i + 1] - 1;
    const float uk = u[k];
    const float vk = v[k];
    const float wk = w[k];
    const float vn = uk * normal_x[i] + vk * normal_y[i] + wk * normal_z[i];
    // Clamp the argument before tanh: in FP32 tanh is evaluated via exp(), and
    // exp() overflows to +inf for arguments beyond ~88, yielding inf/inf = NaN.
    // tanh saturates to +/-1 well within this range, so clamping is exact.
    const float arg = clamp(vn / (uinf * delta), -20.0f, 20.0f);
    const float S0 = 0.5f * (1.0f - tanh(arg));
    x[k] = -0.5f * (uk * uk + vk * vk + wk * wk) * S0;
}
