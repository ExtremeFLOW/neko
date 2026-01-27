#ifndef __SOURCE_TERMS_MAPPING_KERNELS__
#define __SOURCE_TERMS_MAPPING_KERNELS__
/*
 Copyright (c) 2021-2024, The Neko Authors
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
 * Device kernel for smooth step function.
 *
 * @param x  in/out array
 * @param edge0  lower edge
 * @param edge1  upper edge
 * @param n  size of the input array
 */
__kernel void smooth_step_kernel(
    __global real* __restrict__ x, const real edge0, const real edge1,
    const int n) {

    const int idx = get_global_id(0);
    const int str = get_global_size(0);

    for (int i = idx; i < n; i += str) {
        real t = (x[i] - edge0) / (edge1 - edge0);
        t      = min(max(t, 0.0), 1.0);
        x[i]   = pow(t, 3.0) * (t * (t * 6.0 - 15.0) + 10.0);
    }
}

/**
 * Device kernel for step function.
 *
 * @param x  input array
 * @param edge  step location
 * @param left  value to the left of the step
 * @param right  value to the right of the step
 * @param n  size of the input array
 */
__kernel void step_kernel(
    __global real* __restrict__ x, const real edge, const real left,
    const real right, const int n) {

    const int idx = get_global_id(0);
    const int str = get_global_size(0);

    for (int i = idx; i < n; i += str) {
        const real x_i = x[i];
        x[i]           = (x_i < edge) ? left : right;
    }
}

/**
 * Device kernel for the inverse permeability kernel.
 *
 * @param x  input array
 * @param k_0  lower bound of the permeability
 * @param k_1  upper bound of the permeability
 * @param q  parameter
 * @param n  size of the input array
 */
__kernel void permeability_kernel(
    __global real* __restrict__ x, const real k_0, const real k_1, const real q,
    const int n) {

    const int idx = get_global_id(0);
    const int str = get_global_size(0);

    for (int i = idx; i < n; i += str) {
        const real x_i = x[i];
        x[i]           = k_0 + (k_1 - k_0) * x_i * (q + 1.0) / (q + x_i);
    }
}

#endif /* __SOURCE_TERMS_MAPPING_KERNELS__ */