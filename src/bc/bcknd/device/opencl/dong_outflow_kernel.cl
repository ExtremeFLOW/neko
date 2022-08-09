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

/**
 * Device kernel for scalar apply for a dong outflow condition
 */
__kernel
void dong_outflow_apply_scalar_kernel(const int * __restrict__ msk,
				      real * __restrict__ x,
				      const real * normal_x,
				      const real * normal_y,
				      const real * normal_z,
				      const real * u,
				      const real * v,
				      const real * w,
				      const real uinf,
				      const real delta,
				      const int m) {
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < m; i += str) {
    const int k = msk[i + 1] - 1;
    const T uk = u[k];
    const T vk = v[k];
    const T wk = w[k];
    const T vn = uk*normal_x[i] + vk*normal_y[i] + wk*normal_z[i];
    const T S0 = 0.5*(1.0 - tanh(vn/(uinf*delta)));
    x[k] = -0.5*(uk*uk+vk*vk+wk*wk)*S0;
  }
}
