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

#ifndef __COMMON_MAKEOIFS_KERNEL__
#define __COMMON_MAKEOIFS_KERNEL__

template< typename T >
__global__ void makeoifs_kernel(const T * __restrict__ phi_x,
                               const T * __restrict__ phi_y,
                               const T * __restrict__ phi_z,
                               T * __restrict__ bf_x,
                               T * __restrict__ bf_y,
                               T * __restrict__ bf_z,
                               const T rho,
                               const T dt,
                               const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {

    bf_x[i] = bf_x[i] + phi_x[i] * (rho / dt);
    bf_y[i] = bf_y[i] + phi_y[i] * (rho / dt);
    bf_z[i] = bf_z[i] + phi_z[i] * (rho / dt);
  }

}

template< typename T >
__global__ void scalar_makeoifs_kernel(const T * __restrict__ phi_s,
                                      T * __restrict__ bf_s,
                                      const T rho,
                                      const T dt,
                                      const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {

    bf_s[i] = bf_s[i] + phi_s[i] * (rho / dt);
  }

}

#endif // __COMMON_MAKEOIFS_KERNEL__
