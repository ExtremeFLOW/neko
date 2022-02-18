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

template< typename T >
__global__ void sumab_kernel(T * __restrict__ u,
                             T * __restrict__ v,
                             T * __restrict__ w,
                             const T * __restrict__ uu,
                             const T * __restrict__ vv,
                             const T * __restrict__ ww,
                             const T * __restrict__ ulag1,
                             const T * __restrict__ ulag2,
                             const T * __restrict__ vlag1,
                             const T * __restrict__ vlag2,
                             const T * __restrict__ wlag1,
                             const T * __restrict__ wlag2,
                             const real ab1,
                             const real ab2,
                             const real ab3,
                             const int nab,
                             const int n) {
  
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    u[i] = ab1 * uu[i] + ab2 * ulag1[i];
    v[i] = ab1 * vv[i] + ab2 * vlag1[i];
    w[i] = ab1 * ww[i] + ab2 * wlag1[i];
  }

  if (nab == 3) {
    for (int i = idx; i < n; i += str) {
      u[i] = u[i] + ab3 * ulag2[i];
      v[i] = v[i] + ab3 * vlag2[i];
      w[i] = w[i] + ab3 * wlag2[i];
    } 
  }
  
}

