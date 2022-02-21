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
__global__ void makebdf_kernel(T * __restrict__ ta1,
                               T * __restrict__ ta2,
                               T * __restrict__ ta3,
                               T * __restrict__ tb1,
                               T * __restrict__ tb2,
                               T * __restrict__ tb3,
                               const T * __restrict__ ulag1,
                               const T * __restrict__ ulag2,
                               const T * __restrict__ vlag1,
                               const T * __restrict__ vlag2,
                               const T * __restrict__ wlag1,
                               const T * __restrict__ wlag2,
                               T * __restrict__ bfx,
                               T * __restrict__ bfy,
                               T * __restrict__ bfz,
                               const T * __restrict__ u,
                               const T * __restrict__ v,
                               const T * __restrict__ w,
                               const T * __restrict__ B,
                               const T rho,
                               const T dt,
                               const T bd2,
                               const T bd3,
                               const T bd4,
                               const int nbd,
                               const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    tb1[i] = u[i] * B[i] * bd2;
    tb2[i] = v[i] * B[i] * bd2;
    tb3[i] = w[i] * B[i] * bd2;

    ta1[i] = ulag1[i] * B[i] * bd3;
    ta2[i] = vlag1[i] * B[i] * bd3;
    ta3[i] = wlag1[i] * B[i] * bd3;

    tb1[i] += ta1[i];
    tb2[i] += ta2[i];
    tb3[i] += ta3[i];
  }

  if (nbd == 3) {
    for (int i = idx; i < n; i += str) {
      ta1[i] = ulag2[i] * B[i] * bd4;
      ta2[i] = vlag2[i] * B[i] * bd4;
      ta3[i] = wlag2[i] * B[i] * bd4;
      
      tb1[i] += ta1[i];
      tb2[i] += ta2[i];
      tb3[i] += ta3[i];
    }
  }

  for (int i = idx; i < n; i += str) {
    bfx[i] = bfx[i] + tb1[i] * (rho / dt);
    bfy[i] = bfy[i] + tb2[i] * (rho / dt);
    bfz[i] = bfz[i] + tb3[i] * (rho / dt);
  }


}

