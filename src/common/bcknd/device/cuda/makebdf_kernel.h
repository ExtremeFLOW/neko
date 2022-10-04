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
__global__ void makebdf_kernel(const T * __restrict__ ulag1,
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
    T tb1_val = u[i] * B[i] * bd2;
    T tb2_val = v[i] * B[i] * bd2;
    T tb3_val = w[i] * B[i] * bd2;

    T ta1_val = ulag1[i] * B[i] * bd3;
    T ta2_val = vlag1[i] * B[i] * bd3;
    T ta3_val = wlag1[i] * B[i] * bd3;

    tb1_val += ta1_val;
    tb2_val += ta2_val;
    tb3_val += ta3_val;

    if (nbd == 3) {
      tb1_val += ulag2[i] * B[i] * bd4;
      tb2_val += vlag2[i] * B[i] * bd4;
      tb3_val += wlag2[i] * B[i] * bd4;
    }
    
    bfx[i] = bfx[i] + tb1_val * (rho / dt);
    bfy[i] = bfy[i] + tb2_val * (rho / dt);
    bfz[i] = bfz[i] + tb3_val * (rho / dt);
  }
  
}

