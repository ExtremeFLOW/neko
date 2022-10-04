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
__global__ void makeabf_kernel(T * __restrict__ abx1,
                               T * __restrict__ aby1,
                               T * __restrict__ abz1,
                               T * __restrict__ abx2,
                               T * __restrict__ aby2,
                               T * __restrict__ abz2,
                               T * __restrict__ bfx,
                               T * __restrict__ bfy,
                               T * __restrict__ bfz,
                               const T rho,
                               const T ab1,
                               const T ab2,
                               const T ab3,
                               const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    T ta1_val = ab2 * abx1[i] + ab3 * abx2[i];
    T ta2_val = ab2 * aby1[i] + ab3 * aby2[i];
    T ta3_val = ab2 * abz1[i] + ab3 * abz2[i];

    abx2[i] = abx1[i];
    aby2[i] = aby1[i];
    abz2[i] = abz1[i];
    abx1[i] = bfx[i];
    aby1[i] = bfy[i];
    abz1[i] = bfz[i];

    bfx[i] = (ab1 * bfx[i] + ta1_val) * rho;
    bfy[i] = (ab1 * bfy[i] + ta2_val) * rho;
    bfz[i] = (ab1 * bfz[i] + ta3_val) * rho;
  }
  
}

