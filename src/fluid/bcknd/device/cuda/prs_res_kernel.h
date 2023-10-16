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
__global__ void prs_res_part1_kernel(T * __restrict__ ta1,
                                     T * __restrict__ ta2,
                                     T * __restrict__ ta3,
                                     const T * __restrict__ wa1,
                                     const T * __restrict__ wa2,
                                     const T * __restrict__ wa3,
                                     const T * __restrict__ f_u,
                                     const T * __restrict__ f_v,
                                     const T * __restrict__ f_w,
                                     const T * __restrict__ B,
                                     T * __restrict__ h1,
                                     const T mu,
                                     const T rho,
                                     const int n) {
  
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  const T inv_rho = 1.0 / rho;
  
  for (int i = idx; i < n; i += str) {
    h1[i] = inv_rho;
    ta1[i] = (f_u[i] / rho) - ((wa1[i] * (mu / rho)) * B[i]);
    ta2[i] = (f_v[i] / rho) - ((wa2[i] * (mu / rho)) * B[i]);
    ta3[i] = (f_w[i] / rho) - ((wa3[i] * (mu / rho)) * B[i]);
  }

}


template< typename T >
__global__ void prs_res_part2_kernel(T * __restrict__ p_res,
                                     const T * __restrict__ wa1,
                                     const T * __restrict__ wa2,
                                     const T * __restrict__ wa3,
                                     const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  
  for (int i = idx; i < n; i += str) {
    p_res[i] = (-p_res[i]) + (wa1[i] + wa2[i] + wa3[i]);
  }
}

template< typename T >
__global__ void prs_res_part3_kernel(T * __restrict__ p_res,
                                     const T * __restrict__ ta1,
                                     const T * __restrict__ ta2,
                                     const T * __restrict__ ta3,
                                     const T dtbd,
                                     const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  
  for (int i = idx; i < n; i += str) {
    p_res[i] = p_res[i] - (dtbd * (ta1[i] + ta2[i] + ta3[i]));
  }
}
