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

#ifndef __COMMON_MAKEOIFS_KERNEL__
#define __COMMON_MAKEOIFS_KERNEL__

template< typename T >
__global__ void makeoifs_kernel(const T * __restrict__ phix,
                               const T * __restrict__ phiy,
                               const T * __restrict__ phiz,
                               T * __restrict__ bfx,
                               T * __restrict__ bfy,
                               T * __restrict__ bfz,
                               const T rho,
                               const T dt,
                               const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
        
    bfx[i] = bfx[i] + phix[i] * (rho / dt);
    bfy[i] = bfy[i] + phiy[i] * (rho / dt);
    bfz[i] = bfz[i] + phiz[i] * (rho / dt);
  }
  
}

template< typename T >
__global__ void scalar_makeoifs_kernel(const T * __restrict__ phis,
                                      T * __restrict__ bfs,
                                      const T rho,
                                      const T dt,
                                      const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    
    bfs[i] = bfs[i] + phis[i] * (rho / dt);
  }
  
}

#endif // __COMMON_MAKEOIFS_KERNEL__
