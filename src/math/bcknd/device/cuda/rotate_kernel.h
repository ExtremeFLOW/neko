#ifndef __ROTATE_KERNEL_H__
#define __ROTATE_KERNEL_H__
/*
 Copyright (c) 2021-2022, The Neko Authors
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
__global__ void rotate_cyc_kernel(
			   T * __restrict__ vx,
			   T * __restrict__ vy,
			   T * __restrict__ vz,
			   const T * __restrict__ x,
			   const T * __restrict__ y,
			   const T * __restrict__ z,
			   const int * __restrict__ cyc_msk,
			   const T * __restrict__ R11,
			   const T * __restrict__ R12,
			   const int  ncyc,
			   const int  idir) { 
  
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < ncyc-1; i += str) {
    const int j  = cyc_msk[i+1]-1;
    const T vxj  = vx[j];
    const T vyj  = vy[j];
    const T R11i = R11[i];
    const T R12i = R12[i];
    T vnor;
    T vtan;
    if (idir==1){
        vnor =  vxj * R11i + vyj * R12i;
        vtan = -vxj * R12i + vyj * R11i;
    }
    else if (idir==0) {
        vnor =  vxj * R11i - vyj * R12i;
        vtan =  vxj * R12i + vyj * R11i;
    }

    vx[j] = vnor;
    vy[j] = vtan;
  }


}
#endif


