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
 * Kernel for back-substitution of x and update of p
 */
__kernel void gmres_part2_kernel(__global real * w,
                                 __global const real * __restrict__ v,
                                 __global const real * __restrict__ mult,
                                 __global const real * __restrict__ h,
                                 __global real * __restrict__ buf_h1,
                                 const int j,
                                 const int n) {
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  __local real buf1[256];
  real tmp1 = 0.0;

  for (int i = idx; i < n; i+= str) {
    real tmp = 0.0;
    for (int k = 0; k < j; k ++) {
      tmp += -h[k]*v[(k * n) + i];
    }
    w[i] += tmp;
    tmp1 += w[i]*w[i]*mult[i];
  }
  buf1[get_local_id(0)] = tmp1;
  barrier(CLK_LOCAL_MEM_FENCE);

  int i = get_local_size(0)>>1;
  while (i != 0) {
    if (get_local_id(0) < i) {
      buf1[get_local_id(0)] += buf1[get_local_id(0) + i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    i = i>>1;
  }
 
  if (get_local_id(0) == 0) {
    buf_h1[get_group_id(0)] = buf1[0];
  }

}



