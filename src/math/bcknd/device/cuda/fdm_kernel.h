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

template< typename T, int NL >
__global__ void fdm_do_fast_kernel(T  * __restrict__ e,
                                   T * __restrict__ r,
                                   T * __restrict__ s,
                                   T * __restrict__ d) {
  __shared__ T shwork[NL*NL*NL];
  __shared__ T shwork2[NL*NL*NL];
  __shared__ T A[NL*NL];
  __shared__ T Bt[NL*NL];
  __shared__ T Ct[NL*NL];

  const int idx = threadIdx.x;
  const int str = blockDim.x;
  const int el = blockIdx.x;
  if( idx < NL*NL){
     A[idx] = s[idx+NL*NL+el*NL*NL*3*2];
    Bt[idx] = s[idx+2*NL*NL+el*NL*NL*3*2];
    Ct[idx] = s[idx+2*2*NL*NL+el*NL*NL*3*2];
  }    
  __syncthreads();

  for (int ii = idx; ii< NL*NL*NL; ii += str){
    T tmp = 0.0;
    int j = ii/NL;
    int i = ii - j*NL;
    for( int l = 0; l < NL; l++){
      tmp += A[i+l*NL]*r[l+NL*j+el*NL*NL*NL];
    }
    shwork[ii] = tmp;
  }
  __syncthreads();
  for (int ijk = idx; ijk< NL*NL*NL; ijk += str){
    const int jk = ijk / NL;
    const int i = ijk - jk * NL;
    const int k = jk / NL;
    const int j = jk - k * NL;
    T tmp = 0.0;
    const int ik2 = i + k*NL*NL; 
    for( int l = 0; l < NL; l++){
      tmp += Bt[l+j*NL]*shwork[l*NL+ik2];
    }
    shwork2[ijk] = tmp;
  }
  __syncthreads();
  for (int ijk = idx; ijk< NL*NL*NL; ijk += str){
    const int jk = ijk / NL;
    const int i = ijk - jk * NL;
    const int k = jk / NL;
    const int j = jk - k * NL;
    T tmp = 0.0;
    const int ij2 = i + j*NL; 
    for( int l = 0; l < NL; l++){
      tmp += Ct[l+k*NL]*shwork2[ij2 + l*NL*NL];
    }
    r[ijk+el*NL*NL*NL] = tmp*d[ijk+el*NL*NL*NL];
  }
  __syncthreads();
  if( idx < NL*NL){
     A[idx] = s[idx+el*NL*NL*3*2];
    Bt[idx] = s[idx+NL*NL+2*NL*NL+el*NL*NL*3*2];
    Ct[idx] = s[idx+NL*NL+2*2*NL*NL+el*NL*NL*3*2];
  }  
  __syncthreads();

  for (int ii = idx; ii< NL*NL*NL; ii += str){
    T tmp = 0.0;
    int j = ii/NL;
    int i = ii - j*NL;
    for( int l = 0; l < NL; l++){
      tmp += A[i+l*NL]*r[l+NL*j+el*NL*NL*NL];
    }
    shwork[ii] = tmp;
  }
  __syncthreads();
  for (int ijk = idx; ijk< NL*NL*NL; ijk += str){
    const int jk = ijk / NL;
    const int i = ijk - jk * NL;
    const int k = jk / NL;
    const int j = jk - k * NL;
    T tmp = 0.0;
    const int ik2 = i + k*NL*NL; 
    for( int l = 0; l < NL; l++){
      tmp += Bt[l+j*NL]*shwork[l*NL+ik2];
    }
    shwork2[ijk] = tmp;
  }
  __syncthreads();
  for (int ijk = idx; ijk< NL*NL*NL; ijk += str){
    const int jk = ijk / NL;
    const int i = ijk - jk * NL;
    const int k = jk / NL;
    const int j = jk - k * NL;
    T tmp = 0.0;
    const int ij2 = i + j*NL; 
    for( int l = 0; l < NL; l++){
      tmp += Ct[l+k*NL]*shwork2[ij2 + l*NL*NL];
    }
    e[ijk+el*NL*NL*NL] = tmp;
  }

}


