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

/**
 * Device kernel for schwarz extrude
 * We sum the "shell" of a1 that is l1 steps in scaled with f1
 * with the shell of a2 that is l2 steps in and scale with f3.
 * Right now we just throw away all arrays that are not
 * on the first face of dimension (nx-2)(nx-2)
 * It should be noted that a1, a2 are often the same array.
 */
template< typename T, const int NX>
__global__ void schwarz_extrude_kernel(T * a1,
                                       const int l1,
                                       const T f1,
                                       T * a2,
                                       const int l2,
                                       const T f2){

  const int idx = threadIdx.x;
  const int el = blockIdx.x*NX*NX*NX;
  const int x = idx%(NX-2) + 1;
  const int y = idx/(NX-2) + 1;
  int idx1,idx2;

  idx1 = l1 + x*NX + y*NX*NX + el;
  idx2 = l2 + x*NX + y*NX*NX + el;
  a1[idx1] = f1*a1[idx1] + f2*a2[idx2];

  idx1 = (NX-1-l1) + x*NX + y*NX*NX + el;
  idx2 = (NX-1-l2) + x*NX + y*NX*NX + el;
  a1[idx1] = f1*a1[idx1] + f2*a2[idx2];

  __syncthreads();

  idx1 = x + l1*NX + y*NX*NX + el;
  idx2 = x + l2*NX + y*NX*NX + el;
  a1[idx1] = f1*a1[idx1] + f2*a2[idx2];

  idx1 = x + (NX-1-l1)*NX + y*NX*NX + el;
  idx2 = x + (NX-1-l2)*NX + y*NX*NX + el;
  a1[idx1] = f1*a1[idx1] + f2*a2[idx2];

  __syncthreads();

  idx1 = x + y*NX + l1*NX*NX + el;
  idx2 = x + y*NX + l2*NX*NX + el;
  a1[idx1] = f1*a1[idx1] + f2*a2[idx2];

  idx1 = x + y*NX + (NX-1-l1)*NX*NX + el;
  idx2 = x + y*NX + (NX-1-l2)*NX*NX + el;
  a1[idx1] = f1*a1[idx1] + f2*a2[idx2];
}


/**
 * Device kernel for schwarz extrude
 */
template< typename T>
__global__ void schwarz_toext3d_kernel(T * __restrict__ a,
                                       T *__restrict__ b,
                                       const int nx) {

  const int idx = threadIdx.x;
  const int nx2 = nx+2;
  const int el2 = blockIdx.x*nx2*nx2*nx2;
  const int el = blockIdx.x*nx*nx*nx;
  for(int i = idx; i<nx2*nx2*nx2; i+=blockDim.x){
    a[i+el2] = 0.0;
  }
  __syncthreads();
  for(int ijk = idx; ijk<nx*nx*nx; ijk+=blockDim.x){
    const int jk = ijk / nx;
    const int i = ijk - jk * nx;
    const int k = jk / nx;
    const int j = jk - k * nx;
    a[(i+1)+(j+1)*nx2+(k+1)*nx2*nx2+el2] = b[ijk+el];
  }
}

template< typename T>
__global__ void schwarz_toreg3d_kernel(T * __restrict__ b,
                                       T *__restrict__ a,
                                       const int nx) {

  const int idx = threadIdx.x;
  const int nx2 = nx+2;
  const int el2 = blockIdx.x*nx2*nx2*nx2;
  const int el = blockIdx.x*nx*nx*nx;
  for(int ijk = idx; ijk<nx*nx*nx; ijk+=blockDim.x){
    const int jk = ijk / nx;
    const int i = ijk - jk * nx;
    const int k = jk / nx;
    const int j = jk - k * nx;
    b[ijk+el] = a[(i+1)+(j+1)*nx2+(k+1)*nx2*nx2+el2];
  }
}
