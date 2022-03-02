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
 * Device kernel for schwarz extrude
 * This can probably be done i a better way...
 * We sum the "shell" of a1 that is l1 steps in scaled with f1
 * with the shell of a2 that is l2 steps in and scale with f3.
 * Right now we just thorw away all arrays that are not
 * on the first face of dimension (nx-2)(nx-2)
 * It should be noted that a1, a2 are often the same array.
 * If l1,l2 are not the same or if one is not 0 this might lead to a race.
 */
__kernel void schwarz_extrude_kernel(__global real * a1,
                                     const int l1,
                                     const real f1,
                                     __global real * a2,
                                     const int l2,
                                     const real f2,
                                     const int nx) {

  const int idx = get_local_id(0);
  const int el = get_group_id(0) * nx*nx*nx;

  for(int ijk = idx; ijk<nx*nx*nx; ijk+=get_local_size(0)){
     int jk = ijk/nx;
     int i = ijk - nx*jk;
     int k = jk/nx;
     int j = jk -k*nx;
     if(j>0 && j< nx-1 && k > 0 && k < nx -1){
       int idx1 = i + j*nx + k*nx*nx + el;
       if(i == l1){
         int idx2 = l2 + j*nx + k*nx*nx + el;
         a1[idx1] = f1*a1[idx1] + f2*a2[idx2];
       }
       if(i == nx-1-l1){
         int idx2 = nx-1-l2 + j*nx + k*nx*nx + el;
         a1[idx1] = f1*a1[idx1] + f2*a2[idx2];
       }
     }
     if( i > 0 && i < nx-1 && k > 0 && k < nx -1){
       int idx1 = i + j*nx + k*nx*nx + el;
       if(j == l1){
         int idx2 = i + l2*nx + k*nx*nx + el;
         a1[idx1] = f1*a1[idx1] + f2*a2[idx2];
       }
       if(j == nx-1-l1){
         int idx2 = i + (nx-1-l2)*nx + k*nx*nx + el;
         a1[idx1] = f1*a1[idx1] + f2*a2[idx2];
       }
     }
     if( i > 0 && i < nx-1 && j>0 && j< nx-1 ){
       int idx1 = i + j*nx + k*nx*nx + el;
       if(k == l1){
         int idx2 = i + j*nx + l2*nx*nx + el;
         a1[idx1] = f1*a1[idx1] + f2*a2[idx2];
       }
       if(k == nx-1-l1){
         int idx2 = i + j*nx + (nx-l2-1)*nx*nx + el;
         a1[idx1] = f1*a1[idx1] + f2*a2[idx2];
       }
     }
  }    
}

/**
 * Device kernel for schwarz extrude
 */
__kernel void schwarz_toext3d_kernel(__global real * __restrict__ a,
                                     __global real * __restrict__ b,
                                     const int nx) {

  const int idx = get_local_id(0);
  const int nx2 = nx+2;
  const int el2 = get_group_id(0) * nx2*nx2*nx2;
  const int el = get_group_id(0) * nx*nx*nx;
  for(int i = idx; i<nx2*nx2*nx2; i+=get_local_size(0)){
    a[i+el2] = 0.0;
  }
  for(int ijk = idx; ijk<nx*nx*nx; ijk+=get_local_size(0)){
    const int jk = ijk / nx;
    const int i = ijk - jk * nx;
    const int k = jk / nx;
    const int j = jk - k * nx;
    a[(i+1)+(j+1)*nx2+(k+1)*nx2*nx2+el2] = b[ijk+el];
  }
}

__kernel void schwarz_toreg3d_kernel(__global real * __restrict__ b,
                                     __global real * __restrict__ a,
                                     const int nx) {

  const int idx = get_local_id(0);
  const int nx2 = nx+2;
  const int el2 = get_group_id(0) * nx2*nx2*nx2;
  const int el = get_group_id(0) * nx*nx*nx;
  for(int ijk = idx; ijk<nx*nx*nx; ijk+=get_local_size(0)){
    const int jk = ijk / nx;
    const int i = ijk - jk * nx;
    const int k = jk / nx;
    const int j = jk - k * nx;
    b[ijk+el] = a[(i+1)+(j+1)*nx2+(k+1)*nx2*nx2+el2];
  }
}
