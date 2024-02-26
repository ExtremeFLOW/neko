#ifndef __MATH_SCHWARZ_KERNEL_CL__
#define __MATH_SCHWARZ_KERNEL_CL__
/*
 Copyright (c) 2022-2023, The Neko Authors
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
#define DEFINE_SCHWARZ_EXTRUDE_KERNEL(NX)                                      \
__kernel void schwarz_extrude_kernel_nx##NX(__global real * a1,                \
                                     const int l1,                             \
                                     const real f1,                            \
                                     __global real * a2,                       \
                                     const int l2,                             \
                                     const real f2) {                          \
                                                                               \
  __local real x1[NX*NX],x2[NX*NX];                                            \
  __local real y1[NX*NX],y2[NX*NX];                                            \
  __local real z1[NX*NX],z2[NX*NX];                                            \
  const int idx = get_local_id(0);                                             \
  const int el = get_group_id(0) * NX*NX*NX;                                   \
                                                                               \
  if(idx < NX*NX){                                                             \
    int i = idx%NX;                                                            \
    int k = idx/NX;                                                            \
    int idx_x1 = l2 + i*NX + k*NX*NX + el;                                     \
    x1[idx]=a2[idx_x1];                                                        \
    int idx_x2 = NX-1-l2 + i*NX + k*NX*NX + el;                                \
    x2[idx]=a2[idx_x2];                                                        \
                                                                               \
    int idx_y1 = i + l2*NX + k*NX*NX + el;                                     \
    y1[idx]=a2[idx_y1];                                                        \
    int idx_y2 = i + (NX-1-l2)*NX + k*NX*NX + el;                              \
    y2[idx]=a2[idx_y2];                                                        \
                                                                               \
    int idx_z1 = i + k*NX + l2*NX*NX + el;                                     \
    z1[idx]=a2[idx_z1];                                                        \
    int idx_z2 = i + k*NX + (NX-l2-1)*NX*NX + el;                              \
    z2[idx]=a2[idx_z2];                                                        \
  }                                                                            \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for(int ijk = idx; ijk<NX*NX*NX; ijk+=get_local_size(0)){                    \
    int jk = ijk/NX;                                                           \
     int i = ijk - NX*jk;                                                      \
     int k = jk/NX;                                                            \
     int j = jk -k*NX;                                                         \
     if(j>0 && j< NX-1 && k > 0 && k < NX -1){                                 \
       int idx1 = i + j*NX + k*NX*NX + el;                                     \
       if(i == l1){                                                            \
         int idx2 = j + k*NX;                                                  \
         a1[idx1] = f1*a1[idx1] + f2*x1[idx2];                                 \
       }                                                                       \
       if(i == NX-1-l1){                                                       \
         int idx2 = j + k*NX;                                                  \
         a1[idx1] = f1*a1[idx1] + f2*x2[idx2];                                 \
       }                                                                       \
     }                                                                         \
     if( i > 0 && i < NX-1 && k > 0 && k < NX -1){                             \
       int idx1 = i + j*NX + k*NX*NX + el;                                     \
       if(j == l1){                                                            \
         int idx2 = i + k*NX;                                                  \
         a1[idx1] = f1*a1[idx1] + f2*y1[idx2];                                 \
       }                                                                       \
       if(j == NX-1-l1){                                                       \
         int idx2 = i + k*NX;                                                  \
         a1[idx1] = f1*a1[idx1] + f2*y2[idx2];                                 \
       }                                                                       \
     }                                                                         \
     if( i > 0 && i < NX-1 && j>0 && j< NX-1 ){                                \
       int idx1 = i + j*NX + k*NX*NX + el;                                     \
       if(k == l1){                                                            \
         int idx2 = i + j*NX;                                                  \
         a1[idx1] = f1*a1[idx1] + f2*z1[idx2];                                 \
       }                                                                       \
       if(k == NX-1-l1){                                                       \
         int idx2 = i + j*NX;                                                  \
         a1[idx1] = f1*a1[idx1] + f2*z2[idx2];                                 \
       }                                                                       \
     }                                                                         \
  }                                                                            \
}                                                                              

DEFINE_SCHWARZ_EXTRUDE_KERNEL(2);
DEFINE_SCHWARZ_EXTRUDE_KERNEL(3);
DEFINE_SCHWARZ_EXTRUDE_KERNEL(4);
DEFINE_SCHWARZ_EXTRUDE_KERNEL(5);
DEFINE_SCHWARZ_EXTRUDE_KERNEL(6);
DEFINE_SCHWARZ_EXTRUDE_KERNEL(7);
DEFINE_SCHWARZ_EXTRUDE_KERNEL(8);
DEFINE_SCHWARZ_EXTRUDE_KERNEL(9);
DEFINE_SCHWARZ_EXTRUDE_KERNEL(10);
DEFINE_SCHWARZ_EXTRUDE_KERNEL(11);
DEFINE_SCHWARZ_EXTRUDE_KERNEL(12);
DEFINE_SCHWARZ_EXTRUDE_KERNEL(13);
DEFINE_SCHWARZ_EXTRUDE_KERNEL(14);
DEFINE_SCHWARZ_EXTRUDE_KERNEL(15);
DEFINE_SCHWARZ_EXTRUDE_KERNEL(16);

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

  barrier(CLK_LOCAL_MEM_FENCE);
  
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

#endif // __MATH_SCHWARZ_KERNEL_CL__
