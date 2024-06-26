#ifndef __MATH_FDM_KERNEL_CL__
#define __MATH_FDM_KERNEL_CL__
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

#define DEFINE_FDM_DO_FAST_KERNEL(NL)                                          \
__kernel void fdm_do_fast_kernel_nl##NL(__global real * __restrict__ e,        \
                                        __global real * __restrict__ r,        \
                                        __global real * __restrict__ s,        \
                                        __global real * __restrict__ d) {      \
                                                                               \
  __local real shwork[NL*NL*NL];                                               \
  __local real shwork2[NL*NL*NL];                                              \
  __local real A[NL*NL];                                                       \
  __local real Bt[NL*NL];                                                      \
  __local real Ct[NL*NL];                                                      \
                                                                               \
  const int idx = get_local_id(0);                                             \
  const int str = get_local_size(0);                                           \
  const int el = get_group_id(0);                                              \
  if( idx < NL*NL){                                                            \
     A[idx] = s[idx+NL*NL+el*NL*NL*3*2];                                       \
    Bt[idx] = s[idx+2*NL*NL+el*NL*NL*3*2];                                     \
    Ct[idx] = s[idx+2*2*NL*NL+el*NL*NL*3*2];                                   \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for (int ii = idx; ii< NL*NL*NL; ii += str){                                 \
    real tmp = 0.0;                                                            \
    int j = ii/NL;                                                             \
    int i = ii - j*NL;                                                         \
    for( int l = 0; l < NL; l++){                                              \
      tmp += A[i+l*NL]*r[l+NL*j+el*NL*NL*NL];                                  \
    }                                                                          \
    shwork[ii] = tmp;                                                          \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for (int ijk = idx; ijk< NL*NL*NL; ijk += str){                              \
    const int jk = ijk / NL;                                                   \
    const int i = ijk - jk * NL;                                               \
    const int k = jk / NL;                                                     \
    const int j = jk - k * NL;                                                 \
    real tmp = 0.0;                                                            \
    const int ik2 = i + k*NL*NL;                                               \
    for( int l = 0; l < NL; l++){                                              \
      tmp += Bt[l+j*NL]*shwork[l*NL+ik2];                                      \
    }                                                                          \
    shwork2[ijk] = tmp;                                                        \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for (int ijk = idx; ijk< NL*NL*NL; ijk += str){                              \
    const int jk = ijk / NL;                                                   \
    const int i = ijk - jk * NL;                                               \
    const int k = jk / NL;                                                     \
    const int j = jk - k * NL;                                                 \
    real tmp = 0.0;                                                            \
    const int ij2 = i + j*NL;                                                  \
    for( int l = 0; l < NL; l++){                                              \
      tmp += Ct[l+k*NL]*shwork2[ij2 + l*NL*NL];                                \
    }                                                                          \
    r[ijk+el*NL*NL*NL] = tmp*d[ijk+el*NL*NL*NL];                               \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  if( idx < NL*NL){                                                            \
     A[idx] = s[idx+el*NL*NL*3*2];                                             \
    Bt[idx] = s[idx+NL*NL+2*NL*NL+el*NL*NL*3*2];                               \
    Ct[idx] = s[idx+NL*NL+2*2*NL*NL+el*NL*NL*3*2];                             \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for (int ii = idx; ii< NL*NL*NL; ii += str){                                 \
    real tmp = 0.0;                                                            \
    int j = ii/NL;                                                             \
    int i = ii - j*NL;                                                         \
    for( int l = 0; l < NL; l++){                                              \
      tmp += A[i+l*NL]*r[l+NL*j+el*NL*NL*NL];                                  \
    }                                                                          \
    shwork[ii] = tmp;                                                          \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for (int ijk = idx; ijk< NL*NL*NL; ijk += str){                              \
    const int jk = ijk / NL;                                                   \
    const int i = ijk - jk * NL;                                               \
    const int k = jk / NL;                                                     \
    const int j = jk - k * NL;                                                 \
    real tmp = 0.0;                                                            \
    const int ik2 = i + k*NL*NL;                                               \
    for( int l = 0; l < NL; l++){                                              \
      tmp += Bt[l+j*NL]*shwork[l*NL+ik2];                                      \
    }                                                                          \
    shwork2[ijk] = tmp;                                                        \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for (int ijk = idx; ijk< NL*NL*NL; ijk += str){                              \
    const int jk = ijk / NL;                                                   \
    const int i = ijk - jk * NL;                                               \
    const int k = jk / NL;                                                     \
    const int j = jk - k * NL;                                                 \
    real tmp = 0.0;                                                            \
    const int ij2 = i + j*NL;                                                  \
    for( int l = 0; l < NL; l++){                                              \
      tmp += Ct[l+k*NL]*shwork2[ij2 + l*NL*NL];                                \
    }                                                                          \
    e[ijk+el*NL*NL*NL] = tmp;                                                  \
  }                                                                            \
}

DEFINE_FDM_DO_FAST_KERNEL(2)
DEFINE_FDM_DO_FAST_KERNEL(3)
DEFINE_FDM_DO_FAST_KERNEL(4)
DEFINE_FDM_DO_FAST_KERNEL(5)
DEFINE_FDM_DO_FAST_KERNEL(6)
DEFINE_FDM_DO_FAST_KERNEL(7)
DEFINE_FDM_DO_FAST_KERNEL(8)
DEFINE_FDM_DO_FAST_KERNEL(9)
DEFINE_FDM_DO_FAST_KERNEL(10)
DEFINE_FDM_DO_FAST_KERNEL(11)
DEFINE_FDM_DO_FAST_KERNEL(12)
DEFINE_FDM_DO_FAST_KERNEL(13)
DEFINE_FDM_DO_FAST_KERNEL(14)
DEFINE_FDM_DO_FAST_KERNEL(15)


#endif // __MATH_FDM_KERNEL_CL__
