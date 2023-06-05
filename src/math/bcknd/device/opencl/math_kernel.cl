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
 * Device kernel for cmult
 */
__kernel void cmult_kernel(__global real * __restrict__ a,
                           const real c,
                           const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = c * a[i];
  } 
}

/**
 * Device kernel for cmult2
 */
__kernel void cmult2_kernel(__global real * __restrict__ a,
               __global real * __restrict__ b,
                           const real c,
                           const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = c * b[i];
  } 
}

/**
 * Device kernel for cadd
 */
__kernel void cadd_kernel(__global real * __restrict__ a,
                          const real c,
                          const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] + c;
  } 
}

/**
 * Device kernel for cfill
 */
__kernel void cfill_kernel(__global real * __restrict__ a,
                           const real c,
                           const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = c;
  } 
}

/**
 * Device kernel for add2
 */
__kernel void add2_kernel(__global real * __restrict__ a,
                          __global const real * __restrict__ b,
                          const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] + b[i];
  } 
}

/**
 * Device kernel for add2s1
 */
__kernel void add2s1_kernel(__global real * __restrict__ a,
                            __global const real * __restrict__ b,
                            const real c1,
                            const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = c1 * a[i] + b[i];
  } 
}

/**
 * Device kernel for add2s2
 */
__kernel void add2s2_kernel(__global real * __restrict__ a,
                            __global const real * __restrict__ b,
                            const real c1,
                            const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] + c1 * b[i];
  }
}

/**
 * Device kernel for add2s2 many
 */
__kernel void add2s2_many_kernel(__global  real * __restrict__  x,
                                 __global const real * __restrict__  p,
                                 __global const real * __restrict__ alpha,
                                 const int p_cur,
                                 const int n) {

  const int idx = get_global_id(0);
  const int str = get_local_size(0) * get_num_groups(0);

  for (int i = idx; i < n; i+= str) {
    real tmp = 0.0;
    for (int j = 0; j < p_cur; j ++) {
      tmp += p[j * n + i]*alpha[j];
    }
    x[i] += tmp;
  }
}

/**
 * Device kernel for addsqr2s2
 */
__kernel void addsqr2s2_kernel(__global real * __restrict__ a,
                               __global const real * __restrict__ b,
                               const real c1,
                               const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] + c1 * (b[i] * b[i]);
  }
}

/**
 * Device kernel for add3s2
 */
__kernel void add3s2_kernel(__global real * __restrict__ a,
                            __global const real * __restrict__ b,
                            __global const real * __restrict__ c,
                            const real c1,
                            const real c2,
                            const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = c1 * b[i] + c2 * c[i];
  }
}

/**
 * Device kernel for invcol1
 */
__kernel void invcol1_kernel(__global real * __restrict__ a,
                             const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);
  const real one = 1.0;
  
  for (int i = idx; i < n; i += str) {
    a[i] = one / a[i];
  }
}

/** 
 * Device kernel for invcol2
 */
__kernel void invcol2_kernel(__global real * __restrict__ a,
                             __global const real * __restrict__ b,
                             const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] / b[i];
  }  
}

/** 
 * Device kernel for col2
 */
__kernel void col2_kernel(__global real * __restrict__ a,
                          __global const real * __restrict__ b,
                          const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] * b[i];
  }  
}

/** 
 * Device kernel for col3
 */
__kernel void col3_kernel(__global real * __restrict__ a,
                          __global const real * __restrict__ b,
                          __global const real * __restrict__ c,
                          const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = b[i] * c[i];
  }  
}

/** 
 * Device kernel for subcol3
 */
__kernel void subcol3_kernel(__global real * __restrict__ a,
                             __global const real * __restrict__ b,
                             __global const real * __restrict__ c,
                             const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] - b[i] * c[i];
  }  
}

/** 
 * Device kernel for sub2
 */
__kernel void sub2_kernel(__global real * __restrict__ a,
                          __global const real * __restrict__ b,
                          const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] - b[i];
  }  
}

/** 
 * Device kernel for sub3
 */
__kernel void sub3_kernel(__global real * __restrict__ a,
                          __global const real * __restrict__ b,
                          __global const real * __restrict__ c,
                          const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = b[i] - c[i];
  }  
}

/** 
 * Device kernel for addcol3
 */
__kernel void addcol3_kernel(__global real * __restrict__ a,
                             __global const real * __restrict__ b,
                             __global const real * __restrict__ c,
                             const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] + b[i] * c[i];
  }  
}

/** 
 * Device kernel for addcol4
 */
__kernel void addcol4_kernel(__global real * __restrict__ a,
                             __global const real * __restrict__ b,
                             __global const real * __restrict__ c,
                             __global const real * __restrict__ d,
                             const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] + b[i] * c[i] * d[i];
  }  
}

/**
 * Device kernel for vdot3
 */
__kernel void vdot3_kernel(__global real * __restrict__ dot,
                           __global const real * __restrict__ u1,
                           __global const real * __restrict__ u2,
                           __global const real * __restrict__ u3,
                           __global const real * __restrict__ v1,
                           __global const real * __restrict__ v2,
                           __global const real * __restrict__ v3,
                           const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    dot[i] = u1[i] * v1[i]  + u2[i] * v2[i] + u3[i] * v3[i];
  }  
  
}

/**
 * Device kernel for glsc3
 */
__kernel void glsc3_kernel(__global const real * __restrict__ a,
                           __global const real * __restrict__ b,
                           __global const real * __restrict__ c,
                           __global real * __restrict__ buf_h,
                           const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  __local real buf[256]; /* Make this nice...*/
  real tmp = 0.0;

  for (int i = idx; i < n; i+= str) {
    tmp += a[i] * b[i] * c[i];
  }
  buf[get_local_id(0)] = tmp;
  barrier(CLK_LOCAL_MEM_FENCE);

  int i = (get_local_size(0))>>1;
  while (i != 0) {
    if (get_local_id(0) < i) {
      buf[get_local_id(0)] += buf[get_local_id(0) + i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    i = i>>1;
  }
 
  if (get_local_id(0) == 0) {
    buf_h[get_group_id(0)] = buf[0];
  }

}

/**
 * Device kernel for glsc3 many
 */
__kernel void glsc3_many_kernel(__global const real * __restrict__ a,
                                __global const real * __restrict__ b,
                                __global const real * __restrict__ c,
                                __global real * __restrict__ buf_h,
                                const int j,
                                const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);
  const int y = get_local_id(1);

  __local real buf[256]; /* Make this nice...*/
  real tmp = 0;
  if(y < j){
    for (int i = idx; i < n; i+= str) {
      tmp += a[i] * b[get_local_id(1) * n + i] * c[i];
    }
  }

  buf[get_local_id(0) * get_local_size(1) + y] = tmp;
  barrier(CLK_LOCAL_MEM_FENCE);

  int i = get_local_size(0)>>1;
  while (i != 0) {
    if (get_local_id(0) < i) {
      buf[get_local_id(0) * get_local_size(1) + y] +=
        buf[(get_local_id(0) + i) * get_local_size(1) + y];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    i = i>>1;
  }
  
  if (get_local_id(0) == 0) {
    if( y < j) {
      buf_h[j * get_group_id(0) + y] = buf[y];
    }
  }
}

/**
 * Device kernel for glsc2
 */
__kernel void glsc2_kernel(__global const real * __restrict__ a,
                           __global const real * __restrict__ b,
                           __global real * __restrict__ buf_h,
                           const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  __local real buf[256]; /* Make this nice...*/
  real tmp = 0.0;

  for (int i = idx; i < n; i+= str) {
    tmp += a[i] * b[i];
  }
  buf[get_local_id(0)] = tmp;
  barrier(CLK_LOCAL_MEM_FENCE);

  int i = (get_local_size(0))>>1;
  while (i != 0) {
    if (get_local_id(0) < i) {
      buf[get_local_id(0)] += buf[get_local_id(0) + i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    i = i>>1;
  }
 
  if (get_local_id(0) == 0) {
    buf_h[get_group_id(0)] = buf[0];
  }

}

/**
 * Device kernel for glsum
 */
__kernel void glsum_kernel(__global const real * __restrict__ a,
			   __global real * __restrict__ buf_h,
                           const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  __local real buf[256]; /* Make this nice...*/
  real tmp = 0.0;

  for (int i = idx; i < n; i+= str) {
    tmp += a[i];
  }
  buf[get_local_id(0)] = tmp;
  barrier(CLK_LOCAL_MEM_FENCE);

  int i = (get_local_size(0))>>1;
  while (i != 0) {
    if (get_local_id(0) < i) {
      buf[get_local_id(0)] += buf[get_local_id(0) + i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    i = i>>1;
  }
 
  if (get_local_id(0) == 0) {
    buf_h[get_group_id(0)] = buf[0];
  }

}
