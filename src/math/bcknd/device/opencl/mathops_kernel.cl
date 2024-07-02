#ifndef __MATH_MATHOPS_KERNEL_CL__
#define __MATH_MATHOPS_KERNEL_CL__
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
 * Device kernel for opchsign
 */
__kernel void opchsign_kernel(__global real *a1, 
			      __global real *a2, 
			      __global real *a3,
			      const int gdim,
			      const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  if (gdim == 3) {
    for (int i = idx; i < n; i += str) {
      a1[i] = -a1[i];
      a2[i] = -a2[i];
      a3[i] = -a3[i];
    }
  } 
  else {
    for (int i = idx; i < n; i += str) {
      a1[i] = -a1[i];
      a2[i] = -a2[i];
    }
  }

}

/**
 * Device kernel for opcolv
 */
__kernel void opcolv_kernel(__global real *a1, 
			    __global real *a2, 
			    __global real *a3,
			    __global const real *c,
			    const int gdim,
			    const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  if (gdim == 3) {
    for (int i = idx; i < n; i += str) {
      a1[i] = a1[i] * c[i];
      a2[i] = a2[i] * c[i];
      a3[i] = a3[i] * c[i];
    }
  } 
  else {
    for (int i = idx; i < n; i += str) {
      a1[i] = a1[i] * c[i];
      a2[i] = a2[i] * c[i];
    }
  }

}

/**
 * Device kernel for opcolv3c
 */
__kernel void opcolv3c_kernel(__global real *a1, 
			      __global real *a2, 
			      __global real *a3,
			      __global const real *b1, 
			      __global const real *b2, 
			      __global const real *b3,
			      __global const real *c,
			      const real d,
			      const int gdim,
			      const int n) {
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  if (gdim == 3) {
    for (int i = idx; i < n; i += str) {
      a1[i] = b1[i] * c[i] * d;
      a2[i] = b2[i] * c[i] * d;
      a3[i] = b3[i] * c[i] * d;
    }
  } 
  else {
    for (int i = idx; i < n; i += str) {
      a1[i] = b1[i] * c[i] * d;
      a2[i] = b2[i] * c[i] * d;
    }
  }

}

/**
 * Device kernel for opadd2cm
 */
__kernel void opadd2cm_kernel(__global real *a1, 
			      __global real *a2, 
			      __global real *a3,
			      __global const real *b1, 
			      __global const real *b2, 
			      __global const real *b3,
			      const real c,
			      const int gdim,
			      const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);
  
  if (gdim == 3) {
    for (int i = idx; i < n; i += str) {
      a1[i] = a1[i] + b1[i] * c;
      a2[i] = a2[i] + b2[i] * c;
      a3[i] = a3[i] + b3[i] * c;
    }
  } 
  else {
    for (int i = idx; i < n; i += str) {
      a1[i] = a1[i] + b1[i] * c;
      a2[i] = a2[i] + b2[i] * c;
    }
  }

}


/**
 * Device kernel for opadd2col
 */
__kernel void opadd2col_kernel(__global real *a1, 
			       __global real *a2, 
			       __global real *a3,
			       __global const real *b1, 
			       __global const real *b2, 
			       __global const real *b3,
			       __global const real *c,
			       const int gdim,
			       const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);
  
  if (gdim == 3) {
    for (int i = idx; i < n; i += str) {
      a1[i] = a1[i] + b1[i] * c[i];
      a2[i] = a2[i] + b2[i] * c[i];
      a3[i] = a3[i] + b3[i] * c[i];
    }
  } 
  else {
    for (int i = idx; i < n; i += str) {
      a1[i] = a1[i] + b1[i] * c[i];
      a2[i] = a2[i] + b2[i] * c[i];
    }
  }

}


#endif // __MATH_MATHOPS_KERNEL_CL__
