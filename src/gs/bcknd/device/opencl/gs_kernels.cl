/*
 Copyright (c) 2021, The Neko Authors
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
 * Device gather kernel for addition of data
 * \f$ v(dg(i)) = v(dg(i)) + u(gd(i)) \f$
 */
__kernel void gather_kernel_add(__global real * __restrict__ v,
				const int m,
				const int o,
				__global const int * __restrict__ dg,
				__global const real * __restrict__ u,
			        const int n,
				__global const int * __restrict__ gd,
				const int nb,
				__global const int * __restrict__ b,
				__global const int * __restrict__ bo) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < nb; i += str) {
    const int blk_len = b[i];
    const int k = bo[i];
    real tmp = u[gd[k] - 1];
    for (int j = 1; j < blk_len; j++) {
      tmp += u[gd[k + j] - 1];
    }
    v[dg[k] - 1] = tmp;
  }
  
  if (o < 0) {
    for (int i = ((abs(o) - 1) + idx); i < m ; i += str) {
      v[dg[i] - 1] = u[gd[i] - 1];
    }
  }
  else {
    if ((idx%2 == 0)) {
      for (int i = ((o - 1) + idx); i < m ; i += str) {
	real tmp = u[gd[i] - 1] + u[gd[i+1] - 1];
	v[dg[i] - 1] = tmp;
      }
    }
  }
  
}

/**
 * Device gather kernel for multiplication of data
 * \f$ v(dg(i)) = v(dg(i)) \cdot u(gd(i)) \f$
 */
__kernel void gather_kernel_mul(__global real * __restrict__ v,
				const int m,
				const int o,
				__global const int * __restrict__ dg,
				__global const real * __restrict__ u,
				const int n,
				__global const int * __restrict__ gd,
				const int nb,
				__global const int * __restrict__ b,
				__global const int * __restrict__ bo) { 
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < nb; i += str) {
    const int blk_len = b[i];
    const int k = bo[i];
    real tmp = u[gd[k] - 1];
    for (int j = 1; j < blk_len; j++) {
      tmp *= u[gd[k + j] - 1];
    }
    v[dg[k] - 1] = tmp;
  }
  
  if (o < 0) {
    for (int i = ((abs(o) - 1) + idx); i < m ; i += str) {
      v[dg[i] - 1] = u[gd[i] - 1];
    }
  }
  else {
    if ((idx%2 == 0)) {
      for (int i = ((o - 1) + idx); i < m ; i += str) {
	real tmp = u[gd[i] - 1] * u[gd[i+1] - 1];
	v[dg[i] - 1] = tmp;
      }
    }
  }

}

/**
 * Device gather kernel for minimum of data
 * \f$ v(dg(i)) = \min(v(dg(i)), u(gd(i))) \f$
 */
__kernel void gather_kernel_min(__global real * __restrict__ v,
				const int m,
				const int o,
				__global const int * __restrict__ dg,
				__global const real * __restrict__ u,
				const int n,
				__global const int * __restrict__ gd,
				const int nb,
				__global const int * __restrict__ b,
				__global const int * __restrict__ bo) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < nb; i += str) {
    const int blk_len = b[i];
    const int k = bo[i];
    real tmp = u[gd[k] - 1];
    for (int j = 1; j < blk_len; j++) {
      tmp = fmin(u[gd[k + j] - 1], tmp);
    }
    v[dg[k] - 1] = tmp;
  }
  
  if (o < 0) {
    for (int i = ((abs(o) - 1) + idx); i < m ; i += str) {
      v[dg[i] - 1] = u[gd[i] - 1];
    }
  }
  else {
    if ((idx%2 == 0)) {
      for (int i = ((o - 1) + idx); i < m ; i += str) {
	real tmp = fmin(u[gd[i] - 1], u[gd[i+1] - 1]);
	v[dg[i] - 1] = tmp;
      }
    }
  }
  
}

/**
 * Device gather kernel for maximum of data
 * \f$ v(dg(i)) = \max(v(dg(i)), u(gd(i))) \f$
 */
__kernel void gather_kernel_max(__global real * __restrict__ v,
				const int m,
				const int o,
				__global const int * __restrict__ dg,
				__global const real * __restrict__ u,
				const int n,
				__global const int * __restrict__ gd,
				const int nb,
				__global const int *b,
				__global const int *bo) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < nb; i += str) {
    const int blk_len = b[i];
    const int k = bo[i];
    real tmp = u[gd[k] - 1];
    for (int j = 1; j < blk_len; j++) {
      tmp = fmax(u[gd[k + j] - 1], tmp);
    }
    v[dg[k] - 1] = tmp;
  }
  
  if (o < 0) {
    for (int i = ((abs(o) - 1) + idx); i < m ; i += str) {
      v[dg[i] - 1] = u[gd[i] - 1];
    }
  }
  else {
    if ((idx%2 == 0)) {
      for (int i = ((o - 1) + idx); i < m ; i += str) {
	real tmp = fmax(u[gd[i] - 1], u[gd[i+1] - 1]);
	v[dg[i] - 1] = tmp;
      }
    }
  }
  
}

/**
 * Device scatter kernel
 * \f$ u(gd(i) = v(dg(i)) \f$
 */
__kernel void scatter_kernel(__global real * __restrict__ v,
			     const int m,
			     __global const int * __restrict__ dg,
			     __global real *u,
			     const int n,
			     __global const int * __restrict__ gd,
			     const int nb,
			     __global const int * __restrict__ b,
			     __global const int * __restrict__ bo) {
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);
  
  for (int i = idx; i < nb; i += str) {
    const int blk_len = b[i];
    const int k = bo[i];
    real tmp = v[dg[k] - 1];
    for (int j = 0; j < blk_len; j++) {
      u[gd[k + j] - 1] = tmp;
    }      
  }

  const int facet_offset = bo[nb - 1] + b[nb - 1];
  
  for (int i = ((facet_offset - 1) + idx); i < m; i += str) {
    u[gd[i] - 1] = v[dg[i] - 1];
  }
  
}
