/*
 Copyright (c) 2026, The Neko Authors
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

#ifndef __BC_COUPLED_VECTOR_BC_RESOLVER_KERNEL__
#define __BC_COUPLED_VECTOR_BC_RESOLVER_KERNEL__

__kernel void coupled_vector_bc_resolver_apply_kernel(
    __global const int *mixed_msk,
    __global real *x,
    __global real *y,
    __global real *z,
    __global const int *constraint_n,
    __global const int *constraint_t1,
    __global const int *constraint_t2,
    __global const real *n,
    __global const real *t1,
    __global const real *t2,
    const int m) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < m; i += str) {
    const int k = mixed_msk[i];
    const int off = 3 * i;

    real u1 = x[k];
    real u2 = y[k];
    real u3 = z[k];

    real uloc_n = u1 * n[off] + u2 * n[off + 1] + u3 * n[off + 2];
    real uloc_t1 = u1 * t1[off] + u2 * t1[off + 1] + u3 * t1[off + 2];
    real uloc_t2 = u1 * t2[off] + u2 * t2[off + 1] + u3 * t2[off + 2];

    if (constraint_n[i] != 0)
      uloc_n = 0.0;
    if (constraint_t1[i] != 0)
      uloc_t1 = 0.0;
    if (constraint_t2[i] != 0)
      uloc_t2 = 0.0;

    x[k] = uloc_n * n[off] + uloc_t1 * t1[off] + uloc_t2 * t2[off];
    y[k] = uloc_n * n[off + 1] + uloc_t1 * t1[off + 1] +
           uloc_t2 * t2[off + 1];
    z[k] = uloc_n * n[off + 2] + uloc_t1 * t1[off + 2] +
           uloc_t2 * t2[off + 2];
  }
}

#endif
