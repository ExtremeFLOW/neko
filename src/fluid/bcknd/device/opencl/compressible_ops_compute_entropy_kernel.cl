/*
 Copyright (c) 2025, The Neko Authors
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

#ifndef __OPENCL_COMPRESSIBLE_OPS_COMPUTE_ENTROPY_KERNEL__
#define __OPENCL_COMPRESSIBLE_OPS_COMPUTE_ENTROPY_KERNEL__

/**
 * Device kernel for compute_entropy
 */
__kernel void compute_entropy_kernel(__global real * S,
                                     __global const real * p,
                                     __global const real * rho,
                                     const real gamma,
                                     const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    if (i < n) {
      // S = 1/(gamma-1) * rho * (log(p) - gamma * log(rho))
      real log_p = log(p[i]);
      real log_rho = log(rho[i]);
      S[i] = (1.0 / (gamma - 1.0)) * rho[i] * (log_p - gamma * log_rho);
    }
  }

}

#endif // __OPENCL_COMPRESSIBLE_OPS_COMPUTE_ENTROPY_KERNEL__ 
