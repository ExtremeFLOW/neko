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

#ifndef __BC_UTILS_H__
#define __BC_UTILS_H__

/**
 * Computes the linear index for area and normal arrays
 * @note Fortran indexing input, C indexing output
 */
#define coef_normal_area_idx(i, j, k, l, lx, nf) \
  (((i) + (lx) * (((j) - 1) + (lx) * (((k) - 1) + (nf) * (((l) - 1))))) - 1)

/**
 * Device function to compute i,j,k,e indices from a linear index
 * @note Assumes idx is a Fortran index 
 */
void nonlinear_index(const int idx, const int lx, int *index) {
  const int idx2 = idx -1;
  index[3] = idx2/(lx * lx * lx) ;
  index[2] = (idx2 - (lx*lx*lx)*index[3])/(lx * lx);
  index[1] = (idx2 - (lx*lx*lx)*index[3] - (lx*lx) * index[2]) / lx;
  index[0] = (idx2 - (lx*lx*lx)*index[3] - (lx*lx) * index[2]) - lx*index[1];
  index[0]++;
  index[1]++;
  index[2]++;
  index[3]++;
}

#endif 