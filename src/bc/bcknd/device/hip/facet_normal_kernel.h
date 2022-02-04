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
 * Computes the linear index for area and normal arrays
 * @note Fortran indexing input, C indexing output
 */

#define coef_normal_area_idx(i, j, k, l, lx, nf) \
  (((i) + (lx) * (((j) - 1) + (lx) * (((k) - 1) + (nf) * (((l) - 1))))) - 1)

/**
 * Device function to compute i,j,k,e indices from a linear index
 * @note Assumes idx is a Fortran index 
 */
__device__
void nonlinear_index(const int idx, const int lx, int *index) {
  index[3] = idx/(lx * lx * lx);
  index[2] = (idx - (lx*lx*lx)*index[3])/(lx * lx);
  index[1] = (idx - (lx*lx*lx)*index[3] - (lx*lx) * index[2]) / lx;
  index[0] = (idx - (lx*lx*lx)*index[3] - (lx*lx) * index[2]) - lx*index[1];
  index[1]++;
  index[2]++;
  index[3]++;
}


/**
 * Device kernel for vector apply for a symmetry condition
 */
template< typename T >
__global__
void facet_normal_apply_surfvec_kernel(const int * __restrict__ msk,
				       const int * __restrict__ facet,
				       T * __restrict__ x,
				       T * __restrict__ y,
				       T * __restrict__ z,
				       const T * __restrict__ u,
				       const T * __restrict__ v,
				       const T * __restrict__ w,
				       const T * __restrict__ nx,
				       const T * __restrict__ ny,
				       const T * __restrict__ nz,
				       const T * __restrict__ area,
				       const int lx,
				       const int m) {
  int index[4];
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = (idx + 1); i < m; i += str) {
    const int k = (msk[i] - 1);
    const int f = (facet[i]);
    nonlinear_index(msk[i], lx, index);


    switch(f) {
    case 1:
    case 2:
      {
	const int na_idx = coef_normal_area_idx(index[1], index[2],
						f, index[3], lx, 6);
	x[k] = u[k] * nx[na_idx] * area[na_idx];
	y[k] = v[k] * ny[na_idx] * area[na_idx];
	z[k] = w[k] * nz[na_idx] * area[na_idx];
	break;
      }
    case 3:
    case 4:
      {
	const int na_idx = coef_normal_area_idx(index[0], index[2],
						f, index[3], lx, 6);
	x[k] = u[k] * nx[na_idx] * area[na_idx];
	y[k] = v[k] * ny[na_idx] * area[na_idx];
	z[k] = w[k] * nz[na_idx] * area[na_idx];
	break;
      }
    case 5:
    case 6:
      {
	const int na_idx = coef_normal_area_idx(index[0], index[1],
						f, index[3], lx, 6);
	x[k] = u[k] * nx[na_idx] * area[na_idx];
	y[k] = v[k] * ny[na_idx] * area[na_idx];
	z[k] = w[k] * nz[na_idx] * area[na_idx];
	break;
      }    
    }
  }
}

