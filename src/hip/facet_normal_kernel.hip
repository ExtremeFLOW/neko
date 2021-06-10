#include <hip/hip_runtime.h>

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
__global__
void facet_normal_apply_surfvec_kernel(const int * __restrict__ msk,
				       const int * __restrict__ facet,
				       double * __restrict__ x,
				       double * __restrict__ y,
						  double * __restrict__ z,
				       const double * __restrict__ u,
				       const double * __restrict__ v,
				       const double * __restrict__ w,
				       const double * __restrict__ nx,
				       const double * __restrict__ ny,
				       const double * __restrict__ nz,
				       const double * __restrict__ area,
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

extern "C" {

  /** 
   * Fortran wrapper for device facet normal apply surfvec
   */
  void hip_facet_normal_apply_surfvec(void *msk, void *facet,
				      void *x, void *y, void *z,
				      void *u, void *v, void *w,
				      void *nx, void * ny, void *nz,
				      void *area, int *lx, int *m) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m) + 1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(facet_normal_apply_surfvec_kernel,
		       nblcks, nthrds, 0, 0,
		       (int *) msk, (int *) facet,
		       (double *) x, (double *) y, (double *) z,
		       (double *) u, (double *) v, (double *) v,
		       (double *) nx, (double *) ny, (double *) nz,
		       (double *) area, *lx, *m);
  }
 
}
