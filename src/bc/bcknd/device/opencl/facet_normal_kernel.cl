
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
__kernel
void facet_normal_apply_surfvec_kernel(__global const int *msk,
				       __global const int *facet,
				       __global real *x,
				       __global real *y,
				       __global real *z,
				       __global const real *u,
				       __global const real *v,
				       __global const real *w,
				       __global const real *nx,
				       __global const real *ny,
				       __global const real *nz,
				       __global const real *area,
				       const int lx,
				       const int m) {
  int index[4];
  const int idx = get_global_id(0);
  const int str = get_global_size(0);

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

