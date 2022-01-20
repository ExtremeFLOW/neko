/**
 * Device kernel for vector apply for a Dirichlet condition
 */
__kernel void inflow_apply_vector_kernel(__global const int *msk,
					 __global real *x,
					 __global real *y,
					 __global real *z,
					 const real gx,
					 const real gy,
					 const real gz,
					 const int m) {
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < m; i += str) {
    const int k = (msk[i+1] - 1);
    x[k] = gx;
    y[k] = gy;
    z[k] = gz;
  }
}
