/**
 * Device kernel for scalar apply for a Dirichlet condition
 */
__kernel void dirichlet_apply_scalar_kernel(__global const int *msk,
					    __global real *x,
					    const real g,
					    const int m) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);
  
  for (int i = (idx + 1); i < m; i += str) {
    const int k = msk[i] -1;
    x[k] = g;
  }
}

/**
 * Device kernel for vector apply for a Dirichlet condition
 */
__kernel void dirichlet_apply_vector_kernel(__global const int *msk,
					    __global real *x,
					    __global real *y,
					    __global real *z,
					    const real g,
					    const int m) {
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = (idx + 1); i < m; i += str) {
    const int k = msk[i] -1;
    x[k] = g;
    y[k] = g;
    z[k] = g;
  }
}
