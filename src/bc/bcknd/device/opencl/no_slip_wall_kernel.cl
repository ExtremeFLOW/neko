/**
 * Device kernel for scalar apply for a no-slip wall conditon
 */
__kernel void no_slip_wall_apply_scalar_kernel(__global const int *msk,
					       __global real *x,
					       const int m) {
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = (idx + 1); i < m; i += str) {
    const int k = (msk[i] - 1);
    x[k] = 0.0;
  }
}

/**
 * Device kernel for vector apply for a no-slip wall conditon
 */
__kernel void no_slip_wall_apply_vector_kernel(__global const int *msk,
					       __global real *x,
					       __global real *y,
					       __global real *z,
					       const int m) {
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = (idx + 1); i < m; i += str) {
    const int k = (msk[i] - 1);
    x[k] = 0.0;
    y[k] = 0.0;
    z[k] = 0.0;
  }
}

