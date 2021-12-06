/**
 * Device kernel for vector apply for a Blasius profile
 */
__kernel void blasius_apply_vector_kernel(__global const int *msk,
					  __global real *x,
					  __global real *y,
					  __global real *z,
					  __global real *bla_x,
					  __global real *bla_y,
					  __global real *bla_z,
					  const int m) {
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = (idx + 1); i < m; i += str) {
    const int k = msk[i] -1;
    x[k] = bla_x[i];
    y[k] = bla_y[i];
    z[k] = bla_z[i];
  }
}
