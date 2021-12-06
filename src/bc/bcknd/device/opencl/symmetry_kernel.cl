
/**
 * Device kernel for vector apply for a symmetry condition
 */
__kernel void symmetry_apply_vector_kernel(__global const int *xmsk,
					   __global const int *ymsk,
					   __global const int *zmsk,
					   __global real *x,
					   __global real *y,
					   __global real *z,
					   const int m,
					   const int n,
					   const int l) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = (idx + 1); i < m; i += str) {
    const int k = (xmsk[i] - 1);
    x[k] = 0.0;
  }

  for (int i = (idx + 1); i < n; i += str) {
    const int k = (ymsk[i] - 1);
    y[k] = 0.0;
  }

  for (int i = (idx + 1); i < l; i += str) {
    const int k = (zmsk[i] - 1);
    z[k] = 0.0;
  }
  
}
