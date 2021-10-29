/**
 * Device kernel for add2s1
 */
__kernel void add2s1_kernel(__global double *a,
			  __global const double *b,
			  const double c1,
			  const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = c1 * a[i] + b[i];
  } 
}
