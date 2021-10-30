/**
 * Device kernel for add2s1
 */
__kernel void add2s1_kernel(__global real *a,
			    __global const real *b,
			    const real c1,
			    const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = c1 * a[i] + b[i];
  } 
}

/**
 * Device kernel for add2s2
 */
__kernel void add2s2_kernel(__global real *a,
			    __global const real *b,
			    const real c1,
			    const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] + c1 * b[i];
  }
}

/** 
 * Device kernel for col2
 */
__kernel void col2_kernel(__global real *a,
			  __global const real *b,
			  const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] * b[i];
  }  
}
