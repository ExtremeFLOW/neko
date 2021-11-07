/**
 * Device kernel for add2s1
 */
__kernel void add2s1_kernel(__global real * __restrict__ a,
			    __global const real * __restrict__ b,
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
__kernel void add2s2_kernel(__global real * __restrict__ a,
			    __global const real * __restrict__ b,
			    const real c1,
			    const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] + c1 * b[i];
  }
}

/** 
 * Device kernel for invcol2
 */
__kernel void invcol2_kernel(__global real * __restrict__ a,
			     __global const real * __restrict__ b,
			     const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] / b[i];
  }  
}

/** 
 * Device kernel for col2
 */
__kernel void col2_kernel(__global real * __restrict__ a,
			  __global const real * __restrict__ b,
			  const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] * b[i];
  }  
}

/** 
 * Device kernel for col3
 */
__kernel void col3_kernel(__global real * __restrict__ a,
			  __global const real * __restrict__ b,
			  __global const real * __restrict__ c,
			  const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = b[i] * c[i];
  }  
}

/** 
 * Device kernel for sub3
 */
__kernel void sub3_kernel(__global real * __restrict__ a,
			  __global const real * __restrict__ b,
			  __global const real * __restrict__ c,
			  const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = b[i] - c[i];
  }  
}

/** 
 * Device kernel for addcol3
 */
__kernel void addcol3_kernel(__global real * __restrict__ a,
			     __global const real * __restrict__ b,
			     __global const real * __restrict__ c,
			     const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] + b[i] * c[i];
  }  
}

/**
 * Device kernel for glsc3
 */
__kernel void glsc3_kernel(__global const real * __restrict__ a,
			   __global const real * __restrict__ b,
			   __global const real * __restrict__ c,
			   __global real * __restrict__ buf_h,
			   const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  __local real buf[256]; /* Make this nice...*/
  real tmp = 0.0;

  for (int i = idx; i < n; i+= str) {
    tmp += a[i] * b[i] * c[i];
  }
  buf[get_local_id(0)] = tmp;
  barrier(CLK_LOCAL_MEM_FENCE);

  int i = (get_local_size(0))>>1;
  while (i != 0) {
    if (get_local_id(0) < i) {
      buf[get_local_id(0)] += buf[get_local_id(0) + i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    i = i>>1;
  }
 
  if (get_local_id(0) == 0) {
    buf_h[get_group_id(0)] = buf[0];
  }

}

/**
 * Device kernel for glsc2
 */
__kernel void glsc2_kernel(__global const real * __restrict__ a,
			   __global const real * __restrict__ b,
			   __global real * __restrict__ buf_h,
			   const int n) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  __local real buf[256]; /* Make this nice...*/
  real tmp = 0.0;

  for (int i = idx; i < n; i+= str) {
    tmp += a[i] * b[i];
  }
  buf[get_local_id(0)] = tmp;
  barrier(CLK_LOCAL_MEM_FENCE);

  int i = (get_local_size(0))>>1;
  while (i != 0) {
    if (get_local_id(0) < i) {
      buf[get_local_id(0)] += buf[get_local_id(0) + i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    i = i>>1;
  }
 
  if (get_local_id(0) == 0) {
    buf_h[get_group_id(0)] = buf[0];
  }

}
