
/**
 * Device gather kernel for addition of data
 * \f$ v(dg(i)) = v(dg(i)) + u(gd(i)) \f$
 */
__kernel void gather_kernel_add(__global real *v,
				const int m,
				const int o,
				__global const int *dg,
				__global const real *u,
			        const int n,
				__global const int *gd,
				const int nb,
				__global const int *b,
				__global const int *bo) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < nb; i += str) {
    const int blk_len = b[i];
    const int k = bo[i];
    real tmp = u[gd[k] - 1];
    for (int j = 1; j < blk_len; j++) {
      tmp += u[gd[k + j] - 1];
    }
    v[dg[k] - 1] = tmp;
  }
  
  if (o < 0) {
    for (int i = ((abs(o) - 1) + idx); i < m ; i += str) {
      v[dg[i] - 1] = u[gd[i] - 1];
    }
  }
  else {
    if ((idx%2 == 0)) {
      for (int i = ((o - 1) + idx); i < m ; i += str) {
	real tmp = u[gd[i] - 1] + u[gd[i+1] - 1];
	v[dg[i] - 1] = tmp;
      }
    }
  }
  
}

/**
 * Device gather kernel for multiplication of data
 * \f$ v(dg(i)) = v(dg(i)) \cdot u(gd(i)) \f$
 */
__kernel void gather_kernel_mul(__global real *v,
				const int m,
				const int o,
				__global const int *dg,
				__global const real *u,
				const int n,
				__global const int *gd,
				const int nb,
				__global const int *b,
				__global const int *bo) { 
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < nb; i += str) {
    const int blk_len = b[i];
    const int k = bo[i];
    real tmp = u[gd[k] - 1];
    for (int j = 1; j < blk_len; j++) {
      tmp *= u[gd[k + j] - 1];
    }
    v[dg[k] - 1] = tmp;
  }
  
  if (o < 0) {
    for (int i = ((abs(o) - 1) + idx); i < m ; i += str) {
      v[dg[i] - 1] = u[gd[i] - 1];
    }
  }
  else {
    if ((idx%2 == 0)) {
      for (int i = ((o - 1) + idx); i < m ; i += str) {
	real tmp = u[gd[i] - 1] * u[gd[i+1] - 1];
	v[dg[i] - 1] = tmp;
      }
    }
  }

}

/**
 * Device gather kernel for minimum of data
 * \f$ v(dg(i)) = \min(v(dg(i)), u(gd(i))) \f$
 */
__kernel void gather_kernel_min(__global real *v,
				const int m,
				const int o,
				__global const int *dg,
				__global const real *u,
				const int n,
				__global const int *gd,
				const int nb,
				__global const int *b,
				__global const int *bo) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < nb; i += str) {
    const int blk_len = b[i];
    const int k = bo[i];
    real tmp = u[gd[k] - 1];
    for (int j = 1; j < blk_len; j++) {
      tmp = min(u[gd[k + j] - 1], tmp);
    }
    v[dg[k] - 1] = tmp;
  }
  
  if (o < 0) {
    for (int i = ((abs(o) - 1) + idx); i < m ; i += str) {
      v[dg[i] - 1] = u[gd[i] - 1];
    }
  }
  else {
    if ((idx%2 == 0)) {
      for (int i = ((o - 1) + idx); i < m ; i += str) {
	real tmp = min(u[gd[i] - 1], u[gd[i+1] - 1]);
	v[dg[i] - 1] = tmp;
      }
    }
  }
  
}

/**
 * Device gather kernel for maximum of data
 * \f$ v(dg(i)) = \max(v(dg(i)), u(gd(i))) \f$
 */
__kernel void gather_kernel_max(__global real *v,
				const int m,
				const int o,
				__global const int *dg,
				__global const real *u,
				const int n,
				__global const int *gd,
				const int nb,
				__global const int *b,
				__global const int *bo) {

  const int idx = get_global_id(0);
  const int str = get_global_size(0);

  for (int i = idx; i < nb; i += str) {
    const int blk_len = b[i];
    const int k = bo[i];
    real tmp = u[gd[k] - 1];
    for (int j = 1; j < blk_len; j++) {
      tmp = max(u[gd[k + j] - 1], tmp);
    }
    v[dg[k] - 1] = tmp;
  }
  
  if (o < 0) {
    for (int i = ((abs(o) - 1) + idx); i < m ; i += str) {
      v[dg[i] - 1] = u[gd[i] - 1];
    }
  }
  else {
    if ((idx%2 == 0)) {
      for (int i = ((o - 1) + idx); i < m ; i += str) {
	real tmp = max(u[gd[i] - 1], u[gd[i+1] - 1]);
	v[dg[i] - 1] = tmp;
      }
    }
  }
  
}

/**
 * Device scatter kernel
 * \f$ u(gd(i) = v(dg(i)) \f$
 */
__kernel void scatter_kernel(__global real *v,
			     const int m,
			     __global const int *dg,
			     __global real *u,
			     const int n,
			     __global const int *gd,
			     const int nb,
			     __global const int *b,
			     __global const int *bo) {
  
  const int idx = get_global_id(0);
  const int str = get_global_size(0);
  
  for (int i = idx; i < nb; i += str) {
    const int blk_len = b[i];
    const int k = bo[i];
    real tmp = v[dg[k] - 1];
    for (int j = 0; j < blk_len; j++) {
      u[gd[k + j] - 1] = tmp;
    }      
  }

  const int facet_offset = bo[nb - 1] + b[nb - 1];
  
  for (int i = ((facet_offset - 1) + idx); i < m; i += str) {
    u[gd[i] - 1] = v[dg[i] - 1];
  }
  
}
