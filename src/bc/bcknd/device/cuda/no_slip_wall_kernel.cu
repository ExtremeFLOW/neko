/**
 * Device kernel for scalar apply for a no-slip wall conditon
 */
__global__ void no_slip_wall_apply_scalar_kernel(const int * __restrict__ msk,
						 double * __restrict__ x,
						 const int m) {
  
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = (idx + 1); i < m; i += str) {
    const int k = (msk[i] - 1);
    x[k] = 0.0;
  }
}

/**
 * Device kernel for vector apply for a no-slip wall conditon
 */
__global__ void no_slip_wall_apply_vector_kernel(const int * __restrict__ msk,
						 double * __restrict__ x,
						 double * __restrict__ y,
						 double * __restrict__ z,
						 const int m) {
  
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = (idx + 1); i < m; i += str) {
    const int k = (msk[i] - 1);
    x[k] = 0.0;
    y[k] = 0.0;
    z[k] = 0.0;
  }
}

extern "C" {

  /** 
   * Fortran wrapper for device no-slop wall apply vector
   */
  void cuda_no_slip_wall_apply_scalar(void *msk, void *x, int *m) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    no_slip_wall_apply_scalar_kernel<<<nblcks, nthrds>>>((int *) msk,
							 (double *) x, *m);
  }
  
  /** 
   * Fortran wrapper for device no-slop wall apply vector
   */
  void cuda_no_slip_wall_apply_vector(void *msk, void *x, void *y,
				     void *z, int *m) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m)+1024 - 1)/ 1024, 1, 1);

    no_slip_wall_apply_vector_kernel<<<nblcks, nthrds>>>((int *) msk,
							 (double *) x,
							 (double *) y, 
							 (double *) z, *m); 
  }
 
}
