#include <hip/hip_runtime.h>

/**
 * Device kernel for add2s1
 */
__global__ void add2s1_kernel(double * __restrict__ a,
			      const double * __restrict__ b,
			      const double c1,
			      const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = c1 * a[i] + b[i];
  }
}

/**
 * Device kernel for add2s2
 */
__global__ void add2s2_kernel(double * __restrict__ a,
			      const double * __restrict__ b,
			      const double c1,
			      const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = a[i] + c1 * b[i];
  }
}

/** 
 * Device kernel for sub3
 */
__global__ void sub3_kernel(double * __restrict__ a,
			    const double * __restrict__ b,
			    const double * __restrict__ c,
			    const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    a[i] = b[i] - c[i];
  }  
}

/**
 * Device kernel for glsc3
 */
__global__ void glsc3_kernel(const double * a,
			     const double * b,
			     const double * c,
			     double * buf_h,
			     const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  __shared__ double buf[1024];
  double tmp = 0.0;

  for (int i = idx; i < n; i+= str) {
    tmp += a[i] * b[i] * c[i];
  }
  buf[threadIdx.x] = tmp;
  __syncthreads();

  int i = blockDim.x>>1;
  while (i != 0) {
    if (threadIdx.x < i) {
      buf[threadIdx.x] += buf[threadIdx.x + i];
    }
    __syncthreads();
    i = i>>1;
  }
 
  if (threadIdx.x == 0) {
    buf_h[blockIdx.x] = buf[0];
  }
}

extern "C" {

  /** Fortran wrapper for copy
   * Copy a vector \f$ a = b \f$
   */
  void hip_copy(void *a, void *b, int *n) {
    hipMemcpy(a, b, (*n) * sizeof(double), hipMemcpyDeviceToDevice);
  }

  /** Fortran wrapper for rzero
   * Zero a real vector
   */
  void hip_rzero(void *a, int *n) {
    hipMemset(a, 0, (*n) * sizeof(double));
  }

  
  /**
   * Fortran wrapper for add2s1
   * Vector addition with scalar multiplication \f$ a = c_1 a + b \f$
   * (multiplication on first argument) 
   */
  void hip_add2s1(void *a, void *b, double *c1, int *n) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(add2s1_kernel,
		       nblcks, nthrds, 0, 0,
		       (double *) a, (double *) b,
		       *c1, *n);
    
  }

  /**
   * Fortran wrapper for add2s2
   * Vector addition with scalar multiplication \f$ a = a + c_1 b \f$
   * (multiplication on second argument) 
   */
  void hip_add2s2(void *a, void *b, double *c1, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(add2s2_kernel,
		       nblcks, nthrds, 0, 0,
		       (double *) a, (double *) b,
		       *c1, *n);

  }


  /**
   * Fortran wrapper for sub3
   * Vector subtraction \f$ a = b - c \f$
   * (multiplication on second argument) 
   */
  void hip_sub3(void *a, void *b, void *c, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    hipLaunchKernelGGL(sub3_kernel,
		       nblcks, nthrds, 0, 0,
		       (double *) a, (double *) b, (double *) c, *n);
  }


  /**
   * Fortran wrapper glsc3
   * Weighted inner product \f$ a^T b c \f$
   */
  double hip_glsc3(void *a, void *b, void *c, int *n) {
	
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    
    double * buf = (double *) malloc(nb * sizeof(double));
    double * buf_d;

    hipMalloc(&buf_d, nb*sizeof(double));
     
    hipLaunchKernelGGL(glsc3_kernel,
		       nblcks, nthrds, 0, 0,
		       (double *) a, (double *) b,
		       (double *) c, buf_d, *n);

    hipMemcpy(buf, buf_d, nb * sizeof(double), hipMemcpyDeviceToHost);

    double res = 0.0;
    for (int i = 0; i < nb; i++) {
      res += buf[i];
    }

    free(buf);
    hipFree(buf_d);

    return res;
  }


 
}
