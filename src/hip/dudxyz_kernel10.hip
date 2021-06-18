#include <hip/hip_runtime.h>

/**
 * Device kernel for derivative (LX 10)
 */

#define LX 10

__global__ void dudxyz_kernel10(double * __restrict__ du,
				const double * __restrict__ u,
				const double * __restrict__ dr,
				const double * __restrict__ ds,
				const double * __restrict__ dt,
				const double * __restrict__ dx,
				const double * __restrict__ dy,
				const double * __restrict__ dz,
				const double * __restrict__ jacinv,				
				const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

}

