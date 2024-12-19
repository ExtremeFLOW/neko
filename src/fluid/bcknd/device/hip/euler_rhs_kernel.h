#ifndef __FLUID_EULER_RES_KERNEL__
#define __FLUID_EULER_RES_KERNEL__

template< typename T >
__global__ void euler_res_part1_kernel(T * __restrict__ rhs_rho,
                                     const T * __restrict__ Binv,
                                     const T * __restrict__ lap_rho,
                                     const T c_avisc,
                                     const int n) {
  
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  
  for (int i = idx; i < n; i += str) {
    rhs_rho[i] =  rhs_rho[i] + c_avisc * Binv[i] * lap_rho[i];
  }
}

#endif // __FLUID_EULER_RES_KERNEL__
