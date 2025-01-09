#ifndef __FLUID_EULER_RES_KERNEL__
#define __FLUID_EULER_RES_KERNEL__

template< typename T >
__global__ void euler_res_part_visc_kernel(T * __restrict__ rhs_rho,
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

template< typename T >
__global__ void euler_res_part_mx_flux_kernel(T * __restrict__ f_x,
                                     T * __restrict__ f_y,
                                     T * __restrict__ f_z,
                                     const T * __restrict__ m_x,
                                     const T * __restrict__ m_y,
                                     const T * __restrict__ m_z,
                                     const T * __restrict__ rho_field,
                                     const T * __restrict__ p,
                                     const int n) {
  
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  
  for (int i = idx; i < n; i += str) {
    f_x[i] = m_x[i] * m_x[i] / rho_field[i] + p[i];
    f_y[i] = m_x[i] * m_y[i] / rho_field[i];
    f_z[i] = m_x[i] * m_z[i] / rho_field[i];
  }
}

template< typename T >
__global__ void euler_res_part_my_flux_kernel(T * __restrict__ f_x,
                                     T * __restrict__ f_y,
                                     T * __restrict__ f_z,
                                     const T * __restrict__ m_x,
                                     const T * __restrict__ m_y,
                                     const T * __restrict__ m_z,
                                     const T * __restrict__ rho_field,
                                     const T * __restrict__ p,
                                     const int n) {
  
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  
  for (int i = idx; i < n; i += str) {
    f_x[i] = m_y[i] * m_x[i] / rho_field[i];
    f_y[i] = m_y[i] * m_y[i] / rho_field[i] + p[i];
    f_z[i] = m_y[i] * m_z[i] / rho_field[i];
  }
}

template< typename T >
__global__ void euler_res_part_mz_flux_kernel(T * __restrict__ f_x,
                                     T * __restrict__ f_y,
                                     T * __restrict__ f_z,
                                     const T * __restrict__ m_x,
                                     const T * __restrict__ m_y,
                                     const T * __restrict__ m_z,
                                     const T * __restrict__ rho_field,
                                     const T * __restrict__ p,
                                     const int n) {
  
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  
  for (int i = idx; i < n; i += str) {
    f_x[i] = m_z[i] * m_x[i] / rho_field[i];
    f_y[i] = m_z[i] * m_y[i] / rho_field[i];
    f_z[i] = m_z[i] * m_z[i] / rho_field[i] + p[i];
  }
}

template< typename T >
__global__ void euler_res_part_E_flux_kernel(T * __restrict__ f_x,
                                     T * __restrict__ f_y,
                                     T * __restrict__ f_z,
                                     const T * __restrict__ m_x,
                                     const T * __restrict__ m_y,
                                     const T * __restrict__ m_z,
                                     const T * __restrict__ rho_field,
                                     const T * __restrict__ p,
                                     const T * __restrict__ E,
                                     const int n) {
  
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  
  for (int i = idx; i < n; i += str) {
    f_x[i] = (E[i] + p[i]) * m_x[i] / rho_field[i];
    f_y[i] = (E[i] + p[i]) * m_y[i] / rho_field[i];
    f_z[i] = (E[i] + p[i]) * m_z[i] / rho_field[i];
  }
}

template< typename T >
__global__ void euler_res_part_coef_mult_kernel(T * __restrict__ rhs_rho,
                                     T * __restrict__ rhs_m_x,
                                     T * __restrict__ rhs_m_y,
                                     T * __restrict__ rhs_m_z,
                                     T * __restrict__ rhs_E,
                                     const T * __restrict__ mult,
                                     const int n) {
  
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  
  for (int i = idx; i < n; i += str) {
    rhs_rho[i] = rhs_rho[i] * mult[i];
    rhs_m_x[i] = rhs_m_x[i] * mult[i];
    rhs_m_y[i] = rhs_m_y[i] * mult[i];
    rhs_m_z[i] = rhs_m_z[i] * mult[i];
    rhs_E[i] = rhs_E[i] * mult[i];
  }
}

#endif // __FLUID_EULER_RES_KERNEL__
