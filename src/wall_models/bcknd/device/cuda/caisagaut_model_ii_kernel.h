#ifndef __COMMON_CAISAGAUT_MODEL_II_KERNEL_H__
#define __COMMON_CAISAGAUT_MODEL_II_KERNEL_H__

#include <cfloat>
#include <algorithm>
#include <cmath>

/**
 * CUDA kernel for the Caisagaut Model-II wall model.
 * @param u_d The sampled x-velocity field.
 * @param v_d The sampled y-velocity field.
 * @param w_d The sampled z-velocity field.
 * @param ind_r_d The r-index array for sampled GLL points.
 * @param ind_s_d The s-index array for sampled GLL points.
 * @param ind_t_d The t-index array for sampled GLL points.
 * @param ind_e_d The element-index array for sampled GLL points.
 * @param n_x_d The x-component of the wall normals.
 * @param n_y_d The y-component of the wall normals.
 * @param n_z_d The z-component of the wall normals.
 * @param nu_d The sampled kinematic viscosity at wall points.
 * @param rho_w_d The sampled density at wall points.
 * @param h_d The wall-model sampling distances.
 * @param tau_x_d The x-component of the wall shear stress.
 * @param tau_y_d The y-component of the wall shear stress.
 * @param tau_z_d The z-component of the wall shear stress.
 * @param n_nodes The number of wall points.
 * @param lx The number of GLL points per direction.
 * @param kappa The von Karman coefficient.
 * @param B The log-law intercept.
 * @param p The blending exponent.
 * @param s The blending scale.
 */
template<typename T>
__global__ void caisagaut_model_ii_compute(const T * __restrict__ u_d,
                                           const T * __restrict__ v_d,
                                           const T * __restrict__ w_d,
                                           const int * __restrict__ ind_r_d,
                                           const int * __restrict__ ind_s_d,
                                           const int * __restrict__ ind_t_d,
                                           const int * __restrict__ ind_e_d,
                                           const T * __restrict__ n_x_d,
                                           const T * __restrict__ n_y_d,
                                           const T * __restrict__ n_z_d,
                                           const T * __restrict__ nu_d,
                                           const T * __restrict__ rho_w_d,
                                           const T * __restrict__ h_d,
                                           T * __restrict__ tau_x_d,
                                           T * __restrict__ tau_y_d,
                                           T * __restrict__ tau_z_d,
                                           const int n_nodes,
                                           const int lx,
                                           const T kappa,
                                           const T B,
                                           const T p,
                                           const T s) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  const T one = static_cast<T>(1.0);
  const T eps = (sizeof(T) == sizeof(float)) ?
    static_cast<T>(FLT_EPSILON) :
    static_cast<T>(DBL_EPSILON);

  for (int i = idx; i < n_nodes; i += str) {
    const int index = (ind_e_d[i] - 1) * lx * lx * lx +
      (ind_t_d[i] - 1) * lx * lx +
      (ind_s_d[i] - 1) * lx +
      (ind_r_d[i] - 1);

    T ui = u_d[index];
    T vi = v_d[index];
    T wi = w_d[index];
    const T rho = rho_w_d[i];
    const T nx = n_x_d[i];
    const T ny = n_y_d[i];
    const T nz = n_z_d[i];

    const T normu = ui * nx + vi * ny + wi * nz;
    ui -= normu * nx;
    vi -= normu * ny;
    wi -= normu * nz;

    const T magu = sqrt(ui * ui + vi * vi + wi * wi);

    if (magu < eps) {
      tau_x_d[i] = static_cast<T>(0.0);
      tau_y_d[i] = static_cast<T>(0.0);
      tau_z_d[i] = static_cast<T>(0.0);
      continue;
    }

    const T e_const = exp(kappa * B);
    const T rey = magu * h_d[i] / nu_d[i];
    const T blend = exp(-pow(rey / s, p));
    const T warg = kappa * e_const * rey;
    const T a = one /
      (one + static_cast<T>(0.5) * log(one + warg));
    T wlam = log(one + a * warg);
    wlam = wlam / (one + wlam) * (one + log(warg / wlam));
    const T up = blend * sqrt(rey) +
      (one - blend) * wlam / kappa;
    const T utau = magu / (up + eps);

    tau_x_d[i] = -rho * utau * utau * ui / (magu + eps);
    tau_y_d[i] = -rho * utau * utau * vi / (magu + eps);
    tau_z_d[i] = -rho * utau * utau * wi / (magu + eps);
  }
}

#endif
