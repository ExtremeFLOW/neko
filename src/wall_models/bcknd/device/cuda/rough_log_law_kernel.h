#ifndef ROUGH_LOG_LAW_KERNEL_H
#define ROUGH_LOG_LAW_KERNEL_H

#include <cmath>
#include <algorithm>

/**
 * CUDA kernel for the rough log-law wall model.
 */
template<typename T>
__global__ void rough_log_law_compute(const T* __restrict__ u_d,
                                      const T* __restrict__ v_d,
                                      const T* __restrict__ w_d,
                                      const int* __restrict__ ind_r_d,
                                      const int* __restrict__ ind_s_d,
                                      const int* __restrict__ ind_t_d,
                                      const int* __restrict__ ind_e_d,
                                      const T* __restrict__ n_x_d,
                                      const T* __restrict__ n_y_d,
                                      const T* __restrict__ n_z_d,
                                      const T* __restrict__ h_d,
                                      T* __restrict__ tau_x_d,
                                      T* __restrict__ tau_y_d,
                                      T* __restrict__ tau_z_d,
                                      const int n_nodes,
                                      const int lx,
                                      const T kappa,
                                      const T B,
                                      const T z0) {

    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int str = blockDim.x * gridDim.x;
    for (int i = idx; i < n_nodes; i += str) {
        // Sample the velocity
        const int index = (ind_e_d[i] - 1) * lx * lx * lx +
                          (ind_t_d[i] - 1) * lx * lx +
                          (ind_s_d[i] - 1) * lx +
                          (ind_r_d[i] - 1);

        T ui = u_d[index];
        T vi = v_d[index];
        T wi = w_d[index];

        // Load normal vectors and wall shear stress values once
        T nx = n_x_d[i];
        T ny = n_y_d[i];
        T nz = n_z_d[i];
        T h = h_d[i];

        // Project on tangential direction
        T normu = ui * nx + vi * ny + wi * nz;

        ui -= normu * nx;
        vi -= normu * ny;
        wi -= normu * nz;

        T magu = sqrt(ui * ui + vi * vi + wi * wi);

        // Compute the wall shear stress using the rough log-law
        T utau = 0.0;
        if (h > z0) {
            utau = (magu - B) * kappa / log(h / z0);
        }

        // Distribute according to the velocity vector
        tau_x_d[i] = -utau * utau * ui / magu;
        tau_y_d[i] = -utau * utau * vi / magu;
        tau_z_d[i] = -utau * utau * wi / magu;
    }
}

#endif // ROUGH_LOG_LAW_KERNEL_H