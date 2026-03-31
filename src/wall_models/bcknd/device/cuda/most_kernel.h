/*
 Copyright (c) 2025, The Neko Authors
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

   * Redistributions in binary form must reproduce the above
     copyright notice, this list of conditions and the following
     disclaimer in the documentation and/or other materials provided
     with the distribution.

   * Neither the name of the authors nor the names of its
     contributors may be used to endorse or promote products derived
     from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
     LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef MOST_KERNEL_H
#define MOST_KERNEL_H

#include <cmath>
#include <algorithm>

/*
* Similarity laws and corrections for the STABLE regime:
* REFERENCE: Cheng, Y., and W. Brutsaert (2005), Flux-profile relationships
* for wind speed and temperature in the stable atmospheric boundary layer,  
* Bound.-Layer Meteorol., 3, 519-538.
* NOTE: This formulation is chosen for its superior behavior in very stable
* conditions (large z/L), avoiding the numerical decoupling found in older linear (e.g., Webb or Holtslag) functions.
*/

template<typename T>
__device__ T slaw_m_stable(T z, T L_ob, T z0)
{
    return log(z/z0)
           - corr_m_stable<T>(z,L_ob)
           + corr_m_stable<T>(z0,L_ob);
}

template<typename T>
__device__ T slaw_h_stable(T z, T L_ob, T z0h)
{
    return log(z/z0h)
           - corr_h_stable<T>(z,L_ob)
           + corr_h_stable<T>(z0h,L_ob);
}

template<typename T>
__device__ T corr_m_stable(T z, T L_ob)
{
    // Coefficients specific to Cheng & Brutsaert (2005)
    T a = 1.0;
    T b = 2.0/3.0;
    T c = 5.0;
    T d = 0.35;
    T zeta = z / L_ob;

    return -a*zeta - b*(zeta-c/d)*exp(-d*zeta) - b*c/d;
}

template<typename T>
__device__ T corr_h_stable(T z, T L_ob)
{
    // Coefficients specific to Cheng & Brutsaert (2005)
    T a = 1.0;
    T b = 2.0/3.0;
    T c = 5.0;
    T d = 0.35;
    T zeta = z / L_ob;

    return -b*(zeta-c/d)*exp(-d*zeta)
           -pow(1.0 + 2.0/3.0*a*zeta, 1.5)
           -b*c/d + 1.0;
}

template<typename T>
__device__ T f_neumann_stable(T Ri_b, T z, T z0, T z0h, T L_ob)
{
    return Ri_b - z/L_ob / pow(slaw_m_stable<T>(z,L_ob,z0),3);   
}

template<typename T>
__device__ T dfdl_neumann_stable(T l_upper,
                          T l_lower,
                          T z,
                          T z0,
                          T z0h,
                          T fd_h)
{
    T up = -z/l_upper /
           pow(slaw_m_stable<T>(z,l_upper,z0),3);    

    T low =  z/l_lower /
             pow(slaw_m_stable<T>(z,l_lower,z0),3);     

    return (up + low) / (2*fd_h);
}

template<typename T>
__device__ T f_dirichlet_stable(T Ri_b, T z, T z0, T z0h, T L_ob)
{
    return Ri_b - z/L_ob * slaw_h_stable<T>(z,L_ob,z0h) / pow(slaw_m_stable<T>(z,L_ob,z0),2);      
}

template<typename T>
__device__ T dfdl_dirichlet_stable(T l_upper,
                          T l_lower,
                          T z,
                          T z0,
                          T z0h,
                          T fd_h)
{
    T up = -z/l_upper *
           slaw_h_stable<T>(z,l_upper,z0h) / pow(slaw_m_stable<T>(z,l_upper,z0),2);      

    T low =  z/l_lower *
             slaw_h_stable<T>(z,l_lower,z0h) / pow(slaw_m_stable<T>(z,l_lower,z0),2);      

    return (up + low) / (2*fd_h);
}

/*
* Similarity laws and corrections for the UNSTABLE (convective) regime:
* REFERENCE: Dyer, A. J. (1974), A review of flux-profile relationships,
* Bound.-Layer Meteorol., 7, 363-372.
* INTEGRATION: Paulson, C. A. (1970), The mathematical representation
* of wind speed and temperature profiles in the unstable atmospheric 
* surface layer, J. Appl. Meteorol., 9, 857-861.
*/

template<typename T>
__device__ T slaw_m_convective(T z, T L_ob, T z0)
{
    return log(z/z0)
           - corr_m_convective<T>(z,L_ob)
           + corr_m_convective<T>(z0,L_ob);
}

template<typename T>
__device__ T slaw_h_convective(T z, T L_ob, T z0h)
{
    return log(z/z0h)
           - corr_h_convective<T>(z,L_ob)
           + corr_h_convective<T>(z0h,L_ob);
}

template<typename T>
__device__ T corr_m_convective(T z, T L_ob)
{
    T zeta = z / L_ob;
    T pi = 4*atan(1.0);
    // Standard Dyer-Businger coefficient gamma = 16.0
    T xi = sqrt(sqrt((1.0 - 16.0*zeta)));

    return 2*log(0.5*(1 + xi)) + log(0.5*(1 + xi*xi)) - 2*atan(xi) + pi/2;
}

template<typename T>
__device__ T corr_h_convective(T z, T L_ob)
{
    T zeta = z/L_ob;
    T pi = 4*atan(1.0);
    // Standard Dyer-Businger coefficient gamma = 16.0
    T xi = sqrt(sqrt((1.0 - 16.0*zeta)));
    return 2*log(0.5*(1 + xi*xi));
}

template<typename T>
__device__ T f_neumann_convective(T Ri_b, T z, T z0, T z0h, T L_ob)
{
    return Ri_b - z/L_ob / pow(slaw_m_convective<T>(z,L_ob,z0),3);   
}

template<typename T>
__device__ T dfdl_neumann_convective(T l_upper,
                          T l_lower,
                          T z,
                          T z0,
                          T z0h,
                          T fd_h)
{
    T up = -z/l_upper /
           pow(slaw_m_convective<T>(z,l_upper,z0),3);    

    T low =  z/l_lower /
             pow(slaw_m_convective<T>(z,l_lower,z0),3);     

    return (up + low) / (2*fd_h);
}

template<typename T>
__device__ T f_dirichlet_convective(T Ri_b, T z, T z0, T z0h, T L_ob)
{
    return Ri_b - z/L_ob * slaw_h_convective<T>(z,L_ob,z0h) / pow(slaw_m_convective<T>(z,L_ob,z0),2);      
}

template<typename T>
__device__ T dfdl_dirichlet_convective(T l_upper,
                          T l_lower,
                          T z,
                          T z0,
                          T z0h,
                          T fd_h)
{
    T up = -z/l_upper *
           slaw_h_convective<T>(z,l_upper,z0h) / pow(slaw_m_convective<T>(z,l_upper,z0),2);      

    T low =  z/l_lower *
             slaw_h_convective<T>(z,l_lower,z0h) / pow(slaw_m_convective<T>(z,l_lower,z0),2);      

    return (up + low) / (2*fd_h);
}

/*
* Similarity laws and corrections for the NEUTRAL regime:
*/

template<typename T>
__device__ T slaw_m_neutral(T z, T z0)
{
    return log(z/z0);
}

template<typename T>
__device__ T slaw_h_neutral(T z, T z0h)
{
    return log(z/z0h);
}

/*
 * CUDA kernel for the most wall model.   
 */
template<typename T, int BC_TYPE>  
__global__ void most_compute(
    const T* __restrict__ u_d,     
    const T* __restrict__ v_d,
    const T* __restrict__ w_d,
    const T* __restrict__ temp_d,
    const T* __restrict__ h_d,     
    const T* __restrict__ n_x_d,  
    const T* __restrict__ n_y_d,
    const T* __restrict__ n_z_d,
    const int* __restrict__ ind_r_d,
    const int* __restrict__ ind_s_d,
    const int* __restrict__ ind_t_d,
    const int* __restrict__ ind_e_d,
    T* __restrict__ tau_x_d,       
    T* __restrict__ tau_y_d,
    T* __restrict__ tau_z_d,
    int n_nodes,
    int lx,                     
    T kappa,
    T mu,
    T rho, 
    T z0,
    T z0h_in,
    T bc_value
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int str = blockDim.x * gridDim.x;
    if(idx >= n_nodes) return;

    const T g = 9.80665;
    const T Ri_threshold = 1e-4;
    const T tol = 1e-3;
    const T NR_step = 1e-3;
    const int max_iter = 50;

    // Calculate the global 1D offset in the field arrays
    // Layout: e is the slowest, r is the fastest
    // Mapping GLL points from the full field to the boundary thread
    // Neko provides 1-based indices from Fortran, so subtract 1.
    for (int i = idx; i < n_nodes; i += str) {
        const int index = (ind_e_d[i] - 1) * lx * lx * lx +
                          (ind_t_d[i] - 1) * lx * lx +
                          (ind_s_d[i] - 1) * lx +
                          (ind_r_d[i] - 1);  // 'long long' might be needed to avoid
                                             // overflowing 32-bits of 'int'
        T ui = u_d[index];
        T vi = v_d[index];
        T wi = w_d[index];
        T ti = temp_d[index];
        T hi = h_d[i];

        // Extract the local normal vector
        T nx = n_x_d[i];
        T ny = n_y_d[i];
        T nz = n_z_d[i];
        
        // Get the tangential component
        T normu = ui * nx + vi * ny + wi * nz;
        ui -= normu * nx;
        vi -= normu * ny;
        wi -= normu * nz;

        // The magnitude used for MOST is the magnitude of the tangential vector
        T magu = sqrt(ui*ui + vi*vi + wi*wi);
        magu = fmax(magu, (T)1e-6); // Avoid division by zero in extreme cases

        // Initial utau estimate
        T utau = kappa * magu / log(hi/z0);

        // Zilitinkevich 1995 correlation for thermal roughness
        T z0h;
        if (z0h_in < 0.0) {
            const T nu = 1.46e-5; // Kinematic viscosity of air [m^2/s]
            const T Re_z0 = (utau * z0) / nu;
            z0h = z0 * exp(-0.1 * sqrt(Re_z0));
        } else {
            z0h = z0h_in;
        }

        T ts = 0;
        T q  = 0;

        if constexpr (BC_TYPE == 0)      // Neumann
            q = bc_value;
        else                             // Dirichlet
        {
            ts = bc_value;
            q  = kappa*utau*(ts-ti)/log(hi/z0h);
        }

        T Ri_b;

        if constexpr (BC_TYPE == 0)
            Ri_b = -g*hi/ti*q/(magu*magu*magu*kappa*kappa);
        else
            Ri_b =  g*hi/ti*(ti-ts)/(magu*magu);

        T L_ob = 1e10;   // neutral default

        const T L_sign = (Ri_b > 0) ? 1.0 : -1.0;
    
        // Stability branching localy 
        if (fabs(Ri_b) <= Ri_threshold) {
            // NEUTRAL CASe
            utau = kappa * magu / slaw_m_neutral<T>(hi, z0);
            if constexpr (BC_TYPE == 1) q = kappa * utau * (ts - ti) / slaw_h_neutral<T>(hi, z0h);
        } 
        else {
            // STABLE or CONVECTIVE (NR)
            // Initial guess based on stability
            if (Ri_b > 0)
                L_ob = hi / fmax(Ri_b, Ri_threshold);
            else
                L_ob = hi / fmin(Ri_b, -Ri_threshold);
            
            T L_old;
            for (int it = 0; it < max_iter; ++it) {
                L_old = L_ob;
                T f_val, dfdl;
                T fd_h   = NR_step * L_ob;
                T L_upper = L_ob + fd_h;
                T L_lower = L_ob - fd_h;

                // Use the appropriate simlarity law based on stability and b.c. type
                if (Ri_b > 0) { // Stable 
                    if constexpr (BC_TYPE == 0) {
                        f_val = f_neumann_stable<T>(Ri_b, hi, z0, z0h, L_ob);
                        dfdl = dfdl_neumann_stable<T>(L_upper, L_lower, hi, z0, z0h, fd_h);
                    } else {
                        f_val = f_dirichlet_stable<T>(Ri_b, hi, z0, z0h, L_ob);
                        dfdl = dfdl_dirichlet_stable<T>(L_upper, L_lower, hi, z0, z0h, fd_h);
                    }
                } else { // Convective 
                    if constexpr (BC_TYPE == 0) {
                        f_val = f_neumann_convective<T>(Ri_b, hi, z0, z0h, L_ob);
                        dfdl = dfdl_neumann_convective<T>(L_upper, L_lower, hi, z0, z0h, fd_h);  
                    } else {
                        f_val = f_dirichlet_convective<T>(Ri_b, hi, z0, z0h, L_ob);
                        dfdl = dfdl_dirichlet_convective<T>(L_upper, L_lower, hi, z0, z0h, fd_h);
                    }
                }

                if (fabs(dfdl) < (T)1e-12) break; 
                L_ob -= f_val / dfdl;
                if (L_ob * L_sign <= 0) L_ob = 0.5 * L_old;
                L_ob = L_sign * fmax(fmin(fabs(L_ob), (T)1e6), (T)1e-6);
                if (fabs((L_ob - L_old) / L_ob) < tol) break;
            }

            // Final local variables update
            if (Ri_b > 0) {
                utau = kappa * magu / slaw_m_stable<T>(hi, L_ob, z0);
                if constexpr (BC_TYPE == 1) q = kappa * utau * (ts - ti) / slaw_h_stable<T>(hi, L_ob, z0h);
            } else {
                utau = kappa * magu / slaw_m_convective<T>(hi, L_ob, z0);
                if constexpr (BC_TYPE == 1) q = kappa * utau * (ts - ti) / slaw_h_convective<T>(hi, L_ob, z0h);
            }
        }
        tau_x_d[i] = -utau*utau*ui/magu;
        tau_y_d[i] = -utau*utau*vi/magu;
        tau_z_d[i] = -utau*utau*wi/magu;
    }
}

#endif // MOST_KERNEL_H