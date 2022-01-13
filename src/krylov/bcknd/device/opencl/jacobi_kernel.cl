/**
 * Device kernel for jacobi
 */

#define DEFINE_JACOBI_KERNEL(LX)                                               \
__kernel void jacobi_kernel_lx##LX(__global real * __restrict__ du,            \
                                   __global const real * __restrict__ dxt,     \
                                   __global const real * __restrict__ dyt,     \
                                   __global const real * __restrict__ dzt,     \
                                   __global const real * __restrict__ G11,     \
                                   __global const real * __restrict__ G22,     \
                                   __global const real * __restrict__ G33,     \
                                   __global const real * __restrict__ G12,     \
                                   __global const real * __restrict__ G13,     \
                                   __global const real * __restrict__ G23,     \
                                   const int nel) {                            \
                                                                               \
  const int idx = get_global_id(0);                                            \
  const int e = idx / (LX*LX*LX);                                              \
  const int k = idx / (LX*LX) % LX;                                            \
  const int j = idx / LX % LX;                                                 \
  const int i = idx % LX;                                                      \
                                                                               \
  if (e >= nel)                                                                \
    return;                                                                    \
                                                                               \
  real d = 0.0;                                                                \
                                                                               \
  for (int l = 0; l < LX; l++) {                                               \
    real g = G11[l + LX*j + LX*LX*k + LX*LX*LX*e];                             \
    real t = dxt[i + LX*l];                                                    \
    d += g*t*t;                                                                \
  }                                                                            \
                                                                               \
  for (int l = 0; l < LX; l++) {                                               \
    real g = G22[i + LX*l + LX*LX*k + LX*LX*LX*e];                             \
    real t = dyt[j + LX*l];                                                    \
    d += g*t*t;                                                                \
  }                                                                            \
                                                                               \
  for (int l = 0; l < LX; l++) {                                               \
    real g = G33[i + LX*j + LX*LX*l + LX*LX*LX*e];                             \
    real t = dzt[k + LX*l];                                                    \
    d += g*t*t;                                                                \
  }                                                                            \
                                                                               \
  /* Corrections for deformed elements */                                      \
  if (i == 0 || i == LX-1) {                                                   \
    d += G12[i + LX*j + LX*LX*k + LX*LX*LX*e] * dxt[i + LX*i] * dyt[j + LX*j]; \
    d += G13[i + LX*j + LX*LX*k + LX*LX*LX*e] * dxt[i + LX*i] * dzt[k + LX*k]; \
  }                                                                            \
                                                                               \
  if (j == 0 || j == LX-1) {                                                   \
    d += G12[i + LX*j + LX*LX*k + LX*LX*LX*e] * dyt[i + LX*i] * dxt[j + LX*j]; \
    d += G23[i + LX*j + LX*LX*k + LX*LX*LX*e] * dyt[i + LX*i] * dzt[k + LX*k]; \
  }                                                                            \
                                                                               \
  if (k == 0 || k == LX-1) {                                                   \
    d += G13[i + LX*j + LX*LX*k + LX*LX*LX*e] * dzt[i + LX*i] * dxt[j + LX*j]; \
    d += G23[i + LX*j + LX*LX*k + LX*LX*LX*e] * dzt[i + LX*i] * dyt[k + LX*k]; \
  }                                                                            \
                                                                               \
  du[idx] = d;                                                                 \
}

DEFINE_JACOBI_KERNEL(2)
DEFINE_JACOBI_KERNEL(3)
DEFINE_JACOBI_KERNEL(4)
DEFINE_JACOBI_KERNEL(5)
DEFINE_JACOBI_KERNEL(6)
DEFINE_JACOBI_KERNEL(7)
DEFINE_JACOBI_KERNEL(8)
DEFINE_JACOBI_KERNEL(9)
DEFINE_JACOBI_KERNEL(10)
DEFINE_JACOBI_KERNEL(11)
DEFINE_JACOBI_KERNEL(12)
DEFINE_JACOBI_KERNEL(13)
DEFINE_JACOBI_KERNEL(14)
DEFINE_JACOBI_KERNEL(15)
DEFINE_JACOBI_KERNEL(16)

