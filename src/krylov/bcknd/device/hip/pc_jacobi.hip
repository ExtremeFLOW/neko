#include <hip/hip_runtime.h>
#include <device/hip/check.h>
#include <device/device_config.h>

template< typename T, const int LX >
__global__ void jacobi_kernel(T * __restrict__ du,
			      const T * __restrict__ dxt,
			      const T * __restrict__ dyt,
			      const T * __restrict__ dzt,
			      const T * __restrict__ G11,
			      const T * __restrict__ G22,
			      const T * __restrict__ G33,
			      const T * __restrict__ G12,
			      const T * __restrict__ G13,
			      const T * __restrict__ G23,
			      const int nel) {
  const int idx = threadIdx.x + blockIdx.x * blockDim.x;
  const int e = idx / (LX*LX*LX);
  const int ijk = idx - e*LX*LX*LX;
  const int jk = ijk / LX;
  const int i = ijk - jk * LX;
  const int k = jk / LX;
  const int j = jk - k * LX;



  if (e >= nel)
    return;

  T d = 0.0;

  for (int l = 0; l < LX; l++) {
    T g = G11[l + LX*j + LX*LX*k + LX*LX*LX*e];
    T t = dxt[i + LX*l];
    d += g*t*t;
  }

  for (int l = 0; l < LX; l++) {
    T g = G22[i + LX*l + LX*LX*k + LX*LX*LX*e];
    T t = dyt[j + LX*l];
    d += g*t*t;
  }

  for (int l = 0; l < LX; l++) {
    T g = G33[i + LX*j + LX*LX*l + LX*LX*LX*e];
    T t = dzt[k + LX*l];
    d += g*t*t;
  }

  // Corrections for deformed elements
  if (i == 0 || i == LX-1) {
    d += G12[i + LX*j + LX*LX*k + LX*LX*LX*e] * dxt[i + LX*i] * dyt[j + LX*j];
    d += G13[i + LX*j + LX*LX*k + LX*LX*LX*e] * dxt[i + LX*i] * dzt[k + LX*k];
  }

  if (j == 0 || j == LX-1) {
    d += G12[i + LX*j + LX*LX*k + LX*LX*LX*e] * dyt[j + LX*j] * dxt[i + LX*i];
    d += G23[i + LX*j + LX*LX*k + LX*LX*LX*e] * dyt[j + LX*j] * dzt[k + LX*k];
  }

  if (k == 0 || k == LX-1) {
    d += G13[i + LX*j + LX*LX*k + LX*LX*LX*e] * dzt[k + LX*k] * dxt[i + LX*i];
    d += G23[i + LX*j + LX*LX*k + LX*LX*LX*e] * dzt[k + LX*k] * dyt[j + LX*j];
  }



  du[idx] = d;
}

extern "C" {
  void hip_jacobi_update(void *d, 
			 void *dxt, void *dyt, void *dzt,
			 void *G11, void *G22, void *G33,
			 void *G12, void *G13, void *G23,
			 int *nel, int *lxp) {
    
    const int lx = *lxp;
    const int threads = 1024;
    const int blocks = ((*nel * lx*lx*lx) + threads - 1) / threads;

#define CASE(N)\
    case N:\
    hipLaunchKernelGGL(HIP_KERNEL_NAME( jacobi_kernel<real, N> ),\
	blocks, threads, 0, 0,\
	(double *) d,\
	(double *) dxt, (double *) dyt, (double *) dzt,\
	(double *) G11, (double *) G22, (double *) G33,\
	(double *) G12, (double *) G13, (double *) G23,\
	*nel);\
    break

    switch (lx) {
    CASE(1);
    CASE(2);
    CASE(3);
    CASE(4);
    CASE(5);
    CASE(6);
    CASE(7);
    CASE(8);
    CASE(9);
    CASE(10);
    CASE(11);
    CASE(12);
    CASE(13);
    CASE(14);
    CASE(15);
    default:
      fprintf(stderr, __FILE__ ": size not supported: %d\n", lx);
    }

    HIP_CHECK(hipGetLastError());
  }
}
