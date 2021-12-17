#include <stdio.h>
#include <stdlib.h>

template<const int LX>
__global__ void jacobi_kernel(double * __restrict__ du,
			      const double * __restrict__ dxt,
			      const double * __restrict__ dyt,
			      const double * __restrict__ dzt,
			      const double * __restrict__ G11,
			      const double * __restrict__ G22,
			      const double * __restrict__ G33,
			      const int nel) {
  const int idx = threadIdx.x + blockIdx.x * blockDim.x;
  const int e = idx / (LX*LX*LX);
  const int k = idx / (LX*LX) % LX;
  const int j = idx / LX % LX;
  const int i = idx % LX;

  if (e >= nel)
    return;

  double d = 0.0;

  for (int l = 0; l < LX; l++) {
    double g = G11[l + LX*j + LX*LX*k + LX*LX*LX*e];
    double t = dxt[i + LX*l];
    d += g*t*t;
  }

  for (int l = 0; l < LX; l++) {
    double g = G22[i + LX*l + LX*LX*k + LX*LX*LX*e];
    double t = dyt[j + LX*l];
    d += g*t*t;
  }

  for (int l = 0; l < LX; l++) {
    double g = G33[i + LX*j + LX*LX*l + LX*LX*LX*e];
    double t = dzt[k + LX*l];
    d += g*t*t;
  }

  du[idx] = d;
}

extern "C" {
  void cuda_jacobi_update(void *d,
			  void *dxt, void *dyt, void *dzt,
			  void *G11, void *G22, void *G33,
			  void *G12, void *G13, void *G23,
			  int *nel, int *lxp) {

    const int lx = *lxp;
    const int threads = 1024;
    const int blocks = ((*nel * lx*lx*lx) + threads - 1) / threads;

#define CASE(N)\
    case N:\
    jacobi_kernel<N><<<blocks, threads>>>(\
	(double*)d,\
	(double*)dxt, (double*)dyt, (double*)dzt,\
	(double*)G11, (double*)G22, (double*)G33,\
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

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
      fprintf(stderr, __FILE__ ": %s\n", cudaGetErrorString(err));
    }
  }
}
