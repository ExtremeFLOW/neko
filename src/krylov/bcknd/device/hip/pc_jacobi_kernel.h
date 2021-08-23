/**
 * Device kernel for jacobi preconditioner
 */

template< const int LX, const int CHUNKS >
__global__ void jacobi_kernel(double * __restrict__ du,
			      const double * __restrict__ dxt,
			      const double * __restrict__ dyt,
			      const double * __restrict__ dzt,
			      const double * __restrict__ G11,
			      const double * __restrict__ G22,
			      const double * __restrict__ G33,
			      const double * __restrict__ G12,
			      const double * __restrict__ G13,
			      const double * __restrict__ G23) {
  
  
  int i,j,k;
  
  const int e = blockIdx.x;
  const int iii = threadIdx.x;
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;


}

