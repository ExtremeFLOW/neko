/**
 * Device kernel for derivative 
 */

template< const int LX, const int CHUNKS >
__global__ void dudxyz_kernel(double * __restrict__ du,
			      const double * __restrict__ u,
			      const double * __restrict__ dr,
			      const double * __restrict__ ds,
			      const double * __restrict__ dt,
			      const double * __restrict__ dx,
			      const double * __restrict__ dy,
			      const double * __restrict__ dz,
			      const double * __restrict__ jacinv) { 
  
  __shared__ double shu[LX * LX * LX];
  __shared__ double shdr[LX * LX * LX];
  __shared__ double shds[LX * LX * LX];
  __shared__ double shdt[LX * LX * LX];

  __shared__ double shdx[LX * LX];
  __shared__ double shdy[LX * LX];
  __shared__ double shdz[LX * LX];
  
  __shared__ double shjacinv[LX * LX * LX];
  
  int i,j,k;
  
  const int e = blockIdx.x;
  const int iii = threadIdx.x;
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;

  if (iii < (LX * LX)) {
    shdx[iii] = dx[iii];
    shdy[iii] = dy[iii];
    shdz[iii] = dz[iii];
  }

  j = iii;
  while(j < (LX * LX * LX)) {
    shu[j] = u[j + e * LX * LX * LX];
    shdr[j] = dr[j + e * LX * LX * LX];
    shds[j] = ds[j + e * LX * LX * LX];
    shdt[j] = dt[j + e * LX * LX * LX];
    shjacinv[j] = jacinv[j + e * LX * LX * LX];
    j = j + CHUNKS;
  }

  __syncthreads();
  
  for (int n = 0; n < nchunks; n++) {
    const int ijk = iii + n * CHUNKS;
    const int jk = ijk / LX;
    i = ijk - jk * LX;
    k = jk / LX;
    j = jk - k * LX;
    if ( i < LX && j < LX && k < LX) {
      double rtmp = 0.0;
      double stmp = 0.0;
      double ttmp = 0.0;
      for (int l = 0; l < LX; l++) {
	rtmp += shdx[i + l * LX] * shu[l + j * LX + k * LX * LX];
	stmp += shdy[j + l * LX] * shu[i + l * LX + k * LX * LX];
	ttmp += shdz[k + l * LX] * shu[i + j * LX + l * LX * LX];
      }
      du[ijk + e * LX * LX * LX] = ((rtmp * shdr[ijk])
				    + (stmp * shds[ijk])
				    + (ttmp * shdt[ijk]))
	                           * shjacinv[ijk];
      
    }
  }
}

