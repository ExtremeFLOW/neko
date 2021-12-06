/**
 * Device kernel for derivative 
 */

template< typename T, const int LX, const int CHUNKS >
__global__ void dudxyz_kernel(T * __restrict__ du,
			      const T * __restrict__ u,
			      const T * __restrict__ dr,
			      const T * __restrict__ ds,
			      const T * __restrict__ dt,
			      const T * __restrict__ dx,
			      const T * __restrict__ dy,
			      const T * __restrict__ dz,
			      const T * __restrict__ jacinv) { 
  
  __shared__ T shu[LX * LX * LX];
  __shared__ T shdr[LX * LX * LX];
  __shared__ T shds[LX * LX * LX];
  __shared__ T shdt[LX * LX * LX];

  __shared__ T shdx[LX * LX];
  __shared__ T shdy[LX * LX];
  __shared__ T shdz[LX * LX];
  
  __shared__ T shjacinv[LX * LX * LX];
  
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
      T rtmp = 0.0;
      T stmp = 0.0;
      T ttmp = 0.0;
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

