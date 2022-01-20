/**
 * Device kernel for \f$ D^T x \f$
 */

template< typename T, const int LX, const int CHUNKS >
__global__ void cdtp_kernel(T * __restrict__ dtx,
			    const T * __restrict__ x,
			    const T * __restrict__ dr,
			    const T * __restrict__ ds,
			    const T * __restrict__ dt,
			    const T * __restrict__ dxt,
			    const T * __restrict__ dyt,
			    const T * __restrict__ dzt,
			    const T * __restrict__ B,
			    const T * __restrict__ jac) { 
  
  __shared__ T shx[LX * LX * LX];
  __shared__ T shdr[LX * LX * LX];
  __shared__ T shds[LX * LX * LX];
  __shared__ T shdt[LX * LX * LX];

  __shared__ T shdxt[LX * LX];
  __shared__ T shdyt[LX * LX];
  __shared__ T shdzt[LX * LX];
  
  __shared__ T shtar[LX * LX * LX];
  __shared__ T shtas[LX * LX * LX];
  __shared__ T shtat[LX * LX * LX];
  
  int i,j,k;
  
  const int e = blockIdx.x;
  const int iii = threadIdx.x;
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;

  if (iii < (LX * LX)) {
    shdxt[iii] = dxt[iii];
    shdyt[iii] = dyt[iii];
    shdzt[iii] = dzt[iii];
  }

  j = iii;
  while(j < (LX * LX * LX)) {
     // We can probably avoid this division by use of jacinv instead.
    T wx = (x[j + e * LX * LX * LX] * B[j + e * LX * LX * LX])/jac[j + e * LX * LX * LX];

    shtar[j] = wx*dr[j + e * LX * LX * LX];
    shtas[j] = wx*ds[j + e * LX * LX * LX];
    shtat[j] = wx*dt[j + e * LX * LX * LX];

    j = j + CHUNKS;
  }

  __syncthreads();
  for (int n = 0; n < nchunks; n++) {
    const int ijk = iii + n * CHUNKS;
    const int jk = ijk / LX;
    i = ijk - jk * LX;
    k = jk / LX;
    j = jk - k * LX;
    if ( i < LX && j < LX && k < LX && ijk < LX*LX*LX) {
      T rtmp = 0.0;
      T stmp = 0.0;
      T ttmp = 0.0;
      for (int l = 0; l < LX; l++) {
	    rtmp += shdxt[i + l * LX] * shtar[l+j*LX+k*LX*LX];	
	    stmp += shdyt[j + l * LX] * shtas[i+l*LX + k*LX*LX];
	    ttmp += shdzt[k + l * LX] * shtat[i + j*LX + l*LX*LX];
      }
      dtx[ijk + e * LX * LX * LX] = ( rtmp + stmp + ttmp );
      
    }
  }
}

