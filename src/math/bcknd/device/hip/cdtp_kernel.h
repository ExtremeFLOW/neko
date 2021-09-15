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
  
  __shared__ T shjac[LX * LX * LX];
  __shared__ T shB[LX * LX * LX];
  
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
    shx[j] = x[j + e * LX * LX * LX];

    shdr[j] = dr[j + e * LX * LX * LX];
    shds[j] = ds[j + e * LX * LX * LX];
    shdt[j] = dt[j + e * LX * LX * LX];

    shB[j] = B[j + e * LX * LX * LX];
    
    shjac[j] = jac[j + e * LX * LX * LX];

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

	// col3(wx, B, x) + invcol2(wx, jac) + col3(ta1, wx, dr)
	const T ta1_r = (((shB[l + j * LX + k * LX * LX]
				* shx[l + j * LX + k * LX * LX]) /
			       shjac[l + j * LX + k * LX * LX]) *
			      shdr[l + j * LX + k * LX * LX]);

	// col3(wx, B, x) + invcol2(wx, jac) + col3(ta1, wx, ds)
	const T ta1_s = (((shB[i + l * LX + k * LX * LX]
				* shx[i + l * LX + k * LX * LX]) /
			       shjac[i + l * LX + k * LX * LX])
			      * shds[i + l * LX + k * LX * LX]);

	// col3(wx, B, x) + invcol2(wx, jac) + col3(ta1, wx, dt)
	const T ta1_t = (((shB[i + j * LX + l * LX * LX]
				* shx[i + j * LX + l * LX * LX]) /
			       shjac[i + j * LX + l * LX * LX])
			      * shdt[i + j * LX + l * LX * LX]);
		
	rtmp += shdxt[i + l * LX] * ta1_r;	
	stmp += shdyt[j + l * LX] * ta1_s;
	ttmp += shdzt[k + l * LX] * ta1_t;

      }
      
      dtx[ijk + e * LX * LX * LX] = ( rtmp + stmp + ttmp );
      
    }
  }
}

