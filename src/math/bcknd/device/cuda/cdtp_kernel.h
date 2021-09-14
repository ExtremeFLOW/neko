/**
 * Device kernel for \f$ D^T x \f$
 */

template< const int LX, const int CHUNKS >
__global__ void cdtp_kernel(double * __restrict__ dtx,
			    const double * __restrict__ x,
			    const double * __restrict__ dr,
			    const double * __restrict__ ds,
			    const double * __restrict__ dt,
			    const double * __restrict__ dxt,
			    const double * __restrict__ dyt,
			    const double * __restrict__ dzt,
			    const double * __restrict__ B,
			    const double * __restrict__ jac) { 
  
  __shared__ double shx[LX * LX * LX];
  __shared__ double shdr[LX * LX * LX];
  __shared__ double shds[LX * LX * LX];
  __shared__ double shdt[LX * LX * LX];

  __shared__ double shdxt[LX * LX];
  __shared__ double shdyt[LX * LX];
  __shared__ double shdzt[LX * LX];
  
  __shared__ double shjac[LX * LX * LX];
  __shared__ double shB[LX * LX * LX];
  
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
      double rtmp = 0.0;
      double stmp = 0.0;
      double ttmp = 0.0;
      for (int l = 0; l < LX; l++) {

	// col3(wx, B, x) + invcol2(wx, jac) + col3(ta1, wx, dr)
	const double ta1_r = (((shB[l + j * LX + k * LX * LX]
				* shx[l + j * LX + k * LX * LX]) /
			       shjac[l + j * LX + k * LX * LX]) *
			      shdr[l + j * LX + k * LX * LX]);

	// col3(wx, B, x) + invcol2(wx, jac) + col3(ta1, wx, ds)
	const double ta1_s = (((shB[i + l * LX + k * LX * LX]
				* shx[i + l * LX + k * LX * LX]) /
			       shjac[i + l * LX + k * LX * LX])
			      * shds[i + l * LX + k * LX * LX]);

	// col3(wx, B, x) + invcol2(wx, jac) + col3(ta1, wx, dt)
	const double ta1_t = (((shB[i + j * LX + l * LX * LX]
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

