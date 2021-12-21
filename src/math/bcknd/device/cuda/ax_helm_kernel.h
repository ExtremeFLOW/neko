/**
 * Device kernel for axhelm
 */

template< typename T, const int LX >
__global__ void ax_helm_kernel(T * __restrict__ w,
			       const T * __restrict__ u,
			       const T * __restrict__ dx,
			       const T * __restrict__ dy,
			       const T * __restrict__ dz,
			       const T * __restrict__ h1,
			       const T * __restrict__ g11,
			       const T * __restrict__ g22,
			       const T * __restrict__ g33,
			       const T * __restrict__ g12,
			       const T * __restrict__ g13,
			       const T * __restrict__ g23) {

  __shared__ T shdx[LX * LX];
  __shared__ T shdy[LX * LX];
  __shared__ T shdz[LX * LX];
    
  __shared__ T shu[LX * LX];
  __shared__ T shur[LX * LX];
  __shared__ T shus[LX * LX];

  T ru[LX];
  T rw[LX];
  T rut;

  const int e = blockIdx.x;
  const int j = threadIdx.y;
  const int i = threadIdx.x;
  const int ij = i + j*LX;
  const int ele = e*LX*LX*LX;

  shdx[ij] = dx[ij];
  shdy[ij] = dy[ij];
  shdz[ij] = dz[ij];
  
#pragma unroll
  for(int k = 0; k < LX; ++k){
    ru[k] = u[ij + k*LX*LX + ele];
    rw[k] = 0.0;
  }


  __syncthreads();
#pragma unroll
  for (int k = 0; k < LX; ++k){
    const int ijk = ij + k*LX*LX; 
    const T G00 = g11[ijk+ele];
    const T G11 = g22[ijk+ele];
    const T G22 = g33[ijk+ele]; 
    const T G01 = g12[ijk+ele];
    const T G02 = g13[ijk+ele];
    const T G12 = g23[ijk+ele];
    const T H1  = h1[ijk+ele];
    T ttmp = 0.0;
    shu[ij] = ru[k];
    for (int l = 0; l < LX; l++){
      ttmp += shdz[k+l*LX] * ru[l];
    }
    __syncthreads();
    
    T rtmp = 0.0;
    T stmp = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++){
      rtmp += shdx[i+l*LX] * shu[l+j*LX];
      stmp += shdy[j+l*LX] * shu[i+l*LX];
    }
    shur[ij] = H1
	     * (G00 * rtmp
		+ G01 * stmp
		+ G02 * ttmp);
    shus[ij] = H1
	     * (G01 * rtmp
		+ G11 * stmp
		+ G12 * ttmp);
    rut      = H1
	     * (G02 * rtmp
		+ G12 * stmp 
		+ G22 * ttmp);

    __syncthreads();

    T wijke = 0.0;
#pragma unroll
    for (int l = 0; l < LX; l++){
      wijke += shur[l+j*LX] * shdx[l+i*LX];
      rw[l] += rut * shdz[k+l*LX];
      wijke += shus[i+l*LX] * shdy[l + j*LX];
    }
    rw[k] += wijke;
  }
#pragma unroll
  for (int k = 0; k < LX; ++k){
    w[ij + k*LX*LX + ele] = rw[k]; 
  }
}
