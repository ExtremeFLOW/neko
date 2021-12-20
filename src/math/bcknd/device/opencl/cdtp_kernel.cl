/**
 * Device kernel for \f$ D^T x \f$
 */

#define DEFINE_CDTP_KERNEL(LX, CHUNKS)                                         \
__kernel void cdtp_kernel_lx##LX(__global real * __restrict__ dtx,             \
				 __global const real * __restrict__ x,	       \
				 __global const real * __restrict__ dr,	       \
				 __global const real * __restrict__ ds,	       \
				 __global const real * __restrict__ dt,	       \
				 __global const real * __restrict__ dxt,       \
				 __global const real * __restrict__ dyt,       \
				 __global const real * __restrict__ dzt,       \
				 __global const real * __restrict__ B,         \
				 __global const real * __restrict__ jac) {     \
                                                                               \
  __local real shx[LX * LX * LX];                                              \
  __local real shdr[LX * LX * LX];                                             \
  __local real shds[LX * LX * LX];                                             \
  __local real shdt[LX * LX * LX];                                             \
                                                                               \
  __local real shdxt[LX * LX];                                                 \
  __local real shdyt[LX * LX];                                                 \
  __local real shdzt[LX * LX];                                                 \
                                                                               \
  __local real shjac[LX * LX * LX];                                            \
  __local real shB[LX * LX * LX];                                              \
                                                                               \
  int i,j,k;                                                                   \
                                                                               \
  const int e = get_group_id(0);   				               \
  const int iii = get_local_id(0);					       \
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;                         \
                                                                               \
  if (iii < (LX * LX)) {                                                       \
    shdxt[iii] = dxt[iii];                                                     \
    shdyt[iii] = dyt[iii];                                                     \
    shdzt[iii] = dzt[iii];                                                     \
  }                                                                            \
                                                                               \
  j = iii;                                                                     \
  while(j < (LX * LX * LX)) {                                                  \
    shx[j] = x[j + e * LX * LX * LX];                                          \
                                                                               \
    shdr[j] = dr[j + e * LX * LX * LX];                                        \
    shds[j] = ds[j + e * LX * LX * LX];                                        \
    shdt[j] = dt[j + e * LX * LX * LX];                                        \
                                                                               \
    shB[j] = B[j + e * LX * LX * LX];                                          \
                                                                               \
    shjac[j] = jac[j + e * LX * LX * LX];                                      \
                                                                               \
    j = j + CHUNKS;                                                            \
  }                                                                            \
                                                                               \
  barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                               \
  for (int n = 0; n < nchunks; n++) {                                          \
    const int ijk = iii + n * CHUNKS;                                          \
    const int jk = ijk / LX;                                                   \
    i = ijk - jk * LX;                                                         \
    k = jk / LX;                                                               \
    j = jk - k * LX;                                                           \
    if ( i < LX && j < LX && k < LX) {                                         \
      real rtmp = 0.0;                                                         \
      real stmp = 0.0;                                                         \
      real ttmp = 0.0;                                                         \
      for (int l = 0; l < LX; l++) {                                           \
                                                                               \
	const real ta1_r = (((shB[l + j * LX + k * LX * LX]                    \
			      * shx[l + j * LX + k * LX * LX]) /               \
			     shjac[l + j * LX + k * LX * LX])                  \
			    * shdr[l + j * LX + k * LX * LX]);     	       \
                                                                               \
	const real ta1_s = (((shB[i + l * LX + k * LX * LX]                    \
			      * shx[i + l * LX + k * LX * LX]) /               \
			     shjac[i + l * LX + k * LX * LX])                  \
			    * shds[i + l * LX + k * LX * LX]);                 \
                                                                               \
	const real  ta1_t = (((shB[i + j * LX + l * LX * LX]                   \
			       * shx[i + j * LX + l * LX * LX]) /              \
			      shjac[i + j * LX + l * LX * LX])                 \
			     * shdt[i + j * LX + l * LX * LX]);                \
		                                                               \
	rtmp += shdxt[i + l * LX] * ta1_r;	                               \
	stmp += shdyt[j + l * LX] * ta1_s;                                     \
	ttmp += shdzt[k + l * LX] * ta1_t;                                     \
                                                                               \
      }                                                                        \
                                                                               \
      dtx[ijk + e * LX * LX * LX] = ( rtmp + stmp + ttmp );                    \
                                                                               \
    }                                                                          \
  }                                                                            \
}                                                                              

DEFINE_CDTP_KERNEL(2, 256)
DEFINE_CDTP_KERNEL(3, 256)
DEFINE_CDTP_KERNEL(4, 256)
DEFINE_CDTP_KERNEL(5, 256)
DEFINE_CDTP_KERNEL(6, 256)
DEFINE_CDTP_KERNEL(7, 256)
DEFINE_CDTP_KERNEL(8, 256)



