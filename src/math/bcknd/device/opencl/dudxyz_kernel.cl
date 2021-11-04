/**
 * Device kernel for derivative 
 */
                                                              
#define DEFINE_DUDXYZ_KERNEL(LX, CHUNKS)                                       \
__kernel void dudxyz_kernel_lx##LX(__global real *du,                          \
				   __global const real *u,		       \
				   __global const real *dr,		       \
				   __global const real *ds,		       \
				   __global const real *dt,		       \
				   __global const real *dx,		       \
				   __global const real *dy,		       \
				   __global const real *dz,		       \
				   __global const real *jacinv) {	       \
                                                                               \
  __local real shu[LX * LX * LX];                                              \
  __local real shdr[LX * LX * LX];                                             \
  __local real shds[LX * LX * LX];                                             \
  __local real shdt[LX * LX * LX];                                             \
                                                                               \
  __local real shdx[LX * LX];                                                  \
  __local real shdy[LX * LX];                                                  \
  __local real shdz[LX * LX];                                                  \
                                                                               \
  __local real shjacinv[LX * LX * LX];                                         \
                                                                               \
  int i,j,k;                                                                   \
                                                                               \
  const int e = get_group_id(0);                                               \
  const int iii = get_local_id(0);                                             \
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;                         \
                                                                               \
  if (iii < (LX * LX)) {                                                       \
    shdx[iii] = dx[iii];                                                       \
    shdy[iii] = dy[iii];                                                       \
    shdz[iii] = dz[iii];                                                       \
  }                                                                            \
                                                                               \
  j = iii;                                                                     \
  while(j < (LX * LX * LX)) {                                                  \
    shu[j] = u[j + e * LX * LX * LX];                                          \
    shdr[j] = dr[j + e * LX * LX * LX];                                        \
    shds[j] = ds[j + e * LX * LX * LX];                                        \
    shdt[j] = dt[j + e * LX * LX * LX];                                        \
    shjacinv[j] = jacinv[j + e * LX * LX * LX];                                \
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
	rtmp += shdx[i + l * LX] * shu[l + j * LX + k * LX * LX];              \
	stmp += shdy[j + l * LX] * shu[i + l * LX + k * LX * LX];              \
	ttmp += shdz[k + l * LX] * shu[i + j * LX + l * LX * LX];              \
      }                                                                        \
      du[ijk + e * LX * LX * LX] = ((rtmp * shdr[ijk])                         \
				    + (stmp * shds[ijk])                       \
				    + (ttmp * shdt[ijk]))                      \
	                           * shjacinv[ijk];                            \
                                                                               \
    }                                                                          \
  }                                                                            \
}                                                                              \

DEFINE_DUDXYZ_KERNEL(12, 256)
DEFINE_DUDXYZ_KERNEL(11, 256)
DEFINE_DUDXYZ_KERNEL(10, 256)
DEFINE_DUDXYZ_KERNEL(9, 256)
DEFINE_DUDXYZ_KERNEL(8, 256)
DEFINE_DUDXYZ_KERNEL(7, 256)
DEFINE_DUDXYZ_KERNEL(6, 256)
DEFINE_DUDXYZ_KERNEL(5, 256)
DEFINE_DUDXYZ_KERNEL(4, 256)
DEFINE_DUDXYZ_KERNEL(3, 256)
DEFINE_DUDXYZ_KERNEL(2, 256)
