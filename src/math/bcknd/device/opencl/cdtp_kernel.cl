/**
 * Device kernel for \f$ D^T x \f$
 */

#define DEFINE_CDTP_KERNEL(LX, CHUNKS)                                         \
__kernel void cdtp_kernel_lx##LX(__global real * __restrict__ dtx,             \
                                 __global const real * __restrict__ x,         \
                                 __global const real * __restrict__ dr,        \
                                 __global const real * __restrict__ ds,        \
                                 __global const real * __restrict__ dt,        \
                                 __global const real * __restrict__ dxt,       \
                                 __global const real * __restrict__ dyt,       \
                                 __global const real * __restrict__ dzt,       \
                                 __global const real * __restrict__ B,         \
                                 __global const real * __restrict__ jac) {     \
                                                                               \
                                                                               \
  __local real shdxt[LX * LX];                                                 \
  __local real shdyt[LX * LX];                                                 \
  __local real shdzt[LX * LX];                                                 \
                                                                               \
  __local real shtar[LX * LX * LX];                                            \
  __local real shtas[LX * LX * LX];                                            \
  __local real shtat[LX * LX * LX];                                            \
                                                                               \
  int i,j,k;                                                                   \
                                                                               \
  const int e = get_group_id(0);                                               \
  const int iii = get_local_id(0);                                             \
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
                                                                               \
    real wx = (x[j + e * LX * LX * LX] * B[j + e * LX * LX * LX]) /            \
      jac[j + e * LX * LX * LX];                                               \
                                                                               \
    shtar[j] = wx*dr[j + e * LX * LX * LX];                                    \
    shtas[j] = wx*ds[j + e * LX * LX * LX];                                    \
    shtat[j] = wx*dt[j + e * LX * LX * LX];                                    \
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
    if ( i < LX && j < LX && k < LX && ijk < LX*LX*LX) {                       \
      real rtmp = 0.0;                                                         \
      real stmp = 0.0;                                                         \
      real ttmp = 0.0;                                                         \
      for (int l = 0; l < LX; l++) {                                           \
        rtmp += shdxt[i + l * LX] * shtar[l+j*LX+k*LX*LX];                     \
        stmp += shdyt[j + l * LX] * shtas[i+l*LX + k*LX*LX];                   \
        ttmp += shdzt[k + l * LX] * shtat[i + j*LX + l*LX*LX];                 \
      }                                                                        \
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



