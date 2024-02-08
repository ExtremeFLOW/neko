#ifndef __MATH_LAMBDA2_KERNEL_CL__
#define __MATH_LAMBDA2_KERNEL_CL__
/*
 Copyright (c) 2023, The Neko Authors
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

   * Redistributions in binary form must reproduce the above
     copyright notice, this list of conditions and the following
     disclaimer in the documentation and/or other materials provided
     with the distribution.

   * Neither the name of the authors nor the names of its
     contributors may be used to endorse or promote products derived
     from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
*/

/**
 * Device kernel for lambda2
 */

inline real
eigen_val_calc(const real grad11,
               const real grad12,
               const real grad13,
               const real grad21,
               const real grad22,
               const real grad23,
               const real grad31,
               const real grad32,
               const real grad33) {
  
    real s11 = grad11;
    real s22 = grad22;
    real s33 = grad33;
    real s12 = 0.5*(grad12+grad21);
    real s13 = 0.5*(grad13+grad31);
    real s23 = 0.5*(grad23+grad32);

    real o12 = 0.5*(grad12-grad21);
    real o13 = 0.5*(grad13-grad31);
    real o23 = 0.5*(grad23-grad32);

    real a11 = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13;
    real a12 = s11 * s12  +  s12 * s22  +  s13 * s23 - o13 * o23;
    real a13 = s11 * s13  +  s12 * s23  +  s13 * s33 + o12 * o23;
        
    real a22 = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23;
    real a23 = s12 * s13 + s22 * s23 + s23 * s33 - o12 * o13;
    real a33 = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23;
               
               
    real B = -(a11 + a22 + a33);
    real C = -(a12*a12 + a13*a13 + a23*a23 - a11 * a22 - a11 * a33 - a22 * a33);
    real D = -(2.0 * a12 * a13 * a23 - a11 * a23*a23 - a22 * a13*a13
            - a33 * a12*a12  +  a11 * a22 * a33);
                     
                     
    real q = (3.0 * C - B*B) / 9.0;
    real r = (9.0 * C * B - 27.0 * D - 2.0 * B*B*B) / 54.0;

    real theta = acos( r / sqrt(-q*q*q) ); 
    real pi = 4.0*atan(1.0);
    real eigen1 = 2.0 * sqrt(-q) * cos(theta / 3.0) - B / 3.0;
    real eigen2 = 2.0 * sqrt(-q) * cos((theta + 2.0 * pi) / 3.0) - B / 3.0;
    real eigen3 = 2.0 * sqrt(-q) * cos((theta + 4.0 * pi) / 3.0) - B / 3.0;
                                 
    if (eigen1 <= eigen2 && eigen2 <= eigen3)
        return eigen2;
    else if (eigen3 <= eigen2 && eigen2 <= eigen1)
        return eigen2;
    else if (eigen1 <= eigen3 && eigen3 <= eigen2)
        return eigen3;
    else if (eigen2 <= eigen3 && eigen3 <= eigen1)
        return eigen3;
    else if (eigen2 <= eigen1 && eigen1 <= eigen3)
        return eigen1;
    else if (eigen3 <= eigen1 && eigen1 <= eigen2)
        return eigen1;
    return 0.0;
}


#define DEFINE_LAMBDA2_KERNEL(LX, CHUNKS)                                      \
__kernel void lambda2_kernel_lx##LX(__global real * __restrict__ lambda2,      \
                                    __global const real * __restrict__ u,      \
                                    __global const real * __restrict__ v,      \
                                    __global const real * __restrict__ w,      \
                                    __global const real * __restrict__ dx,     \
                                    __global const real * __restrict__ dy,     \
                                    __global const real * __restrict__ dz,     \
                                    __global const real * __restrict__ drdx,   \
                                    __global const real * __restrict__ dsdx,   \
                                    __global const real * __restrict__ dtdx,   \
                                    __global const real * __restrict__ drdy,   \
                                    __global const real * __restrict__ dsdy,   \
                                    __global const real * __restrict__ dtdy,   \
                                    __global const real * __restrict__ drdz,   \
                                    __global const real * __restrict__ dsdz,   \
                                    __global const real * __restrict__ dtdz,   \
                                    __global const real * __restrict__ jacinv){\
                                                                               \
                                                                               \
  __local real shu[LX * LX * LX];                                              \
  __local real shv[LX * LX * LX];                                              \
  __local real shw[LX * LX * LX];                                              \
                                                                               \
  __local real shdx[LX * LX];                                                  \
  __local real shdy[LX * LX];                                                  \
  __local real shdz[LX * LX];                                                  \
                                                                               \
  int i,j,k;                                                                   \
                                                                               \
  const int e = get_group_id(0);                                               \
  const int ele = get_group_id(0) * LX * LX * LX;                              \
  const int iii = get_local_id(0);                                             \
  const int nchunks = (LX * LX * LX - 1)/CHUNKS + 1;                           \
                                                                               \
                                                                               \
  if (iii < (LX * LX)) {                                                       \
    shdx[iii] = dx[iii];                                                       \
    shdy[iii] = dy[iii];                                                       \
    shdz[iii] = dz[iii];                                                       \
  }                                                                            \
                                                                               \
  j = iii;                                                                     \
  while(j < (LX * LX * LX)) {                                                  \
    shu[j] = u[j + ele];                                                       \
    shv[j] = v[j + ele];                                                       \
    shw[j] = w[j + ele];                                                       \
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
    if ( i < LX && j < LX && k < LX ) {                                        \
                                                                               \
      real rtmpu = 0.0;                                                        \
      real stmpu = 0.0;                                                        \
      real ttmpu = 0.0;                                                        \
                                                                               \
      real rtmpv = 0.0;                                                        \
      real stmpv = 0.0;                                                        \
      real ttmpv = 0.0;                                                        \
                                                                               \
      real rtmpw = 0.0;                                                        \
      real stmpw = 0.0;                                                        \
      real ttmpw = 0.0;                                                        \
                                                                               \
      for (int l = 0; l < LX; l++) {		                               \
        rtmpu += shdx[i + l * LX] * shu[l + j * LX + k * LX * LX];	       \
        stmpu += shdy[j + l * LX] * shu[i + l * LX + k * LX * LX];             \
        ttmpu += shdz[k + l * LX] * shu[i + j * LX + l * LX * LX];             \
	                                                                       \
        rtmpv += shdx[i + l * LX] * shv[l + j * LX + k * LX * LX];	       \
        stmpv += shdy[j + l * LX] * shv[i + l * LX + k * LX * LX];             \
        ttmpv += shdz[k + l * LX] * shv[i + j * LX + l * LX * LX];             \
	                                                                       \
        rtmpw += shdx[i + l * LX] * shw[l + j * LX + k * LX * LX];	       \
        stmpw += shdy[j + l * LX] * shw[i + l * LX + k * LX * LX];             \
        ttmpw += shdz[k + l * LX] * shw[i + j * LX + l * LX * LX];             \
      }                                                                        \
                                                                               \
      real jinv = jacinv[ijk + ele];                                           \
                                                                               \
      real grad11 = jinv                                                       \
	  * (drdx[ijk + ele] * rtmpu                                           \
	   + dsdx[ijk + ele] * stmpu                                           \
	   + dtdx[ijk + ele] * ttmpu);                                         \
                                                                               \
      real grad12 = jinv                                                       \
	  * (drdy[ijk + ele] * rtmpu                                           \
	   + dsdy[ijk + ele] * stmpu                                           \
	   + dtdy[ijk + ele] * ttmpu);                                         \
                                                                               \
      real grad13 = jinv                                                       \
	  * (drdz[ijk + ele] * rtmpu                                           \
	   + dsdz[ijk + ele] * stmpu                                           \
	   + dtdz[ijk + ele] * ttmpu);                                         \
                                                                               \
      real grad21 = jinv                                                       \
	  * (drdx[ijk + ele] * rtmpv                                           \
	   + dsdx[ijk + ele] * stmpv                                           \
	   + dtdx[ijk + ele] * ttmpv);                                         \
                                                                               \
      real grad22 = jinv                                                       \
	  * (drdy[ijk + ele] * rtmpv                                           \
	   + dsdy[ijk + ele] * stmpv                                           \
	   + dtdy[ijk + ele] * ttmpv);                                         \
                                                                               \
      real grad23 = jinv                                                       \
	  * (drdz[ijk + ele] * rtmpv                                           \
	   + dsdz[ijk + ele] * stmpv                                           \
	   + dtdz[ijk + ele] * ttmpv);                                         \
                                                                               \
      real grad31 = jinv                                                       \
	  * (drdx[ijk + ele] * rtmpw                                           \
	   + dsdx[ijk + ele] * stmpw                                           \
	   + dtdx[ijk + ele] * ttmpw);                                         \
                                                                               \
      real grad32 = jinv                                                       \
	  * (drdy[ijk + ele] * rtmpw                                           \
	   + dsdy[ijk + ele] * stmpw                                           \
	   + dtdy[ijk + ele] * ttmpw);                                         \
                                                                               \
      real grad33 = jinv                                                       \
	  * (drdz[ijk + ele] * rtmpw                                           \
	   + dsdz[ijk + ele] * stmpw                                           \
	   + dtdz[ijk + ele] * ttmpw);                                         \
      lambda2[ijk + e*LX*LX*LX] = eigen_val_calc( grad11, grad12, grad13,      \
                                                  grad21, grad22, grad23,      \
                                                  grad31, grad32, grad33);     \
    }                                                                          \
  }                                                                            \
}                              
                                                                               
DEFINE_LAMBDA2_KERNEL(1, 256)
DEFINE_LAMBDA2_KERNEL(2, 256)
DEFINE_LAMBDA2_KERNEL(3, 256)
DEFINE_LAMBDA2_KERNEL(4, 256)
DEFINE_LAMBDA2_KERNEL(5, 256)
DEFINE_LAMBDA2_KERNEL(6, 256)
DEFINE_LAMBDA2_KERNEL(7, 256)
DEFINE_LAMBDA2_KERNEL(8, 256)
DEFINE_LAMBDA2_KERNEL(9, 256)
DEFINE_LAMBDA2_KERNEL(10, 256)
DEFINE_LAMBDA2_KERNEL(11, 256)
DEFINE_LAMBDA2_KERNEL(12, 256)


#endif // __MATH_LAMBDA2_KERNEL_CL__
