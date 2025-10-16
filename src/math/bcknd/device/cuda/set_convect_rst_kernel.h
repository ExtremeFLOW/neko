#ifndef __MATH_SET_CONVECT_RST_KERNEL_H__
#define __MATH_SET_CONVECT_RST_KERNEL_H__
/*
 Copyright (c) 2021-2023, The Neko Authors
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
 * Device kernel for convective terms
 */

template< typename T, const int LX, const int CHUNKS >
__global__ void set_convect_rst_kernel_1d(T * __restrict__ cr,
                                 T * __restrict__ cs,
                                 T * __restrict__ ct,
                                 const T * __restrict__ cx,
                                 const T * __restrict__ cy,
                                 const T * __restrict__ cz,
                                 const T * __restrict__ drdx,
                                 const T * __restrict__ dsdx,
                                 const T * __restrict__ dtdx,
                                 const T * __restrict__ drdy,
                                 const T * __restrict__ dsdy,
                                 const T * __restrict__ dtdy,
                                 const T * __restrict__ drdz,
                                 const T * __restrict__ dsdz,
                                 const T * __restrict__ dtdz,
                                 const T * __restrict__ w3) {

  int i,j,k;

  const int e = blockIdx.x;
  const int iii = threadIdx.x;
  const int nchunks = (LX * LX * LX - 1) / CHUNKS + 1;

  for (int n = 0; n < nchunks; n++) {
    const int ijk = iii + n * CHUNKS;
    const int jk = ijk / LX;
    i = ijk - jk * LX;
    k = jk / LX;
    j = jk - k * LX;
    if ( i < LX && j < LX && k < LX ) {

      cr[ijk + e * LX * LX * LX] = w3[ijk]
	* (drdx[ijk + e * LX * LX * LX] * cx[ijk + e * LX * LX * LX]
	   + drdy[ijk + e * LX * LX * LX] * cy[ijk + e * LX * LX * LX]
	   + drdz[ijk + e * LX * LX * LX] * cz[ijk + e * LX * LX * LX]);

      cs[ijk + e * LX * LX * LX] = w3[ijk]
	* (dsdx[ijk + e * LX * LX * LX] * cx[ijk + e * LX * LX * LX]
	   + dsdy[ijk + e * LX * LX * LX] * cy[ijk + e * LX * LX * LX]
	   + dsdz[ijk + e * LX * LX * LX] * cz[ijk + e * LX * LX * LX]);

      ct[ijk + e * LX * LX * LX] = w3[ijk]
	* (dtdx[ijk + e * LX * LX * LX] * cx[ijk + e * LX * LX * LX]
	   + dtdy[ijk + e * LX * LX * LX] * cy[ijk + e * LX * LX * LX]
	   + dtdz[ijk + e * LX * LX * LX] * cz[ijk + e * LX * LX * LX]);

    }
  }

}

template< typename T, const int LX >
__global__ void __launch_bounds__(LX*LX,3)
set_convect_rst_kernel_kstep(T * __restrict__ cr,
                    T * __restrict__ cs,
                    T * __restrict__ ct,
                    const T * __restrict__ cx,
                    const T * __restrict__ cy,
                    const T * __restrict__ cz,
                    const T * __restrict__ drdx,
                    const T * __restrict__ dsdx,
                    const T * __restrict__ dtdx,
                    const T * __restrict__ drdy,
                    const T * __restrict__ dsdy,
                    const T * __restrict__ dtdy,
                    const T * __restrict__ drdz,
                    const T * __restrict__ dsdz,
                    const T * __restrict__ dtdz,
                    const T * __restrict__ w3) {

  const int e = blockIdx.x;
  const int j = threadIdx.y;
  const int i = threadIdx.x;
  const int ij = i + j * LX;
  const int ele = e*LX*LX*LX;

  #pragma unroll
  for (int k = 0; k < LX; ++k) {
    const int ijk = ij + k*LX*LX;
    const T W3 = w3[ijk];

    cr[ijk + ele] = W3 * (drdx[ijk + ele] * cx[ijk + ele]
                          + drdy[ijk + ele] * cy[ijk + ele]
                          + drdz[ijk + ele] * cz[ijk + ele]);

    cs[ijk + ele] = W3 * (dsdx[ijk + ele] * cx[ijk + ele]
                          + dsdy[ijk + ele] * cy[ijk + ele]
                          + dsdz[ijk + ele] * cz[ijk + ele]);

    ct[ijk + ele] = W3 * (dtdx[ijk + ele] * cx[ijk + ele]
                          + dtdy[ijk + ele] * cy[ijk + ele]
                          + dtdz[ijk + ele] * cz[ijk + ele]);
  }
}


#endif // __MATH_SET_CONVECT_RST_KERNEL_H__
