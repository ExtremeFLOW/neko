/*
 Copyright (c) 2026, The Neko Authors
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

#ifndef __LPT_PERIODIC_PARAMS_H__
#define __LPT_PERIODIC_PARAMS_H__

#include <device/device_config.h>

/**
 * Periodic direction data for the translational wrapping kernels, packed
 * into one kernel argument so that the OpenCL and Metal backends stay
 * within their kernel argument limits.
 *
 * Each triple is the (x, y, z) part of periodic direction @a j, so
 * @a dir[3*j], @a dir[3*j+1] and @a dir[3*j+2] are direction @a j.
 *
 * @note The layout is mirrored by lpt_periodic_params in
 * opencl/lpt_periodic_bc_kernel.cl and metal/lpt_periodic_bc_kernel.metal.
 */
typedef struct {
  real dir[9];
  real pmin[3];
  real pmax[3];
  real shift[9];
  real len[3];
} lpt_periodic_params;

/** Pack the flat Fortran arguments into an lpt_periodic_params */
static inline void lpt_periodic_params_pack(
    lpt_periodic_params *prm,
    const real *dir_x1, const real *dir_y1, const real *dir_z1,
    const real *dir_x2, const real *dir_y2, const real *dir_z2,
    const real *dir_x3, const real *dir_y3, const real *dir_z3,
    const real *min1, const real *min2, const real *min3,
    const real *max1, const real *max2, const real *max3,
    const real *shift_x1, const real *shift_y1, const real *shift_z1,
    const real *shift_x2, const real *shift_y2, const real *shift_z2,
    const real *shift_x3, const real *shift_y3, const real *shift_z3,
    const real *len1, const real *len2, const real *len3) {
  prm->dir[0] = *dir_x1;
  prm->dir[1] = *dir_y1;
  prm->dir[2] = *dir_z1;
  prm->dir[3] = *dir_x2;
  prm->dir[4] = *dir_y2;
  prm->dir[5] = *dir_z2;
  prm->dir[6] = *dir_x3;
  prm->dir[7] = *dir_y3;
  prm->dir[8] = *dir_z3;

  prm->pmin[0] = *min1;
  prm->pmin[1] = *min2;
  prm->pmin[2] = *min3;

  prm->pmax[0] = *max1;
  prm->pmax[1] = *max2;
  prm->pmax[2] = *max3;

  prm->shift[0] = *shift_x1;
  prm->shift[1] = *shift_y1;
  prm->shift[2] = *shift_z1;
  prm->shift[3] = *shift_x2;
  prm->shift[4] = *shift_y2;
  prm->shift[5] = *shift_z2;
  prm->shift[6] = *shift_x3;
  prm->shift[7] = *shift_y3;
  prm->shift[8] = *shift_z3;

  prm->len[0] = *len1;
  prm->len[1] = *len2;
  prm->len[2] = *len3;
}

#endif /* __LPT_PERIODIC_PARAMS_H__ */
