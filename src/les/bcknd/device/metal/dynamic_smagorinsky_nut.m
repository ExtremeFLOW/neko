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

/**
 * Metal host-side dispatch for the dynamic Smagorinsky model.
 *
 * @note Apple GPUs do not support FP64. This backend operates in FP32.
 */

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include <device/device_config.h>
#include <device/metal/kernel_utils.h>

/** Fortran wrapper for the Metal strain rate magnitude kernel */
void metal_s_abs_compute(void *s_abs, void *s11, void *s22, void *s33,
                         void *s12, void *s13, void *s23, int *n) {
  if (*n < 1) return;

  const int n_r = *n;

  neko_metal_dispatch_1d(neko_metal_pipeline(@"s_abs_compute_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) s_abs offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) s11 offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) s22 offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>) s33 offset:0 atIndex:3];
      [enc setBuffer:(__bridge id<MTLBuffer>) s12 offset:0 atIndex:4];
      [enc setBuffer:(__bridge id<MTLBuffer>) s13 offset:0 atIndex:5];
      [enc setBuffer:(__bridge id<MTLBuffer>) s23 offset:0 atIndex:6];
      [enc setBytes:&n_r length:sizeof(int) atIndex:7];
    }, (NSUInteger) n_r);
}

/** Fortran wrapper for the Metal Leonard stress kernel, part 1 */
void metal_lij_compute_part1(void *l11, void *l22, void *l33,
                             void *l12, void *l13, void *l23,
                             void *u, void *v, void *w,
                             void *fu, void *fv, void *fw,
                             void *fuu, void *fvv, void *fww,
                             void *fuv, void *fuw, void *fvw, int *n) {
  if (*n < 1) return;

  const int n_r = *n;

  neko_metal_dispatch_1d(neko_metal_pipeline(@"lij_compute_part1_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) l11 offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) l22 offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) l33 offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>) l12 offset:0 atIndex:3];
      [enc setBuffer:(__bridge id<MTLBuffer>) l13 offset:0 atIndex:4];
      [enc setBuffer:(__bridge id<MTLBuffer>) l23 offset:0 atIndex:5];
      [enc setBuffer:(__bridge id<MTLBuffer>) u offset:0 atIndex:6];
      [enc setBuffer:(__bridge id<MTLBuffer>) v offset:0 atIndex:7];
      [enc setBuffer:(__bridge id<MTLBuffer>) w offset:0 atIndex:8];
      [enc setBuffer:(__bridge id<MTLBuffer>) fu offset:0 atIndex:9];
      [enc setBuffer:(__bridge id<MTLBuffer>) fv offset:0 atIndex:10];
      [enc setBuffer:(__bridge id<MTLBuffer>) fw offset:0 atIndex:11];
      [enc setBuffer:(__bridge id<MTLBuffer>) fuu offset:0 atIndex:12];
      [enc setBuffer:(__bridge id<MTLBuffer>) fvv offset:0 atIndex:13];
      [enc setBuffer:(__bridge id<MTLBuffer>) fww offset:0 atIndex:14];
      [enc setBuffer:(__bridge id<MTLBuffer>) fuv offset:0 atIndex:15];
      [enc setBuffer:(__bridge id<MTLBuffer>) fuw offset:0 atIndex:16];
      [enc setBuffer:(__bridge id<MTLBuffer>) fvw offset:0 atIndex:17];
      [enc setBytes:&n_r length:sizeof(int) atIndex:18];
    }, (NSUInteger) n_r);
}

/** Fortran wrapper for the Metal Leonard stress kernel, part 2 */
void metal_lij_compute_part2(void *l11, void *l22, void *l33,
                             void *l12, void *l13, void *l23,
                             void *fuu, void *fvv, void *fww,
                             void *fuv, void *fuw, void *fvw, int *n) {
  if (*n < 1) return;

  const int n_r = *n;

  neko_metal_dispatch_1d(neko_metal_pipeline(@"lij_compute_part2_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) l11 offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) l22 offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) l33 offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>) l12 offset:0 atIndex:3];
      [enc setBuffer:(__bridge id<MTLBuffer>) l13 offset:0 atIndex:4];
      [enc setBuffer:(__bridge id<MTLBuffer>) l23 offset:0 atIndex:5];
      [enc setBuffer:(__bridge id<MTLBuffer>) fuu offset:0 atIndex:6];
      [enc setBuffer:(__bridge id<MTLBuffer>) fvv offset:0 atIndex:7];
      [enc setBuffer:(__bridge id<MTLBuffer>) fww offset:0 atIndex:8];
      [enc setBuffer:(__bridge id<MTLBuffer>) fuv offset:0 atIndex:9];
      [enc setBuffer:(__bridge id<MTLBuffer>) fuw offset:0 atIndex:10];
      [enc setBuffer:(__bridge id<MTLBuffer>) fvw offset:0 atIndex:11];
      [enc setBytes:&n_r length:sizeof(int) atIndex:12];
    }, (NSUInteger) n_r);
}

/** Fortran wrapper for the Metal model tensor kernel, part 1 */
void metal_mij_compute_part1(void *m11, void *m22, void *m33,
                             void *m12, void *m13, void *m23,
                             void *s_abs, void *s11, void *s22, void *s33,
                             void *s12, void *s13, void *s23,
                             void *fs_abs, void *fs11, void *fs22, void *fs33,
                             void *fs12, void *fs13, void *fs23,
                             void *fsabss11, void *fsabss22, void *fsabss33,
                             void *fsabss12, void *fsabss13, void *fsabss23,
                             real *delta_ratio2, int *n) {
  if (*n < 1) return;

  const real delta_ratio2_r = *delta_ratio2;
  const int n_r = *n;

  neko_metal_dispatch_1d(neko_metal_pipeline(@"mij_compute_part1_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) m11 offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) m22 offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) m33 offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>) m12 offset:0 atIndex:3];
      [enc setBuffer:(__bridge id<MTLBuffer>) m13 offset:0 atIndex:4];
      [enc setBuffer:(__bridge id<MTLBuffer>) m23 offset:0 atIndex:5];
      [enc setBuffer:(__bridge id<MTLBuffer>) s_abs offset:0 atIndex:6];
      [enc setBuffer:(__bridge id<MTLBuffer>) s11 offset:0 atIndex:7];
      [enc setBuffer:(__bridge id<MTLBuffer>) s22 offset:0 atIndex:8];
      [enc setBuffer:(__bridge id<MTLBuffer>) s33 offset:0 atIndex:9];
      [enc setBuffer:(__bridge id<MTLBuffer>) s12 offset:0 atIndex:10];
      [enc setBuffer:(__bridge id<MTLBuffer>) s13 offset:0 atIndex:11];
      [enc setBuffer:(__bridge id<MTLBuffer>) s23 offset:0 atIndex:12];
      [enc setBuffer:(__bridge id<MTLBuffer>) fs_abs offset:0 atIndex:13];
      [enc setBuffer:(__bridge id<MTLBuffer>) fs11 offset:0 atIndex:14];
      [enc setBuffer:(__bridge id<MTLBuffer>) fs22 offset:0 atIndex:15];
      [enc setBuffer:(__bridge id<MTLBuffer>) fs33 offset:0 atIndex:16];
      [enc setBuffer:(__bridge id<MTLBuffer>) fs12 offset:0 atIndex:17];
      [enc setBuffer:(__bridge id<MTLBuffer>) fs13 offset:0 atIndex:18];
      [enc setBuffer:(__bridge id<MTLBuffer>) fs23 offset:0 atIndex:19];
      [enc setBuffer:(__bridge id<MTLBuffer>) fsabss11 offset:0 atIndex:20];
      [enc setBuffer:(__bridge id<MTLBuffer>) fsabss22 offset:0 atIndex:21];
      [enc setBuffer:(__bridge id<MTLBuffer>) fsabss33 offset:0 atIndex:22];
      [enc setBuffer:(__bridge id<MTLBuffer>) fsabss12 offset:0 atIndex:23];
      [enc setBuffer:(__bridge id<MTLBuffer>) fsabss13 offset:0 atIndex:24];
      [enc setBuffer:(__bridge id<MTLBuffer>) fsabss23 offset:0 atIndex:25];
      [enc setBytes:&delta_ratio2_r length:sizeof(real) atIndex:26];
      [enc setBytes:&n_r length:sizeof(int) atIndex:27];
    }, (NSUInteger) n_r);
}

/** Fortran wrapper for the Metal model tensor and eddy viscosity kernel */
void metal_mij_nut_compute_part2(void *m11, void *m22, void *m33,
                                 void *m12, void *m13, void *m23,
                                 void *l11, void *l22, void *l33,
                                 void *l12, void *l13, void *l23,
                                 void *fsabss11, void *fsabss22,
                                 void *fsabss33, void *fsabss12,
                                 void *fsabss13, void *fsabss23,
                                 void *num, void *den, void *c_dyn,
                                 void *delta, void *s_abs, void *nut,
                                 real *alpha, int *n) {
  if (*n < 1) return;

  const real alpha_r = *alpha;
  const int n_r = *n;

  neko_metal_dispatch_1d(neko_metal_pipeline(@"mij_nut_compute_part2_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) m11 offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) m22 offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) m33 offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>) m12 offset:0 atIndex:3];
      [enc setBuffer:(__bridge id<MTLBuffer>) m13 offset:0 atIndex:4];
      [enc setBuffer:(__bridge id<MTLBuffer>) m23 offset:0 atIndex:5];
      [enc setBuffer:(__bridge id<MTLBuffer>) l11 offset:0 atIndex:6];
      [enc setBuffer:(__bridge id<MTLBuffer>) l22 offset:0 atIndex:7];
      [enc setBuffer:(__bridge id<MTLBuffer>) l33 offset:0 atIndex:8];
      [enc setBuffer:(__bridge id<MTLBuffer>) l12 offset:0 atIndex:9];
      [enc setBuffer:(__bridge id<MTLBuffer>) l13 offset:0 atIndex:10];
      [enc setBuffer:(__bridge id<MTLBuffer>) l23 offset:0 atIndex:11];
      [enc setBuffer:(__bridge id<MTLBuffer>) fsabss11 offset:0 atIndex:12];
      [enc setBuffer:(__bridge id<MTLBuffer>) fsabss22 offset:0 atIndex:13];
      [enc setBuffer:(__bridge id<MTLBuffer>) fsabss33 offset:0 atIndex:14];
      [enc setBuffer:(__bridge id<MTLBuffer>) fsabss12 offset:0 atIndex:15];
      [enc setBuffer:(__bridge id<MTLBuffer>) fsabss13 offset:0 atIndex:16];
      [enc setBuffer:(__bridge id<MTLBuffer>) fsabss23 offset:0 atIndex:17];
      [enc setBuffer:(__bridge id<MTLBuffer>) num offset:0 atIndex:18];
      [enc setBuffer:(__bridge id<MTLBuffer>) den offset:0 atIndex:19];
      [enc setBuffer:(__bridge id<MTLBuffer>) c_dyn offset:0 atIndex:20];
      [enc setBuffer:(__bridge id<MTLBuffer>) delta offset:0 atIndex:21];
      [enc setBuffer:(__bridge id<MTLBuffer>) s_abs offset:0 atIndex:22];
      [enc setBuffer:(__bridge id<MTLBuffer>) nut offset:0 atIndex:23];
      [enc setBytes:&alpha_r length:sizeof(real) atIndex:24];
      [enc setBytes:&n_r length:sizeof(int) atIndex:25];
    }, (NSUInteger) n_r);
}
