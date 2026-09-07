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
 * Metal host-side dispatch for the Richardson wall model.
 *
 * @note Apple GPUs do not support FP64. This backend operates in FP32.
 */

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include <device/device_config.h>
#include <device/metal/kernel_utils.h>

/**
 * Fortran wrapper for the Metal richardson wall model kernel
 * @note @a tstep is unused; it is kept for interface compatibility
 * with the other wall models.
 */
void metal_richardson_compute(void *u_d, void *v_d, void *w_d, void *temp_d,
                        void *temp_w_d,
                        void *n_x_d, void *n_y_d, void *n_z_d, void *h_d,
                        void *tau_x_d, void *tau_y_d, void *tau_z_d,
                        int *n_nodes, real *kappa, void *mu_w_d,
                        void *rho_w_d, real *g, real *Pr, real *z0,
                        real *z0h_in, int *bc_type, real *bc_value,
                        int *tstep,
                        void *Ri_b_diagn, void *L_ob_diagn, void *utau_diagn,
                        void *magu_diagn, void *ti_diagn, void *ts_diagn,
                        void *q_diagn) {
  if (*n_nodes < 1) return;
  if (*bc_type != 0 && *bc_type != 1) return;

  const int n_r = *n_nodes;
  const real kappa_r = *kappa;
  const real g1 = g[0];
  const real g2 = g[1];
  const real g3 = g[2];
  const real Pr_r = *Pr;
  const real z0_r = *z0;
  const real z0h_in_r = *z0h_in;
  const real bc_value_r = *bc_value;

  NSString *name = (*bc_type == 0) ? @"richardson_compute_neumann_kernel"
                                   : @"richardson_compute_dirichlet_kernel";

  neko_metal_dispatch_1d(neko_metal_pipeline(name),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) u_d offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) v_d offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) w_d offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>) temp_d offset:0 atIndex:3];
      [enc setBuffer:(__bridge id<MTLBuffer>) temp_w_d offset:0 atIndex:4];
      [enc setBuffer:(__bridge id<MTLBuffer>) h_d offset:0 atIndex:5];
      [enc setBuffer:(__bridge id<MTLBuffer>) n_x_d offset:0 atIndex:6];
      [enc setBuffer:(__bridge id<MTLBuffer>) n_y_d offset:0 atIndex:7];
      [enc setBuffer:(__bridge id<MTLBuffer>) n_z_d offset:0 atIndex:8];
      [enc setBuffer:(__bridge id<MTLBuffer>) tau_x_d offset:0 atIndex:9];
      [enc setBuffer:(__bridge id<MTLBuffer>) tau_y_d offset:0 atIndex:10];
      [enc setBuffer:(__bridge id<MTLBuffer>) tau_z_d offset:0 atIndex:11];
      [enc setBytes:&n_r length:sizeof(int) atIndex:12];
      [enc setBytes:&kappa_r length:sizeof(real) atIndex:13];
      [enc setBuffer:(__bridge id<MTLBuffer>) mu_w_d offset:0 atIndex:14];
      [enc setBuffer:(__bridge id<MTLBuffer>) rho_w_d offset:0 atIndex:15];
      [enc setBytes:&g1 length:sizeof(real) atIndex:16];
      [enc setBytes:&g2 length:sizeof(real) atIndex:17];
      [enc setBytes:&g3 length:sizeof(real) atIndex:18];
      [enc setBytes:&Pr_r length:sizeof(real) atIndex:19];
      [enc setBytes:&z0_r length:sizeof(real) atIndex:20];
      [enc setBytes:&z0h_in_r length:sizeof(real) atIndex:21];
      [enc setBytes:&bc_value_r length:sizeof(real) atIndex:22];
      [enc setBuffer:(__bridge id<MTLBuffer>) Ri_b_diagn offset:0 atIndex:23];
      [enc setBuffer:(__bridge id<MTLBuffer>) L_ob_diagn offset:0 atIndex:24];
      [enc setBuffer:(__bridge id<MTLBuffer>) utau_diagn offset:0 atIndex:25];
      [enc setBuffer:(__bridge id<MTLBuffer>) magu_diagn offset:0 atIndex:26];
      [enc setBuffer:(__bridge id<MTLBuffer>) ti_diagn offset:0 atIndex:27];
      [enc setBuffer:(__bridge id<MTLBuffer>) ts_diagn offset:0 atIndex:28];
      [enc setBuffer:(__bridge id<MTLBuffer>) q_diagn offset:0 atIndex:29];
    }, (NSUInteger) n_r);
}
