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
 * Metal host-side dispatch for Reichardt's wall model.
 *
 * @note Apple GPUs do not support FP64. This backend operates in FP32.
 */

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include <device/device_config.h>
#include <device/metal/kernel_utils.h>

/**
 * Fortran wrapper for the Metal Reichardt wall model kernel.
 * @param u_d The sampled x-velocity.
 * @param v_d The sampled y-velocity.
 * @param w_d The sampled z-velocity.
 * @param n_x_d The x-component of the wall normals.
 * @param n_y_d The y-component of the wall normals.
 * @param n_z_d The z-component of the wall normals.
 * @param nu_d The kinematic viscosity at wall points.
 * @param rho_w_d The density at wall points.
 * @param h_d The wall-model sampling distances.
 * @param tau_x_d The x-component of the wall shear stress.
 * @param tau_y_d The y-component of the wall shear stress.
 * @param tau_z_d The z-component of the wall shear stress.
 * @param n_nodes The number of wall points.
 * @param kappa The von Karman coefficient.
 * @param C The amplitude of the exponential correction.
 * @param B1 The damping length scale, in wall units.
 * @param B2 The decay length scale of the second exponential term.
 * @param tstep The current time-step.
 */
void metal_reichardt_compute(void *u_d, void *v_d, void *w_d,
                             void *n_x_d, void *n_y_d, void *n_z_d,
                             void *nu_d, void *rho_w_d, void *h_d,
                             void *tau_x_d, void *tau_y_d, void *tau_z_d,
                             int *n_nodes, real *kappa, real *C, real *B1,
                             real *B2, int *tstep) {
  if (*n_nodes < 1) return;

  const int n_r = *n_nodes;
  const real kappa_r = *kappa;
  const real C_r = *C;
  const real B1_r = *B1;
  const real B2_r = *B2;
  const int tstep_r = *tstep;

  neko_metal_dispatch_1d(neko_metal_pipeline(@"reichardt_compute_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) u_d offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) v_d offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) w_d offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>) n_x_d offset:0 atIndex:3];
      [enc setBuffer:(__bridge id<MTLBuffer>) n_y_d offset:0 atIndex:4];
      [enc setBuffer:(__bridge id<MTLBuffer>) n_z_d offset:0 atIndex:5];
      [enc setBuffer:(__bridge id<MTLBuffer>) nu_d offset:0 atIndex:6];
      [enc setBuffer:(__bridge id<MTLBuffer>) rho_w_d offset:0 atIndex:7];
      [enc setBuffer:(__bridge id<MTLBuffer>) h_d offset:0 atIndex:8];
      [enc setBuffer:(__bridge id<MTLBuffer>) tau_x_d offset:0 atIndex:9];
      [enc setBuffer:(__bridge id<MTLBuffer>) tau_y_d offset:0 atIndex:10];
      [enc setBuffer:(__bridge id<MTLBuffer>) tau_z_d offset:0 atIndex:11];
      [enc setBytes:&n_r length:sizeof(int) atIndex:12];
      [enc setBytes:&kappa_r length:sizeof(real) atIndex:13];
      [enc setBytes:&C_r length:sizeof(real) atIndex:14];
      [enc setBytes:&B1_r length:sizeof(real) atIndex:15];
      [enc setBytes:&B2_r length:sizeof(real) atIndex:16];
      [enc setBytes:&tstep_r length:sizeof(int) atIndex:17];
    }, (NSUInteger) n_r);
}
