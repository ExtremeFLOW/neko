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
 * Metal host-side dispatch for ALE mesh kinematics.
 *
 * @note Apple GPUs do not support FP64. This backend operates in FP32.
 */

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include <device/device_config.h>
#include <device/metal/kernel_utils.h>

/**
 * Rigid body kinematics, passed by value.
 * @note The layout must match kinematics_params_t in the device kernels
 * and the bind(c) type in ale_routines_device.F90.
 */
typedef struct {
  real cx, cy, cz;
  real vtx, vty, vtz;
  real vax, vay, vaz;
  real px, py, pz;
  real r11, r12, r13;
  real r21, r22, r23;
  real r31, r32, r33;
} kinematics_params_t;

/** Fortran wrapper for the Metal mesh velocity kinematics kernel */
void add_kinematics_to_mesh_velocity_metal(void *wx, void *wy, void *wz,
                                           void *x_ref, void *y_ref,
                                           void *z_ref, void *phi,
                                           void *x, void *y, void *z,
                                           kinematics_params_t kin_params,
                                           int n) {
  if (n < 1) return;

  neko_metal_dispatch_1d(neko_metal_pipeline(@"ale_add_kinematics_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBytes:&n length:sizeof(int) atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) wx offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) wy offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>) wz offset:0 atIndex:3];
      [enc setBuffer:(__bridge id<MTLBuffer>) x_ref offset:0 atIndex:4];
      [enc setBuffer:(__bridge id<MTLBuffer>) y_ref offset:0 atIndex:5];
      [enc setBuffer:(__bridge id<MTLBuffer>) z_ref offset:0 atIndex:6];
      [enc setBuffer:(__bridge id<MTLBuffer>) phi offset:0 atIndex:7];
      [enc setBuffer:(__bridge id<MTLBuffer>) x offset:0 atIndex:8];
      [enc setBuffer:(__bridge id<MTLBuffer>) y offset:0 atIndex:9];
      [enc setBuffer:(__bridge id<MTLBuffer>) z offset:0 atIndex:10];
      [enc setBytes:&kin_params length:sizeof(kinematics_params_t)
            atIndex:11];
    }, (NSUInteger) n);
}

/** Fortran wrapper for the Metal cheap distance kernel */
void compute_cheap_dist_metal(void *d_d, void *x_d, void *y_d, void *z_d,
                              int lx, int ly, int lz, int nel,
                              int local_iters, void *nchange_d) {
  if (nel < 1) return;

  neko_metal_dispatch_1d(neko_metal_pipeline(@"compute_cheap_dist_kernel"),
    ^(id<MTLComputeCommandEncoder> enc) {
      [enc setBuffer:(__bridge id<MTLBuffer>) d_d offset:0 atIndex:0];
      [enc setBuffer:(__bridge id<MTLBuffer>) x_d offset:0 atIndex:1];
      [enc setBuffer:(__bridge id<MTLBuffer>) y_d offset:0 atIndex:2];
      [enc setBuffer:(__bridge id<MTLBuffer>) z_d offset:0 atIndex:3];
      [enc setBytes:&lx length:sizeof(int) atIndex:4];
      [enc setBytes:&ly length:sizeof(int) atIndex:5];
      [enc setBytes:&lz length:sizeof(int) atIndex:6];
      [enc setBytes:&nel length:sizeof(int) atIndex:7];
      [enc setBytes:&local_iters length:sizeof(int) atIndex:8];
      [enc setBuffer:(__bridge id<MTLBuffer>) nchange_d offset:0 atIndex:9];
    }, (NSUInteger) nel);
}
