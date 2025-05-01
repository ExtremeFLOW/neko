/*
 Copyright (c) 2025, The Neko Authors
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

   * Redistributions in binary form must reproduce the above
     copyright notice, this list of conditions and/or other materials
     provided with the distribution.

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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "wall_model_kernel.h"

extern "C" {
  void hip_wall_model_compute_mag_field(void *tau_x_d, void *tau_y_d, void *tau_z_d,
                                         void *tau_field_d, void *msk_d, int *m) {
    // Define hip grid and block dimensions
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m) + 1024 - 1) / 1024, 1, 1);
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;

    // Launch the hip kernel
    wall_model_compute_mag_field<real>
    <<<nblcks, nthrds, 0, stream>>>((real *) tau_x_d, (real *) tau_y_d, (real *) tau_z_d,
                                    (real *) tau_field_d, (int *) msk_d, *m);

    // Check for hip errors
    HIP_CHECK(hipGetLastError());
  }
}