#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "ale_kinematics_kernel.h"

extern "C" {
  void add_kinematics_to_mesh_velocity_hip(
      void *wx, void *wy, void *wz,
      void *x_ref, void *y_ref, void *z_ref,
      void *phi, void *x, void *y, void *z,
      void *center, void *vel_trans, void *vel_ang,
      void *pivot, void *rot_mat, int *n)
  {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    real cx = ((real*)center)[0], cy = ((real*)center)[1], cz = ((real*)center)[2];
    real vtx = ((real*)vel_trans)[0], vty = ((real*)vel_trans)[1], vtz = ((real*)vel_trans)[2];
    real vax = ((real*)vel_ang)[0], vay = ((real*)vel_ang)[1], vaz = ((real*)vel_ang)[2];
    real px = ((real*)pivot)[0], py = ((real*)pivot)[1], pz = ((real*)pivot)[2];

    real *r = (real*)rot_mat;
    real r11 = r[0], r21 = r[1], r31 = r[2];
    real r12 = r[3], r22 = r[4], r32 = r[5];
    real r13 = r[6], r23 = r[7], r33 = r[8];

    hipLaunchKernelGGL(
        HIP_KERNEL_NAME(ale_add_kinematics_kernel<real>),
        nblcks, nthrds, 0, (hipStream_t)glb_cmd_queue,
        *n,
        (real*)wx, (real*)wy, (real*)wz,
        (real*)x_ref, (real*)y_ref, (real*)z_ref,
        (real*)phi, (real*)x, (real*)y, (real*)z,
        cx, cy, cz, vtx, vty, vtz, vax, vay, vaz, px, py, pz,
        r11, r12, r13, r21, r22, r23, r31, r32, r33
    );

    HIP_CHECK(hipGetLastError());
  }
}