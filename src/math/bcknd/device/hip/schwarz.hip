#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "schwarz_kernel.h"

extern "C" {

  /** 
   * Fortran wrapper for device extrude
   */
  void hip_schwarz_extrude(void *arr1, int * l1, real * f1,
                            void *arr2, int * l2, real * f2,
                            int * nx, int * nel) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);
    hipLaunchKernelGGL(HIP_KERNEL_NAME(schwarz_extrude_kernel<real>),
                       nblcks, nthrds, 0, 0,
                       (real *) arr1,* l1, * f1, 
                       (real *) arr2, *l2, *f2, *nx);
    HIP_CHECK(hipGetLastError());

  } 

  void hip_schwarz_toext3d(void *a, void *b,int * nx, int * nel){
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);
                                \
    hipLaunchKernelGGL(HIP_KERNEL_NAME(schwarz_toext3d_kernel<real>),
                                       nblcks, nthrds, 0, 0,
                                       (real *) a,(real *) b, * nx);
    HIP_CHECK(hipGetLastError());
  } 

  void hip_schwarz_toreg3d(void *b, void *a,int * nx, int * nel){
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);
    hipLaunchKernelGGL(HIP_KERNEL_NAME(schwarz_toreg3d_kernel<real>),
                                       nblcks, nthrds, 0, 0,
                                       (real *) b,(real *) a, * nx);
    HIP_CHECK(hipGetLastError());
  } 

}
