#include "schwarz_kernel.h"
#include <device/device_config.h>
#include <device/cuda/check.h>
#include <stdio.h>

extern "C" {

  /** 
   * Fortran wrapper for device extrude
   */
  void cuda_schwarz_extrude(void *arr1, int * l1, real * f1,
                            void *arr2, int * l2, real * f2,
                            int * nx, int * nel) {
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);

    schwarz_extrude_kernel<real>
    <<<nblcks, nthrds>>>((real *) arr1,* l1, * f1, 
                         (real *) arr2, *l2, *f2, *nx);  
    CUDA_CHECK(cudaGetLastError());

  } 

  void cuda_schwarz_toext3d(void *a, void *b,int * nx, int * nel){
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);
                                \
    schwarz_toext3d_kernel<real>
    <<<nblcks, nthrds>>>((real *) a,(real *) b, * nx);  
    CUDA_CHECK(cudaGetLastError());
  } 

  void cuda_schwarz_toreg3d(void *b, void *a,int * nx, int * nel){
    
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((*nel), 1, 1);
                                \
    schwarz_toreg3d_kernel<real>
    <<<nblcks, nthrds>>>((real *) b,(real *) a, * nx);  
    CUDA_CHECK(cudaGetLastError());
  } 

}
