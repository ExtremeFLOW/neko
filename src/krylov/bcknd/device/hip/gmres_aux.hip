#include <stdio.h>
#include <hip/hip_runtime.h>
#include "gmres_kernel.h"
#include <device/device_config.h>
#include <device/hip/check.h>

extern "C" {
  real * gmres_bf1;
  real * gmres_bfd1;
  
  real hip_gmres_part2(void *w, void *v, void *h, void * mult, int *j, int *n) {
	
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    
    if (!gmres_bf1){
      gmres_bf1 = (real *) malloc(nb * sizeof(real));
      HIP_CHECK(hipMalloc(&gmres_bfd1, nb*sizeof(real)));
    }
     
    hipLaunchKernelGGL(HIP_KERNEL_NAME(gmres_part2_kernel<real>),
		       nblcks, nthrds, 0, 0, (real *) w, (real **) v,
		       (real *) mult, (real *) h, gmres_bfd1, *j, *n);
    HIP_CHECK(hipGetLastError());

    HIP_CHECK(hipMemcpy(gmres_bf1, gmres_bfd1, nb * sizeof(real),
			hipMemcpyDeviceToHost));

    real res1 = 0.0;
    for (int i = 0; i < nb; i++) {
      res1 += gmres_bf1[i];
    }

    return res1;
  }
}
