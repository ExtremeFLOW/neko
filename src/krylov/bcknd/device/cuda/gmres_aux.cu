#include <stdio.h>
#include "gmres_kernel.h"
#include <device/device_config.h>
extern "C" {
real * gmres_bf1;
real * gmres_bfd1;
  
  real cuda_gmres_part2(void *w, void *v, void *h, void * mult, int *j, int *n) {
	
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    
    if (!gmres_bf1){
      gmres_bf1 = (real *) malloc(nb * sizeof(real));
      cudaMalloc(&gmres_bfd1, nb*sizeof(real));
    }
     
    gmres_part2_kernel<real><<<nblcks, nthrds>>>((real *) w, (real **) v,
                                                 (real *) mult,(real *) h,
                                                 gmres_bfd1, *j, *n);

    cudaMemcpy(gmres_bf1, gmres_bfd1, nb * sizeof(real), cudaMemcpyDeviceToHost);

    real res1 = 0.0;
    for (int i = 0; i < nb; i++) {
      res1 += gmres_bf1[i];
    }

    return res1;
  }
}
