#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include "pipecg_kernel.h"
extern "C" {
  
  void hip_cg_update_xp(void *x, void *p, void *u, void *alpha, void *beta, int *p_cur, int *p_space, int *n) {
	
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    
    cg_update_xp_kernel<real><<<nblcks, nthrds>>>((real *) x, (real *) p,(real **) u, (real *) alpha, (real *) beta, *p_cur, *p_space, *n);

  }

  real * buf1;
  real * buf2;
  real * buf3;
  real * buf_d1;
  real * buf_d2;
  real * buf_d3;
  void hip_pipecg_vecops(void *p, void *q, void *r, void *s, void *u1,
                            void *u2, void *w, void *z, 
                            void *ni, void *mi, real *alpha, 
                            real *beta, void *mult, 
                            real *reduction, int *n) {
	
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    
    if (!buf1){
      buf1 = (real *) malloc(nb * sizeof(real));
      buf2 = (real *) malloc(nb * sizeof(real));
      buf3 = (real *) malloc(nb * sizeof(real));
      hipMalloc(&buf_d1, nb*sizeof(real));
      hipMalloc(&buf_d2, nb*sizeof(real));
      hipMalloc(&buf_d3, nb*sizeof(real));
    }
     
    pipecg_vecops_kernel<real><<<nblcks, nthrds>>>((real *) p, (real *) q,
			                    	    (real *) r, (real *) s,
                                        (real *) u1, (real *) u2,
                                        (real *) w, (real *) z,
                                        (real *) ni, (real *) mi, 
                                        *alpha, *beta, (real *)mult, 
                                        buf_d1, buf_d2, buf_d3, *n);

    hipMemcpy(buf1, buf_d1, nb * sizeof(real), hipMemcpyDeviceToHost);
    hipMemcpy(buf2, buf_d2, nb * sizeof(real), hipMemcpyDeviceToHost);
    hipMemcpy(buf3, buf_d3, nb * sizeof(real), hipMemcpyDeviceToHost);

    real res1 = 0.0;
    real res2 = 0.0;
    real res3 = 0.0;
    for (int i = 0; i < nb; i++) {
      res1 += buf1[i];
      res2 += buf2[i];
      res3 += buf3[i];
    }

    reduction[0] = res1;
    reduction[1] = res2;
    reduction[2] = res3;
  }
}
