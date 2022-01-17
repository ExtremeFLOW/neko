#include <hip/hip_runtime.h>
#include <device/device_config.h>
#include <device/hip/check.h>
#include "pipecg_kernel.h"

/**
 * @todo Make sure that this gets deleted at some point...
 */
real *buf1 = NULL;
real *buf2 = NULL;
real *buf3 = NULL;
real *buf_d1 = NULL;
real *buf_d2 = NULL;
real *buf_d3 = NULL;

extern "C" {
  
  void hip_cg_update_xp(void *x, void *p, void *u, void *alpha, void *beta,
                        int *p_cur, int *p_space, int *n) {
        
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1) / 1024;
    
    hipLaunchKernelGGL(HIP_KERNEL_NAME(cg_update_xp_kernel<real>),
                       nblcks, nthrds, 0, 0,
                       (real *) x, (real *) p, (real **) u, (real *) alpha,
                       (real *) beta, *p_cur, *p_space, *n);
    HIP_CHECK(hipGetLastError());
  }

  void hip_pipecg_vecops(void *p, void *q, void *r, void *s, void *u1,
                            void *u2, void *w, void *z, 
                            void *ni, void *mi, real *alpha, 
                            real *beta, void *mult, 
                            real *reduction, int *n) {
        
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    
    if (!buf1){
      buf1 = (real *) malloc(nb * sizeof(real));
      buf2 = (real *) malloc(nb * sizeof(real));
      buf3 = (real *) malloc(nb * sizeof(real));
      HIP_CHECK(hipMalloc(&buf_d1, nb*sizeof(real)));
      HIP_CHECK(hipMalloc(&buf_d2, nb*sizeof(real)));
      HIP_CHECK(hipMalloc(&buf_d3, nb*sizeof(real)));
    }
     
    hipLaunchKernelGGL(HIP_KERNEL_NAME(pipecg_vecops_kernel<real>),
                       nblcks, nthrds, 0, 0,
                       (real *) p, (real *) q, (real *) r, (real *) s,
                       (real *) u1, (real *) u2, (real *) w, (real *) z,
                       (real *) ni, (real *) mi, *alpha, *beta,
                       (real *) mult, buf_d1, buf_d2, buf_d3, *n);
    HIP_CHECK(hipGetLastError());

    HIP_CHECK(hipMemcpy(buf1, buf_d1, nb * sizeof(real),
                        hipMemcpyDeviceToHost));
    HIP_CHECK(hipMemcpy(buf2, buf_d2, nb * sizeof(real),
                        hipMemcpyDeviceToHost));
    HIP_CHECK(hipMemcpy(buf3, buf_d3, nb * sizeof(real),
                        hipMemcpyDeviceToHost));

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
