#include "facet_normal_kernel.h"

extern "C" {

  /** 
   * Fortran wrapper for device facet normal apply surfvec
   */
  void cuda_facet_normal_apply_surfvec(void *msk, void *facet,
				       void *x, void *y, void *z,
				       void *u, void *v, void *w,
				       void *nx, void * ny, void *nz,
				       void *area, int *lx, int *m) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m) + 1024 - 1)/ 1024, 1, 1);

    facet_normal_apply_surfvec_kernel<double>
      <<<nblcks, nthrds>>>((int *) msk, (int *) facet,
			   (double *) x, (double *) y, (double *) z,
			   (double *) u, (double *) v, (double *) v,
			   (double *) nx,(double *) ny, (double *) nz,
			   (double *) area, *lx, *m);
  }
  
}
