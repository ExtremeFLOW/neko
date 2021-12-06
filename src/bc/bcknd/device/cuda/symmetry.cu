#include <algorithm>
#include <device/device_config.h>

#include "symmetry_kernel.h"


extern "C" {

  /** 
   * Fortran wrapper for device symmetry apply vector
   */
  void cuda_symmetry_apply_vector(void *xmsk, void *ymsk, void *zmsk,
				 void *x, void *y, void *z,
				 int *m, int *n, int *l) {

    const int max_len = std::max(std::max(*m, *n), *l);
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((max_len) + 1024 - 1)/ 1024, 1, 1);

    symmetry_apply_vector_kernel<real>
      <<<nblcks, nthrds>>>((int *) xmsk, (int *) ymsk, (int *) zmsk,
			   (real *) x, (real *) y, (real *) z, *m, *n, *l);
  }
 
}
