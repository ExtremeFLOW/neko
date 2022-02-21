#ifndef __PRGM_LIB_H
#define __PRGM_LIB_H

/**
 * OpenCL program library (see prgm_lib.F90)
 */

/** Device math kernels */
extern void *math_program;

/** Device mathops kernels */
extern void *mathops_program;

/** Device Dirichlet kernels */
extern void *dirichlet_program;

/** Device Inflow kernels */
extern void *inflow_program;

/** Device No-slip wall kernels */
extern void *no_slip_wall_program;

/** Device Symmetry kernels */
extern void *symmetry_program;

/** Device Facet normal kernels */
extern void *facet_normal_program;

/** Device Blasius profile kernel */
extern void *blasius_program;

/** Device Derivative kernels */
extern void *dudxyz_program;

/** Device \f$ D^T X \f$ kernels */
extern void *cdtp_program;

/** Device convective kernels */
extern void *conv1_program;

/** Device Velocity gradient kernels */
extern void *opgrad_program;

/** Device Gather-Scatter kernels */
extern void *gs_program;

/** Device Ax helm kernels */
extern void *ax_helm_program;

/** Device jacobi kernels */
extern void *jacobi_program;

/** Device abbdf kernels */
extern void *abbdf_program;

/** Device pnpn residual kernels */
extern void *pnpn_res_program;


#endif
