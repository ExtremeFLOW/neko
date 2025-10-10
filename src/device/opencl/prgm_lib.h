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

/** Device zero dirichlet kernels */
extern void *zero_dirichlet_program;

/** Device Symmetry kernels */
extern void *symmetry_program;

/** Device Facet normal kernels */
extern void *facet_normal_program;

/** Device Neumann kernels */
extern void *neumann_program;

/** Device Blasius profile kernel */
extern void *inhom_dirichlet_program;

/** Device Derivative kernels */
extern void *dudxyz_program;

/** Device \f$ D^T X \f$ kernels */
extern void *cdtp_program;

/** Device convective kernels */
extern void *conv1_program;

/** Device convect_scalar kernels */
extern void *convect_scalar_program;

/** Device set_convect_rst kernels */
extern void *set_convect_rst_program;

/** Device CFL kernels */
extern void *cfl_program;

/** Device Velocity gradient kernels */
extern void *opgrad_program;

/** Device Gather-Scatter kernels */
extern void *gs_program;

/** Device Ax helm kernels */
extern void *ax_helm_program;

/** Device Ax helm full kernels */
extern void *ax_helm_full_program;

/** Device jacobi kernels */
extern void *jacobi_program;

/** Device rhs_maker kernels */
extern void *rhs_maker_program;

/** Device pnpn residual kernels */
extern void *pnpn_res_program;

/** Device pnpn residual kernels (stress formulation) */
extern void *pnpn_stress_res_program;

/** Device euler residual kernels */
extern void *euler_res_program;

/** Device compressible ops kernels */
extern void *compressible_ops_compute_max_wave_speed_program;
extern void *compressible_ops_compute_entropy_program;

/** Device fdm kernels */
extern void *fdm_program;

/** Device tensor kernels */
extern void *tensor_program;

/** Device schwarz kernels */
extern void *schwarz_program;

/** Device dong kernels */
extern void *dong_program;

/** Device coef kernels */
extern void *coef_program;

/** Device scalar residual kernels */
extern void *scalar_residual_program;

/** Device lambda2 kernel */
extern void *lambda2_program;

/** Device compute_max_wave_speed kernel */
extern void *compute_max_wave_speed_program;

/** Device filter kernel */
extern void *mapping_program;

/** Device find rst kernel */
extern void *find_rst_legendre_program;

#endif
