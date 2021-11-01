/**
 * OpenCL program library (see prgm_lib.F90)
 */

/** Device math kernels */
extern void *math_program;

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
