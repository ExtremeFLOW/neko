#ifndef __NEKO_H

#include <stdint.h>

/**
 * Floating point precision in the Neko API
 * @note Set to the same C/C++ equivalent type as @a rp
 */
typedef double neko_real;

/* Initialise Neko */
void neko_init();

/* Finalize Neko */
void neko_finalize();

/* Display job information */
void neko_job_info();

/* Initialise a Neko device layer */
void neko_device_init();

/* Destroy a Neko device layer */
void neko_device_finalize();

/* Initialise a Neko field registry */
void neko_field_registry_init();

/* Destroy a Neko field registry */
void neko_field_registry_free();

/* Allocate memory for a Neko case */
void neko_case_allocate(int **case_iptr);

/* Initialise a Neko case */
void neko_case_init(const char **case_json, int case_len, int **case_iptr);

/* Destroy a Neko case */
void neko_case_free(int **case_iptr);

/* Retrive the current time of a case */
double neko_case_time(int **case_iptr);

/* Retrive the end time of a case */
double neko_case_end_time(int **case_iptr);

/* Retrive the current time-step of a case */
int neko_case_tstep(int **case_iptr);

/* Solve a Neko case */
void neko_solve(int **case_iptr);

/* Compute a time-step for a neko case */
void neko_step(int **case_iptr);

/* Retrive a pointer to a flow field */
neko_real *neko_field(char *field_name);

/* Retrive the order of a field */
int neko_field_order(char *field_name);

/* Retrive the number of elements in a field */
int neko_field_nelements(char *field_name);

/* Retrive the total number of degrees of freedom of a field */
int neko_field_size(char *field_name);

/* Retrive the dofmap associated with a field */
void neko_field_dofmap(char *field_name, int64_t **dof, neko_real **x,
                       neko_real **y, neko_real **z);

/* Retrive the dofmapassociated with a case's fluid solver */
void neko_case_fluid_dofmap(int **case_iptr, int64_t **dof, neko_real **x,
                            neko_real **y, neko_real **z);

/* Retrive the space associated with a field */
void neko_field_space(char *field_name, int *lx, neko_real **zg,
                      neko_real **dr, neko_real **ds, neko_real **dt,
                      neko_real **wx, neko_real **wy, neko_real **wz,
                      neko_real **dx, neko_real **dy, neko_real **dz);

/* Retrive the space associated with a field */
void neko_case_fluid_space(int **case_iptr, int *lx, neko_real **zg,
                           neko_real **dr, neko_real **ds, neko_real **dt,
                           neko_real **wx, neko_real **wy, neko_real **wz,
                           neko_real **dx, neko_real **dy, neko_real **dz);

/* Retrive the coefficient associated with a case's fluid solver */
void neko_case_fluid_coef(int **case_iptr, neko_real **G11, neko_real **G22,
                          neko_real **G33, neko_real **G12, neko_real **G13,
                          neko_real **G23, neko_real **mult, neko_real **dxdr,
                          neko_real **dydr, neko_real **dzdr, neko_real **dxds,
                          neko_real **dyds, neko_real **dzds, neko_real **dxdt,
                          neko_real **dydt, neko_real **dzdt, neko_real **drdx,
                          neko_real **drdy, neko_real **drdz, neko_real **dsdx,
                          neko_real **dsdy, neko_real **dsdz, neko_real **dtdx,
                          neko_real **dtdy, neko_real **dtdz, neko_real **jac,
                          neko_real **B, neko_real **area, neko_real **nx,
                          neko_real **ny, neko_real **nz);
#endif
