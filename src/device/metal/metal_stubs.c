// Auto-generated Metal stub implementations
// 59 functions that need real Metal implementations
//
// Each stub calls abort() with an error message.
// Replace with real implementations as Metal kernels are developed.
//
#include <stdio.h>
#include <stdlib.h>

#define METAL_STUB(name) do { \
    fprintf(stderr, "Metal kernel not yet implemented: %s\n", name); \
    abort(); \
} while(0)

// From math/bcknd/device/ax_helm_full_device.F90
void metal_ax_helm_stress_vector(void *au_d, void *av_d, void *aw_d, void *u_d, void *v_d, void *w_d, void *dx_d, void *dy_d, void *dz_d, void *dxt_d, void *dyt_d, void *dzt_d, void *h1_d, void *drdx_d, void *drdy_d, void *drdz_d, void *dsdx_d, void *dsdy_d, void *dsdz_d, void *dtdx_d, void *dtdy_d, void *dtdz_d, void *jacinv_d, void *weight3_d, int *nelv, int *lx) {
    METAL_STUB("metal_ax_helm_stress_vector");
}

// From math/bcknd/device/ax_helm_full_device.F90
void metal_ax_helm_stress_vector_part2(void *au_d, void *av_d, void *aw_d, void *u_d, void *v_d, void *w_d, void *h2_d, void *B_d, int *n) {
    METAL_STUB("metal_ax_helm_stress_vector_part2");
}

// From math/bcknd/device/opr_device.F90
void metal_cdtp(void *dtx_d, void *x_d, void *dr_d, void *ds_d, void *dt_d, void *dxt_d, void *dyt_d, void *dzt_d, void *w3_d, int *nel, int *lx) {
    METAL_STUB("metal_cdtp");
}

// From krylov/bcknd/device/pipecg_device.F90
void metal_cg_update_xp(void *x_d, void *p_d, void *u_d_d, void *alpha, void *beta, void *p_cur, int *p_space, int *n) {
    METAL_STUB("metal_cg_update_xp");
}

// From krylov/bcknd/device/cheby_device.F90
void metal_cheby_part1(void *d_d, void *x_d, double *inv_tha, int *n, void *strm) {
    METAL_STUB("metal_cheby_part1");
}

// From krylov/bcknd/device/cheby_device.F90
void metal_cheby_part2(void *d_d, void *w_d, void *x_d, double *tmp1, double *tmp2, int *n, void *strm) {
    METAL_STUB("metal_cheby_part2");
}

// From math/bcknd/device/opr_device.F90
void metal_conv1(void *du_d, void *u_d, void *vx_d, void *vy_d, void *vz_d, void *dx_d, void *dy_d, void *dz_d, void *drdx_d, void *dsdx_d, void *dtdx_d, void *drdy_d, void *dsdy_d, void *dtdy_d, void *drdz_d, void *dsdz_d, void *dtdz_d, void *jacinv_d, int *nel, int *gdim, int *lx) {
    METAL_STUB("metal_conv1");
}

// From math/bcknd/device/opr_device.F90
void metal_convect_scalar(void *du_d, void *u_d, void *cr_d, void *cs_d, void *ct_d, void *dx_d, void *dy_d, void *dz_d, int *nel, int *lx) {
    METAL_STUB("metal_convect_scalar");
}





// From math/bcknd/device/fdm_device.F90
void metal_fdm_do_fast(void *e_d, void *r_d, void *s_d, void *d_d, int *nl, int *nelv, void *stream) {
    METAL_STUB("metal_fdm_do_fast");
}

// From sem/bcknd/device/device_local_interpolation.F90
void metal_find_rst_legendre(void *rst, void *pt_x, void *pt_y, void *pt_z, void *x_hat, void *y_hat, void *z_hat, void *resx, void *resy, void *resz, void *lx, void *el_ids, int *n_pt, double *tol, void *conv_pts) {
    METAL_STUB("metal_find_rst_legendre");
}

// From krylov/bcknd/device/fusedcg_cpld_device.F90
void metal_fusedcg_cpld_part1(void *a1_d, void *a2_d, void *a3_d, void *b1_d, void *b2_d, void *b3_d, void *tmp_d, int *n) {
    METAL_STUB("metal_fusedcg_cpld_part1");
}

// From krylov/bcknd/device/fusedcg_cpld_device.F90
double metal_fusedcg_cpld_part2(void *a1_d, void *a2_d, void *a3_d, void *b_d, void *c1_d, void *c2_d, void *c3_d, void *alpha_d, double *alpha, int *p_cur, int *n) {
    METAL_STUB("metal_fusedcg_cpld_part2");
    return 0.0;
}

// From krylov/bcknd/device/fusedcg_cpld_device.F90
void metal_fusedcg_cpld_update_p(void *p1_d, void *p2_d, void *p3_d, void *z1_d, void *z2_d, void *z3_d, void *po1_d, void *po2_d, void *po3_d, double *beta, int *n) {
    METAL_STUB("metal_fusedcg_cpld_update_p");
}

// From krylov/bcknd/device/fusedcg_cpld_device.F90
void metal_fusedcg_cpld_update_x(void *x1_d, void *x2_d, void *x3_d, void *p1_d, void *p2_d, void *p3_d, void *alpha, int *p_cur, int *n) {
    METAL_STUB("metal_fusedcg_cpld_update_x");
}

// From krylov/bcknd/device/fusedcg_device.F90
double metal_fusedcg_part2(void *a_d, void *b_d, void *c_d, void *alpha_d, double *alpha, void *p_cur, int *n) {
    METAL_STUB("metal_fusedcg_part2");
    return 0.0;
}

// From krylov/bcknd/device/fusedcg_device.F90
void metal_fusedcg_update_p(void *p_d, void *z_d, void *po_d, double *beta, int *n) {
    METAL_STUB("metal_fusedcg_update_p");
}

// From krylov/bcknd/device/fusedcg_device.F90
void metal_fusedcg_update_x(void *x_d, void *p_d, void *alpha, int *p_cur, int *n) {
    METAL_STUB("metal_fusedcg_update_x");
}

// From krylov/bcknd/device/gmres_device.F90
double metal_gmres_part2(void *w_d, void *v_d_d, void *h_d, void *mult_d, int *j, int *n) {
    METAL_STUB("metal_gmres_part2");
    return 0.0;
}



// From math/bcknd/device/opr_device.F90
void metal_lambda2(void *lambda2_d, void *u_d, void *v_d, void *w_d, void *dx_d, void *dy_d, void *dz_d, void *drdx_d, void *dsdx_d, void *dtdx_d, void *drdy_d, void *dsdy_d, void *dtdy_d, void *drdz_d, void *dsdz_d, void *dtdz_d, void *jacinv_d, int *nel, int *lx) {
    METAL_STUB("metal_lambda2");
}




// From math/bcknd/device/opr_device.F90
void metal_opgrad(void *ux_d, void *uy_d, void *uz_d, void *u_d, void *dx_d, void *dy_d, void *dz_d, void *drdx_d, void *dsdx_d, void *dtdx_d, void *drdy_d, void *dsdy_d, void *dtdy_d, void *drdz_d, void *dsdz_d, void *dtdz_d, void *w3_d, int *nel, int *lx) {
    METAL_STUB("metal_opgrad");
}

// From krylov/bcknd/device/pipecg_device.F90
void metal_pipecg_vecops(void *p_d, void *q_d, void *r_d, void *s_d, void *u_d1, void *u_d2, void *w_d, void *z_d, void *ni_d, void *mi_d, double *alpha, double *beta, void *mult_d, void *reduction, int *n) {
    METAL_STUB("metal_pipecg_vecops");
}

// From common/bcknd/device/device_projection.F90
void metal_project_on(void *a_d, void *b_d, void *x_d_d, void *b_d_d, void *mult_d, void *x_d, int *j, int *n) {
    METAL_STUB("metal_project_on");
}

// From common/bcknd/device/device_projection.F90
void metal_project_ortho(void *a_d, void *b_d, void *x_d_d, void *b_d_d, void *w_d, void *xm_d, int *j, int *n, double *nrm) {
    METAL_STUB("metal_project_ortho");
}

// From math/bcknd/device/device_schwarz.F90
void metal_schwarz_extrude(void *arr1_d, int *l1, double *f1, void *arr2_d, int *l2, double *f2, void *nx, int *nelv, void *stream) {
    METAL_STUB("metal_schwarz_extrude");
}

// From math/bcknd/device/device_schwarz.F90
void metal_schwarz_toext3d(void *a_d, void *b_d, int *nx, int *nelv, void *stream) {
    METAL_STUB("metal_schwarz_toext3d");
}

// From math/bcknd/device/device_schwarz.F90
void metal_schwarz_toreg3d(void *b_d, void *a_d, int *nx, int *nelv, void *stream) {
    METAL_STUB("metal_schwarz_toreg3d");
}

// From math/bcknd/device/opr_device.F90
void metal_set_convect_rst(void *cr_d, void *cs_d, void *ct_d, void *cx_d, void *cy_d, void *cz_d, void *drdx_d, void *dsdx_d, void *dtdx_d, void *drdy_d, void *dsdy_d, void *dtdy_d, void *drdz_d, void *dsdz_d, void *dtdz_d, void *w3_d, int *nel, int *lx) {
    METAL_STUB("metal_set_convect_rst");
}





// From fluid/bcknd/device/pnpn_res_device.F90
void pnpn_prs_res_part1_metal(void *ta1_d, void *ta2_d, void *ta3_d, void *wa1_d, void *wa2_d, void *wa3_d, void *f_u_d, void *f_v_d, void *f_w_d, void *B_d, void *h1_d, double *mu, double *rho, int *n) {
    METAL_STUB("pnpn_prs_res_part1_metal");
}

// From fluid/stress_formulation/bcknd/device/pnpn_res_stress_device.F90
void pnpn_prs_res_part2_metal(void *p_res_d, void *wa1_d, void *wa2_d, void *wa3_d, int *n) {
    METAL_STUB("pnpn_prs_res_part2_metal");
}

// From fluid/bcknd/device/pnpn_res_device.F90
void pnpn_prs_res_part3_metal(void *p_res_d, void *ta1_d, void *ta2_d, void *ta3_d, double *dtbd, void *n) {
    METAL_STUB("pnpn_prs_res_part3_metal");
}

// From fluid/stress_formulation/bcknd/device/pnpn_res_stress_device.F90
void pnpn_prs_stress_res_part1_metal(void *ta1_d, void *ta2_d, void *ta3_d, void *wa1_d, void *wa2_d, void *wa3_d, void *s11_d, void *s22_d, void *s33_d, void *s12_d, void *s13_d, void *s23_d, void *f_u_d, void *f_v_d, void *f_w_d, void *B_d, void *h1_d, void *rho_d, int *n) {
    METAL_STUB("pnpn_prs_stress_res_part1_metal");
}

// From fluid/stress_formulation/bcknd/device/pnpn_res_stress_device.F90
void pnpn_prs_stress_res_part3_metal(void *p_res_d, void *ta1_d, void *ta2_d, void *ta3_d, void *wa1_d, void *wa2_d, void *wa3_d, double *dtbd, int *n) {
    METAL_STUB("pnpn_prs_stress_res_part3_metal");
}

// From fluid/stress_formulation/bcknd/device/pnpn_res_stress_device.F90
void pnpn_vel_res_update_metal(void *u_res_d, void *v_res_d, void *w_res_d, void *ta1_d, void *ta2_d, void *ta3_d, void *f_u_d, void *f_v_d, void *f_w_d, int *n) {
    METAL_STUB("pnpn_vel_res_update_metal");
}

// From fluid/bcknd/device/fluid_abbdf_device.F90
void rhs_maker_bdf_metal(void *ulag1_d, void *ulag2_d, void *vlag1_d, void *vlag2_d, void *wlag1_d, void *wlag2_d, void *bfx_d, void *bfy_d, void *bfz_d, void *u_d, void *v_d, void *w_d, void *B_d, void *rho, void *dt, void *bd2, void *bd3, void *bd4, int *nbd, int *n) {
    METAL_STUB("rhs_maker_bdf_metal");
}

// From fluid/bcknd/device/fluid_abbdf_device.F90
void rhs_maker_ext_metal(void *abx1_d, void *aby1_d, void *abz1_d, void *abx2_d, void *aby2_d, void *abz2_d, void *bfx_d, void *bfy_d, void *bfz_d, void *rho, double *ab1, double *ab2, double *ab3, int *n) {
    METAL_STUB("rhs_maker_ext_metal");
}

// From common/bcknd/device/rhs_maker_device.F90
void rhs_maker_oifs_metal(void *phi_x_d, void *phi_y_d, void *phi_z_d, void *bf_x_d, void *bf_y_d, void *bf_z_d, void *rho, void *dt, int *n) {
    METAL_STUB("rhs_maker_oifs_metal");
}

// From fluid/bcknd/device/fluid_abbdf_device.F90
void rhs_maker_sumab_metal(void *u_d, void *v_d, void *w_d, void *uu_d, void *vv_d, void *ww_d, void *uulag1, void *uulag2, void *vvlag1, void *vvlag2, void *wwlag1, void *wwlag2, double *ab1, double *ab2, double *ab3, int *nab, int *n) {
    METAL_STUB("rhs_maker_sumab_metal");
}

// From scalar/bcknd/device/scalar_residual_device.F90
void scalar_residual_update_metal(void *s_res_d, void *f_s_d, int *n) {
    METAL_STUB("scalar_residual_update_metal");
}

// From common/bcknd/device/rhs_maker_device.F90
void scalar_rhs_maker_bdf_metal(void *s_lag_d, void *s_laglag_d, void *fs_d, void *s_d, void *B_d, void *rho, void *dt, void *bd2, void *bd3, void *bd4, int *nbd, int *n) {
    METAL_STUB("scalar_rhs_maker_bdf_metal");
}

// From common/bcknd/device/rhs_maker_device.F90
void scalar_rhs_maker_ext_metal(void *fs_lag_d, void *fs_laglag_d, void *fs_d, double *rho, void *ext1, double *ext2, double *ext3, int *n) {
    METAL_STUB("scalar_rhs_maker_ext_metal");
}

// From common/bcknd/device/rhs_maker_device.F90
void scalar_rhs_maker_oifs_metal(void *phi_s_d, void *bf_s_d, void *rho, void *dt, int *n) {
    METAL_STUB("scalar_rhs_maker_oifs_metal");
}


// From math/bcknd/device/opr_device.F90
void metal_rotate_cyc(void *vx_d, void *vy_d, void *vz_d, void *x_d, void *y_d, void *z_d, void *cyc_msk_d, void *R11_d, void *R12_d, int *ncyc, int *idir) {
    METAL_STUB("metal_rotate_cyc");
}

// From math/bcknd/device/device_math.F90
void metal_masked_scatter_copy_aligned(void *a_d, void *b_d, void *mask_d, int *n, int *n_mask, void *strm) {
    METAL_STUB("metal_masked_scatter_copy_aligned");
}

// From math/bcknd/device/device_math.F90
void metal_face_masked_gather_copy(void *a_d, void *b_d, void *mask_d, void *facet_d, int *n1, int *n2, int *lx, int *ly, int *lz, int *n_mask, void *strm) {
    METAL_STUB("metal_face_masked_gather_copy");
}

// From math/bcknd/device/device_math.F90
void metal_cwrap(void *a_d, float *min_val, float *max_val, int *n, void *strm) {
    METAL_STUB("metal_cwrap");
}

// From filter/bcknd/device/mappings_device.F90
void metal_smooth_step(void *x, double *edge0, double *edge1, int *n) {
    METAL_STUB("metal_smooth_step");
}

// From filter/bcknd/device/mappings_device.F90
void metal_step_function(void *x, double *edge, double *left, double *right, int *n) {
    METAL_STUB("metal_step_function");
}

// From filter/bcknd/device/mappings_device.F90
void metal_permeability(void *x, double *k_0, double *k_1, double *q, int *n) {
    METAL_STUB("metal_permeability");
}
