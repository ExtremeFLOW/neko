/*
 Copyright (c) 2026, The Neko Authors
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

   * Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
   * Neither the name of the authors nor the names of its contributors may be
     used to endorse or promote products derived from this software without
     specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES ARE DISCLAIMED.
*/

/* Fortran entry points and the ordered sequence of launches are shared by
   all device backends. Buffer arguments carry their element type for CUDA/HIP,
   while OpenCL and Metal bind the corresponding device memory object. */
static inline int idp_blocks(int n) { return n / 256 + (n % 256 != 0); }

void gpu_euler_idp_primitives(void *rho, void *mx, void *my, void *mz, void *energy,
                              void *u, void *v, void *w, void *p, void *sound,
                              void *internal, void *status, real *gamma, real *floor,
                              int *n) {
  IDP_LAUNCH(idp_primitives_kernel, idp_blocks(*n), 256, RBUF(rho), RBUF(mx), RBUF(my),
             RBUF(mz), RBUF(energy), RBUF(u), RBUF(v), RBUF(w), RBUF(p), RBUF(sound),
             RBUF(internal), RBUF(status), RVAL(*gamma), RVAL(*floor), IVAL(*n));
}

void gpu_euler_idp_graph_viscosity(void *rho, void *mx, void *my, void *mz,
                                   void *energy, void *user_speed, int *has_user,
                                   void *left, void *right, void *coef, void *edge_visc,
                                   void *vis_sum, void *edge_speed, real *gamma, int *n,
                                   int *nedge) {
  IDP_ZERO(vis_sum, (*n));
  IDP_LAUNCH(idp_graph_viscosity_kernel, idp_blocks(*nedge), 256, RBUF(rho), RBUF(mx),
             RBUF(my), RBUF(mz), RBUF(energy), RBUF(user_speed), IVAL(*has_user),
             IBUF(left), IBUF(right), RBUF(coef), RBUF(edge_visc), RBUF(vis_sum),
             RBUF(edge_speed), RVAL(*gamma), IVAL(*nedge));
}

void gpu_euler_idp_graph_rate(void *vis, void *mass, void *rate, int *n) {
  IDP_LAUNCH(idp_graph_rate_kernel, idp_blocks(*n), 256, RBUF(vis), RBUF(mass),
             RBUF(rate), IVAL(*n));
}
void gpu_euler_idp_nodal_wave(void *u, void *v, void *w, void *sound, void *wave,
                              int *n) {
  IDP_LAUNCH(idp_nodal_wave_kernel, idp_blocks(*n), 256, RBUF(u), RBUF(v), RBUF(w),
             RBUF(sound), RBUF(wave), IVAL(*n));
}

void gpu_euler_idp_low_residual(void *rho, void *mx, void *my, void *mz, void *energy,
                                void *diag, void *left, void *right, void *coef,
                                void *visc, void *r0, void *r1, void *r2, void *r3,
                                void *r4, real *gamma, int *periodic, int *n,
                                int *nedge) {
  IDP_LAUNCH(idp_low_init_kernel, idp_blocks(*n), 256, RBUF(rho), RBUF(mx), RBUF(my),
             RBUF(mz), RBUF(energy), RBUF(diag), RBUF(r0), RBUF(r1), RBUF(r2), RBUF(r3),
             RBUF(r4), RVAL(*gamma), IVAL(*periodic), IVAL(*n));
  IDP_LAUNCH(idp_low_edge_kernel, idp_blocks(*nedge), 256, RBUF(rho), RBUF(mx),
             RBUF(my), RBUF(mz), RBUF(energy), IBUF(left), IBUF(right), RBUF(coef),
             RBUF(visc), RBUF(r0), RBUF(r1), RBUF(r2), RBUF(r3), RBUF(r4), RVAL(*gamma),
             IVAL(*periodic), IVAL(*nedge));
}

void gpu_euler_idp_scale_residual(void *r0, void *r1, void *r2, void *r3, void *r4,
                                  void *binv, int *n) {
  IDP_LAUNCH(idp_scale_residual_kernel, idp_blocks(*n), 256, RBUF(r0), RBUF(r1),
             RBUF(r2), RBUF(r3), RBUF(r4), RBUF(binv), IVAL(*n));
}

void gpu_euler_idp_reconstruct(
    void *rho, void *mx, void *my, void *mz, void *energy, void *diag, void *left,
    void *right, void *coef, void *const residual[5], void *correction, void *error,
    void *B, void *jacinv, void *const dr[3], void *const ds[3], void *const dt[3],
    void *const derivative[3], void *const weight[3], real *gamma, int *lx, int *ly,
    int *lz, int *nelv, int *edges_per_element, int *affine, int *n, int *nedge) {
  IDP_ZERO(correction, 5 * (*nedge));
  IDP_ZERO(error, 5);
  if (*affine) {
    const int nx = (*nelv) * (*ly) * (*lz), ny = (*nelv) * (*lx) * (*lz),
              nz = (*nelv) * (*lx) * (*ly);
    IDP_LAUNCH(idp_reconstruct_direction_kernel, idp_blocks(nx), 256, IVAL(0),
               RBUF(rho), RBUF(mx), RBUF(my), RBUF(mz), RBUF(energy), RBUF(coef),
               RBUF(B), RBUF(jacinv), RBUF(dr[0]), RBUF(dr[1]), RBUF(dr[2]),
               RBUF(derivative[0]), RBUF(correction), RBUF(error), RVAL(*gamma),
               IVAL(*lx), IVAL(*ly), IVAL(*lz), IVAL(*edges_per_element), IVAL(nx));
    IDP_LAUNCH(idp_reconstruct_direction_kernel, idp_blocks(ny), 256, IVAL(1),
               RBUF(rho), RBUF(mx), RBUF(my), RBUF(mz), RBUF(energy), RBUF(coef),
               RBUF(B), RBUF(jacinv), RBUF(ds[0]), RBUF(ds[1]), RBUF(ds[2]),
               RBUF(derivative[1]), RBUF(correction), RBUF(error), RVAL(*gamma),
               IVAL(*lx), IVAL(*ly), IVAL(*lz), IVAL(*edges_per_element), IVAL(ny));
    if (*lz > 1)
      IDP_LAUNCH(idp_reconstruct_direction_kernel, idp_blocks(nz), 256, IVAL(2),
                 RBUF(rho), RBUF(mx), RBUF(my), RBUF(mz), RBUF(energy), RBUF(coef),
                 RBUF(B), RBUF(jacinv), RBUF(dt[0]), RBUF(dt[1]), RBUF(dt[2]),
                 RBUF(derivative[2]), RBUF(correction), RBUF(error), RVAL(*gamma),
                 IVAL(*lx), IVAL(*ly), IVAL(*lz), IVAL(*edges_per_element), IVAL(nz));
  } else {
    IDP_LAUNCH(idp_correction_residual_init_kernel, idp_blocks(*n), 256, RBUF(rho),
               RBUF(mx), RBUF(my), RBUF(mz), RBUF(energy), RBUF(diag), RBUF(dr[0]),
               RBUF(dr[1]), RBUF(dr[2]), RBUF(ds[0]), RBUF(ds[1]), RBUF(ds[2]),
               RBUF(dt[0]), RBUF(dt[1]), RBUF(dt[2]), RBUF(derivative[0]),
               RBUF(derivative[1]), RBUF(derivative[2]), RBUF(weight[0]),
               RBUF(weight[1]), RBUF(weight[2]), RBUF(residual[0]), RBUF(residual[1]),
               RBUF(residual[2]), RBUF(residual[3]), RBUF(residual[4]), RVAL(*gamma),
               IVAL(*lx), IVAL(*ly), IVAL(*lz), IVAL(*n));
    IDP_LAUNCH(idp_correction_residual_edge_kernel, idp_blocks(*nedge), 256, RBUF(rho),
               RBUF(mx), RBUF(my), RBUF(mz), RBUF(energy), IBUF(left), IBUF(right),
               RBUF(coef), RBUF(residual[0]), RBUF(residual[1]), RBUF(residual[2]),
               RBUF(residual[3]), RBUF(residual[4]), RVAL(*gamma), IVAL(*nedge));
    IDP_LAUNCH(idp_reconstruct_general_kernel, idp_blocks(*nelv), 256,
               RBUF(residual[0]), RBUF(residual[1]), RBUF(residual[2]),
               RBUF(residual[3]), RBUF(residual[4]), RBUF(weight[0]), RBUF(weight[1]),
               RBUF(weight[2]), RBUF(correction), RBUF(error), IVAL(*lx), IVAL(*ly),
               IVAL(*lz), IVAL(*edges_per_element), IVAL(*nelv));
  }
}

void gpu_euler_idp_floor_timestep(void *rho, void *mx, void *my, void *mz, void *energy,
                                  void *r0, void *r1, void *r2, void *r3, void *r4,
                                  void *limit, real *floor, real *upper, int *n) {
  IDP_LAUNCH(idp_floor_timestep_kernel, idp_blocks(*n), 256, RBUF(rho), RBUF(mx),
             RBUF(my), RBUF(mz), RBUF(energy), RBUF(r0), RBUF(r1), RBUF(r2), RBUF(r3),
             RBUF(r4), RBUF(limit), RVAL(*floor), RVAL(*upper), IVAL(*n));
}
void gpu_euler_idp_low_update(void *rho, void *mx, void *my, void *mz, void *energy,
                              void *r0, void *r1, void *r2, void *r3, void *r4,
                              real *dt, int *scalar_mode, int *n) {
  IDP_LAUNCH(idp_low_update_kernel, idp_blocks(*n), 256, RBUF(rho), RBUF(mx), RBUF(my),
             RBUF(mz), RBUF(energy), RBUF(r0), RBUF(r1), RBUF(r2), RBUF(r3), RBUF(r4),
             RVAL(*dt), IVAL(*scalar_mode), IVAL(*n));
}

void gpu_euler_idp_bounds_init(void *rho, void *mx, void *my, void *mz, void *energy,
                               void *lower, void *upper, void *entropy, void *work,
                               real *gamma, int *use_entropy, int *n) {
  IDP_LAUNCH(idp_bounds_init_kernel, idp_blocks(*n), 256, RBUF(rho), RBUF(mx), RBUF(my),
             RBUF(mz), RBUF(energy), RBUF(lower), RBUF(upper), RBUF(entropy),
             RBUF(work), RVAL(*gamma), IVAL(*use_entropy), IVAL(*n));
}
void gpu_euler_idp_bounds_edges(void *rho, void *mx, void *my, void *mz, void *left,
                                void *right, void *coef, void *visc,
                                void *second_weight_left, void *second_weight_right,
                                void *lower, void *upper, void *entropy_stage,
                                void *entropy, void *work, int *relax, int *use_entropy,
                                int *nedge) {
  IDP_LAUNCH(idp_bounds_edge_kernel, idp_blocks(*nedge), 256, RBUF(rho), RBUF(mx),
             RBUF(my), RBUF(mz), IBUF(left), IBUF(right), RBUF(coef), RBUF(visc),
             RBUF(second_weight_left), RBUF(second_weight_right), RBUF(lower),
             RBUF(upper), RBUF(entropy_stage), RBUF(entropy), RBUF(work), IVAL(*relax),
             IVAL(*use_entropy), IVAL(*nedge));
}
void gpu_euler_idp_relax_edges(void *left, void *right, void *direction, void *degree0,
                               void *degree1, void *degree2, void *first, void *second,
                               int *n, int *nedge) {
  IDP_ZERO(second, (*n));
  IDP_LAUNCH(idp_relax_edge_kernel, idp_blocks(*nedge), 256, IBUF(left), IBUF(right),
             IBUF(direction), RBUF(degree0), RBUF(degree1), RBUF(degree2), RBUF(first),
             RBUF(second), IVAL(*nedge));
}
void gpu_euler_idp_relax_finalize(void *lower, void *upper, void *second, real *factor,
                                  real *nodal_mass, real *volume, int *dimensions,
                                  int *n) {
  IDP_LAUNCH(idp_relax_finalize_kernel, idp_blocks(*n), 256, RBUF(lower), RBUF(upper),
             RBUF(second), RVAL(*factor), RVAL(*nodal_mass), RVAL(*volume),
             IVAL(*dimensions), IVAL(*n));
}

void gpu_euler_idp_entropy_relax_edges(void *rho, void *mx, void *my, void *mz,
                                       void *energy, void *left, void *right,
                                       void *midpoint_max, real *gamma, int *nedge) {
  IDP_LAUNCH(idp_entropy_relax_edge_kernel, idp_blocks(*nedge), 256, RBUF(rho),
             RBUF(mx), RBUF(my), RBUF(mz), RBUF(energy), IBUF(left), IBUF(right),
             RBUF(midpoint_max), RVAL(*gamma), IVAL(*nedge));
}

void gpu_euler_idp_entropy_relax_finalize(void *lower, void *midpoint_max, real *factor,
                                          real *nodal_mass, real *volume,
                                          int *dimensions, int *n) {
  IDP_LAUNCH(idp_entropy_relax_finalize_kernel, idp_blocks(*n), 256, RBUF(lower),
             RBUF(midpoint_max), RVAL(*factor), RVAL(*nodal_mass), RVAL(*volume),
             IVAL(*dimensions), IVAL(*n));
}

void gpu_euler_idp_blend(void *rho, void *mx, void *my, void *mz, void *energy,
                         void *left, void *right, void *visc, void *ev,
                         void *correction, real *dt, int *has_ev, int *low_only,
                         int *scalar_mode, int *nedge) {
  IDP_LAUNCH(idp_blend_kernel, idp_blocks(*nedge), 256, RBUF(rho), RBUF(mx), RBUF(my),
             RBUF(mz), RBUF(energy), IBUF(left), IBUF(right), RBUF(visc), RBUF(ev),
             RBUF(correction), RVAL(*dt), IVAL(*has_ev), IVAL(*low_only),
             IVAL(*scalar_mode), IVAL(*nedge));
}
void gpu_euler_idp_limiter(void *q0, void *q1, void *q2, void *q3, void *q4,
                           void *lower, void *upper, void *entropy, void *mass,
                           void *degree0, void *degree1, void *degree2, void *left,
                           void *right, void *direction, void *correction,
                           void *edge_limit, void *limited, void *density_flag,
                           void *energy_flag, void *entropy_flag, real *gamma,
                           real *floor, int *enforce_energy, int *enforce_entropy,
                           int *check_base, int *diagnostics_level, int *dimensions,
                           int *nedge) {
  IDP_LAUNCH(idp_limiter_kernel, idp_blocks(*nedge), 256, RBUF(q0), RBUF(q1), RBUF(q2),
             RBUF(q3), RBUF(q4), RBUF(lower), RBUF(upper), RBUF(entropy), RBUF(mass),
             RBUF(degree0), RBUF(degree1), RBUF(degree2), IBUF(left), IBUF(right),
             IBUF(direction), RBUF(correction), RBUF(edge_limit), RBUF(limited),
             RBUF(density_flag), RBUF(energy_flag), RBUF(entropy_flag), RVAL(*gamma),
             RVAL(*floor), IVAL(*enforce_energy), IVAL(*enforce_entropy),
             IVAL(*check_base), IVAL(*diagnostics_level), IVAL(*dimensions),
             IVAL(*nedge));
}
void gpu_euler_idp_incidence(void *left, void *right, void *correction, void *q0,
                             void *q1, void *q2, void *q3, void *q4, void *r0, void *r1,
                             void *r2, void *r3, void *r4, int *n, int *nedge) {
  void *r[5] = {r0, r1, r2, r3, r4};
  for (int c = 0; c < 5; ++c)
    IDP_ZERO(r[c], (*n));
  IDP_LAUNCH(idp_incidence_kernel, idp_blocks(*nedge), 256, IBUF(left), IBUF(right),
             RBUF(correction), RBUF(r0), RBUF(r1), RBUF(r2), RBUF(r3), RBUF(r4),
             IVAL(*nedge));
}
void gpu_euler_idp_correction_update(void *q0, void *q1, void *q2, void *q3, void *q4,
                                     void *r0, void *r1, void *r2, void *r3, void *r4,
                                     void *mass, int *n) {
  IDP_LAUNCH(idp_correction_update_kernel, idp_blocks(*n), 256, RBUF(q0), RBUF(q1),
             RBUF(q2), RBUF(q3), RBUF(q4), RBUF(r0), RBUF(r1), RBUF(r2), RBUF(r3),
             RBUF(r4), RBUF(mass), IVAL(*n));
}

void gpu_euler_idp_validate(void *q0, void *q1, void *q2, void *q3, void *q4,
                            void *lower, void *upper, void *entropy_lower,
                            void *lower_v, void *upper_v, void *entropy_v,
                            void *entropy_value, real *gamma, int *use_entropy,
                            int *n) {
  IDP_LAUNCH(idp_validation_kernel, idp_blocks(*n), 256, RBUF(q0), RBUF(q1), RBUF(q2),
             RBUF(q3), RBUF(q4), RBUF(lower), RBUF(upper), RBUF(entropy_lower),
             RBUF(lower_v), RBUF(upper_v), RBUF(entropy_v), RBUF(entropy_value),
             RVAL(*gamma), IVAL(*use_entropy), IVAL(*n));
}

void gpu_euler_idp_full_diagnostics(void *rho, void *mx, void *my, void *mz,
                                    void *energy, real *gamma, void *entropy_fraction,
                                    int *has_entropy, void *edge_limit, void *limited,
                                    void *density_flag, void *energy_flag,
                                    void *entropy_flag, void *correction, void *r0,
                                    void *r1, void *r2, void *r3, void *r4,
                                    void *directional_error, void *summary,
                                    void *partial, int *n, int *nedge) {
  const int node_blocks = idp_blocks(*n), edge_blocks = idp_blocks(*nedge);
  IDP_LAUNCH(idp_full_diagnostics_node_partial_kernel, node_blocks, 256, RBUF(rho),
             RBUF(mx), RBUF(my), RBUF(mz), RBUF(energy), RVAL(*gamma),
             RBUF(entropy_fraction), IVAL(*has_entropy), RBUF(r0), RBUF(r1), RBUF(r2),
             RBUF(r3), RBUF(r4), RBUF(partial), IVAL(*n));
  IDP_LAUNCH(idp_full_diagnostics_edge_partial_kernel, edge_blocks, 256,
             RBUF(edge_limit), RBUF(limited), RBUF(density_flag), RBUF(energy_flag),
             RBUF(entropy_flag), RBUF(correction), RBUF(partial), IVAL(*nedge),
             IVAL(node_blocks));
  IDP_LAUNCH(idp_full_diagnostics_finalize_kernel, 1, 32, RBUF(partial),
             RBUF(directional_error), RBUF(summary), IVAL(node_blocks), IVAL(edge_blocks));
}

void gpu_euler_idp_validation_summary(void *q0, void *q1, void *q2, void *q3, void *q4,
                                      void *lower, void *upper, void *entropy_lower,
                                      void *summary, real *gamma, int *use_entropy,
                                      int *n) {
  IDP_LAUNCH(idp_diagnostic_summary_init_kernel, 1, 1, RBUF(summary), IVAL(9),
             IVAL(validation_summary));
  IDP_LAUNCH(idp_validation_summary_kernel, idp_blocks(*n), 256, RBUF(q0), RBUF(q1),
             RBUF(q2), RBUF(q3), RBUF(q4), RBUF(lower), RBUF(upper),
             RBUF(entropy_lower), RBUF(summary), RVAL(*gamma), IVAL(*use_entropy),
             IVAL(*n));
}

void gpu_euler_idp_observation_summary(void *rho, void *internal, void *pressure,
                                       void *u, void *v, void *w, void *sound,
                                       void *summary, int *n) {
  IDP_LAUNCH(idp_diagnostic_summary_init_kernel, 1, 1, RBUF(summary), IVAL(4),
             IVAL(observation_summary));
  IDP_LAUNCH(idp_observation_summary_kernel, idp_blocks(*n), 256, RBUF(rho),
             RBUF(internal), RBUF(pressure), RBUF(u), RBUF(v), RBUF(w), RBUF(sound),
             RBUF(summary), IVAL(*n));
}

void gpu_euler_idp_limiter_status(void *edge_limit, void *summary, int *nedge) {
  IDP_LAUNCH(idp_diagnostic_summary_init_kernel, 1, 1, RBUF(summary), IVAL(3),
             IVAL(limiter_status_summary));
  IDP_LAUNCH(idp_limiter_status_kernel, idp_blocks(*nedge), 256, RBUF(edge_limit),
             RBUF(summary), IVAL(*nedge));
}

void gpu_euler_idp_update_uvw(void *u, void *v, void *w, void *mx, void *my, void *mz,
                              void *rho, int *n) {
  IDP_LAUNCH(idp_update_uvw_kernel_idp, idp_blocks(*n), 256, RBUF(u), RBUF(v), RBUF(w),
             RBUF(mx), RBUF(my), RBUF(mz), RBUF(rho), IVAL(*n));
}
void gpu_euler_idp_update_momentum_pressure(void *mx, void *my, void *mz, void *p,
                                            void *kin, void *u, void *v, void *w,
                                            void *energy, void *rho, real *gamma,
                                            int *n) {
  IDP_LAUNCH(idp_update_momentum_pressure_kernel, idp_blocks(*n), 256, RBUF(mx),
             RBUF(my), RBUF(mz), RBUF(p), RBUF(kin), RBUF(u), RBUF(v), RBUF(w),
             RBUF(energy), RBUF(rho), RVAL(*gamma), IVAL(*n));
}
void gpu_euler_idp_update_energy(void *energy, void *p, void *kin, real *gamma,
                                 real *floor, int *n) {
  IDP_LAUNCH(idp_update_energy_kernel, idp_blocks(*n), 256, RBUF(energy), RBUF(p),
             RBUF(kin), RVAL(*gamma), RVAL(*floor), IVAL(*n));
}
