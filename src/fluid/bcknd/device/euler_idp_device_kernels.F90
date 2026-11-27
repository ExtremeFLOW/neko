! Copyright (c) 2026, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
!   * Redistributions of source code must retain the above copyright notice,
!     this list of conditions and the following disclaimer.
!   * Redistributions in binary form must reproduce the above copyright
!     notice, this list of conditions and the following disclaimer in the
!     documentation and/or other materials provided with the distribution.
!   * Neither the name of the authors nor the names of its contributors may be
!     used to endorse or promote products derived from this software without
!     specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES ARE DISCLAIMED.
!
!> GPU launch interfaces used by the Euler IDP device backend.
module euler_idp_device_kernels
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int
  use num_types, only : c_rp
  implicit none
  private

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  public :: gpu_euler_idp_primitives, gpu_euler_idp_graph_viscosity
  public :: gpu_euler_idp_graph_rate, gpu_euler_idp_nodal_wave
  public :: gpu_euler_idp_low_residual, gpu_euler_idp_scale_residual
  public :: gpu_euler_idp_reconstruct, gpu_euler_idp_floor_timestep
  public :: gpu_euler_idp_low_update, gpu_euler_idp_bounds_init
  public :: gpu_euler_idp_bounds_edges, gpu_euler_idp_relax_edges
  public :: gpu_euler_idp_relax_finalize, gpu_euler_idp_blend
  public :: gpu_euler_idp_entropy_relax_edges
  public :: gpu_euler_idp_entropy_relax_finalize
  public :: gpu_euler_idp_limiter, gpu_euler_idp_incidence
  public :: gpu_euler_idp_correction_update, gpu_euler_idp_validate
  public :: gpu_euler_idp_full_diagnostics
  public :: gpu_euler_idp_validation_summary
  public :: gpu_euler_idp_observation_summary
  public :: gpu_euler_idp_limiter_status
  public :: gpu_euler_idp_update_uvw
  public :: gpu_euler_idp_update_momentum_pressure
  public :: gpu_euler_idp_update_energy

  interface
     subroutine gpu_euler_idp_primitives(rho, mx, my, mz, energy, u, v, w, &
          p, sound, internal, status, gamma, floor, n) &
          bind(c, name = 'gpu_euler_idp_primitives')
       import c_ptr, c_int, c_rp
       type(c_ptr), value :: rho, mx, my, mz, energy, u, v, w, p
       type(c_ptr), value :: sound, internal, status
       real(c_rp) :: gamma, floor
       integer(c_int) :: n
     end subroutine gpu_euler_idp_primitives

     subroutine gpu_euler_idp_graph_viscosity(rho, mx, my, mz, energy, &
          user_speed, has_user, left, right, coefficient, edge_viscosity, &
          viscosity_sum, edge_speed, gamma, n, n_edges) &
          bind(c, name = 'gpu_euler_idp_graph_viscosity')
       import c_ptr, c_int, c_rp
       type(c_ptr), value :: rho, mx, my, mz, energy, user_speed
       type(c_ptr), value :: left, right, coefficient, edge_viscosity
       type(c_ptr), value :: viscosity_sum, edge_speed
       integer(c_int) :: has_user, n, n_edges
       real(c_rp) :: gamma
     end subroutine gpu_euler_idp_graph_viscosity

     subroutine gpu_euler_idp_graph_rate(viscosity_sum, mass, rate, n) &
          bind(c, name = 'gpu_euler_idp_graph_rate')
       import c_ptr, c_int
       type(c_ptr), value :: viscosity_sum, mass, rate
       integer(c_int) :: n
     end subroutine gpu_euler_idp_graph_rate

     subroutine gpu_euler_idp_nodal_wave(u, v, w, sound, wave, n) &
          bind(c, name = 'gpu_euler_idp_nodal_wave')
       import c_ptr, c_int
       type(c_ptr), value :: u, v, w, sound, wave
       integer(c_int) :: n
     end subroutine gpu_euler_idp_nodal_wave

     subroutine gpu_euler_idp_low_residual(rho, mx, my, mz, energy, &
          diagonal, left, right, coefficient, edge_viscosity, r0, r1, r2, &
          r3, r4, gamma, periodic, n, n_edges) &
          bind(c, name = 'gpu_euler_idp_low_residual')
       import c_ptr, c_int, c_rp
       type(c_ptr), value :: rho, mx, my, mz, energy, diagonal
       type(c_ptr), value :: left, right, coefficient, edge_viscosity
       type(c_ptr), value :: r0, r1, r2, r3, r4
       real(c_rp) :: gamma
       integer(c_int) :: periodic, n, n_edges
     end subroutine gpu_euler_idp_low_residual

     subroutine gpu_euler_idp_scale_residual(r0, r1, r2, r3, r4, binv, n) &
          bind(c, name = 'gpu_euler_idp_scale_residual')
       import c_ptr, c_int
       type(c_ptr), value :: r0, r1, r2, r3, r4, binv
       integer(c_int) :: n
     end subroutine gpu_euler_idp_scale_residual

     subroutine gpu_euler_idp_reconstruct(rho, mx, my, mz, energy, &
          diagonal, left, right, coefficient, residual, correction, error, &
          b, jacinv, dr, ds, dt, derivative, weight, gamma, lx, ly, lz, &
          nelv, edges_per_element, affine, n, n_edges) &
          bind(c, name = 'gpu_euler_idp_reconstruct')
       import c_ptr, c_int, c_rp
       type(c_ptr), value :: rho, mx, my, mz, energy, diagonal
       type(c_ptr), value :: left, right, coefficient, correction, error
       type(c_ptr), value :: b, jacinv
       type(c_ptr), intent(in) :: residual(5), dr(3), ds(3), dt(3)
       type(c_ptr), intent(in) :: derivative(3), weight(3)
       real(c_rp) :: gamma
       integer(c_int) :: lx, ly, lz, nelv, edges_per_element
       integer(c_int) :: affine, n, n_edges
     end subroutine gpu_euler_idp_reconstruct

     subroutine gpu_euler_idp_floor_timestep(rho, mx, my, mz, energy, r0, &
          r1, r2, r3, r4, limit, floor, upper, n) &
          bind(c, name = 'gpu_euler_idp_floor_timestep')
       import c_ptr, c_int, c_rp
       type(c_ptr), value :: rho, mx, my, mz, energy, r0, r1, r2, r3, r4
       type(c_ptr), value :: limit
       real(c_rp) :: floor, upper
       integer(c_int) :: n
     end subroutine gpu_euler_idp_floor_timestep

     subroutine gpu_euler_idp_low_update(rho, mx, my, mz, energy, r0, r1, &
          r2, r3, r4, dt, scalar_mode, n) &
          bind(c, name = 'gpu_euler_idp_low_update')
       import c_ptr, c_int, c_rp
       type(c_ptr), value :: rho, mx, my, mz, energy, r0, r1, r2, r3, r4
       real(c_rp) :: dt
       integer(c_int) :: scalar_mode, n
     end subroutine gpu_euler_idp_low_update

     subroutine gpu_euler_idp_bounds_init(rho, mx, my, mz, energy, lower, &
          upper, entropy, work, gamma, use_entropy, n) &
          bind(c, name = 'gpu_euler_idp_bounds_init')
       import c_ptr, c_int, c_rp
       type(c_ptr), value :: rho, mx, my, mz, energy
       type(c_ptr), value :: lower, upper, entropy, work
       real(c_rp) :: gamma
       integer(c_int) :: use_entropy, n
     end subroutine gpu_euler_idp_bounds_init

     subroutine gpu_euler_idp_bounds_edges(rho, mx, my, mz, left, right, &
          coefficient, edge_viscosity, second_weight_left, &
          second_weight_right, lower, upper, entropy_stage, entropy, work, &
          relax, use_entropy, n_edges) &
          bind(c, name = 'gpu_euler_idp_bounds_edges')
       import c_ptr, c_int
       type(c_ptr), value :: rho, mx, my, mz, left, right
       type(c_ptr), value :: coefficient, edge_viscosity
       type(c_ptr), value :: second_weight_left, second_weight_right
       type(c_ptr), value :: lower, upper, entropy_stage, entropy, work
       integer(c_int) :: relax, use_entropy, n_edges
     end subroutine gpu_euler_idp_bounds_edges

     subroutine gpu_euler_idp_relax_edges(left, right, direction, degree0, &
          degree1, degree2, first, second, n, n_edges) &
          bind(c, name = 'gpu_euler_idp_relax_edges')
       import c_ptr, c_int
       type(c_ptr), value :: left, right, direction, degree0, degree1, degree2
       type(c_ptr), value :: first, second
       integer(c_int) :: n, n_edges
     end subroutine gpu_euler_idp_relax_edges

     subroutine gpu_euler_idp_relax_finalize(lower, upper, second, factor, &
          nodal_mass, volume, dimensions, n) &
          bind(c, name = 'gpu_euler_idp_relax_finalize')
       import c_ptr, c_int, c_rp
       type(c_ptr), value :: lower, upper, second
       real(c_rp) :: factor, nodal_mass, volume
       integer(c_int) :: dimensions, n
     end subroutine gpu_euler_idp_relax_finalize

     subroutine gpu_euler_idp_entropy_relax_edges(rho, mx, my, mz, energy, &
          left, right, midpoint_max, gamma, n_edges) &
          bind(c, name = 'gpu_euler_idp_entropy_relax_edges')
       import c_ptr, c_int, c_rp
       type(c_ptr), value :: rho, mx, my, mz, energy, left, right
       type(c_ptr), value :: midpoint_max
       real(c_rp) :: gamma
       integer(c_int) :: n_edges
     end subroutine gpu_euler_idp_entropy_relax_edges

     subroutine gpu_euler_idp_entropy_relax_finalize(lower, midpoint_max, &
          factor, nodal_mass, volume, dimensions, n) &
          bind(c, name = 'gpu_euler_idp_entropy_relax_finalize')
       import c_ptr, c_int, c_rp
       type(c_ptr), value :: lower, midpoint_max
       real(c_rp) :: factor, nodal_mass, volume
       integer(c_int) :: dimensions, n
     end subroutine gpu_euler_idp_entropy_relax_finalize

     subroutine gpu_euler_idp_blend(rho, mx, my, mz, energy, left, right, &
          edge_viscosity, entropy_fraction, correction, dt, has_entropy, &
          low_only, scalar_mode, n_edges) bind(c, name = 'gpu_euler_idp_blend')
       import c_ptr, c_int, c_rp
       type(c_ptr), value :: rho, mx, my, mz, energy, left, right
       type(c_ptr), value :: edge_viscosity, entropy_fraction, correction
       real(c_rp) :: dt
       integer(c_int) :: has_entropy, low_only, scalar_mode, n_edges
     end subroutine gpu_euler_idp_blend

     subroutine gpu_euler_idp_limiter(q0, q1, q2, q3, q4, lower, upper, &
          entropy, mass, degree0, degree1, degree2, left, right, direction, &
          correction, edge_limit, limited, density_flag, energy_flag, &
          entropy_flag, gamma, floor, enforce_energy, enforce_entropy, &
          check_base, diagnostics_level, dimensions, n_edges) bind(c, name = &
          'gpu_euler_idp_limiter')
       import c_ptr, c_int, c_rp
       type(c_ptr), value :: q0, q1, q2, q3, q4, lower, upper, entropy, mass
       type(c_ptr), value :: degree0, degree1, degree2, left, right, direction
       type(c_ptr), value :: correction, edge_limit, limited
       type(c_ptr), value :: density_flag, energy_flag, entropy_flag
       real(c_rp) :: gamma, floor
       integer(c_int) :: enforce_energy, enforce_entropy, check_base
       integer(c_int) :: diagnostics_level, dimensions, n_edges
     end subroutine gpu_euler_idp_limiter

     subroutine gpu_euler_idp_incidence(left, right, correction, q0, q1, &
          q2, q3, q4, r0, r1, r2, r3, r4, n, n_edges) &
          bind(c, name = 'gpu_euler_idp_incidence')
       import c_ptr, c_int
       type(c_ptr), value :: left, right, correction, q0, q1, q2, q3, q4
       type(c_ptr), value :: r0, r1, r2, r3, r4
       integer(c_int) :: n, n_edges
     end subroutine gpu_euler_idp_incidence

     subroutine gpu_euler_idp_correction_update(q0, q1, q2, q3, q4, r0, &
          r1, r2, r3, r4, mass, n) &
          bind(c, name = 'gpu_euler_idp_correction_update')
       import c_ptr, c_int
       type(c_ptr), value :: q0, q1, q2, q3, q4, r0, r1, r2, r3, r4, mass
       integer(c_int) :: n
     end subroutine gpu_euler_idp_correction_update

     subroutine gpu_euler_idp_validate(q0, q1, q2, q3, q4, lower, upper, &
          entropy_lower, lower_violation, upper_violation, entropy_violation, &
          entropy_value, gamma, use_entropy, n) &
          bind(c, name = 'gpu_euler_idp_validate')
       import c_ptr, c_int, c_rp
       type(c_ptr), value :: q0, q1, q2, q3, q4, lower, upper, entropy_lower
       type(c_ptr), value :: lower_violation, upper_violation
       type(c_ptr), value :: entropy_violation, entropy_value
       real(c_rp) :: gamma
       integer(c_int) :: use_entropy, n
     end subroutine gpu_euler_idp_validate

     subroutine gpu_euler_idp_full_diagnostics(rho, mx, my, mz, energy, &
          gamma, entropy_fraction, has_entropy, edge_limit, limited, &
          density_flag, energy_flag, entropy_flag, correction, r0, r1, r2, &
          r3, r4, directional_error, summary, partial, n, n_edges) &
          bind(c, name = 'gpu_euler_idp_full_diagnostics')
       import c_ptr, c_int, c_rp
       type(c_ptr), value :: rho, mx, my, mz, energy, entropy_fraction
       type(c_ptr), value :: edge_limit, limited, density_flag, energy_flag
       type(c_ptr), value :: entropy_flag, correction
       type(c_ptr), value :: r0, r1, r2, r3, r4
       type(c_ptr), value :: directional_error, summary, partial
       real(c_rp) :: gamma
       integer(c_int) :: has_entropy, n, n_edges
     end subroutine gpu_euler_idp_full_diagnostics

     subroutine gpu_euler_idp_validation_summary(q0, q1, q2, q3, q4, &
          lower, upper, entropy_lower, summary, gamma, use_entropy, n) &
          bind(c, name = 'gpu_euler_idp_validation_summary')
       import c_ptr, c_int, c_rp
       type(c_ptr), value :: q0, q1, q2, q3, q4, lower, upper
       type(c_ptr), value :: entropy_lower, summary
       real(c_rp) :: gamma
       integer(c_int) :: use_entropy, n
     end subroutine gpu_euler_idp_validation_summary

     subroutine gpu_euler_idp_observation_summary(rho, internal, pressure, &
          u, v, w, sound, summary, n) bind(c, name = &
          'gpu_euler_idp_observation_summary')
       import c_ptr, c_int
       type(c_ptr), value :: rho, internal, pressure, u, v, w, sound, summary
       integer(c_int) :: n
     end subroutine gpu_euler_idp_observation_summary

     subroutine gpu_euler_idp_limiter_status(edge_limit, summary, n_edges) &
          bind(c, name = 'gpu_euler_idp_limiter_status')
       import c_ptr, c_int
       type(c_ptr), value :: edge_limit, summary
       integer(c_int) :: n_edges
     end subroutine gpu_euler_idp_limiter_status

     subroutine gpu_euler_idp_update_uvw(u, v, w, mx, my, mz, rho, n) &
          bind(c, name = 'gpu_euler_idp_update_uvw')
       import c_ptr, c_int
       type(c_ptr), value :: u, v, w, mx, my, mz, rho
       integer(c_int) :: n
     end subroutine gpu_euler_idp_update_uvw

     subroutine gpu_euler_idp_update_momentum_pressure(mx, my, mz, p, &
          kinetic, u, v, w, energy, rho, gamma, n) &
          bind(c, name = 'gpu_euler_idp_update_momentum_pressure')
       import c_ptr, c_int, c_rp
       type(c_ptr), value :: mx, my, mz, p, kinetic, u, v, w, energy, rho
       real(c_rp) :: gamma
       integer(c_int) :: n
     end subroutine gpu_euler_idp_update_momentum_pressure

     subroutine gpu_euler_idp_update_energy(energy, p, kinetic, gamma, &
          floor, n) bind(c, name = 'gpu_euler_idp_update_energy')
       import c_ptr, c_int, c_rp
       type(c_ptr), value :: energy, p, kinetic
       real(c_rp) :: gamma, floor
       integer(c_int) :: n
     end subroutine gpu_euler_idp_update_energy
  end interface
#endif

end module euler_idp_device_kernels
