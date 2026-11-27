! Copyright (c) 2026, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above copyright
!     notice, this list of conditions and the following disclaimer in the
!     documentation and/or other materials provided with the distribution.
!
!   * Neither the name of the authors nor the names of its contributors may
!     be used to endorse or promote products derived from this software
!     without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Backend contract and common data for the Euler GLL-IDP solver.
module euler_idp_backend
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use num_types, only : rp, dp
  use field, only : field_t
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t
  use bc_list, only : bc_list_t
  use time_state, only : time_state_t
  implicit none
  private

  integer, public, parameter :: EULER_IDP_NCOMP = 5

  integer, public, parameter :: EULER_IDP_DIAGNOSTICS_OFF = 0
  integer, public, parameter :: EULER_IDP_DIAGNOSTICS_SAFETY = 1
  integer, public, parameter :: EULER_IDP_DIAGNOSTICS_FULL = 2
  public :: euler_idp_stage_time
  public :: euler_idp_maximum_wave_speed, euler_idp_flux_dot_vector
  public :: euler_idp_bar_state, euler_idp_internal_energy
  public :: euler_idp_internal_energy_timestep
  public :: euler_idp_state_is_admissible, euler_idp_specific_entropy
  public :: euler_idp_entropy_tolerance, euler_idp_entropy_is_admissible
  public :: euler_idp_local_entropy_bounds, euler_idp_relax_density_bounds
  public :: euler_idp_relax_entropy_bound
  public :: euler_idp_limit_endpoint, euler_idp_limit_edge

  !> State and graph data immediately around a strong boundary map.
  type, public :: euler_idp_state_observation_t
     real(kind=rp) :: time = 0.0_rp
     real(kind=rp) :: min_density = huge(1.0_rp)
     real(kind=rp) :: min_internal_energy = huge(1.0_rp)
     real(kind=rp) :: min_pressure = huge(1.0_rp)
     real(kind=rp) :: max_nodal_wave_speed = 0.0_rp
     real(kind=rp) :: max_graph_wave_speed = 0.0_rp
     real(kind=rp) :: max_graph_rate = 0.0_rp
  end type euler_idp_state_observation_t

  !> Diagnostics returned by one limited Forward Euler map.
  type, public :: euler_idp_diagnostics_t
     integer :: stage = 0
     real(kind=rp) :: stage_time = 0.0_rp
     real(kind=rp) :: min_density = huge(1.0_rp)
     real(kind=rp) :: min_internal_energy = huge(1.0_rp)
     real(kind=rp) :: min_pressure = huge(1.0_rp)
     real(kind=rp) :: min_specific_entropy = huge(1.0_rp)
     real(kind=rp) :: max_graph_cfl = 0.0_rp
     real(kind=rp) :: max_nodal_wave_speed = 0.0_rp
     real(kind=rp) :: max_graph_wave_speed = 0.0_rp
     real(kind=rp) :: max_graph_rate = 0.0_rp
     real(kind=rp) :: min_convex_weight = 1.0_rp
     real(kind=rp) :: maximum_floor_timestep = huge(1.0_rp)
     logical :: entropy_viscosity_enabled = .false.
     real(kind=rp) :: max_entropy_viscosity = 0.0_rp
     real(kind=rp) :: mean_entropy_viscosity = 0.0_rp
     real(kind=rp) :: entropy_viscosity_conservation(EULER_IDP_NCOMP) = &
          0.0_rp
     real(kind=rp) :: entropy_viscosity_element_compatibility( &
          EULER_IDP_NCOMP) = 0.0_rp
     real(kind=rp) :: high_order_conservation(EULER_IDP_NCOMP) = 0.0_rp
     real(kind=rp) :: low_order_conservation(EULER_IDP_NCOMP) = 0.0_rp
     real(kind=rp) :: limited_conservation(EULER_IDP_NCOMP) = 0.0_rp
     real(kind=rp) :: correction_assembly_error(EULER_IDP_NCOMP) = 0.0_rp
     real(kind=rp) :: correction_global_compatibility(EULER_IDP_NCOMP) = &
          0.0_rp
     real(kind=rp) :: correction_element_compatibility(EULER_IDP_NCOMP) = &
          0.0_rp
     real(kind=rp) :: correction_face_mismatch(EULER_IDP_NCOMP) = 0.0_rp
     real(kind=rp) :: element_compatibility(EULER_IDP_NCOMP) = 0.0_rp
     real(kind=rp) :: graph_compatibility(EULER_IDP_NCOMP) = 0.0_rp
     real(kind=rp) :: reconstruction_residual(EULER_IDP_NCOMP) = 0.0_rp
     real(kind=rp) :: forward_euler_time = 0.0_rp
     real(kind=rp) :: reconstruction_time = 0.0_rp
     real(kind=rp) :: graph_gs_time = 0.0_rp
     real(kind=rp) :: graph_reduction_time = 0.0_rp
     real(kind=rp) :: max_correction_flux = 0.0_rp
     real(kind=rp) :: rms_correction_flux = 0.0_rp
     real(kind=rp) :: min_limiter = 1.0_rp
     real(kind=rp) :: mean_limiter = 1.0_rp
     real(kind=rp) :: max_limiter = 1.0_rp
     real(kind=rp) :: limited_edge_fraction = 0.0_rp
     real(kind=rp) :: limiter_weight_error = 0.0_rp
     real(kind=rp) :: min_density_lower_bound = huge(1.0_rp)
     real(kind=rp) :: max_density_upper_bound = -huge(1.0_rp)
     real(kind=rp) :: max_density_bound_relaxation = 0.0_rp
     real(kind=rp) :: max_density_lower_violation = 0.0_rp
     real(kind=rp) :: max_density_upper_violation = 0.0_rp
     real(kind=rp) :: max_entropy_lower_violation = 0.0_rp
     integer :: density_limited_edges = 0
     integer :: internal_energy_limited_edges = 0
     integer :: entropy_limited_edges = 0
     type(euler_idp_state_observation_t) :: before_boundary
     type(euler_idp_state_observation_t) :: after_boundary
   contains
     procedure, pass(this) :: reset => euler_idp_diagnostics_reset
  end type euler_idp_diagnostics_t

  !> Large-grain operations implemented by an execution backend.
  type, public, abstract :: euler_idp_backend_t
     logical :: initialized = .false.
     integer :: diagnostics_level = EULER_IDP_DIAGNOSTICS_FULL
     real(kind=rp) :: max_graph_rate = 0.0_rp
     real(kind=rp) :: max_graph_wave_speed = 0.0_rp
     real(kind=rp) :: maximum_graph_timestep = huge(1.0_rp)
   contains
     procedure(euler_idp_backend_init_intrf), pass(this), deferred :: init
     procedure(euler_idp_backend_init_graph_intrf), pass(this), deferred :: &
          init_graph
     procedure(euler_idp_backend_free_intrf), pass(this), deferred :: free
     procedure(euler_idp_backend_boundary_intrf), pass(this), deferred :: &
          apply_boundary_conditions
     procedure(euler_idp_backend_prepare_intrf), pass(this), deferred :: &
          prepare_stage
     procedure(euler_idp_backend_forward_euler_intrf), pass(this), deferred :: &
          forward_euler
     procedure(euler_idp_backend_candidate_boundary_intrf), pass(this), &
          deferred :: apply_candidate_boundary
     procedure(euler_idp_backend_candidate_observe_intrf), pass(this), &
          deferred :: observe_candidate
     procedure(euler_idp_backend_candidate_validate_intrf), pass(this), &
          deferred :: validate_candidate
     procedure(euler_idp_backend_save_state_intrf), pass(this), deferred :: &
          save_state
     procedure(euler_idp_backend_combine_stage_intrf), pass(this), deferred :: &
          combine_stage
     procedure(euler_idp_backend_copy_primitives_intrf), pass(this), &
          deferred :: copy_primitives
     procedure, pass(this) :: graph_cfl => euler_idp_backend_graph_cfl
  end type euler_idp_backend_t

  abstract interface
     subroutine euler_idp_backend_init_intrf(this, dof)
       import :: dofmap_t, euler_idp_backend_t
       class(euler_idp_backend_t), intent(inout) :: this
       type(dofmap_t), target, intent(in) :: dof
     end subroutine euler_idp_backend_init_intrf

     subroutine euler_idp_backend_init_graph_intrf(this, coef, gs, &
          relax_density_bounds, relax_entropy_bounds, low_order_only, &
          limit_internal_energy, limit_entropy, &
          density_bound_relaxation_factor, entropy_bound_relaxation_factor, &
          time_order, diagnostics_level, correction_tolerance)
       import :: coef_t, gs_t, rp, euler_idp_backend_t
       class(euler_idp_backend_t), intent(inout) :: this
       type(coef_t), target, intent(in) :: coef
       type(gs_t), intent(inout) :: gs
       logical, intent(in) :: relax_density_bounds, relax_entropy_bounds
       logical, intent(in) :: low_order_only
       logical, intent(in) :: limit_internal_energy, limit_entropy
       real(kind=rp), intent(in) :: density_bound_relaxation_factor
       real(kind=rp), intent(in) :: entropy_bound_relaxation_factor
       real(kind=rp), intent(in) :: correction_tolerance
       integer, intent(in) :: time_order, diagnostics_level
     end subroutine euler_idp_backend_init_graph_intrf

     subroutine euler_idp_backend_free_intrf(this)
       import :: euler_idp_backend_t
       class(euler_idp_backend_t), intent(inout) :: this
     end subroutine euler_idp_backend_free_intrf

     subroutine euler_idp_backend_boundary_intrf(this, rho, m_x, m_y, m_z, &
          energy, density_bcs, velocity_bcs, pressure_bcs, gamma, &
          internal_energy_floor, time, label, refresh)
       import :: bc_list_t, euler_idp_backend_t, field_t, rp, time_state_t
       class(euler_idp_backend_t), intent(inout) :: this
       type(field_t), intent(inout) :: rho, m_x, m_y, m_z, energy
       type(bc_list_t), intent(inout) :: density_bcs, velocity_bcs
       type(bc_list_t), intent(inout) :: pressure_bcs
       real(kind=rp), intent(in) :: gamma, internal_energy_floor
       type(time_state_t), intent(in), optional :: time
       character(len=*), intent(in) :: label
       logical, intent(in), optional :: refresh
     end subroutine euler_idp_backend_boundary_intrf

     subroutine euler_idp_backend_prepare_intrf(this, rho, m_x, m_y, m_z, &
          energy, gs, gamma, internal_energy_floor, primitives_valid, &
          graph_wave_speed)
       import :: euler_idp_backend_t, field_t, gs_t, rp
       class(euler_idp_backend_t), intent(inout) :: this
       type(field_t), intent(in) :: rho, m_x, m_y, m_z, energy
       type(gs_t), intent(inout) :: gs
       real(kind=rp), intent(in) :: gamma, internal_energy_floor
       logical, intent(in) :: primitives_valid
       type(field_t), intent(in), optional :: graph_wave_speed
     end subroutine euler_idp_backend_prepare_intrf

     subroutine euler_idp_backend_forward_euler_intrf(this, rho, m_x, m_y, &
          m_z, energy, coef, gs, gamma, internal_energy_floor, dt, time, &
          diagnostics, stage, &
          entropy_viscosity_fraction, graph_wave_speed)
       import :: coef_t, euler_idp_backend_t
       import :: euler_idp_diagnostics_t, field_t, gs_t, rp, time_state_t
       class(euler_idp_backend_t), intent(inout) :: this
       type(field_t), intent(in) :: rho, m_x, m_y, m_z, energy
       type(coef_t), intent(inout) :: coef
       type(gs_t), intent(inout) :: gs
       real(kind=rp), intent(in) :: gamma, internal_energy_floor, dt
       type(time_state_t), intent(in) :: time
       type(euler_idp_diagnostics_t), intent(inout) :: diagnostics
       integer, intent(in) :: stage
       type(field_t), intent(in), optional :: entropy_viscosity_fraction
       type(field_t), intent(in), optional :: graph_wave_speed
     end subroutine euler_idp_backend_forward_euler_intrf

     subroutine euler_idp_backend_candidate_boundary_intrf(this, density_bcs, &
          velocity_bcs, pressure_bcs, gamma, internal_energy_floor, time, &
          label)
       import :: bc_list_t, euler_idp_backend_t, rp, time_state_t
       class(euler_idp_backend_t), intent(inout) :: this
       type(bc_list_t), intent(inout) :: density_bcs, velocity_bcs
       type(bc_list_t), intent(inout) :: pressure_bcs
       real(kind=rp), intent(in) :: gamma, internal_energy_floor
       type(time_state_t), intent(in) :: time
       character(len=*), intent(in) :: label
     end subroutine euler_idp_backend_candidate_boundary_intrf

     subroutine euler_idp_backend_candidate_observe_intrf(this, gs, gamma, &
          internal_energy_floor, time, label, primitives_valid, observation, &
          graph_wave_speed)
       import :: euler_idp_backend_t, euler_idp_state_observation_t
       import :: field_t, gs_t, rp, time_state_t
       class(euler_idp_backend_t), intent(inout) :: this
       type(gs_t), intent(inout) :: gs
       real(kind=rp), intent(in) :: gamma, internal_energy_floor
       type(time_state_t), intent(in) :: time
       character(len=*), intent(in) :: label
       logical, intent(in) :: primitives_valid
       type(euler_idp_state_observation_t), intent(out) :: observation
       type(field_t), intent(in), optional :: graph_wave_speed
     end subroutine euler_idp_backend_candidate_observe_intrf

     subroutine euler_idp_backend_candidate_validate_intrf(this, gamma, &
          diagnostics, label)
       import :: euler_idp_backend_t, euler_idp_diagnostics_t, rp
       class(euler_idp_backend_t), intent(inout) :: this
       real(kind=rp), intent(in) :: gamma
       type(euler_idp_diagnostics_t), intent(inout) :: diagnostics
       character(len=*), intent(in) :: label
     end subroutine euler_idp_backend_candidate_validate_intrf

     subroutine euler_idp_backend_save_state_intrf(this, rho, m_x, m_y, m_z, &
          energy)
       import :: euler_idp_backend_t, field_t
       class(euler_idp_backend_t), intent(inout) :: this
       type(field_t), intent(in) :: rho, m_x, m_y, m_z, energy
     end subroutine euler_idp_backend_save_state_intrf

     subroutine euler_idp_backend_combine_stage_intrf(this, rho, m_x, m_y, &
          m_z, energy, saved_weight, candidate_weight)
       import :: euler_idp_backend_t, field_t, rp
       class(euler_idp_backend_t), intent(in) :: this
       type(field_t), intent(inout) :: rho, m_x, m_y, m_z, energy
       real(kind=rp), intent(in) :: saved_weight, candidate_weight
     end subroutine euler_idp_backend_combine_stage_intrf

     subroutine euler_idp_backend_copy_primitives_intrf(this, u, v, w, p)
       import :: euler_idp_backend_t, field_t
       class(euler_idp_backend_t), intent(in) :: this
       type(field_t), intent(inout) :: u, v, w, p
     end subroutine euler_idp_backend_copy_primitives_intrf
  end interface

contains

  !> Return the physical time of a Forward Euler evaluation. The incoming
  !! time is Neko's end-of-step time.
  pure real(kind=dp) function euler_idp_stage_time(time, dt, order, stage) &
       result(value)
    real(kind=dp), intent(in) :: time
    real(kind=rp), intent(in) :: dt
    integer, intent(in) :: order, stage

    value = time - dt
    if (order .eq. 3) then
       select case (stage)
       case (1)
          value = time - dt
       case (2)
          value = time
       case (3)
          value = time - 0.5_rp * dt
       end select
    end if
  end function euler_idp_stage_time

  !> Return the graph CFL from already prepared graph data.
  pure real(kind=rp) function euler_idp_backend_graph_cfl(this, dt) &
       result(cfl)
    class(euler_idp_backend_t), intent(in) :: this
    real(kind=rp), intent(in) :: dt

    cfl = dt * this%max_graph_rate
  end function euler_idp_backend_graph_cfl

  !> Reset diagnostics before a Forward Euler map or SSPRK stage.
  subroutine euler_idp_diagnostics_reset(this)
    class(euler_idp_diagnostics_t), intent(inout) :: this

    this%stage = 0
    this%stage_time = 0.0_rp
    this%min_density = huge(1.0_rp)
    this%min_internal_energy = huge(1.0_rp)
    this%min_pressure = huge(1.0_rp)
    this%min_specific_entropy = huge(1.0_rp)
    this%max_graph_cfl = 0.0_rp
    this%max_nodal_wave_speed = 0.0_rp
    this%max_graph_wave_speed = 0.0_rp
    this%max_graph_rate = 0.0_rp
    this%min_convex_weight = 1.0_rp
    this%maximum_floor_timestep = huge(1.0_rp)
    this%entropy_viscosity_enabled = .false.
    this%max_entropy_viscosity = 0.0_rp
    this%mean_entropy_viscosity = 0.0_rp
    this%entropy_viscosity_conservation = 0.0_rp
    this%entropy_viscosity_element_compatibility = 0.0_rp
    this%high_order_conservation = 0.0_rp
    this%low_order_conservation = 0.0_rp
    this%limited_conservation = 0.0_rp
    this%correction_assembly_error = 0.0_rp
    this%correction_global_compatibility = 0.0_rp
    this%correction_element_compatibility = 0.0_rp
    this%correction_face_mismatch = 0.0_rp
    this%element_compatibility = 0.0_rp
    this%graph_compatibility = 0.0_rp
    this%reconstruction_residual = 0.0_rp
    this%forward_euler_time = 0.0_rp
    this%reconstruction_time = 0.0_rp
    this%graph_gs_time = 0.0_rp
    this%graph_reduction_time = 0.0_rp
    this%max_correction_flux = 0.0_rp
    this%rms_correction_flux = 0.0_rp
    this%min_limiter = 1.0_rp
    this%mean_limiter = 1.0_rp
    this%max_limiter = 1.0_rp
    this%limited_edge_fraction = 0.0_rp
    this%limiter_weight_error = 0.0_rp
    this%min_density_lower_bound = huge(1.0_rp)
    this%max_density_upper_bound = -huge(1.0_rp)
    this%max_density_bound_relaxation = 0.0_rp
    this%max_density_lower_violation = 0.0_rp
    this%max_density_upper_violation = 0.0_rp
    this%max_entropy_lower_violation = 0.0_rp
    this%density_limited_edges = 0
    this%internal_energy_limited_edges = 0
    this%entropy_limited_edges = 0
    this%before_boundary = euler_idp_state_observation_t()
    this%after_boundary = euler_idp_state_observation_t()
  end subroutine euler_idp_diagnostics_reset


  !> Guaranteed upper bound for the two-state Euler Riemann fan speed.
  pure real(kind=rp) function euler_idp_maximum_wave_speed(left, right, &
       normal, gamma) result(speed)
    real(kind=rp), intent(in) :: left(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: right(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: normal(3), gamma
    real(kind=rp) :: reverse_normal(3)

    reverse_normal = -normal
    speed = max(euler_idp_ordered_wave_speed(left, right, normal, gamma), &
         euler_idp_ordered_wave_speed(right, left, reverse_normal, gamma))
  end function euler_idp_maximum_wave_speed

  !> Euler flux contracted with a Cartesian vector.
  pure subroutine euler_idp_flux_dot_vector(state, vector, gamma, flux)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: vector(3), gamma
    real(kind=rp), intent(out) :: flux(EULER_IDP_NCOMP)
    real(kind=rp) :: pressure, velocity(3), velocity_dot_vector

    velocity = state(2:4) / state(1)
    pressure = (gamma - 1.0_rp) * euler_idp_internal_energy(state)
    velocity_dot_vector = dot_product(velocity, vector)
    flux(1) = dot_product(state(2:4), vector)
    flux(2:4) = state(2:4) * velocity_dot_vector + pressure * vector
    flux(5) = (state(5) + pressure) * velocity_dot_vector
  end subroutine euler_idp_flux_dot_vector

  !> Construct the common bar state associated with one undirected edge.
  pure subroutine euler_idp_bar_state(left, right, flux_difference, &
       viscosity, bar_state)
    real(kind=rp), intent(in) :: left(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: right(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: flux_difference(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: viscosity
    real(kind=rp), intent(out) :: bar_state(EULER_IDP_NCOMP)

    if (abs(viscosity) .le. tiny(1.0_rp)) then
       bar_state = 0.5_rp * (left + right)
    else
       bar_state = 0.5_rp * (left + right) - &
            flux_difference / (2.0_rp * viscosity)
    end if
  end subroutine euler_idp_bar_state

  !> Internal-energy density of one conserved Euler state.
  pure real(kind=rp) function euler_idp_internal_energy(state) result(value)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)

    value = state(5) - 0.5_rp * dot_product(state(2:4), state(2:4)) / &
         state(1)
  end function euler_idp_internal_energy

  !> Maximum timestep along one update direction respecting admissibility.
  pure real(kind=rp) function euler_idp_internal_energy_timestep(state, &
       residual, internal_energy_floor, upper_bound) result(dt_limit)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: residual(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: internal_energy_floor
    real(kind=rp), intent(in) :: upper_bound
    real(kind=rp) :: coefficient_a, coefficient_b, coefficient_c
    real(kind=rp) :: scaled_a, scaled_b, scaled_c, coefficient_scale
    real(kind=rp) :: discriminant, square_root, q, guarded_floor
    real(kind=rp) :: first_root, second_root, root
    real(kind=rp) :: density_limit, linear_tolerance, discriminant_tolerance
    real(kind=rp) :: lower, upper, midpoint
    real(kind=rp) :: trial(EULER_IDP_NCOMP)
    integer :: iteration
    logical :: analytic_root

    dt_limit = max(0.0_rp, upper_bound)
    if (dt_limit .le. 0.0_rp) return

    ! Keep the density strictly positive. The downward margin also absorbs
    ! rounding in the division that locates the density root.
    if (residual(1) .gt. 0.0_rp) then
       density_limit = state(1) / residual(1)
       if (density_limit .le. dt_limit) then
          dt_limit = max(0.0_rp, density_limit * &
               (1.0_rp - 64.0_rp * epsilon(1.0_rp)))
       end if
    end if
    if (dt_limit .le. 0.0_rp) return

    trial = state - dt_limit * residual
    if (euler_idp_state_is_admissible(trial, internal_energy_floor)) return
    if (.not. euler_idp_state_is_admissible(state, &
         internal_energy_floor)) then
       dt_limit = 0.0_rp
       return
    end if
    upper = dt_limit
    ! The caller forms the candidate again, possibly with different rounding.
    guarded_floor = internal_energy_floor + 16.0_rp * epsilon(1.0_rp) * &
         max(abs(state(5)), dt_limit * abs(residual(5)))

    ! For U(t) = U - t R, multiply the internal-energy constraint by rho(t):
    !   a*t**2 + b*t + c >= 0.
    ! Scaling all coefficients equally makes the stable quadratic formula
    ! insensitive to the absolute magnitude of the conserved state.
    coefficient_a = residual(1) * residual(5) - &
         0.5_rp * dot_product(residual(2:4), residual(2:4))
    coefficient_b = -(state(1) * residual(5) + residual(1) * &
         (state(5) - internal_energy_floor)) + &
         dot_product(state(2:4), residual(2:4))
    coefficient_c = state(1) * (state(5) - internal_energy_floor) - &
         0.5_rp * dot_product(state(2:4), state(2:4))
    coefficient_scale = max(abs(coefficient_a), abs(coefficient_b), &
         abs(coefficient_c), tiny(1.0_rp))
    analytic_root = all(ieee_is_finite([coefficient_a, coefficient_b, &
         coefficient_c, coefficient_scale]))
    root = huge(1.0_rp)

    if (analytic_root) then
       scaled_a = coefficient_a / coefficient_scale
       scaled_b = coefficient_b / coefficient_scale
       scaled_c = coefficient_c / coefficient_scale
       linear_tolerance = 64.0_rp * epsilon(1.0_rp) * max( &
            abs(scaled_b) / dt_limit, &
            abs(scaled_c) / max(dt_limit * dt_limit, tiny(1.0_rp)), &
            tiny(1.0_rp))

       if (abs(scaled_a) .le. linear_tolerance) then
          if (scaled_b .lt. 0.0_rp) then
             root = -scaled_c / scaled_b
          else
             analytic_root = .false.
          end if
       else
          discriminant = scaled_b * scaled_b - &
               4.0_rp * scaled_a * scaled_c
          discriminant_tolerance = 64.0_rp * epsilon(1.0_rp) * max( &
               scaled_b * scaled_b, &
               abs(4.0_rp * scaled_a * scaled_c), tiny(1.0_rp))
          if (discriminant .lt. -discriminant_tolerance) then
             analytic_root = .false.
          else
             discriminant = max(0.0_rp, discriminant)
             square_root = sqrt(discriminant)
             q = -0.5_rp * (scaled_b + sign(square_root, scaled_b))
             first_root = huge(1.0_rp)
             second_root = huge(1.0_rp)
             if (abs(scaled_a) .gt. tiny(1.0_rp)) then
                first_root = q / scaled_a
             end if
             if (abs(q) .gt. tiny(1.0_rp)) then
                second_root = scaled_c / q
             end if
             if (first_root .gt. 0.0_rp .and. &
                  first_root .le. dt_limit * (1.0_rp + &
                  64.0_rp * epsilon(1.0_rp))) root = first_root
             if (second_root .gt. 0.0_rp .and. &
                  second_root .le. dt_limit * (1.0_rp + &
                  64.0_rp * epsilon(1.0_rp))) then
                root = min(root, second_root)
             end if
             if (coefficient_c .le. 0.0_rp .and. &
                  (coefficient_b .lt. 0.0_rp .or. &
                  (coefficient_b .eq. 0.0_rp .and. &
                  coefficient_a .lt. 0.0_rp))) root = 0.0_rp
             if (.not. ieee_is_finite(root) .or. &
                  root .eq. huge(1.0_rp)) analytic_root = .false.
          end if
       end if
    end if

    if (analytic_root) then
       dt_limit = max(0.0_rp, min(dt_limit, root) * &
            (1.0_rp - 64.0_rp * epsilon(1.0_rp)))
       trial = state - dt_limit * residual
       if (euler_idp_state_is_admissible(trial, &
            guarded_floor)) return
    end if

    ! Degenerate or numerically ambiguous polynomials retain a robust fallback.
    lower = 0.0_rp
    do iteration = 1, 64
       midpoint = 0.5_rp * (lower + upper)
       trial = state - midpoint * residual
       if (euler_idp_state_is_admissible(trial, &
            guarded_floor)) then
          lower = midpoint
       else
          upper = midpoint
       end if
    end do
    dt_limit = lower
  end function euler_idp_internal_energy_timestep

  !> Non-iterative Guermond-Popov wave-speed estimate for ordered Riemann data.
  !! The estimate is a guaranteed upper bound for 1 < gamma <= 5/3.
  pure real(kind=rp) function euler_idp_ordered_wave_speed(left, right, &
       normal, gamma) result(speed)
    real(kind=rp), intent(in) :: left(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: right(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: normal(3), gamma
    real(kind=rp) :: a_left, a_right, a_min, a_max
    real(kind=rp) :: pressure_left, pressure_right, pressure_min, pressure_max
    real(kind=rp) :: density_min
    real(kind=rp) :: velocity_left, velocity_right
    real(kind=rp) :: exponent, pressure_ratio, phi_min, phi_max
    real(kind=rp) :: pressure_two_rarefaction, pressure_upper
    real(kind=rp) :: coefficient_a, coefficient_b, numerator
    real(kind=rp) :: left_speed, right_speed

    pressure_left = (gamma - 1.0_rp) * &
         euler_idp_internal_energy(left)
    pressure_right = (gamma - 1.0_rp) * &
         euler_idp_internal_energy(right)
    a_left = sqrt(gamma * pressure_left / left(1))
    a_right = sqrt(gamma * pressure_right / right(1))
    velocity_left = dot_product(left(2:4), normal) / left(1)
    velocity_right = dot_product(right(2:4), normal) / right(1)

    if (pressure_left .le. pressure_right) then
       pressure_min = pressure_left
       pressure_max = pressure_right
       density_min = left(1)
       a_min = a_left
       a_max = a_right
    else
       pressure_min = pressure_right
       pressure_max = pressure_left
       density_min = right(1)
       a_min = a_right
       a_max = a_left
    end if

    exponent = (gamma - 1.0_rp) / (2.0_rp * gamma)
    pressure_ratio = (pressure_min / pressure_max)**exponent
    phi_min = 2.0_rp * a_max * (pressure_ratio - 1.0_rp) / &
         (gamma - 1.0_rp) + velocity_right - velocity_left
    if (phi_min .ge. 0.0_rp) then
       speed = max(max(-(velocity_left - a_left), 0.0_rp), &
            max(velocity_right + a_right, 0.0_rp))
       return
    end if

    coefficient_a = 2.0_rp / ((gamma + 1.0_rp) * density_min)
    coefficient_b = pressure_min * (gamma - 1.0_rp) / &
         (gamma + 1.0_rp)
    phi_max = (pressure_max - pressure_min) * &
         sqrt(coefficient_a / (pressure_max + coefficient_b)) + &
         velocity_right - velocity_left
    numerator = a_min + a_max - 0.5_rp * (gamma - 1.0_rp) * &
         (velocity_right - velocity_left)
    pressure_two_rarefaction = pressure_min * &
         (numerator / (a_min + a_max * pressure_ratio))** &
         (2.0_rp * gamma / (gamma - 1.0_rp))
    if (phi_max .lt. 0.0_rp) then
       pressure_upper = pressure_two_rarefaction
    else
       pressure_upper = min(pressure_max, pressure_two_rarefaction)
    end if

    left_speed = velocity_left - a_left * sqrt(1.0_rp + &
         max((pressure_upper - pressure_left) / pressure_left, 0.0_rp) * &
         (gamma + 1.0_rp) / (2.0_rp * gamma))
    right_speed = velocity_right + a_right * sqrt(1.0_rp + &
         max((pressure_upper - pressure_right) / pressure_right, 0.0_rp) * &
         (gamma + 1.0_rp) / (2.0_rp * gamma))
    speed = max(max(-left_speed, 0.0_rp), max(right_speed, 0.0_rp))
  end function euler_idp_ordered_wave_speed

  !> Test positive density and the internal-energy-density threshold.
  pure logical function euler_idp_state_is_admissible(state, &
       internal_energy_floor) result(admissible)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: internal_energy_floor
    real(kind=rp) :: energy_margin

    admissible = .false.
    if (.not. all(ieee_is_finite(state))) return
    if (state(1) .le. 0.0_rp) return
    energy_margin = state(1) * (state(5) - internal_energy_floor) - &
         0.5_rp * dot_product(state(2:4), state(2:4))
    admissible = ieee_is_finite(energy_margin) .and. &
         energy_margin .ge. 0.0_rp
  end function euler_idp_state_is_admissible

  !> Return log(p) - gamma*log(rho) for a gamma-law gas.
  pure real(kind=rp) function euler_idp_specific_entropy(state, gamma) &
       result(entropy)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: gamma
    real(kind=rp) :: internal_energy, pressure

    entropy = -huge(1.0_rp)
    if (.not. all(ieee_is_finite(state))) return
    if (state(1) .le. 0.0_rp) return
    internal_energy = state(5) - &
         0.5_rp * dot_product(state(2:4), state(2:4)) / state(1)
    pressure = (gamma - 1.0_rp) * internal_energy
    if (.not. ieee_is_finite(pressure) .or. pressure .le. 0.0_rp) return
    entropy = log(pressure) - gamma * log(state(1))
  end function euler_idp_specific_entropy

  !> Roundoff tolerance for entropy recovered from a conserved state.
  pure real(kind=rp) function euler_idp_entropy_tolerance(state, &
       entropy_lower) result(tolerance)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: entropy_lower
    real(kind=rp) :: conditioning, internal_energy, kinetic_energy

    conditioning = 1.0_rp
    if (all(ieee_is_finite(state)) .and. state(1) .gt. 0.0_rp) then
       kinetic_energy = 0.5_rp * dot_product(state(2:4), state(2:4)) / &
            state(1)
       internal_energy = state(5) - kinetic_energy
       if (ieee_is_finite(internal_energy) .and. &
            internal_energy .gt. tiny(1.0_rp)) then
          conditioning = (abs(state(5)) + abs(kinetic_energy)) / &
               internal_energy
       end if
    end if
    ! Cap the conditioning correction so it cannot mask a physical violation.
    tolerance = min(sqrt(epsilon(1.0_rp)), &
         64.0_rp * epsilon(1.0_rp) * &
         max(1.0_rp, abs(entropy_lower), conditioning))
  end function euler_idp_entropy_tolerance

  !> Test the Euler admissible set and a local minimum entropy constraint.
  pure logical function euler_idp_entropy_is_admissible(state, gamma, &
       entropy_lower, internal_energy_floor) result(admissible)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: gamma, entropy_lower
    real(kind=rp), intent(in) :: internal_energy_floor
    real(kind=rp) :: entropy, tolerance

    admissible = euler_idp_state_is_admissible(state, &
         internal_energy_floor)
    if (.not. admissible) return
    entropy = euler_idp_specific_entropy(state, gamma)
    tolerance = euler_idp_entropy_tolerance(state, entropy_lower)
    admissible = ieee_is_finite(entropy) .and. &
         entropy .ge. entropy_lower - tolerance
  end function euler_idp_entropy_is_admissible

  !> Compute one-ring entropy minima from an immutable stage field.
  pure subroutine euler_idp_local_entropy_bounds(stage_entropy, left, right, &
       entropy_lower_bound)
    real(kind=rp), intent(in) :: stage_entropy(:,:,:,:)
    integer, intent(in) :: left(:,:), right(:,:)
    real(kind=rp), intent(out) :: entropy_lower_bound(:,:,:,:)
    integer :: edge

    entropy_lower_bound = stage_entropy
    do edge = 1, size(left, 2)
       associate(a => left(:, edge), b => right(:, edge))
         entropy_lower_bound(a(1),a(2),a(3),a(4)) = min( &
              entropy_lower_bound(a(1),a(2),a(3),a(4)), &
              stage_entropy(b(1),b(2),b(3),b(4)))
         entropy_lower_bound(b(1),b(2),b(3),b(4)) = min( &
              entropy_lower_bound(b(1),b(2),b(3),b(4)), &
              stage_entropy(a(1),a(2),a(3),a(4)))
       end associate
    end do
  end subroutine euler_idp_local_entropy_bounds

  !> Apply the averaging relaxation of Guermond et al. to density bounds.
  pure subroutine euler_idp_relax_density_bounds(strict_lower, strict_upper, &
       second_difference, relaxation_factor, nodal_mass, domain_volume, &
       dimension, relaxed_lower, relaxed_upper)
    real(kind=rp), intent(in) :: strict_lower, strict_upper
    real(kind=rp), intent(in) :: second_difference, relaxation_factor
    real(kind=rp), intent(in) :: nodal_mass, domain_volume
    integer, intent(in) :: dimension
    real(kind=rp), intent(out) :: relaxed_lower, relaxed_upper
    real(kind=rp) :: relaxation, r_h

    r_h = (nodal_mass / domain_volume)**(1.5_rp / real(dimension, rp))
    relaxation = relaxation_factor * abs(second_difference)
    relaxed_lower = max((1.0_rp - r_h) * strict_lower, &
         strict_lower - relaxation)
    relaxed_upper = strict_upper + relaxation
  end subroutine euler_idp_relax_density_bounds

  !> Relax a strict logarithmic specific-entropy lower bound.
  !!
  !! The positive entropy variable is chi = exp(s).  The midpoint maximum
  !! supplies the graph-linear variation used by the entropy relaxation of
  !! Guermond et al.  Its cap r_h = (m_star / |D|)^(3/(2d)) follows the mesh
  !! scale instead of remaining a fixed input.  Evaluation in logarithmic
  !! variables avoids overflow and makes the result invariant under an
  !! additive shift of the entropy.
  pure subroutine euler_idp_relax_entropy_bound(strict_lower, midpoint_max, &
       relaxation_factor, nodal_mass, domain_volume, dimension, relaxed_lower)
    real(kind=rp), intent(in) :: strict_lower, midpoint_max
    real(kind=rp), intent(in) :: relaxation_factor
    real(kind=rp), intent(in) :: nodal_mass, domain_volume
    integer, intent(in) :: dimension
    real(kind=rp), intent(out) :: relaxed_lower
    real(kind=rp) :: cutoff, entropy_difference, relative_variation
    real(kind=rp) :: relaxation, r_h

    r_h = (nodal_mass / domain_volume)**(1.5_rp / real(dimension, rp))
    cutoff = log(1.0_rp + r_h / relaxation_factor)
    entropy_difference = max(0.0_rp, midpoint_max - strict_lower)
    relative_variation = exp(min(entropy_difference, cutoff)) - 1.0_rp
    relaxation = min(r_h, &
         relaxation_factor * max(0.0_rp, relative_variation))
    relaxed_lower = strict_lower + log(1.0_rp - relaxation)
  end subroutine euler_idp_relax_entropy_bound

  !> Limit one directed auxiliary correction to the Euler admissible set.
  pure subroutine euler_idp_limit_endpoint(base, correction, density_lower, &
       density_upper, entropy_lower, gamma, internal_energy_floor, limit, &
       density_limited, energy_limited, entropy_limited, &
       enforce_internal_energy, enforce_entropy, check_base_constraints)
    real(kind=rp), intent(in) :: base(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: correction(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: density_lower, density_upper, entropy_lower
    real(kind=rp), intent(in) :: gamma, internal_energy_floor
    real(kind=rp), intent(out) :: limit
    logical, intent(out) :: density_limited, energy_limited, entropy_limited
    logical, intent(in), optional :: enforce_internal_energy, enforce_entropy
    logical, intent(in), optional :: check_base_constraints
    real(kind=rp) :: trial(EULER_IDP_NCOMP)
    logical :: enforce_energy_constraint, enforce_entropy_constraint
    logical :: validate_base

    enforce_energy_constraint = .true.
    enforce_entropy_constraint = .true.
    if (present(enforce_internal_energy)) then
       enforce_energy_constraint = enforce_internal_energy
    end if
    if (present(enforce_entropy)) then
       enforce_entropy_constraint = enforce_entropy
    end if
    validate_base = .true.
    if (present(check_base_constraints)) then
       validate_base = check_base_constraints
    end if
    limit = 1.0_rp
    density_limited = .false.
    energy_limited = .false.
    entropy_limited = .false.
    if (validate_base) then
       if (.not. all(ieee_is_finite(base))) then
          limit = 0.0_rp
          energy_limited = enforce_energy_constraint
          entropy_limited = enforce_entropy_constraint
          return
       end if
       if (base(1) .le. 0.0_rp .or. base(1) .lt. density_lower .or. &
            base(1) .gt. density_upper) then
          limit = 0.0_rp
          density_limited = .true.
          return
       end if
       if (enforce_energy_constraint .and. &
            .not. euler_idp_state_is_admissible(base, &
            internal_energy_floor)) then
          limit = 0.0_rp
          energy_limited = .true.
          return
       end if
       if (enforce_entropy_constraint .and. &
            .not. euler_idp_entropy_is_admissible(base, gamma, &
            entropy_lower, internal_energy_floor)) then
          limit = 0.0_rp
          entropy_limited = .true.
          return
       end if
    end if

    if (correction(1) .lt. 0.0_rp) then
       limit = min(limit, (base(1) - density_lower) / (-correction(1)))
    else if (correction(1) .gt. 0.0_rp) then
       limit = min(limit, (density_upper - base(1)) / correction(1))
    end if
    if (correction(1) .ne. 0.0_rp) then
       limit = max(0.0_rp, limit)
       density_limited = limit .lt. 1.0_rp
       if (density_limited .and. limit .gt. 0.0_rp) then
          limit = limit * (1.0_rp - 32.0_rp * epsilon(1.0_rp))
       end if
    end if

    trial = base + limit * correction
    energy_limited = enforce_energy_constraint .and. &
         .not. euler_idp_state_is_admissible(trial, internal_energy_floor)
    entropy_limited = enforce_entropy_constraint .and. &
         .not. euler_idp_entropy_is_admissible(trial, gamma, entropy_lower, &
         internal_energy_floor)
    if (.not. energy_limited .and. .not. entropy_limited) return

    if (energy_limited) then
       call limit_internal_energy(base, correction, internal_energy_floor, &
            limit)
    end if

    trial = base + limit * correction
    if (enforce_entropy_constraint .and. &
         .not. euler_idp_entropy_is_admissible(trial, gamma, entropy_lower, &
         internal_energy_floor)) then
       call limit_entropy(base, correction, entropy_lower, gamma, limit)
       entropy_limited = .true.
    end if
    if (limit .gt. 0.0_rp) then
       limit = limit * (1.0_rp - 32.0_rp * epsilon(1.0_rp))
    end if
    if (enforce_entropy_constraint .and. limit .gt. 0.0_rp) then
       trial = base + limit * correction
       if (.not. euler_idp_entropy_is_admissible(trial, gamma, &
            entropy_lower, internal_energy_floor)) then
          call limit_entropy_bisection(base, correction, entropy_lower, &
               gamma, internal_energy_floor, limit)
          entropy_limited = .true.
       end if
    end if
  end subroutine euler_idp_limit_endpoint

  !> Limit a segment with a safeguarded Newton--secant line search.
  pure subroutine limit_internal_energy(base, correction, energy_floor, &
       limit)
    real(kind=rp), intent(in) :: base(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: correction(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: energy_floor
    real(kind=rp), intent(inout) :: limit
    real(kind=rp) :: left, right, next_left, next_right
    real(kind=rp) :: value_left, value_right, value_next, value_tolerance
    real(kind=rp) :: slope
    real(kind=rp) :: trial(EULER_IDP_NCOMP)
    integer :: iteration

    left = 0.0_rp
    right = limit
    value_left = internal_energy_margin(base, energy_floor)
    trial = base + right * correction
    value_right = internal_energy_margin(trial, energy_floor)
    value_tolerance = 256.0_rp * epsilon(1.0_rp) * &
         max(1.0_rp, abs(value_left), abs(value_right))
    if (value_left .le. 0.0_rp) then
       limit = 0.0_rp
       return
    end if

    do iteration = 1, 16
       if (right - left .le. 256.0_rp * epsilon(1.0_rp) * &
            max(1.0_rp, right)) exit
       if (value_left .le. value_right) exit

       slope = (value_right - value_left) / (right - left)
       next_left = left - value_left / slope
       if (.not. ieee_is_finite(next_left) .or. next_left .le. left .or. &
            next_left .ge. right) exit
       trial = base + next_left * correction
       value_next = internal_energy_margin(trial, energy_floor)
       if (.not. ieee_is_finite(value_next) .or. &
            value_next .lt. -value_tolerance) exit
       left = next_left
       if (abs(value_next) .le. value_tolerance) exit
       value_left = value_next

       trial = base + right * correction
       slope = internal_energy_margin_derivative(trial, correction)
       if (.not. ieee_is_finite(slope) .or. slope .ge. 0.0_rp) exit
       next_right = right - value_right / slope
       if (.not. ieee_is_finite(next_right) .or. next_right .le. left .or. &
            next_right .ge. right) exit
       trial = base + next_right * correction
       value_next = internal_energy_margin(trial, energy_floor)
       if (.not. ieee_is_finite(value_next) .or. &
            value_next .gt. value_tolerance) exit
       if (value_next .ge. 0.0_rp) then
          left = next_right
          exit
       end if
       right = next_right
       value_right = value_next
    end do
    limit = left
  end subroutine limit_internal_energy

  !> Limit the concave gamma-law entropy constraint by Newton--secant.
  pure subroutine limit_entropy(base, correction, entropy_lower, gamma, &
       limit)
    real(kind=rp), intent(in) :: base(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: correction(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: entropy_lower, gamma
    real(kind=rp), intent(inout) :: limit
    real(kind=rp) :: left, right, next_left, next_right
    real(kind=rp) :: value_left, value_right, value_next, value_tolerance
    real(kind=rp) :: slope
    real(kind=rp) :: trial(EULER_IDP_NCOMP)
    integer :: iteration

    left = 0.0_rp
    right = limit
    value_left = entropy_margin(base, entropy_lower, gamma)
    trial = base + right * correction
    value_right = entropy_margin(trial, entropy_lower, gamma)
    value_tolerance = 256.0_rp * epsilon(1.0_rp) * &
         max(1.0_rp, abs(value_left), abs(value_right))
    if (value_left .le. 0.0_rp) then
       limit = 0.0_rp
       return
    end if

    do iteration = 1, 16
       if (right - left .le. 256.0_rp * epsilon(1.0_rp) * &
            max(1.0_rp, right)) exit
       if (value_left .le. value_right) exit

       slope = (value_right - value_left) / (right - left)
       next_left = left - value_left / slope
       if (.not. ieee_is_finite(next_left) .or. next_left .le. left .or. &
            next_left .ge. right) exit
       trial = base + next_left * correction
       value_next = entropy_margin(trial, entropy_lower, gamma)
       if (.not. ieee_is_finite(value_next) .or. &
            value_next .lt. -value_tolerance) exit
       left = next_left
       if (abs(value_next) .le. value_tolerance) exit
       value_left = value_next

       trial = base + right * correction
       slope = entropy_margin_derivative(trial, correction, entropy_lower, &
            gamma)
       if (.not. ieee_is_finite(slope) .or. slope .ge. 0.0_rp) exit
       next_right = right - value_right / slope
       if (.not. ieee_is_finite(next_right) .or. next_right .le. left .or. &
            next_right .ge. right) exit
       trial = base + next_right * correction
       value_next = entropy_margin(trial, entropy_lower, gamma)
       if (.not. ieee_is_finite(value_next) .or. &
            value_next .gt. value_tolerance) exit
       if (value_next .ge. 0.0_rp) then
          left = next_right
          exit
       end if
       right = next_right
       value_right = value_next
    end do
    limit = left
  end subroutine limit_entropy

  !> Restore the public entropy predicate when the fast line solve ends in
  !! the roundoff band around the admissible-set boundary.
  pure subroutine limit_entropy_bisection(base, correction, entropy_lower, &
       gamma, internal_energy_floor, limit)
    real(kind=rp), intent(in) :: base(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: correction(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: entropy_lower, gamma
    real(kind=rp), intent(in) :: internal_energy_floor
    real(kind=rp), intent(inout) :: limit
    real(kind=rp) :: left, right, midpoint
    real(kind=rp) :: trial(EULER_IDP_NCOMP)
    integer :: iteration

    trial = base + limit * correction
    if (euler_idp_entropy_is_admissible(trial, gamma, entropy_lower, &
         internal_energy_floor)) return
    if (.not. euler_idp_entropy_is_admissible(base, gamma, entropy_lower, &
         internal_energy_floor)) then
       limit = 0.0_rp
       return
    end if

    left = 0.0_rp
    right = limit
    do iteration = 1, 64
       midpoint = left + 0.5_rp * (right - left)
       if (midpoint .le. left .or. midpoint .ge. right) exit
       trial = base + midpoint * correction
       if (euler_idp_entropy_is_admissible(trial, gamma, entropy_lower, &
            internal_energy_floor)) then
          left = midpoint
       else
          right = midpoint
       end if
    end do
    limit = left
  end subroutine limit_entropy_bisection

  !> Internal-energy-density margin along a positive-density segment.
  pure real(kind=rp) function internal_energy_margin(state, energy_floor) &
       result(margin)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: energy_floor

    margin = state(5) - energy_floor - &
         0.5_rp * dot_product(state(2:4), state(2:4)) / state(1)
  end function internal_energy_margin

  !> Directional derivative of the internal-energy-density margin.
  pure real(kind=rp) function internal_energy_margin_derivative(state, &
       correction) result(derivative)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: correction(EULER_IDP_NCOMP)
    real(kind=rp) :: momentum_squared

    momentum_squared = dot_product(state(2:4), state(2:4))
    derivative = correction(5) - &
         dot_product(state(2:4), correction(2:4)) / state(1) + &
         0.5_rp * momentum_squared * correction(1) / state(1)**2
  end function internal_energy_margin_derivative

  !> Concave form of the gamma-law minimum-entropy constraint.
  pure real(kind=rp) function entropy_margin(state, entropy_lower, gamma) &
       result(margin)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: entropy_lower, gamma
    real(kind=rp) :: pressure

    pressure = (gamma - 1.0_rp) * internal_energy_margin(state, 0.0_rp)
    margin = pressure - exp(entropy_lower) * state(1)**gamma
  end function entropy_margin

  !> Directional derivative of the concave entropy margin.
  pure real(kind=rp) function entropy_margin_derivative(state, correction, &
       entropy_lower, gamma) result(derivative)
    real(kind=rp), intent(in) :: state(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: correction(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: entropy_lower, gamma

    derivative = (gamma - 1.0_rp) * &
         internal_energy_margin_derivative(state, correction) - &
         exp(entropy_lower) * gamma * state(1)**(gamma - 1.0_rp) * &
         correction(1)
  end function entropy_margin_derivative

  !> Return one symmetric coefficient for a complete vector edge correction.
  pure subroutine euler_idp_limit_edge(left_base, right_base, &
       left_correction, right_correction, left_density_lower, &
       left_density_upper, right_density_lower, right_density_upper, &
       left_entropy_lower, right_entropy_lower, gamma, &
       internal_energy_floor, limit, density_limited, energy_limited, &
       entropy_limited, enforce_internal_energy, enforce_entropy, &
       check_base_constraints)
    real(kind=rp), intent(in) :: left_base(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: right_base(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: left_correction(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: right_correction(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: left_density_lower, left_density_upper
    real(kind=rp), intent(in) :: right_density_lower, right_density_upper
    real(kind=rp), intent(in) :: left_entropy_lower, right_entropy_lower
    real(kind=rp), intent(in) :: gamma, internal_energy_floor
    real(kind=rp), intent(out) :: limit
    logical, intent(out) :: density_limited, energy_limited, entropy_limited
    logical, intent(in), optional :: enforce_internal_energy, enforce_entropy
    logical, intent(in), optional :: check_base_constraints
    real(kind=rp) :: left_limit, right_limit
    logical :: left_density_limited, right_density_limited
    logical :: left_energy_limited, right_energy_limited
    logical :: left_entropy_limited, right_entropy_limited
    logical :: enforce_entropy_constraint
    real(kind=rp) :: left_trial(EULER_IDP_NCOMP)
    real(kind=rp) :: right_trial(EULER_IDP_NCOMP)

    enforce_entropy_constraint = .true.
    if (present(enforce_entropy)) then
       enforce_entropy_constraint = enforce_entropy
    end if

    call euler_idp_limit_endpoint(left_base, left_correction, &
         left_density_lower, left_density_upper, left_entropy_lower, gamma, &
         internal_energy_floor, left_limit, &
         left_density_limited, left_energy_limited, left_entropy_limited, &
         enforce_internal_energy, enforce_entropy, check_base_constraints)
    call euler_idp_limit_endpoint(right_base, right_correction, &
         right_density_lower, right_density_upper, right_entropy_lower, &
         gamma, internal_energy_floor, right_limit, &
         right_density_limited, right_energy_limited, &
         right_entropy_limited, enforce_internal_energy, enforce_entropy, &
         check_base_constraints)
    limit = min(left_limit, right_limit)
    density_limited = left_density_limited .or. right_density_limited
    energy_limited = left_energy_limited .or. right_energy_limited
    entropy_limited = left_entropy_limited .or. right_entropy_limited
    if (enforce_entropy_constraint .and. limit .gt. 0.0_rp .and. &
         limit .lt. 1.0_rp) then
       left_trial = left_base + limit * left_correction
       right_trial = right_base + limit * right_correction
       if (.not. euler_idp_entropy_is_admissible(left_trial, gamma, &
            left_entropy_lower, internal_energy_floor) .or. &
            .not. euler_idp_entropy_is_admissible(right_trial, gamma, &
            right_entropy_lower, internal_energy_floor)) then
          call limit_edge_entropy_bisection(left_base, right_base, &
               left_correction, right_correction, left_entropy_lower, &
               right_entropy_lower, gamma, internal_energy_floor, limit)
          entropy_limited = .true.
       end if
    end if
  end subroutine euler_idp_limit_edge

  !> Restore admissibility of both endpoints after taking their symmetric
  !! minimum coefficient in the entropy roundoff band.
  pure subroutine limit_edge_entropy_bisection(left_base, right_base, &
       left_correction, right_correction, left_entropy_lower, &
       right_entropy_lower, gamma, internal_energy_floor, limit)
    real(kind=rp), intent(in) :: left_base(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: right_base(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: left_correction(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: right_correction(EULER_IDP_NCOMP)
    real(kind=rp), intent(in) :: left_entropy_lower, right_entropy_lower
    real(kind=rp), intent(in) :: gamma, internal_energy_floor
    real(kind=rp), intent(inout) :: limit
    real(kind=rp) :: left, right, midpoint
    real(kind=rp) :: left_trial(EULER_IDP_NCOMP)
    real(kind=rp) :: right_trial(EULER_IDP_NCOMP)
    logical :: left_admissible, right_admissible
    integer :: iteration

    left = 0.0_rp
    right = limit
    do iteration = 1, 64
       midpoint = left + 0.5_rp * (right - left)
       if (midpoint .le. left .or. midpoint .ge. right) exit
       left_trial = left_base + midpoint * left_correction
       right_trial = right_base + midpoint * right_correction
       left_admissible = euler_idp_entropy_is_admissible(left_trial, &
            gamma, left_entropy_lower, internal_energy_floor)
       right_admissible = euler_idp_entropy_is_admissible(right_trial, &
            gamma, right_entropy_lower, internal_energy_floor)
       if (left_admissible .and. right_admissible) then
          left = midpoint
       else
          right = midpoint
       end if
    end do
    limit = left
  end subroutine limit_edge_entropy_bisection

end module euler_idp_backend
