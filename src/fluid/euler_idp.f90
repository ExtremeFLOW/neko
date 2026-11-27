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
!> Configuration and orchestration for the Euler GLL-IDP solver.
module euler_idp
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use json_module, only : json_file
  use json_utils, only : json_get_or_default
  use num_types, only : rp, dp
  use field, only : field_t
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t
  use bc_list, only : bc_list_t
  use time_state, only : time_state_t
  use euler_idp_backend, only : EULER_IDP_NCOMP, &
       EULER_IDP_DIAGNOSTICS_OFF, EULER_IDP_DIAGNOSTICS_SAFETY, &
       EULER_IDP_DIAGNOSTICS_FULL, euler_idp_backend_t, &
       euler_idp_state_observation_t, euler_idp_diagnostics_t, &
       euler_idp_stage_time
  use euler_idp_cpu, only : euler_idp_cpu_t
  use euler_idp_device, only : euler_idp_device_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  implicit none
  private

  public :: EULER_IDP_NCOMP
  public :: EULER_IDP_DIAGNOSTICS_OFF, EULER_IDP_DIAGNOSTICS_SAFETY
  public :: EULER_IDP_DIAGNOSTICS_FULL
  public :: euler_idp_state_observation_t, euler_idp_diagnostics_t
  public :: euler_idp_stage_time

  !> Configuration of the opt-in Euler GLL-IDP path.
  type, public :: euler_idp_config_t
     logical :: enabled = .false.
     logical :: low_order_only = .false.
     logical :: relax_density_bounds = .false.
     logical :: relax_entropy_bounds = .false.
     logical :: limit_internal_energy = .true.
     logical :: limit_entropy = .true.
     real(kind=rp) :: density_bound_relaxation_factor = 1.0_rp
     real(kind=rp) :: entropy_bound_relaxation_factor = 1.0_rp
     real(kind=rp) :: internal_energy_floor = 1.0e-12_rp
     real(kind=rp) :: correction_tolerance = 1.0e-10_rp
     integer :: diagnostics_level = EULER_IDP_DIAGNOSTICS_FULL
     integer :: diagnostics_interval = 1
   contains
     procedure, pass(this) :: init => euler_idp_config_init
     procedure, pass(this) :: valid => euler_idp_config_valid
  end type euler_idp_config_t

  !> Validity state for cached data derived from the current conserved state.
  type, public :: euler_idp_stage_t
     integer :: state_epoch = 0
     integer :: prepared_epoch = -1
     logical :: state_valid = .false.
     logical :: boundary_valid = .false.
     logical :: primitive_valid = .false.
     logical :: graph_valid = .false.
   contains
     procedure, pass(this) :: invalidate => euler_idp_stage_invalidate
     procedure, pass(this) :: state_changed => euler_idp_stage_state_changed
  end type euler_idp_stage_t

  !> Boundary metadata kept separately from graph topology.
  type, public :: euler_idp_boundary_t
     logical :: valid = .false.
     real(kind=dp) :: time = -huge(1.0_dp)
  end type euler_idp_boundary_t

  !> Backend-independent Euler GLL-IDP lifecycle and stage cache.
  type, public :: euler_idp_t
     logical :: initialized = .false.
     integer :: time_order = 1
     type(euler_idp_config_t) :: config
     type(euler_idp_stage_t) :: stage
     type(euler_idp_boundary_t) :: boundary
     type(euler_idp_diagnostics_t) :: diagnostics
     type(euler_idp_diagnostics_t) :: stage_diagnostics(3)
     class(euler_idp_backend_t), allocatable :: backend
   contains
     procedure, pass(this) :: init => euler_idp_init
     procedure, pass(this) :: validate_setup => euler_idp_validate_setup
     procedure, pass(this) :: prepare => euler_idp_prepare
     procedure, pass(this) :: cfl => euler_idp_cfl
     procedure, pass(this) :: advance => euler_idp_advance
     procedure, pass(this) :: apply_boundary => euler_idp_apply_boundary
     procedure, pass(this) :: invalidate => euler_idp_invalidate
     procedure, pass(this) :: copy_primitives => euler_idp_copy_primitives
     procedure, pass(this) :: free => euler_idp_free
  end type euler_idp_t

contains

  !> Read the Euler GLL-IDP configuration from a case file.
  subroutine euler_idp_config_init(this, params)
    class(euler_idp_config_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
    character(len=*), parameter :: root = 'euler_idp.'
    character(len=:), allocatable :: diagnostics_level

    this%enabled = .false.
    this%low_order_only = .false.
    this%relax_density_bounds = .false.
    this%relax_entropy_bounds = .false.
    this%limit_internal_energy = .true.
    this%limit_entropy = .true.
    this%density_bound_relaxation_factor = 1.0_rp
    this%entropy_bound_relaxation_factor = 1.0_rp
    this%internal_energy_floor = 1.0e-12_rp
    this%correction_tolerance = 1.0e-10_rp
    this%diagnostics_level = EULER_IDP_DIAGNOSTICS_FULL
    this%diagnostics_interval = 1

    if (.not. params%valid_path('euler_idp')) return

    call json_get_or_default(params, root // 'enabled', this%enabled, .false.)
    call json_get_or_default(params, root // 'low_order_only', &
         this%low_order_only, .false.)
    call json_get_or_default(params, root // 'relax_density_bounds', &
         this%relax_density_bounds, .false.)
    call json_get_or_default(params, root // 'relax_entropy_bounds', &
         this%relax_entropy_bounds, .false.)
    call json_get_or_default(params, root // 'limit_internal_energy', &
         this%limit_internal_energy, .true.)
    call json_get_or_default(params, root // 'limit_entropy', &
         this%limit_entropy, .true.)
    call json_get_or_default(params, &
         root // 'density_bound_relaxation_factor', &
         this%density_bound_relaxation_factor, 1.0_rp)
    call json_get_or_default(params, &
         root // 'entropy_bound_relaxation_factor', &
         this%entropy_bound_relaxation_factor, 1.0_rp)
    call json_get_or_default(params, root // 'internal_energy_floor', &
         this%internal_energy_floor, 1.0e-12_rp)
    call json_get_or_default(params, root // 'correction_tolerance', &
         this%correction_tolerance, 1.0e-10_rp)
    call json_get_or_default(params, root // 'diagnostics_level', &
         diagnostics_level, 'full')
    call json_get_or_default(params, root // 'diagnostics_interval', &
         this%diagnostics_interval, 1)

    select case (trim(diagnostics_level))
    case ('off')
       this%diagnostics_level = EULER_IDP_DIAGNOSTICS_OFF
    case ('safety')
       this%diagnostics_level = EULER_IDP_DIAGNOSTICS_SAFETY
    case ('full')
       this%diagnostics_level = EULER_IDP_DIAGNOSTICS_FULL
    case default
       this%diagnostics_level = -1
    end select
  end subroutine euler_idp_config_init

  !> Validate scalar Euler GLL-IDP configuration values.
  logical function euler_idp_config_valid(this, message) result(valid)
    class(euler_idp_config_t), intent(in) :: this
    character(len=*), intent(out) :: message

    valid = .false.
    message = ''

    if (.not. ieee_is_finite(this%internal_energy_floor) .or. &
         this%internal_energy_floor .le. 0.0_rp) then
       message = 'internal_energy_floor must be finite and positive'
       return
    end if
    if (.not. ieee_is_finite(this%density_bound_relaxation_factor) .or. &
         this%density_bound_relaxation_factor .le. 0.0_rp) then
       message = 'density_bound_relaxation_factor must be finite and positive'
       return
    end if
    if (.not. ieee_is_finite(this%entropy_bound_relaxation_factor) .or. &
         this%entropy_bound_relaxation_factor .le. 0.0_rp) then
       message = 'entropy_bound_relaxation_factor must be finite and positive'
       return
    end if
    if (.not. ieee_is_finite(this%correction_tolerance) .or. &
         this%correction_tolerance .le. 0.0_rp) then
       message = 'correction_tolerance must be finite and positive'
       return
    end if
    if (this%diagnostics_level .lt. EULER_IDP_DIAGNOSTICS_OFF .or. &
         this%diagnostics_level .gt. EULER_IDP_DIAGNOSTICS_FULL) then
       message = 'diagnostics_level must be off, safety, or full'
       return
    end if
    if (this%diagnostics_interval .lt. 1) then
       message = 'diagnostics_interval must be positive'
       return
    end if
    valid = .true.
  end function euler_idp_config_valid

  !> Initialize static topology, backend storage, and cache state.
  subroutine euler_idp_init(this, dof, coef, gs, config, time_order)
    class(euler_idp_t), intent(inout) :: this
    type(dofmap_t), target, intent(in) :: dof
    type(coef_t), target, intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    type(euler_idp_config_t), intent(in) :: config
    integer, intent(in) :: time_order
    integer :: i

    call this%free()
    call this%validate_setup(config, time_order)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       allocate(euler_idp_device_t :: this%backend)
    else
       allocate(euler_idp_cpu_t :: this%backend)
    end if
    call this%backend%init(dof)
    call this%backend%init_graph(coef, gs, config%relax_density_bounds, &
         config%relax_entropy_bounds, config%low_order_only, &
         config%limit_internal_energy, config%limit_entropy, &
         config%density_bound_relaxation_factor, &
         config%entropy_bound_relaxation_factor, time_order, &
         config%diagnostics_level, config%correction_tolerance)
    this%config = config
    this%time_order = time_order
    this%stage = euler_idp_stage_t()
    call this%stage%state_changed()
    this%boundary = euler_idp_boundary_t()
    call this%diagnostics%reset()
    do i = 1, 3
       call this%stage_diagnostics(i)%reset()
    end do
    this%initialized = .true.
  end subroutine euler_idp_init

  !> Validate backend-independent setup before allocating graph storage.
  subroutine euler_idp_validate_setup(this, config, time_order)
    class(euler_idp_t), intent(in) :: this
    type(euler_idp_config_t), intent(in) :: config
    integer, intent(in) :: time_order
    character(len=256) :: message

    if (.not. config%valid(message)) then
       call neko_error('Invalid Euler IDP configuration: ' // trim(message))
    end if
    if (time_order .ne. 1 .and. time_order .ne. 3) then
       call neko_error('Euler IDP requires Forward Euler or SSPRK3')
    end if
  end subroutine euler_idp_validate_setup

  !> Prepare the boundary-valid primitive and graph data for a state.
  subroutine euler_idp_prepare(this, rho, m_x, m_y, m_z, energy, gs, &
       density_bcs, velocity_bcs, pressure_bcs, gamma, time, &
       graph_wave_speed)
    class(euler_idp_t), intent(inout) :: this
    type(field_t), intent(inout) :: rho, m_x, m_y, m_z, energy
    type(gs_t), intent(inout) :: gs
    type(bc_list_t), intent(inout) :: density_bcs, velocity_bcs, pressure_bcs
    real(kind=rp), intent(in) :: gamma
    type(time_state_t), intent(in), optional :: time
    type(field_t), intent(in), optional :: graph_wave_speed
    logical :: needs_boundary

    call euler_idp_assert_initialized(this)
    needs_boundary = .not. this%stage%boundary_valid
    if (present(time)) then
       needs_boundary = needs_boundary .or. &
            .not. euler_idp_same_time(this%boundary%time, time%t)
    end if
    if (needs_boundary) then
       call this%apply_boundary(rho, m_x, m_y, m_z, energy, density_bcs, &
            velocity_bcs, pressure_bcs, gamma, time, &
            'prepared state after boundary conditions')
    end if
    if (.not. this%stage%graph_valid .or. &
         this%stage%prepared_epoch .ne. this%stage%state_epoch) then
       call this%backend%prepare_stage(rho, m_x, m_y, m_z, energy, gs, &
            gamma, this%config%internal_energy_floor, &
            this%stage%primitive_valid, graph_wave_speed)
       this%stage%primitive_valid = .true.
       this%stage%graph_valid = .true.
       this%stage%prepared_epoch = this%stage%state_epoch
    end if
  end subroutine euler_idp_prepare

  !> Return the graph CFL from the prepared stage cache.
  real(kind=rp) function euler_idp_cfl(this, dt) result(cfl)
    class(euler_idp_t), intent(in) :: this
    real(kind=rp), intent(in) :: dt

    call euler_idp_assert_initialized(this)
    if (.not. this%stage%graph_valid .or. &
         this%stage%prepared_epoch .ne. this%stage%state_epoch) then
       call neko_error('Euler IDP CFL requested before stage preparation')
    end if
    cfl = this%backend%graph_cfl(dt)
  end function euler_idp_cfl

  !> Advance one complete limited timestep through the selected backend.
  subroutine euler_idp_advance(this, rho, m_x, m_y, m_z, energy, coef, gs, &
       density_bcs, velocity_bcs, pressure_bcs, gamma, dt, time, &
       diagnostics, entropy_viscosity_fraction, graph_wave_speed)
    class(euler_idp_t), intent(inout) :: this
    type(field_t), intent(inout) :: rho, m_x, m_y, m_z, energy
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    type(bc_list_t), intent(inout) :: density_bcs, velocity_bcs, pressure_bcs
    real(kind=rp), intent(in) :: gamma, dt
    type(time_state_t), intent(in) :: time
    type(euler_idp_diagnostics_t), intent(out) :: diagnostics
    type(field_t), intent(in), optional :: entropy_viscosity_fraction
    type(field_t), intent(in), optional :: graph_wave_speed
    type(time_state_t) :: stage_time
    integer :: i

    call euler_idp_assert_initialized(this)
    do i = 1, 3
       call this%stage_diagnostics(i)%reset()
    end do

    stage_time = time
    stage_time%t = euler_idp_stage_time(time%t, dt, this%time_order, 1)
    call euler_idp_prepare_stage(this, rho, m_x, m_y, m_z, energy, gs, &
         density_bcs, velocity_bcs, pressure_bcs, gamma, stage_time, &
         graph_wave_speed, 'stage input after boundary conditions')
    if (this%time_order .eq. 3) then
       call this%backend%save_state(rho, m_x, m_y, m_z, energy)
    end if
    call euler_idp_forward_euler_stage(this, rho, m_x, m_y, m_z, energy, &
         coef, gs, density_bcs, velocity_bcs, pressure_bcs, gamma, dt, &
         stage_time, 1, this%stage_diagnostics(1), &
         entropy_viscosity_fraction, graph_wave_speed)

    if (this%time_order .eq. 1) then
       call this%backend%combine_stage(rho, m_x, m_y, m_z, energy, &
            0.0_rp, 1.0_rp)
       call this%stage%state_changed()
       this%boundary = euler_idp_boundary_t()
       call this%apply_boundary(rho, m_x, m_y, m_z, energy, density_bcs, &
            velocity_bcs, pressure_bcs, gamma, time, &
            'Forward Euler final state after boundary conditions')
       diagnostics = this%stage_diagnostics(1)
       this%diagnostics = diagnostics
       return
    end if

    ! U(1) = FE(U(n)).
    call this%backend%combine_stage(rho, m_x, m_y, m_z, energy, &
         0.0_rp, 1.0_rp)
    call this%stage%state_changed()
    this%boundary = euler_idp_boundary_t()

    stage_time = time
    stage_time%t = euler_idp_stage_time(time%t, dt, this%time_order, 2)
    call euler_idp_prepare_stage(this, rho, m_x, m_y, m_z, energy, gs, &
         density_bcs, velocity_bcs, pressure_bcs, gamma, stage_time, &
         graph_wave_speed, 'SSPRK3 stage 2 input after boundary conditions')
    call euler_idp_forward_euler_stage(this, rho, m_x, m_y, m_z, energy, &
         coef, gs, density_bcs, velocity_bcs, pressure_bcs, gamma, dt, &
         stage_time, 2, this%stage_diagnostics(2), &
         entropy_viscosity_fraction, graph_wave_speed)

    ! U(2) = 3/4 U(n) + 1/4 FE(U(1)).
    call this%backend%combine_stage(rho, m_x, m_y, m_z, energy, &
         0.75_rp, 0.25_rp)
    call this%stage%state_changed()
    this%boundary = euler_idp_boundary_t()

    stage_time = time
    stage_time%t = euler_idp_stage_time(time%t, dt, this%time_order, 3)
    call euler_idp_prepare_stage(this, rho, m_x, m_y, m_z, energy, gs, &
         density_bcs, velocity_bcs, pressure_bcs, gamma, stage_time, &
         graph_wave_speed, 'SSPRK3 stage 3 input after boundary conditions')
    call euler_idp_forward_euler_stage(this, rho, m_x, m_y, m_z, energy, &
         coef, gs, density_bcs, velocity_bcs, pressure_bcs, gamma, dt, &
         stage_time, 3, this%stage_diagnostics(3), &
         entropy_viscosity_fraction, graph_wave_speed)

    ! U(n+1) = 1/3 U(n) + 2/3 FE(U(2)).
    call this%backend%combine_stage(rho, m_x, m_y, m_z, energy, &
         1.0_rp / 3.0_rp, 2.0_rp / 3.0_rp)
    call this%stage%state_changed()
    this%boundary = euler_idp_boundary_t()
    call this%apply_boundary(rho, m_x, m_y, m_z, energy, density_bcs, &
         velocity_bcs, pressure_bcs, gamma, time, &
         'SSPRK3 final state after boundary conditions')
    diagnostics = this%stage_diagnostics(3)
    this%diagnostics = diagnostics
  end subroutine euler_idp_advance

  !> Prepare one boundary-valid state, reusing matching cached graph data.
  subroutine euler_idp_prepare_stage(this, rho, m_x, m_y, m_z, energy, gs, &
       density_bcs, velocity_bcs, pressure_bcs, gamma, time, graph_wave_speed, &
       label)
    class(euler_idp_t), intent(inout) :: this
    type(field_t), intent(inout) :: rho, m_x, m_y, m_z, energy
    type(gs_t), intent(inout) :: gs
    type(bc_list_t), intent(inout) :: density_bcs, velocity_bcs, pressure_bcs
    real(kind=rp), intent(in) :: gamma
    type(time_state_t), intent(in) :: time
    type(field_t), intent(in), optional :: graph_wave_speed
    character(len=*), intent(in) :: label

    if (.not. this%stage%boundary_valid .or. &
         .not. euler_idp_same_time(this%boundary%time, time%t)) then
       call this%apply_boundary(rho, m_x, m_y, m_z, energy, density_bcs, &
            velocity_bcs, pressure_bcs, gamma, time, label)
    end if
    if (.not. this%stage%graph_valid .or. &
         this%stage%prepared_epoch .ne. this%stage%state_epoch) then
       call this%backend%prepare_stage(rho, m_x, m_y, m_z, energy, gs, &
            gamma, this%config%internal_energy_floor, &
            this%stage%primitive_valid, graph_wave_speed)
       this%stage%primitive_valid = .true.
       this%stage%graph_valid = .true.
       this%stage%prepared_epoch = this%stage%state_epoch
    end if
  end subroutine euler_idp_prepare_stage

  !> Execute one complete limited map with boundary composition explicit.
  subroutine euler_idp_forward_euler_stage(this, rho, m_x, m_y, m_z, &
       energy, coef, gs, density_bcs, velocity_bcs, pressure_bcs, gamma, dt, &
       time, stage, diagnostics, entropy_viscosity_fraction, graph_wave_speed)
    class(euler_idp_t), intent(inout) :: this
    type(field_t), intent(in) :: rho, m_x, m_y, m_z, energy
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    type(bc_list_t), intent(inout) :: density_bcs, velocity_bcs, pressure_bcs
    real(kind=rp), intent(in) :: gamma, dt
    type(time_state_t), intent(in) :: time
    integer, intent(in) :: stage
    type(euler_idp_diagnostics_t), intent(inout) :: diagnostics
    type(field_t), intent(in), optional :: entropy_viscosity_fraction
    type(field_t), intent(in), optional :: graph_wave_speed

    call this%backend%forward_euler(rho, m_x, m_y, m_z, energy, coef, gs, &
         gamma, this%config%internal_energy_floor, dt, time, diagnostics, &
         stage, entropy_viscosity_fraction, graph_wave_speed)
    if (this%config%diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL) then
       call this%backend%observe_candidate(gs, gamma, &
            this%config%internal_energy_floor, time, &
            'candidate before boundary conditions', &
            .false., diagnostics%before_boundary, graph_wave_speed)
       call this%backend%validate_candidate(gamma, diagnostics, &
            'limited candidate before boundary conditions')
    end if
    call this%backend%apply_candidate_boundary(density_bcs, velocity_bcs, &
         pressure_bcs, gamma, this%config%internal_energy_floor, time, &
         'limited candidate after boundary conditions')
    if (this%config%diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL) then
       call this%backend%observe_candidate(gs, gamma, &
            this%config%internal_energy_floor, time, &
            'candidate after boundary conditions', .true., &
            diagnostics%after_boundary, graph_wave_speed)
    end if
    call this%backend%validate_candidate(gamma, diagnostics, &
         'limited candidate after boundary conditions')
  end subroutine euler_idp_forward_euler_stage

  !> Apply the strong boundary map without changing graph topology.
  subroutine euler_idp_apply_boundary(this, rho, m_x, m_y, m_z, energy, &
       density_bcs, velocity_bcs, pressure_bcs, gamma, time, label)
    class(euler_idp_t), intent(inout) :: this
    type(field_t), intent(inout) :: rho, m_x, m_y, m_z, energy
    type(bc_list_t), intent(inout) :: density_bcs, velocity_bcs, pressure_bcs
    real(kind=rp), intent(in) :: gamma
    type(time_state_t), intent(in), optional :: time
    character(len=*), intent(in) :: label

    call euler_idp_assert_initialized(this)
    call this%backend%apply_boundary_conditions(rho, m_x, m_y, m_z, energy, &
         density_bcs, velocity_bcs, pressure_bcs, gamma, &
         this%config%internal_energy_floor, time, label, refresh = .true.)
    this%stage%boundary_valid = .true.
    this%stage%primitive_valid = .true.
    this%stage%graph_valid = .false.
    this%boundary%valid = .true.
    if (present(time)) this%boundary%time = time%t
  end subroutine euler_idp_apply_boundary

  !> Invalidate data derived from state changed outside this object.
  subroutine euler_idp_invalidate(this)
    class(euler_idp_t), intent(inout) :: this

    call this%stage%state_changed()
    this%boundary = euler_idp_boundary_t()
  end subroutine euler_idp_invalidate

  !> Copy the backend primitive cache into scheme-owned fields.
  subroutine euler_idp_copy_primitives(this, u, v, w, p)
    class(euler_idp_t), intent(in) :: this
    type(field_t), intent(inout) :: u, v, w, p

    call euler_idp_assert_initialized(this)
    if (.not. this%stage%primitive_valid) then
       call neko_error('Euler IDP primitive cache is not valid')
    end if
    call this%backend%copy_primitives(u, v, w, p)
  end subroutine euler_idp_copy_primitives

  !> Release backend resources; safe after incomplete initialization.
  subroutine euler_idp_free(this)
    class(euler_idp_t), intent(inout) :: this

    if (allocated(this%backend)) then
       call this%backend%free()
       deallocate(this%backend)
    end if
    this%initialized = .false.
    this%time_order = 1
    this%stage = euler_idp_stage_t()
    this%boundary = euler_idp_boundary_t()
  end subroutine euler_idp_free

  !> Invalidate every cache without advancing the state epoch.
  subroutine euler_idp_stage_invalidate(this)
    class(euler_idp_stage_t), intent(inout) :: this

    this%prepared_epoch = -1
    this%state_valid = .false.
    this%boundary_valid = .false.
    this%primitive_valid = .false.
    this%graph_valid = .false.
  end subroutine euler_idp_stage_invalidate

  !> Record a mutation of the conserved state.
  subroutine euler_idp_stage_state_changed(this)
    class(euler_idp_stage_t), intent(inout) :: this

    this%state_epoch = this%state_epoch + 1
    call this%invalidate()
    this%state_valid = .true.
  end subroutine euler_idp_stage_state_changed

  !> Require a fully initialized solver object.
  subroutine euler_idp_assert_initialized(this)
    class(euler_idp_t), intent(in) :: this

    if (.not. this%initialized .or. .not. allocated(this%backend)) then
       call neko_error('Euler IDP object is not initialised')
    end if
  end subroutine euler_idp_assert_initialized

  !> Compare boundary times with a tolerance scaled to the physical time.
  pure logical function euler_idp_same_time(left, right) result(same)
    real(kind=dp), intent(in) :: left, right
    real(kind=dp) :: scale

    scale = max(1.0_dp, abs(left), abs(right))
    same = abs(left - right) .le. 64.0_dp * epsilon(1.0_dp) * scale
  end function euler_idp_same_time

end module euler_idp
