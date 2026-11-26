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
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!> Implements `lpt_t`. (Lagrangian Particle Tracking)
module lpt
  use num_types, only : rp, dp
  use json_module, only : json_file
  use registry, only : neko_registry
  use field, only : field_t
  use case, only : case_t
  use mesh, only : mesh_t
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use json_utils, only : json_get, json_get_or_default, &
       json_get_subdict_or_empty
  use time_state, only : time_state_t
  use time_based_controller, only : time_based_controller_t
  use global_interpolation, only : global_interpolation_t
  use logger, only : neko_log, LOG_SIZE
  use utils, only : neko_error
  use file, only : file_t
  use matrix, only : matrix_t
  use vector_math, only : vector_add2s2, vector_cfill, vector_cmult2, &
       vector_col2, vector_col3, vector_invcol2, vector_sqrt_inplace, &
       vector_power, vector_sub3, vector_vdot3, vector_cmult, &
       vector_cadd
  use ab_time_scheme, only : ab_time_scheme_t
  use lpt_periodic_bc, only : lpt_periodic_bc_t
  use lpt_migrate, only : lpt_migrate_t, LPT_MIGRATE_TO_OWNER, &
       LPT_MIGRATE_NONE
  use lpt_output, only : lpt_output_t
  use comm, only : pe_rank, NEKO_COMM
  use mpi_f08, only : MPI_Allreduce, MPI_INTEGER, MPI_SUM
  use csv_file, only : csv_file_t
  use vector, only : vector_t, vector_ptr_t
  use particles, only : particles_t
  use device, only : DEVICE_TO_HOST
  use checkpoint_payload, only : checkpoint_payload_t, checkpoint_hook_t, &
       checkpoint_array_t, checkpoint_array_ptr_t
  use profiler, only : profiler_start_region, profiler_end_region
  use scratch_registry, only : neko_scratch_registry
  implicit none
  private

  !> Number of entries in the replicated LPT checkpoint state array:
  !! initialised flag, history length, time step, t, dt, tlag(10), dtlag(10).
  integer, parameter :: LPT_CHKP_STATE_SIZE = 25

  !> Pre-write checkpoint hook packing the particle state of an LPT instance.
  type, extends(checkpoint_hook_t) :: lpt_checkpoint_hook_t
     !> LPT instance owning the hook.
     type(lpt_t), pointer :: lpt => null()
   contains
     !> Pack the particle state into the checkpoint buffers.
     procedure, pass(this) :: prepare_write => lpt_checkpoint_hook_prepare_write
  end type lpt_checkpoint_hook_t

  !> Passive Lagrangian particle tracking.
  type, public :: lpt_t
     !> LPT instance name.
     character(:), allocatable :: name
     ! Fields based on the fluid solution space
     type(field_t), pointer :: u_field => null()
     type(field_t), pointer :: v_field => null()
     type(field_t), pointer :: w_field => null()
     type(field_t), pointer :: mu_fluid => null()
     type(field_t), pointer :: rho_fluid => null()
     type(mesh_t), pointer :: msh => null()
     type(dofmap_t), pointer :: dm_Xh => null()
     type(coef_t), pointer :: coef => null()
     integer :: time_order, lag_len
     integer :: history_len = 0
     logical :: inertia = .false.
     real(kind=rp) :: nonlinear_coefficient, nonlinear_exponent
     logical :: elastic_wall_enabled = .false.
     integer, allocatable :: wall_zone_indices(:)
     logical, allocatable :: wall_facet_mask(:, :)
     type(time_state_t) :: lpt_time
     logical :: lpt_time_initialized = .false.
     type(global_interpolation_t) :: global_interp
     type(lpt_periodic_bc_t) :: periodic_bc
     type(lpt_migrate_t) :: migration
     type(particles_t) :: particles
     type(lpt_output_t) :: output
     type(time_based_controller_t) :: output_controller
     logical :: output_enabled = .false.
     logical :: log = .true.
     real(kind=rp) :: start_time = -huge(0.0_rp)
     !> Hook packing the particle state before a checkpoint is written.
     type(lpt_checkpoint_hook_t) :: chkp_hook
     !> Checkpoint payload holding the particle state.
     type(checkpoint_payload_t), pointer :: chkp_payload => null()
     !> Particle vectors stored in the checkpoint.
     type(vector_ptr_t), allocatable :: chkp_vectors(:)
     !> Checkpoint arrays paired with `chkp_vectors`.
     type(checkpoint_array_ptr_t), allocatable :: chkp_arrays(:)
     !> Checkpoint array holding the global particle ids.
     type(checkpoint_array_t), pointer :: chkp_ids => null()
     !> Replicated LPT time-integration state stored in the checkpoint.
     real(kind=dp) :: chkp_state(LPT_CHKP_STATE_SIZE) = 0.0_dp
     !> Serialised LPT configuration kept on restart, where seeding is
     !! deferred until it is known whether the checkpoint holds particles.
     character(len=:), allocatable :: seed_params
     !> Offset added to the solver time step, which starts over at a
     !! restart, so that the LPT step count continues across restarts.
     integer :: tstep_offset = 0
   contains
     procedure, pass(this) :: init => lpt_init_from_json
     procedure, pass(this) :: free => lpt_free
     procedure, pass(this) :: preprocess => lpt_preprocess
     procedure, pass(this) :: compute => lpt_compute
     procedure, pass(this) :: restart => lpt_restart
     procedure, pass(this) :: pack_checkpoint => lpt_pack_checkpoint
     procedure, private, pass(this) :: register_checkpoint
     procedure, private, pass(this) :: seed_particles
     procedure, private, pass(this) :: read_particles_json
     procedure, private, pass(this) :: read_particles_csv
     procedure, private, pass(this) :: evaluate_velocity
     procedure, private, pass(this) :: evaluate_acceleration
     procedure, private, pass(this) :: sync_time_controller
     procedure, private, pass(this) :: ODE_integrate_ab_3c
     procedure, private, pass(this) :: update_current_rhs
     procedure, private, pass(this) :: write_output
     procedure, private, pass(this) :: log_status
  end type lpt_t
  private :: update_lags

  interface
     !> Build a mask of elastic wall facets from configured mesh zone ids.
     module subroutine lpt_init_wall_facet_mask(wall_facet_mask, msh, &
          wall_zone_indices)
       logical, allocatable, intent(inout) :: wall_facet_mask(:, :)
       type(mesh_t), intent(in) :: msh
       integer, intent(in) :: wall_zone_indices(:)
     end subroutine lpt_init_wall_facet_mask

     !> Reflect particles and histories after elastic wall contact.
     module subroutine lpt_handle_elastic_wall_collisions(this, x_old, y_old, &
          z_old, u_old, v_old, w_old)
       class(lpt_t), intent(inout) :: this
       type(vector_t), intent(in) :: x_old, y_old, z_old
       type(vector_t), intent(inout) :: u_old, v_old, w_old
     end subroutine lpt_handle_elastic_wall_collisions
  end interface

contains

  !> Construct from JSON.
  !! Supported particle input is either a flat `coordinates` array
  !! `[x1,y1,z1,x2,y2,z2,...]` or a three-column CSV `points_file`.
  subroutine lpt_init_from_json(this, json, case)
    class(lpt_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    type(json_file) :: interp_subdict
    character(len=:), allocatable :: name
    character(len=:), allocatable :: migration_strategy
    character(len=:), allocatable :: output_filename
    character(len=:), allocatable :: output_format
    character(len=:), allocatable :: snapshots_per_file_str
    character(len=:), allocatable :: output_path
    integer :: migration_strategy_id
    integer :: snapshots_per_file
    integer :: snapshots_per_file_type
    logical :: snapshots_per_file_found
    character(len=:), allocatable :: restart_file
    logical :: restarting

    call this%free()

    call json_get_or_default(json, "name", name, "lpt")
    call json_get_or_default(json, "log", this%log, .true.)
    call json_get_or_default(json, "start_time", this%start_time, &
         -huge(0.0_rp))

    this%name = name
    this%time_order = case%fluid%ext_bdf%advection_time_order
    this%msh => case%fluid%msh
    this%dm_Xh => case%fluid%dm_Xh
    this%coef => case%fluid%c_Xh

    this%lag_len = this%time_order - 1
    call json_get_or_default(json, "migration_strategy", migration_strategy, &
         "owner")
    select case (trim(migration_strategy))
    case ("owner")
       migration_strategy_id = LPT_MIGRATE_TO_OWNER
    case ("none")
       migration_strategy_id = LPT_MIGRATE_NONE
    case default
       call neko_error("lpt migration_strategy must be 'owner' or 'none'")
    end select
    call this%migration%init(this%lag_len, migration_strategy_id)

    call json_get(json, "inertia", this%inertia)

    if (this%inertia) then
       call json_get_or_default(json, "nonlinear_coefficient", &
            this%nonlinear_coefficient, 0.15_rp)
       call json_get_or_default(json, "nonlinear_exponent", &
            this%nonlinear_exponent, 0.687_rp)
       this%mu_fluid => neko_registry%get_field_by_name( &
            case%fluid%name // "_mu")
       this%rho_fluid => neko_registry%get_field_by_name( &
            case%fluid%name // "_rho")
    end if
    this%u_field => neko_registry%get_field_by_name("u")
    this%v_field => neko_registry%get_field_by_name("v")
    this%w_field => neko_registry%get_field_by_name("w")

    if (case%params%valid_path("case.fluid.wall_zone_indices")) then
       call json_get(case%params, "case.fluid.wall_zone_indices", &
            this%wall_zone_indices)
       if (.not. this%inertia .and. size(this%wall_zone_indices) .gt. 0) then
          call neko_error("lpt wall_zone_indices requires inertia = true")
       end if
       this%elastic_wall_enabled = size(this%wall_zone_indices) .gt. 0
       if (this%elastic_wall_enabled .and. &
            migration_strategy_id .eq. LPT_MIGRATE_NONE) then
          call neko_error("lpt migration_strategy = none is not " // &
               "compatible with elastic wall collisions")
       end if
       if (this%elastic_wall_enabled) then
          call lpt_init_wall_facet_mask(this%wall_facet_mask, this%msh, &
               this%wall_zone_indices)
       end if
    end if

    ! On restart the particles come from the checkpoint, so the initial
    ! condition is only read if the checkpoint turns out to hold none.
    restarting = .false.
    if (case%params%valid_path("case.restart_file")) then
       call json_get(case%params, "case.restart_file", restart_file)
       restarting = len_trim(restart_file) .gt. 0
    end if
    if (restarting) then
       call json%print_to_string(this%seed_params)
    else
       call this%seed_particles(json)
    end if

    call json_get_subdict_or_empty(json, "interpolation", interp_subdict)
    call this%global_interp%init(case%fluid%dm_Xh, &
         params_subdict = interp_subdict)
    call this%periodic_bc%init(case%fluid%msh, case%fluid%dm_Xh, &
         case%fluid%c_Xh)
    call this%migration%migrate_particles(this%global_interp, &
         this%periodic_bc, this%inertia, this%particles)
    call this%sync_time_controller(case%time)
    call this%update_current_rhs()
    call this%register_checkpoint(case)

    call json_get_or_default(json, "output_filename", output_filename, &
         trim(this%name))
    call json_get_or_default(json, "output_format", output_format, "csv")

    call json%info("snapshots_per_file", found = snapshots_per_file_found, &
         var_type = snapshots_per_file_type)
    if (snapshots_per_file_found) then
       select case (snapshots_per_file_type)
       case (5)
          call json_get(json, "snapshots_per_file", snapshots_per_file)
          if (snapshots_per_file .lt. 1) then
             call neko_error("lpt snapshots_per_file must be a positive " // &
                  "integer or 'all'")
          end if
       case (7)
          call json_get(json, "snapshots_per_file", snapshots_per_file_str)
          if (trim(snapshots_per_file_str) .eq. "all") then
             snapshots_per_file = 0
          else
             call neko_error("lpt snapshots_per_file must be a positive " // &
                  "integer or 'all'")
          end if
       case default
          call neko_error("lpt snapshots_per_file must be a positive " // &
               "integer or 'all'")
       end select
    else
       snapshots_per_file = 0
       call json%add("snapshots_per_file", "all")
    end if
    output_path = case%output_directory // trim(output_filename) // "." // &
         trim(output_format)
    ! On restart the output files are opened in `lpt_restart`, which
    ! continues those of the previous run.
    call this%output%init(output_path, this%inertia, snapshots_per_file, &
         defer_open = restarting)

    ! Output at the initialisation. On restart the restored particles
    ! continue the previous trajectory, and seeded ones are written in
    ! `lpt_restart` at the restart time.
    this%output_enabled = .true.
    if (.not. restarting) call this%write_output(case%time)

    call this%log_status()
  end subroutine lpt_init_from_json

  !> Seed the particles from the initial condition in the configuration.
  !! @param json LPT configuration.
  subroutine seed_particles(this, json)
    class(lpt_t), intent(inout) :: this
    type(json_file), intent(inout) :: json

    call this%read_particles_json(json)
    call this%migration%initialize_particle_distribution(this%inertia, &
         this%particles)
  end subroutine seed_particles

  !> Read particle data from JSON and initialise particles on rank 0.
  !! @param json LPT configuration containing either `coordinates` or
  !! `points_file`.
  subroutine read_particles_json(this, json)
    class(lpt_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    real(kind=rp), allocatable :: coords(:)
    real(kind=rp), allocatable :: vels(:)
    real(kind=rp), allocatable :: diams(:)
    real(kind=rp), allocatable :: densities(:)
    real(kind=rp), pointer, dimension(:) :: x, y, z, u, v, w
    integer :: n_particles, ind(6)

    if (pe_rank .eq. 0) then
       if (json%valid_path("coordinates")) then
          call json_get(json, "coordinates", coords)
          if (mod(size(coords), 3) .ne. 0) then
             call neko_error("lpt coordinates must contain 3 values per " // &
                  "particle")
          end if
          n_particles = size(coords) / 3
          if (this%inertia) then
             call json_get(json, "velocities", vels)
             if (mod(size(vels), 3) .ne. 0) then
                call neko_error("lpt velocities must contain 3 values per " // &
                     "particle")
             end if
             call json_get(json, "diameters", diams)
             call json_get(json, "densities", densities)
             if (size(vels) / 3 .ne. n_particles .or. &
                  size(diams) .ne. n_particles .or. &
                  size(densities) .ne. n_particles) then
                call neko_error("lpt coordinates, velocities, diameters " // &
                     "and densities must describe the same number of " // &
                     "particles")
             end if
          else
             allocate(vels(size(coords)))
             allocate(diams(n_particles))
             allocate(densities(n_particles))
             vels = 0.0_rp
             diams = 0.0_rp
             densities = 0.0_rp
          end if
          call neko_scratch_registry%request(x, ind(1), n_particles, .false.)
          call neko_scratch_registry%request(y, ind(2), n_particles, .false.)
          call neko_scratch_registry%request(z, ind(3), n_particles, .false.)
          call neko_scratch_registry%request(u, ind(4), n_particles, .false.)
          call neko_scratch_registry%request(v, ind(5), n_particles, .false.)
          call neko_scratch_registry%request(w, ind(6), n_particles, .false.)
          x = coords(1::3)
          y = coords(2::3)
          z = coords(3::3)
          u = vels(1::3)
          v = vels(2::3)
          w = vels(3::3)
          call this%particles%init(x, y, z, this%time_order, u, v, &
               w, diams, densities)
          deallocate(coords)
          deallocate(vels)
          deallocate(diams)
          deallocate(densities)
          call neko_scratch_registry%relinquish(ind)

       else if (json%valid_path("points_file")) then
          call this%read_particles_csv(json)
       else
          call neko_error("lpt requires either coordinates or points_file")
       end if
    else
       return
    end if
  end subroutine read_particles_json

  !> Read particle data from a CSV file and initialise particles on rank 0.
  !! @param json LPT configuration containing `points_file`.
  subroutine read_particles_csv(this, json)
    class(lpt_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: points_file
    type(file_t) :: file_in
    type(matrix_t) :: mat_in
    real(kind=rp), pointer, dimension(:) :: x, y, z, u, v, w
    real(kind=rp), allocatable :: diams(:)
    real(kind=rp), allocatable :: densities(:)
    integer :: n_particles, ind_basic(3), ind_inertia(6)

    if (pe_rank .ne. 0) return

    call json_get(json, "points_file", points_file)
    call file_in%init(trim(points_file))

    select type (ft => file_in%file_type)
    type is (csv_file_t)
       if (this%inertia) then
          call mat_in%init(ft%count_lines(), 8)
          call ft%read(mat_in)
          n_particles = mat_in%get_nrows()
          call neko_scratch_registry%request(x, ind_inertia(1), n_particles, &
               .false.)
          call neko_scratch_registry%request(y, ind_inertia(2), n_particles, &
               .false.)
          call neko_scratch_registry%request(z, ind_inertia(3), n_particles, &
               .false.)
          call neko_scratch_registry%request(u, ind_inertia(4), n_particles, &
               .false.)
          call neko_scratch_registry%request(v, ind_inertia(5), n_particles, &
               .false.)
          call neko_scratch_registry%request(w, ind_inertia(6), n_particles, &
               .false.)
          x = mat_in%x(:, 1)
          y = mat_in%x(:, 2)
          z = mat_in%x(:, 3)
          u = mat_in%x(:, 4)
          v = mat_in%x(:, 5)
          w = mat_in%x(:, 6)
          diams = mat_in%x(:, 7)
          densities = mat_in%x(:, 8)
          call this%particles%init(x, y, z, this%time_order, u, v, &
               w, diams, densities)
          deallocate(diams)
          deallocate(densities)
          call neko_scratch_registry%relinquish(ind_inertia)
       else
          call mat_in%init(ft%count_lines(), 3)
          call ft%read(mat_in)
          n_particles = mat_in%get_nrows()
          call neko_scratch_registry%request(x, ind_basic(1), n_particles, .false.)
          call neko_scratch_registry%request(y, ind_basic(2), n_particles, .false.)
          call neko_scratch_registry%request(z, ind_basic(3), n_particles, .false.)
          x = mat_in%x(:, 1)
          y = mat_in%x(:, 2)
          z = mat_in%x(:, 3)
          call this%particles%init(x, y, z, this%time_order)
          call neko_scratch_registry%relinquish(ind_basic)
       end if
    class default
       call neko_error("lpt points_file must be a csv file")
    end select
    call mat_in%free()
    call file_in%free()
  end subroutine read_particles_csv

  !> Interpolate the carrier velocity at the local particles.
  !! @param u_fluid Interpolated x velocity at particle locations.
  !! @param v_fluid Interpolated y velocity at particle locations.
  !! @param w_fluid Interpolated z velocity at particle locations.
  subroutine evaluate_velocity(this, u_fluid, v_fluid, w_fluid)
    class(lpt_t), intent(inout) :: this
    type(vector_t), intent(inout) :: u_fluid, v_fluid, w_fluid
    logical :: do_interp_on_host

    if (this%particles%n .eq. 0) return

    do_interp_on_host = .false.
    call this%global_interp%evaluate(u_fluid%x, this%u_field%x, &
         do_interp_on_host)
    call this%global_interp%evaluate(v_fluid%x, this%v_field%x, &
         do_interp_on_host)
    call this%global_interp%evaluate(w_fluid%x, this%w_field%x, &
         do_interp_on_host)

  end subroutine evaluate_velocity

  !> Estimate particle acceleration from local carrier-fluid velocity.
  !! @param acc_x Particle acceleration in x.
  !! @param acc_y Particle acceleration in y.
  !! @param acc_z Particle acceleration in z.
  !! @param u_fluid Interpolated x velocity at particle locations.
  !! @param v_fluid Interpolated y velocity at particle locations.
  !! @param w_fluid Interpolated z velocity at particle locations.
  subroutine evaluate_acceleration(this, acc_x, acc_y, acc_z, &
       u_fluid, v_fluid, w_fluid)
    class(lpt_t), intent(inout) :: this
    type(vector_t), intent(in) :: u_fluid, v_fluid, w_fluid
    type(vector_t), intent(inout) :: acc_x, acc_y, acc_z
    type(vector_t), pointer :: tau_p, Re_p, f, rho_fluid_local
    type(vector_t), pointer :: mu_fluid_local, nu_fluid_local
    integer :: ind(5)

    integer :: n
    logical :: do_interp_on_host

    if (this%particles%n .eq. 0) return
    n = this%particles%n

    ! Request the scratch storage used throughout the computation.
    call neko_scratch_registry%request(rho_fluid_local, &
         ind(1), n, .false.)
    call neko_scratch_registry%request(mu_fluid_local, &
         ind(2), n, .false.)
    call neko_scratch_registry%request(tau_p, &
         ind(3), n, .false.)
    call neko_scratch_registry%request(Re_p, &
         ind(4), n, .false.)
    call neko_scratch_registry%request(f, &
         ind(5), n, .false.)

    ! Compute the local fluid properties and particle time scale.
    do_interp_on_host = .false.
    call this%global_interp%evaluate(mu_fluid_local%x, this%mu_fluid%x, &
         do_interp_on_host)
    call this%global_interp%evaluate(rho_fluid_local%x, this%rho_fluid%x, &
         do_interp_on_host)

    ! compute the time scale
    call vector_cfill(tau_p, 1.0_rp/18.0_rp)
    call vector_col2(tau_p, this%particles%rho)
    call vector_invcol2(tau_p, mu_fluid_local)
    call vector_col2(tau_p, this%particles%d)
    call vector_col2(tau_p, this%particles%d)

    ! The dynamic viscosity is no longer needed, so reuse its storage for the
    ! kinematic viscosity.
    call vector_invcol2(mu_fluid_local, rho_fluid_local)
    ! reuse tempoerary array but with a different name for clarity
    nu_fluid_local => mu_fluid_local

    ! Compute the relative velocity and particle Reynolds number.
    ! now use acc_xyz as work arrays
    call vector_sub3(acc_x, u_fluid, this%particles%u)
    call vector_sub3(acc_y, v_fluid, this%particles%v)
    call vector_sub3(acc_z, w_fluid, this%particles%w)
    call vector_vdot3(Re_p, acc_x, acc_y, acc_z, acc_x, acc_y, acc_z)
    call vector_sqrt_inplace(Re_p)
    call vector_col2(Re_p, this%particles%d)
    call vector_invcol2(Re_p, nu_fluid_local)

    ! Compute the nonlinear drag correction.
    call vector_power(f, Re_p, this%nonlinear_exponent)
    call vector_cmult(f, this%nonlinear_coefficient)
    call vector_cadd(f, 1.0_rp)

    ! Assemble the particle acceleration.
    call vector_col2(acc_x, f)
    call vector_col2(acc_y, f)
    call vector_col2(acc_z, f)
    call vector_invcol2(acc_x, tau_p)
    call vector_invcol2(acc_y, tau_p)
    call vector_invcol2(acc_z, tau_p)

    call neko_scratch_registry%relinquish(ind)

  end subroutine evaluate_acceleration

  !> Refresh particle RHS values using the current fluid solution.
  subroutine update_current_rhs(this)
    class(lpt_t), intent(inout) :: this
    type(vector_t), pointer :: u_fluid, v_fluid, w_fluid
    integer :: ind(3)

    call profiler_start_region('LPT_migrate_interp')

    call this%migration%migrate_particles(this%global_interp, &
         this%periodic_bc, this%inertia, this%particles)

    call neko_scratch_registry%request(u_fluid, ind(1), &
         this%particles%n, .false.)
    call neko_scratch_registry%request(v_fluid, ind(2), &
         this%particles%n, .false.)
    call neko_scratch_registry%request(w_fluid, ind(3), &
         this%particles%n, .false.)

    call this%evaluate_velocity(u_fluid, v_fluid, w_fluid)

    if (this%inertia) then
       call this%evaluate_acceleration(this%particles%acc_x, &
            this%particles%acc_y, this%particles%acc_z, &
            u_fluid, v_fluid, w_fluid)
    else
       this%particles%u = u_fluid
       this%particles%v = v_fluid
       this%particles%w = w_fluid
    end if

    call neko_scratch_registry%relinquish(ind)

    call profiler_end_region('LPT_migrate_interp')
  end subroutine update_current_rhs

  !> Shift one particle history level and store new current values.
  !! @param lag Previous-step values.
  !! @param laglag Two-steps-back values.
  !! @param new_values Values to store in `lag`.
  subroutine update_lags(lag, laglag, new_values)
    type(vector_t), intent(inout) :: lag, laglag
    type(vector_t), intent(in) :: new_values

    laglag = lag
    lag = new_values

  end subroutine update_lags

  !> Advance particle positions and, for inertial particles, velocities.
  !! @param time Current simulation time used to build the LPT local step.
  subroutine lpt_preprocess(this, time)
    class(lpt_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(vector_t), pointer :: x_old, y_old, z_old, u_old, v_old, w_old
    integer :: ind(6)

    associate(x => this%particles%x, y => this%particles%y, &
         z => this%particles%z, u => this%particles%u, &
         v => this%particles%v, w => this%particles%w, &
         acc_x => this%particles%acc_x, &
         acc_y => this%particles%acc_y, &
         acc_z => this%particles%acc_z, &
         u_lag => this%particles%u_lag, &
         v_lag => this%particles%v_lag, &
         w_lag => this%particles%w_lag, &
         u_laglag => this%particles%u_laglag, &
         v_laglag => this%particles%v_laglag, &
         w_laglag => this%particles%w_laglag, &
         acc_xlag => this%particles%acc_xlag, &
         acc_ylag => this%particles%acc_ylag, &
         acc_zlag => this%particles%acc_zlag, &
         acc_xlaglag => this%particles%acc_xlaglag, &
         acc_ylaglag => this%particles%acc_ylaglag, &
         acc_zlaglag => this%particles%acc_zlaglag, &
         n => this%particles%n)
      if (time%t .lt. this%start_time) return
      call this%sync_time_controller(time)
      if (abs(this%lpt_time%dt) .le. epsilon(1.0_rp)) return

      call profiler_start_region('LPT_time_integration')

      call neko_scratch_registry%request(x_old, ind(1), n, .false.)
      call neko_scratch_registry%request(y_old, ind(2), n, .false.)
      call neko_scratch_registry%request(z_old, ind(3), n, .false.)
      call neko_scratch_registry%request(u_old, ind(4), n, .false.)
      call neko_scratch_registry%request(v_old, ind(5), n, .false.)
      call neko_scratch_registry%request(w_old, ind(6), n, .false.)

      x_old = x
      y_old = y
      z_old = z
      u_old = u
      v_old = v
      w_old = w

      ! Advance the particle state from the previously stored RHS.
      if (this%inertia) then
         call this%ODE_integrate_ab_3c(u, v, w, acc_x, acc_y, acc_z, &
              acc_xlag, acc_ylag, acc_zlag, acc_xlaglag, acc_ylaglag, &
              acc_zlaglag, n)
      end if

      ! Advance the coordinates using the velocity history available at step
      ! entry, before the fluid solve refreshes the current RHS.
      call this%ODE_integrate_ab_3c(x, y, z, u_old, v_old, w_old, &
           u_lag, v_lag, w_lag, u_laglag, v_laglag, w_laglag, n)

      ! Handle the wall collisions with the pre-step RHS.
      if (this%inertia .and. this%elastic_wall_enabled) then
         call lpt_handle_elastic_wall_collisions(this, x_old, y_old, z_old, &
              u_old, v_old, w_old)
      end if

      ! Update lag histories for the next Adams-Bashforth step.
      if (this%lag_len .gt. 0) then
         call update_lags(u_lag, u_laglag, u_old)
         call update_lags(v_lag, v_laglag, v_old)
         call update_lags(w_lag, w_laglag, w_old)
         if (this%inertia) then
            call update_lags(acc_xlag, acc_xlaglag, acc_x)
            call update_lags(acc_ylag, acc_ylaglag, acc_y)
            call update_lags(acc_zlag, acc_zlaglag, acc_z)
         end if
         this%history_len = min(this%history_len + 1, this%lag_len)
      end if

      call neko_scratch_registry%relinquish(ind)

    end associate

    call profiler_end_region('LPT_time_integration')

  end subroutine lpt_preprocess

  !> Refresh particle/fluid coupling after the fluid step and emit output.
  !! @param time Current simulation time used for output scheduling.
  subroutine lpt_compute(this, time)
    class(lpt_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    if (time%t .lt. this%start_time) return

    call this%update_current_rhs()

    if (this%output_enabled) then
       if (this%output_controller%check(time)) then
          call this%write_output(time)
          call this%output_controller%register_execution(time)
       end if
    end if
  end subroutine lpt_compute

  !> Build an LPT-local time-step history from the times at which LPT runs.
  !! @param time Current simulation time from the main solver.
  subroutine sync_time_controller(this, time)
    class(lpt_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp) :: dt_local
    real(kind=rp) :: t_ref
    integer :: i

    if (.not. this%lpt_time_initialized) then
       this%lpt_time = time
       t_ref = time%t
       if (this%start_time .gt. time%t) t_ref = this%start_time
       this%lpt_time%t = t_ref
       this%lpt_time%tlag = t_ref
       this%lpt_time%dt = 0.0_rp
       this%lpt_time%dtlag = 0.0_rp
       this%lpt_time_initialized = .true.
       return
    end if

    dt_local = time%t - this%lpt_time%t
    if (abs(dt_local) .le. epsilon(1.0_rp)) then
       this%lpt_time%t = time%t
       this%lpt_time%tstep = time%tstep + this%tstep_offset
       this%lpt_time%dt = 0.0_rp
       return
    end if

    do i = size(this%lpt_time%dtlag), 2, -1
       this%lpt_time%dtlag(i) = this%lpt_time%dtlag(i - 1)
       this%lpt_time%tlag(i) = this%lpt_time%tlag(i - 1)
    end do
    this%lpt_time%dtlag(1) = this%lpt_time%dt
    this%lpt_time%tlag(1) = this%lpt_time%t
    this%lpt_time%dt = dt_local
    this%lpt_time%t = time%t
    this%lpt_time%tstep = time%tstep + this%tstep_offset
  end subroutine sync_time_controller

  !> Advance a three-component state with variable-step Adams-Bashforth.
  !! @param sol_x State component in x.
  !! @param sol_y State component in y.
  !! @param sol_z State component in z.
  !! @param n Number of local particles to advance.
  subroutine ODE_integrate_ab_3c(this, sol_x, sol_y, sol_z, &
       rhs_x, rhs_y, rhs_z, rhs_xlag, rhs_ylag, rhs_zlag, &
       rhs_xlaglag, rhs_ylaglag, rhs_zlaglag, n)
    class(lpt_t), intent(inout) :: this
    type(vector_t), intent(inout) :: sol_x, sol_y, sol_z
    type(vector_t), intent(in) :: rhs_x, rhs_y, rhs_z
    type(vector_t), intent(in) :: rhs_xlag, rhs_ylag, rhs_zlag
    type(vector_t), intent(in) :: rhs_xlaglag, rhs_ylaglag
    type(vector_t), intent(in) :: rhs_zlaglag
    integer, intent(in) :: n
    type(ab_time_scheme_t) :: ab_scheme
    real(kind=rp) :: ab_coeffs(4), dt_history(10)
    real(kind=rp) :: dtc
    integer :: i
    integer :: nadv

    if (n .eq. 0) return

    ! set up AB coefficients based on the history length available
    nadv = this%time_order
    nadv = min(nadv, this%history_len + 1)

    dt_history = 0.0_rp
    dt_history(1) = this%lpt_time%dt
    dt_history(2) = this%lpt_time%dtlag(1)
    dt_history(3) = this%lpt_time%dtlag(2)
    call ab_scheme%compute_coeffs(ab_coeffs, dt_history, nadv)

    ! contribution from the current velocity
    dtc = this%lpt_time%dt * ab_coeffs(1)

    call vector_add2s2(sol_x, rhs_x, dtc, n)
    call vector_add2s2(sol_y, rhs_y, dtc, n)
    call vector_add2s2(sol_z, rhs_z, dtc, n)

    if (nadv .ge. 2) then
       dtc = this%lpt_time%dt * ab_coeffs(2)
       call vector_add2s2(sol_x, rhs_xlag, dtc, n)
       call vector_add2s2(sol_y, rhs_ylag, dtc, n)
       call vector_add2s2(sol_z, rhs_zlag, dtc, n)
    end if

    if (nadv .ge. 3) then
       dtc = this%lpt_time%dt * ab_coeffs(3)
       call vector_add2s2(sol_x, rhs_xlaglag, dtc, n)
       call vector_add2s2(sol_y, rhs_ylaglag, dtc, n)
       call vector_add2s2(sol_z, rhs_zlaglag, dtc, n)
    end if

  end subroutine ODE_integrate_ab_3c

  !> Write one trajectory snapshot.
  !! @param time Current simulation time stored with particle data.
  subroutine write_output(this, time)
    class(lpt_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp), allocatable :: local_data(:,:)
    integer :: n_local
    integer :: i
    integer :: n_data

    n_local = this%particles%n

    call this%particles%device_sync(DEVICE_TO_HOST)

    if (this%inertia) then
       n_data = 11
    else
       n_data = 9
    end if
    allocate(local_data(n_data, n_local))
    do i = 1, n_local
       local_data(1,i) = real(time%tstep + this%tstep_offset, rp)
       local_data(2,i) = time%t
       local_data(3,i) = real(this%particles%ids(i), rp)
       local_data(4,i) = this%particles%x%x(i)
       local_data(5,i) = this%particles%y%x(i)
       local_data(6,i) = this%particles%z%x(i)
       local_data(7,i) = this%particles%u%x(i)
       local_data(8,i) = this%particles%v%x(i)
       local_data(9,i) = this%particles%w%x(i)
       if (this%inertia) then
          local_data(10,i) = this%particles%d%x(i)
          local_data(11,i) = this%particles%rho%x(i)
       end if
    end do

    call this%output%write(local_data, n_local)
    deallocate(local_data)
  end subroutine write_output

  !> Free all LPT-owned state and reset pointers/flags.
  subroutine lpt_free(this)
    class(lpt_t), intent(inout) :: this

    call this%particles%free()
    call this%global_interp%free()
    call this%periodic_bc%free()
    call this%migration%free()
    call this%output%free()
    call this%output_controller%free()

    this%u_field => null()
    this%v_field => null()
    this%w_field => null()
    this%msh => null()
    this%dm_Xh => null()
    this%coef => null()
    if (allocated(this%wall_zone_indices)) deallocate(this%wall_zone_indices)
    if (allocated(this%wall_facet_mask)) deallocate(this%wall_facet_mask)
    this%elastic_wall_enabled = .false.
    this%output_enabled = .false.
    this%log = .true.
    this%start_time = -huge(0.0_rp)
    this%history_len = 0
    if (allocated(this%name)) deallocate(this%name)
    call this%lpt_time%reset()
    this%lpt_time_initialized = .false.

    ! The payload and its descriptors belong to the case checkpoint, which may
    ! already be freed, so only drop the references here.
    nullify(this%chkp_payload)
    nullify(this%chkp_ids)
    nullify(this%chkp_hook%lpt)
    if (allocated(this%chkp_vectors)) deallocate(this%chkp_vectors)
    if (allocated(this%chkp_arrays)) deallocate(this%chkp_arrays)
    this%chkp_state = 0.0_dp
    this%tstep_offset = 0
    if (allocated(this%seed_params)) deallocate(this%seed_params)
  end subroutine lpt_free

  !> Register the particle state in the case checkpoint.
  !! The payload is optional, so restarting from a checkpoint without particle
  !! data keeps the seeded particles. Both the HDF5 and `.chkp` formats store
  !! the payload.
  !! @param case Case owning the checkpoint.
  subroutine register_checkpoint(this, case)
    class(lpt_t), intent(inout), target :: this
    class(case_t), intent(inout), target :: case
    integer :: n_vectors

    this%chkp_payload => case%chkp%add_payload("lpt/" // trim(this%name))
    if (associated(this%chkp_payload%hook)) then
       call neko_error("lpt name '" // trim(this%name) // "' is not unique")
    end if
    this%chkp_payload%is_optional = .true.
    this%chkp_hook%lpt => this
    this%chkp_payload%hook => this%chkp_hook

    if (this%inertia) then
       n_vectors = 20
    else
       n_vectors = 12
    end if
    allocate(this%chkp_vectors(n_vectors))
    allocate(this%chkp_arrays(n_vectors))
    n_vectors = 0

    ! Keep x, y, z, u, v, w first and d, rho next: restart relies on it.
    call add_vector("x", this%particles%x)
    call add_vector("y", this%particles%y)
    call add_vector("z", this%particles%z)
    call add_vector("u", this%particles%u)
    call add_vector("v", this%particles%v)
    call add_vector("w", this%particles%w)
    if (this%inertia) then
       call add_vector("d", this%particles%d)
       call add_vector("rho", this%particles%rho)
    end if
    call add_vector("u_lag", this%particles%u_lag)
    call add_vector("v_lag", this%particles%v_lag)
    call add_vector("w_lag", this%particles%w_lag)
    call add_vector("u_laglag", this%particles%u_laglag)
    call add_vector("v_laglag", this%particles%v_laglag)
    call add_vector("w_laglag", this%particles%w_laglag)
    if (this%inertia) then
       call add_vector("acc_xlag", this%particles%acc_xlag)
       call add_vector("acc_ylag", this%particles%acc_ylag)
       call add_vector("acc_zlag", this%particles%acc_zlag)
       call add_vector("acc_xlaglag", this%particles%acc_xlaglag)
       call add_vector("acc_ylaglag", this%particles%acc_ylaglag)
       call add_vector("acc_zlaglag", this%particles%acc_zlaglag)
    end if

    this%chkp_ids => this%chkp_payload%add_resizable_array("ids")
    call this%chkp_payload%add_array_dp("state", this%chkp_state, &
         replicated = .true.)

  contains

    subroutine add_vector(name, vec)
      character(len=*), intent(in) :: name
      type(vector_t), target, intent(inout) :: vec

      n_vectors = n_vectors + 1
      this%chkp_vectors(n_vectors)%ptr => vec
      this%chkp_arrays(n_vectors)%ptr => &
           this%chkp_payload%add_resizable_array(name)
    end subroutine add_vector

  end subroutine register_checkpoint

  !> Pack the local particles and the LPT time state into the checkpoint
  !! buffers. Called by the checkpoint before every write.
  subroutine lpt_pack_checkpoint(this)
    class(lpt_t), intent(inout) :: this
    integer :: i, k, n

    if (.not. associated(this%chkp_payload)) return

    n = this%particles%n
    call this%particles%device_sync(DEVICE_TO_HOST)

    do k = 1, size(this%chkp_vectors)
       if (allocated(this%chkp_arrays(k)%ptr%buf)) then
          deallocate(this%chkp_arrays(k)%ptr%buf)
       end if
       allocate(this%chkp_arrays(k)%ptr%buf(n))
       do i = 1, n
          this%chkp_arrays(k)%ptr%buf(i) = &
               real(this%chkp_vectors(k)%ptr%x(i), dp)
       end do
    end do

    if (allocated(this%chkp_ids%buf)) deallocate(this%chkp_ids%buf)
    allocate(this%chkp_ids%buf(n))
    do i = 1, n
       this%chkp_ids%buf(i) = real(this%particles%ids(i), dp)
    end do

    this%chkp_state(1) = merge(1.0_dp, 0.0_dp, this%lpt_time_initialized)
    this%chkp_state(2) = real(this%history_len, dp)
    this%chkp_state(3) = real(this%lpt_time%tstep, dp)
    this%chkp_state(4) = this%lpt_time%t
    this%chkp_state(5) = this%lpt_time%dt
    this%chkp_state(6:15) = this%lpt_time%tlag
    this%chkp_state(16:25) = this%lpt_time%dtlag

  end subroutine lpt_pack_checkpoint

  !> Restore the particles after the case checkpoint has been read.
  !! Particles are read evenly distributed over the ranks and then migrated
  !! to their owners. Without particle data in the checkpoint, the particles
  !! are seeded from the initial condition at the restart time and the LPT
  !! time history restarts.
  !! @param time Restored simulation time.
  subroutine lpt_restart(this, time)
    class(lpt_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    character(len=LOG_SIZE) :: log_buf
    real(kind=rp), allocatable :: d(:), rho(:)
    type(json_file) :: seed_json
    integer :: i, k, n, ierr

    if (.not. associated(this%chkp_payload)) return

    if (.not. this%chkp_payload%loaded) then
       if (allocated(this%seed_params)) then
          if (this%log) then
             call neko_log%message("lpt '" // trim(this%name) // "': no " // &
                  "particle data in checkpoint, seeding from the initial " // &
                  "condition")
          end if
          call seed_json%load_from_string(this%seed_params)
          call this%seed_particles(seed_json)
          call seed_json%destroy()
          deallocate(this%seed_params)
       else if (this%log) then
          call neko_log%message("lpt '" // trim(this%name) // "': no " // &
               "particle data in checkpoint, keeping the seeded particles")
       end if
       this%lpt_time_initialized = .false.
       this%history_len = 0
       call this%sync_time_controller(time)
       call this%update_current_rhs()
       if (this%output_enabled) then
          call this%output%open_new()
          call this%write_output(time)
       end if
       call this%output_controller%set_counter(time)
       return
    end if
    if (allocated(this%seed_params)) deallocate(this%seed_params)

    n = size(this%chkp_ids%buf)
    call this%particles%free()
    if (n .gt. 0) then
       associate(a => this%chkp_arrays)
         if (this%inertia) then
            d = real(a(7)%ptr%buf, rp)
            rho = real(a(8)%ptr%buf, rp)
            call this%particles%init(real(a(1)%ptr%buf, rp), &
                 real(a(2)%ptr%buf, rp), real(a(3)%ptr%buf, rp), &
                 this%time_order, real(a(4)%ptr%buf, rp), &
                 real(a(5)%ptr%buf, rp), real(a(6)%ptr%buf, rp), d, rho)
         else
            call this%particles%init(real(a(1)%ptr%buf, rp), &
                 real(a(2)%ptr%buf, rp), real(a(3)%ptr%buf, rp), &
                 this%time_order)
         end if
       end associate
       do k = 1, size(this%chkp_vectors)
          this%chkp_vectors(k)%ptr = real(this%chkp_arrays(k)%ptr%buf, rp)
       end do
       do i = 1, n
          this%particles%ids(i) = nint(this%chkp_ids%buf(i))
       end do
    end if
    call MPI_Allreduce(n, this%particles%n_global, 1, MPI_INTEGER, MPI_SUM, &
         NEKO_COMM, ierr)

    this%lpt_time_initialized = this%chkp_state(1) .gt. 0.5_dp
    this%history_len = min(nint(this%chkp_state(2)), this%lag_len)
    this%lpt_time%tstep = nint(this%chkp_state(3))
    this%tstep_offset = this%lpt_time%tstep - time%tstep
    this%lpt_time%t = this%chkp_state(4)
    this%lpt_time%dt = this%chkp_state(5)
    this%lpt_time%tlag = this%chkp_state(6:15)
    this%lpt_time%dtlag = this%chkp_state(16:25)

    ! Migrate the particles to their owners and refresh the current RHS.
    call this%update_current_rhs()
    call this%output_controller%set_counter(time)
    if (this%output_enabled) call this%output%resume(time%t)

    if (this%log) then
       write(log_buf, '(A,A,A,I0,A)') "lpt '", trim(this%name), &
            "': restored ", this%particles%n_global, " particles"
       call neko_log%message(log_buf)
    end if

  end subroutine lpt_restart

  !> Pack the particle state of the owning LPT instance.
  subroutine lpt_checkpoint_hook_prepare_write(this)
    class(lpt_checkpoint_hook_t), intent(inout) :: this

    if (associated(this%lpt)) call this%lpt%pack_checkpoint()
  end subroutine lpt_checkpoint_hook_prepare_write

  !> Emit a setup summary for the configured LPT instance.
  subroutine log_status(this)
    class(lpt_t), intent(in) :: this
    character(len=LOG_SIZE) :: log_buf

    if (.not. this%log) return

    call neko_log%section("Lagrangian particle tracking")
    write(log_buf, '(A,A)') "Name: ", trim(this%name)
    call neko_log%message(log_buf)
    if (allocated(this%seed_params)) then
       call neko_log%message("Particles: restored from the checkpoint")
    else
       write(log_buf, '(A,I0)') "Global seeded particles: ", &
            this%particles%n_global
       call neko_log%message(log_buf)
    end if
    if (this%periodic_bc%periodic_enabled) then
       write(log_buf, '(A,I0)') "Periodic wrap directions: ", &
            this%periodic_bc%n_periodic_dirs
       call neko_log%message(log_buf)
    end if
    if (this%periodic_bc%rotational_periodic_enabled) then
       write(log_buf, '(A,3(ES13.5,A),ES13.5)') &
            "Rotational periodic sector: theta_min = ", &
            this%periodic_bc%rotational_theta_min, ", theta_max = ", &
            this%periodic_bc%rotational_theta_max, ", theta_len = ", &
            this%periodic_bc%rotational_theta_len, ""
       call neko_log%message(log_buf)
    end if
    if (this%elastic_wall_enabled) then
       write(log_buf, '(A,I0)') "Elastic wall zones configured: ", &
            size(this%wall_zone_indices)
       call neko_log%message(log_buf)
    end if
    write(log_buf, '(A,I0)') "Local particles on rank 0 at init: ", &
         this%particles%n
    if (pe_rank .eq. 0 .and. .not. allocated(this%seed_params)) then
       call neko_log%message(log_buf)
    end if
    call neko_log%end_section()
  end subroutine log_status

end module lpt
