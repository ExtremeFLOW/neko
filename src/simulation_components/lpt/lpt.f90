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
module lagrangian_particle_tracking
  use num_types, only : rp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use registry, only : neko_registry
  use field, only : field_t
  use case, only : case_t
  use mesh, only : mesh_t
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use json_utils, only : json_get, json_get_or_default, &
       json_get_subdict_or_empty
  use time_state, only : time_state_t
  use global_interpolation, only : global_interpolation_t
  use logger, only : neko_log, LOG_SIZE
  use utils, only : neko_error
  use file, only : file_t
  use matrix, only : matrix_t
  use vector_math, only : vector_add2s2, vector_cfill, vector_cmult2, &
      vector_col2, vector_col3, vector_invcol2, vector_sqrt_inplace, &
      vector_power, vector_sub3, vector_vdot3, vector_cmult, &
      vector_cadd2, vector_invcol3
  use ab_time_scheme, only : ab_time_scheme_t
  use lpt_periodic_bc, only : lpt_periodic_bc_t
  use lpt_migrate, only : lpt_migrate_t, LPT_MIGRATE_TO_OWNER, &
       LPT_MIGRATE_NONE
  use lpt_wall_collision, only : lpt_handle_elastic_wall_collisions, &
       lpt_init_wall_facet_mask
  use lpt_output, only : lpt_output_t
  use comm, only : pe_rank
  use csv_file, only : csv_file_t
  use vector, only : vector_t
  use device, only : DEVICE_TO_HOST
  implicit none
  private

  type, private :: particles_t
     ! general particle info
     type(vector_t) :: x, y, z
     type(vector_t) :: u, v, w
     integer, allocatable :: ids(:)
     type(vector_t) :: u_lag, v_lag, w_lag
     type(vector_t) :: u_laglag, v_laglag, w_laglag
     integer :: time_order, lag_len
     integer :: n = 0
     integer :: n_global = 0
     ! inertia-specific info
     logical :: inertia = .false.
     type(vector_t) :: acc_x, acc_y, acc_z
     type(vector_t) :: d ! diameter of particles
     type(vector_t) :: rho ! density of particles
     type(vector_t) :: acc_xlag, acc_ylag, acc_zlag
     type(vector_t) :: acc_xlaglag, acc_ylaglag, acc_zlaglag
   contains
     procedure, pass(this) :: init => particles_init
     procedure, pass(this) :: free => particles_free
     procedure, pass(this) :: device_sync => particles_device_sync
  end type particles_t

  !> A simulation component for passive Lagrangian particle tracking.
  type, public, extends(simulation_component_t) :: lpt_t
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
     logical :: output_enabled = .false.
     logical :: log = .true.
     real(kind=rp) :: start_time = -huge(0.0_rp)
   contains
     procedure, pass(this) :: init => lpt_init_from_json
     procedure, pass(this) :: free => lpt_free
     procedure, pass(this) :: preprocess_ => lpt_preprocess
     procedure, pass(this) :: compute_ => lpt_compute
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

contains

  subroutine particles_init(this, x, y, z, time_order, u, v, w, diameter, &
       density)
    class(particles_t), intent(inout) :: this
    real(kind=rp), intent(in) :: x(:), y(:), z(:)
    integer, intent(in) :: time_order
    real(kind=rp), intent(in), optional :: u(:), v(:), w(:)
    real(kind=rp), intent(in), optional :: diameter(:)
    real(kind=rp), intent(in), optional :: density(:)

    integer :: i

    if (size(y) .ne. size(x) .or. size(z) .ne. size(x)) then
       call neko_error("particle coordinate arrays must have the same size")
    end if

    call this%free()

    this%n = size(x)
    this%n_global = this%n
    this%time_order = time_order
    this%lag_len = time_order - 1
    call this%x%init(this%n)
    call this%y%init(this%n)
    call this%z%init(this%n)
    call this%u%init(this%n)
    call this%v%init(this%n)
    call this%w%init(this%n)
    call this%acc_x%init(this%n)
    call this%acc_y%init(this%n)
    call this%acc_z%init(this%n)
    allocate(this%ids(this%n))
    call this%u_lag%init(this%n)
    call this%v_lag%init(this%n)
    call this%w_lag%init(this%n)
    call this%u_laglag%init(this%n)
    call this%v_laglag%init(this%n)
    call this%w_laglag%init(this%n)
    call this%acc_xlag%init(this%n)
    call this%acc_ylag%init(this%n)
    call this%acc_zlag%init(this%n)
    call this%acc_xlaglag%init(this%n)
    call this%acc_ylaglag%init(this%n)
    call this%acc_zlaglag%init(this%n)
    call this%d%init(this%n)
    call this%rho%init(this%n)
    this%x = x
    this%y = y
    this%z = z
    if (present(u) .and. present(v) .and. present(w)) then
       this%u = u
       this%v = v
       this%w = w
    else
       this%u = 0.0_rp
       this%v = 0.0_rp
       this%w = 0.0_rp
    end if
    if (present(diameter)) then
       this%d = diameter
    else
       this%d = 0.0_rp
    end if
    if (present(density)) then
       this%rho = density
    else
       this%rho = 0.0_rp
    end if
    do i = 1, this%n
       this%ids(i) = i
    end do
  end subroutine particles_init

  subroutine particles_free(this)
    class(particles_t), intent(inout) :: this

    call this%x%free()
    call this%y%free()
    call this%z%free()
    if (allocated(this%ids)) deallocate(this%ids)
    call this%u%free()
    call this%v%free()
    call this%w%free()
    call this%acc_x%free()
    call this%acc_y%free()
    call this%acc_z%free()
    call this%u_lag%free()
    call this%v_lag%free()
    call this%w_lag%free()
    call this%u_laglag%free()
    call this%v_laglag%free()
    call this%w_laglag%free()
    call this%acc_xlag%free()
    call this%acc_ylag%free()
    call this%acc_zlag%free()
    call this%acc_xlaglag%free()
    call this%acc_ylaglag%free()
    call this%acc_zlaglag%free()
    call this%d%free()
    call this%rho%free()
    this%n = 0
    this%n_global = 0
  end subroutine particles_free

  !> Synchronise the host and device data for particles
  subroutine particles_device_sync(this, memdir)
    class(particles_t), intent(inout) :: this
    integer, intent(in) :: memdir

    call this%x%copy_from(memdir, .false.)
    call this%y%copy_from(memdir, .false.)
    call this%z%copy_from(memdir, .false.)
    call this%u%copy_from(memdir, .false.)
    call this%v%copy_from(memdir, .false.)
    call this%w%copy_from(memdir, .true.)
    if (this%inertia) then
       call this%d%copy_from(memdir, .false.)
       call this%rho%copy_from(memdir, .true.)
    end if

  end subroutine particles_device_sync

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
    character(len=:), allocatable :: output_path
    integer :: migration_strategy_id

    call this%free()

    call json_get_or_default(json, "name", name, "lpt")
    call json_get_or_default(json, "log", this%log, .true.)
    call json_get_or_default(json, "start_time", this%start_time, &
         -huge(0.0_rp))
    call this%init_base(json, case)
    this%preprocess_controller = this%compute_controller

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

    if (json%valid_path("wall_zone_indices")) then
       call json_get(json, "wall_zone_indices", this%wall_zone_indices)
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

    call this%read_particles_json(json)
    call this%migration%initialize_particle_distribution(this%inertia, &
         this%particles%x, this%particles%y, this%particles%z, &
         this%particles%ids, &
         this%particles%u_lag, this%particles%v_lag, this%particles%w_lag, &
         this%particles%u_laglag, this%particles%v_laglag, &
         this%particles%w_laglag, &
         this%particles%u, this%particles%v, this%particles%w, &
         this%particles%acc_x, this%particles%acc_y, this%particles%acc_z, &
         this%particles%d, this%particles%rho, &
         this%particles%acc_xlag, this%particles%acc_ylag, &
         this%particles%acc_zlag, this%particles%acc_xlaglag, &
         this%particles%acc_ylaglag, this%particles%acc_zlaglag, &
         this%particles%n)

    call json_get_subdict_or_empty(json, "interpolation", interp_subdict)
    call this%global_interp%init(case%fluid%dm_Xh, &
         params_subdict = interp_subdict)
    call this%periodic_bc%init(case%fluid%msh, case%fluid%dm_Xh, &
         case%fluid%c_Xh)
    call this%migration%migrate_particles(this%global_interp, &
         this%periodic_bc, this%inertia, &
         this%particles%x, this%particles%y, this%particles%z, &
         this%particles%ids, &
         this%particles%u_lag, this%particles%v_lag, this%particles%w_lag, &
         this%particles%u_laglag, this%particles%v_laglag, &
         this%particles%w_laglag, &
         this%particles%u, this%particles%v, this%particles%w, &
         this%particles%acc_x, this%particles%acc_y, this%particles%acc_z, &
         this%particles%d, this%particles%rho, &
         this%particles%acc_xlag, this%particles%acc_ylag, &
         this%particles%acc_zlag, this%particles%acc_xlaglag, &
         this%particles%acc_ylaglag, this%particles%acc_zlaglag, &
         this%particles%n, this%particles%n_global)
    call this%sync_time_controller(case%time)
    call this%update_current_rhs()

    call json_get_or_default(json, "output_filename", output_filename, &
         trim(this%name) // ".csv")
    output_path = case%output_directory // trim(output_filename)
    call this%output%init(output_path, this%inertia)

    ! output at the initialisation
    this%output_enabled = .true.
    call this%write_output(case%time)

    call this%log_status()
  end subroutine lpt_init_from_json

  !> Read particle coordinates from JSON.
  subroutine read_particles_json(this, json)
    class(lpt_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    real(kind=rp), allocatable :: coords(:)
    real(kind=rp), allocatable :: vels(:)
    real(kind=rp), allocatable :: diams(:)
    real(kind=rp), allocatable :: densities(:)
    real(kind=rp), allocatable :: x(:), y(:), z(:)
    real(kind=rp), allocatable :: u(:), v(:), w(:)

    if (pe_rank .eq. 0) then
       if (json%valid_path("coordinates")) then
          call json_get(json, "coordinates", coords)
          if (mod(size(coords), 3) .ne. 0) then
             call neko_error("lpt coordinates must contain 3 values per " // &
                  "particle")
          end if
          if (this%inertia) then
             call json_get(json, "velocities", vels)
             if (mod(size(vels), 3) .ne. 0) then
                call neko_error("lpt velocities must contain 3 values per " // &
                     "particle")
             end if
             call json_get(json, "diameters", diams)
             call json_get(json, "densities", densities)
          else
             allocate(vels(size(coords)))
             allocate(diams(size(coords)/3))
             allocate(densities(size(coords)/3))
             vels = 0.0_rp
             diams = 0.0_rp
             densities = 0.0_rp
          end if
          allocate(x(size(coords) / 3))
          allocate(y(size(coords) / 3))
          allocate(z(size(coords) / 3))
          allocate(u(size(vels) / 3))
          allocate(v(size(vels) / 3))
          allocate(w(size(vels) / 3))
          x = coords(1::3)
          y = coords(2::3)
          z = coords(3::3)
          u = vels(1::3)
          v = vels(2::3)
          w = vels(3::3)
          call this%particles%init(x, y, z, this%time_order, u, v, w, &
               diams, densities)
          deallocate(coords)
          deallocate(vels)
          deallocate(diams)
          deallocate(densities)
          deallocate(x)
          deallocate(y)
          deallocate(z)
          deallocate(u)
          deallocate(v)
          deallocate(w)
       else if (json%valid_path("points_file")) then
          call this%read_particles_csv(json)
       else
          call neko_error("lpt requires either coordinates or points_file")
       end if
    else
       return
    end if
  end subroutine read_particles_json

  !> Read particle coordinates from a three-column CSV file.
  subroutine read_particles_csv(this, json)
    class(lpt_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: points_file
    type(file_t) :: file_in
    type(matrix_t) :: mat_in
    real(kind=rp), allocatable :: x(:), y(:), z(:)
    real(kind=rp), allocatable :: u(:), v(:), w(:)
    real(kind=rp), allocatable :: diams(:)
    real(kind=rp), allocatable :: densities(:)

    if (pe_rank .ne. 0) return

    call json_get(json, "points_file", points_file)
    call file_in%init(trim(points_file))

    select type (ft => file_in%file_type)
    type is (csv_file_t)
       if (this%inertia) then
          call mat_in%init(ft%count_lines(), 8)
          call ft%read(mat_in)
          x = mat_in%x(:, 1)
          y = mat_in%x(:, 2)
          z = mat_in%x(:, 3)
          u = mat_in%x(:, 4)
          v = mat_in%x(:, 5)
          w = mat_in%x(:, 6)
          diams = mat_in%x(:, 7)
          densities = mat_in%x(:, 8)
          call this%particles%init(x, y, z, this%time_order, u, v, w, &
               diams, densities)
          deallocate(x)
          deallocate(y)
          deallocate(z)
          deallocate(u)
          deallocate(v)
          deallocate(w)
          deallocate(diams)
          deallocate(densities)
       else
          call mat_in%init(ft%count_lines(), 3)
          call ft%read(mat_in)
          x = mat_in%x(:, 1)
          y = mat_in%x(:, 2)
          z = mat_in%x(:, 3)
          call this%particles%init(x, y, z, this%time_order)
          deallocate(x)
          deallocate(y)
          deallocate(z)
       end if
    class default
       call neko_error("lpt points_file must be a csv file")
    end select

    call mat_in%free()
    call file_in%free()
  end subroutine read_particles_csv

  !> Interpolate the carrier velocity at the local particles.
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

  !> Estimate the local particile acceleration
  subroutine evaluate_acceleration(this, acc_x, acc_y, acc_z, &
             u_fluid, v_fluid, w_fluid)
    class(lpt_t), intent(inout) :: this
    type(vector_t), intent(in) :: u_fluid, v_fluid, w_fluid
    type(vector_t), intent(inout) :: acc_x, acc_y, acc_z
    type(vector_t) :: tau_p, Re_p, f, rho_fluid_local
    type(vector_t) :: mu_fluid_local,  nu_fluid_local
    type(vector_t) :: u_rel, v_rel, w_rel, vel_rel_mag
    type(vector_t) :: wa

    integer :: i
    logical :: do_interp_on_host

    if (this%particles%n .eq. 0) return

    call tau_p%init(this%particles%n)
    call Re_p%init(this%particles%n)
    call f%init(this%particles%n)
    call rho_fluid_local%init(this%particles%n)
    call mu_fluid_local%init(this%particles%n)
    call nu_fluid_local%init(this%particles%n)
    call u_rel%init(this%particles%n)
    call v_rel%init(this%particles%n)
    call w_rel%init(this%particles%n)
    call vel_rel_mag%init(this%particles%n)
    call wa%init(this%particles%n)

    do_interp_on_host = .false.
    call this%global_interp%evaluate(mu_fluid_local%x, this%mu_fluid%x, &
         do_interp_on_host)
    call this%global_interp%evaluate(rho_fluid_local%x, this%rho_fluid%x, &
         do_interp_on_host)
    call vector_invcol3(nu_fluid_local, mu_fluid_local, rho_fluid_local)

    ! compute the time scale
    call vector_cfill(tau_p, 1.0_rp/18.0_rp)
    call vector_col2(tau_p, this%particles%rho)
    call vector_invcol2(tau_p, mu_fluid_local)
    call vector_col2(tau_p, this%particles%d)
    call vector_col2(tau_p, this%particles%d)

    ! compute the particle Reynolds number
    call vector_sub3(u_rel, u_fluid, this%particles%u)
    call vector_sub3(v_rel, v_fluid, this%particles%v)
    call vector_sub3(w_rel, w_fluid, this%particles%w)
    call vector_vdot3(vel_rel_mag, u_rel, v_rel, w_rel, u_rel, v_rel, w_rel)
    call vector_sqrt_inplace(vel_rel_mag)
    call vector_col3(Re_p, vel_rel_mag, this%particles%d)
    call vector_invcol2(Re_p, nu_fluid_local)
    
    ! compute f
    call vector_power(wa, Re_p, this%nonlinear_exponent)
    call vector_cmult(wa, this%nonlinear_coefficient)
    call vector_cadd2(f, wa, 1.0_rp)

    ! assemble compute the acceleration
    call vector_col3(acc_x, u_rel, f)
    call vector_col3(acc_y, v_rel, f)
    call vector_col3(acc_z, w_rel, f)
    call vector_invcol2(acc_x, tau_p)
    call vector_invcol2(acc_y, tau_p)
    call vector_invcol2(acc_z, tau_p)

    call tau_p%free()
    call Re_p%free()
    call f%free()
    call rho_fluid_local%free()
    call mu_fluid_local%free()
    call nu_fluid_local%free()
    call u_rel%free()
    call v_rel%free()
    call w_rel%free()
    call vel_rel_mag%free()
    call wa%free()

  end subroutine evaluate_acceleration

  !> Refresh the particle RHS using the current fluid solution.
  subroutine update_current_rhs(this)
    class(lpt_t), intent(inout) :: this
    type(vector_t) :: u_fluid, v_fluid, w_fluid

    call this%migration%migrate_particles(this%global_interp, &
         this%periodic_bc, this%inertia, &
         this%particles%x, this%particles%y, this%particles%z, &
         this%particles%ids, &
         this%particles%u_lag, this%particles%v_lag, this%particles%w_lag, &
         this%particles%u_laglag, this%particles%v_laglag, &
         this%particles%w_laglag, &
         this%particles%u, this%particles%v, this%particles%w, &
         this%particles%acc_x, this%particles%acc_y, this%particles%acc_z, &
         this%particles%d, this%particles%rho, &
         this%particles%acc_xlag, this%particles%acc_ylag, &
         this%particles%acc_zlag, this%particles%acc_xlaglag, &
         this%particles%acc_ylaglag, this%particles%acc_zlaglag, &
         this%particles%n, this%particles%n_global)

    call u_fluid%init(this%particles%n)
    call v_fluid%init(this%particles%n)
    call w_fluid%init(this%particles%n)

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

    call u_fluid%free()
    call v_fluid%free()
    call w_fluid%free()
  end subroutine update_current_rhs

  !> Refresh the particle lags.
  subroutine update_lags(lag, laglag, new_values)
    type(vector_t), intent(inout) :: lag, laglag
    type(vector_t), intent(in) :: new_values

    laglag = lag
    lag = new_values

  end subroutine update_lags

  !> Advance particles with local Adams-Bashforth coefficients only.
  subroutine lpt_preprocess(this, time)
    class(lpt_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(vector_t) :: x_old, y_old, z_old, u_old, v_old, w_old
   !  ! pieces to test the accuracy: force exact solution at the first 2 steps
   !  real(kind=rp) :: v0
   !  real(kind=rp) :: x0
   !  real(kind=rp) :: tau_p

   !  v0 = 1.0_rp
   !  x0 = 0.0_rp
! write(*,*) "rank: ", pe_rank, ", particle ids: ", this%particles%ids

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
      !  ! pieces to test the accuracy: force exact solution at the first 2 steps
      !  if (time%tstep .le. 2) then
      !     tau_p = 1.0_rp/18.0_rp * this%particles%rho%x(1) * this%particles%d%x(1)**2
      !     u%x(1) = v0 * exp(-time%t/tau_p)
      !  end if
    end if

    ! Advance the coordinates using the velocity history available at step
    ! entry, before the fluid solve refreshes the current RHS.
    call this%ODE_integrate_ab_3c(x, y, z, u_old, v_old, w_old, &
         u_lag, v_lag, w_lag, u_laglag, v_laglag, w_laglag, n)
   !  ! pieces to test the accuracy: force exact solution at the first 2 steps
   !  if (time%tstep .le. 2) then
   !     x%x(1) = x0 + tau_p * v0 * (1.0_rp - exp(-time%t/tau_p))
   !  end if

    ! Handle the wall collisions with the pre-step RHS.
    if (this%inertia .and. this%elastic_wall_enabled) then
       call lpt_handle_elastic_wall_collisions(this%global_interp, this%msh, &
            this%dm_Xh, this%coef, this%wall_facet_mask, x_old%x, y_old%x, z_old%x, &
            x%x, y%x, z%x, &
            this%particles%d%x, &
            u%x, v%x, w%x, &
            u_lag%x, v_lag%x, w_lag%x, u_laglag%x, v_laglag%x, w_laglag%x, &
            acc_xlag%x, acc_ylag%x, acc_zlag%x, acc_xlaglag%x, acc_ylaglag%x, &
            acc_zlaglag%x, &
            u_old%x, v_old%x, w_old%x, acc_x%x, acc_y%x, acc_z%x, this%lag_len)
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

    call x_old%free()
    call y_old%free()
    call z_old%free()
    call u_old%free()
    call v_old%free()
    call w_old%free()
    end associate

  end subroutine lpt_preprocess

  !> Refresh particle/fluid coupling after the fluid step and emit output.
  subroutine lpt_compute(this, time)
    class(lpt_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    if (time%t .lt. this%start_time) return

    call this%update_current_rhs()

    if (this%output_enabled) then
       if (this%output_controller%check(time)) then
          call this%write_output(time)
          call this%output_controller%register_execution()
       end if
    end if
  end subroutine lpt_compute

  !> Build an LPT-local time-step history from the times at which LPT runs.
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
       this%lpt_time%tstep = time%tstep
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
    this%lpt_time%tstep = time%tstep
  end subroutine sync_time_controller

  !> Performing ODE integration by Adam-Bashforth scheme
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
       local_data(1,i) = real(time%tstep, rp)
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

  !> Free the component.
  subroutine lpt_free(this)
    class(lpt_t), intent(inout) :: this

    call this%particles%free()
    call this%global_interp%free()
    call this%periodic_bc%free()
    call this%migration%free()
    call this%output%free()

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
    call this%lpt_time%reset()
    this%lpt_time_initialized = .false.
    call this%free_base()
  end subroutine lpt_free

  !> Emit a setup summary.
  subroutine log_status(this)
    class(lpt_t), intent(in) :: this
    character(len=LOG_SIZE) :: log_buf

    if (.not. this%log) return

    call neko_log%section("Lagrangian particle tracking")
    write(log_buf, '(A,A)') "Name: ", trim(this%name)
    call neko_log%message(log_buf)
    write(log_buf, '(A,I0)') "Global seeded particles: ", &
         this%particles%n_global
    call neko_log%message(log_buf)
    if (this%periodic_bc%periodic_enabled) then
       write(log_buf, '(A,I0)') "Periodic wrap directions: ", &
            this%periodic_bc%n_periodic_dirs
       call neko_log%message(log_buf)
    end if
    if (this%periodic_bc%rotational_periodic_enabled) then
       write(log_buf, '(A,3(ES13.5,A),ES13.5)') &
            "Rotational periodic sector: theta_min=", &
            this%periodic_bc%rotational_theta_min, ", theta_max=", &
            this%periodic_bc%rotational_theta_max, ", theta_len=", &
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
    if (pe_rank .eq. 0) call neko_log%message(log_buf)
    call neko_log%end_section()
  end subroutine log_status

end module lagrangian_particle_tracking
