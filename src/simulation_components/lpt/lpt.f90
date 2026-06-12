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
  use vector, only : vector_t
  use point, only : point_t
  use point_interpolator, only : point_interpolator_t
  use math, only : add2s2, cfill, cmult2, col2, col3, invcol2, sqrt_inplace, &
                   power, sub3, vdot3, cmult, cadd2, invcol3
  use ab_time_scheme, only : ab_time_scheme_t
  use tensor, only : trsp
  use lpt_periodic_bc, only : lpt_periodic_bc_t
  use lpt_redistribute, only : lpt_redistribute_t
  use comm, only : pe_rank, pe_size, NEKO_COMM, MPI_REAL_PRECISION
  use mpi_f08, only : MPI_Gather, MPI_Gatherv, MPI_INTEGER
  use csv_file, only : csv_file_t
  implicit none
  private

  type, private :: particles_t
     ! general particle info
     real(kind=rp), allocatable :: xyz(:,:)
     integer, allocatable :: ids(:)
     real(kind=rp), allocatable :: vel_lag(:,:,:) ! velocity lags
     integer :: time_order, lag_len
     integer :: n = 0
     integer :: n_global = 0
     ! inertia-specific info
     logical :: inertia = .false.
     real(kind=rp), allocatable :: vel(:,:)
     real(kind=rp), allocatable :: acc(:,:)
     real(kind=rp), allocatable :: d(:) ! diameter of particles
     real(kind=rp), allocatable :: rho(:) ! density of particles
     real(kind=rp), allocatable :: acc_lag(:,:,:) ! accelaration lags
   contains
     procedure, pass(this) :: init => particles_init
     procedure, pass(this) :: free => particles_free
  end type particles_t

  !> A simulation component for passive Lagrangian particle tracking.
  !! Particles are redistributed to the rank that owns their current location.
  type, public, extends(simulation_component_t) :: lpt_t
     type(field_t), pointer :: u => null()
     type(field_t), pointer :: v => null()
     type(field_t), pointer :: w => null()
     type(field_t), pointer :: mu_fluid => null()
     type(field_t), pointer :: rho_fluid => null()
     type(mesh_t), pointer :: msh => null()
     type(dofmap_t), pointer :: dm_Xh => null()
     type(coef_t), pointer :: coef => null()
     integer :: time_order, lag_len
     integer :: history_len = 0
     logical :: inertia = .false.
     logical :: elastic_wall_enabled = .false.
     integer, allocatable :: wall_zone_indices(:)
     logical, allocatable :: wall_facet_mask(:, :)
     type(time_state_t) :: lpt_time
     logical :: lpt_time_initialized = .false.
     type(global_interpolation_t) :: global_interp
     type(lpt_periodic_bc_t) :: periodic_bc
     type(lpt_redistribute_t) :: redistributor
     type(particles_t) :: particles
     type(file_t) :: output_file
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
     procedure, private, pass(this) :: handle_elastic_wall_collisions
     procedure, private, pass(this) :: update_current_rhs
     procedure, private, pass(this) :: init_wall_facet_mask
     procedure, private, pass(this) :: identify_wall_facet
     procedure, private, pass(this) :: wall_facet_normal
     procedure, private, pass(this) :: reflect_position
     procedure, private, nopass :: reflect_vector
     procedure, private, pass(this) :: write_output
     procedure, private, pass(this) :: log_status
  end type lpt_t

contains

  subroutine particles_init(this, xyz, time_order, vel, diameter, density)
    class(particles_t), intent(inout) :: this
    real(kind=rp), intent(in) :: xyz(:,:)
    integer, intent(in) :: time_order
    real(kind=rp), intent(in), optional :: vel(:,:)
    real(kind=rp), intent(in), optional :: diameter(:)
    real(kind=rp), intent(in), optional :: density(:)

    integer :: i

    if (size(xyz, 1) .ne. 3) then
       call neko_error("particles coordinates must have size 3 in dim 1")
    end if

    call this%free()

    this%n = size(xyz, 2)
    this%n_global = this%n
    this%time_order = time_order
    this%lag_len = time_order - 1
    allocate(this%xyz(3, this%n))
    allocate(this%vel(3, this%n))
    allocate(this%acc(3, this%n))
    allocate(this%ids(this%n))
    allocate(this%vel_lag(3, this%lag_len, this%n))
    allocate(this%acc_lag(3, this%lag_len, this%n))
    allocate(this%d(this%n))
    allocate(this%rho(this%n))
    this%xyz = xyz
    if (present(vel)) then
       this%vel = vel
    else
       this%vel = 0.0_rp
    end if
    this%acc = 0.0_rp
    this%vel_lag = 0.0_rp
    this%acc_lag = 0.0_rp
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

    if (allocated(this%xyz)) deallocate(this%xyz)
    if (allocated(this%ids)) deallocate(this%ids)
    if (allocated(this%vel)) deallocate(this%vel)
    if (allocated(this%acc)) deallocate(this%acc)
    if (allocated(this%vel_lag)) deallocate(this%vel_lag)
    if (allocated(this%acc_lag)) deallocate(this%acc_lag)
    if (allocated(this%d)) deallocate(this%d)
    if (allocated(this%rho)) deallocate(this%rho)
    this%n = 0
    this%n_global = 0
  end subroutine particles_free

  !> Construct from JSON.
  !! Supported particle input is either a flat `coordinates` array
  !! `[x1,y1,z1,x2,y2,z2,...]` or a three-column CSV `points_file`.
  subroutine lpt_init_from_json(this, json, case)
    class(lpt_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    type(json_file) :: interp_subdict
    character(len=:), allocatable :: name
    character(len=:), allocatable :: output_filename

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
    call this%redistributor%init(this%lag_len)

    call json_get(json, "inertia", this%inertia)

    if (this%inertia) then
       this%mu_fluid => neko_registry%get_field_by_name( &
                        case%fluid%name // "_mu")
       this%rho_fluid => neko_registry%get_field_by_name( &
                        case%fluid%name // "_rho")
    end if
    this%u => neko_registry%get_field_by_name("u")
    this%v => neko_registry%get_field_by_name("v")
    this%w => neko_registry%get_field_by_name("w")

    if (json%valid_path("wall_zone_indices")) then
       call json_get(json, "wall_zone_indices", this%wall_zone_indices)
       if (.not. this%inertia .and. size(this%wall_zone_indices) .gt. 0) then
          call neko_error("lpt wall_zone_indices requires inertia = true")
       end if
       this%elastic_wall_enabled = size(this%wall_zone_indices) .gt. 0
       if (this%elastic_wall_enabled) call this%init_wall_facet_mask()
    end if

    call this%read_particles_json(json)

    call json_get_subdict_or_empty(json, "interpolation", interp_subdict)
    call this%global_interp%init(case%fluid%dm_Xh, &
         params_subdict = interp_subdict)
    call this%periodic_bc%init(case%fluid%msh, case%fluid%dm_Xh, &
         case%fluid%c_Xh)
    call this%redistributor%redistribute_particles(this%global_interp, &
         this%periodic_bc, this%inertia, this%particles%xyz, &
         this%particles%ids, this%particles%vel_lag, this%particles%vel, &
         this%particles%acc, this%particles%d, this%particles%rho, &
         this%particles%acc_lag, this%particles%n, this%particles%n_global)
    call this%sync_time_controller(case%time)
    call this%update_current_rhs()

    call json_get_or_default(json, "output_filename", output_filename, &
         trim(this%name) // ".csv")
    call this%output_file%init(case%output_directory // &
         trim(output_filename), &
         header = "tstep,time,particle_id,x,y,z,u,v,w", overwrite = .true.)

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
    real(kind=rp), allocatable :: empty_xyz(:,:)
    real(kind=rp), allocatable :: empty_vel(:,:)
    real(kind=rp), allocatable :: empty_diams(:)
    real(kind=rp), allocatable :: empty_densities(:)

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
          call this%particles%init(reshape(coords, [3, size(coords) / 3]), &
                                   this%time_order, &
                                   reshape(vels, [3, size(vels) / 3]), &
                                   diams, densities)
          deallocate(coords)
          deallocate(vels)
          deallocate(diams)
          deallocate(densities)
       else if (json%valid_path("points_file")) then
          call this%read_particles_csv(json)
       else
          call neko_error("lpt requires either coordinates or points_file")
       end if
    else
       allocate(empty_xyz(3, 0))
       allocate(empty_vel(3, 0))
       allocate(empty_diams(0))
       allocate(empty_densities(0))
       call this%particles%init(empty_xyz, this%time_order, empty_vel, empty_diams, &
                                empty_densities)
       deallocate(empty_xyz)
       deallocate(empty_vel)
       deallocate(empty_diams)
       deallocate(empty_densities)
    end if
  end subroutine read_particles_json

  !> Read particle coordinates from a three-column CSV file.
  subroutine read_particles_csv(this, json)
    class(lpt_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: points_file
    type(file_t) :: file_in
    type(matrix_t) :: mat_in
    real(kind=rp), allocatable :: coords(:,:)
    real(kind=rp), allocatable :: vels(:,:)
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
          coords = mat_in%x(:, 1:3)
          vels = mat_in%x(:, 4:6)
          diams = mat_in%x(:, 7)
          densities = mat_in%x(:, 8)
          call this%particles%init(transpose(coords), this%time_order, &
                                   transpose(vels), &
                                   diams, densities)
          deallocate(coords)
          deallocate(vels)
          deallocate(diams)
          deallocate(densities)
       else
          call mat_in%init(ft%count_lines(), 3)
          call ft%read(mat_in)
          coords = mat_in%x(:, 1:3)
          call this%particles%init(transpose(coords), this%time_order)
          deallocate(coords)
       end if
    class default
       call neko_error("lpt points_file must be a csv file")
    end select

    call mat_in%free()
    call file_in%free()
  end subroutine read_particles_csv

  !> Interpolate the carrier velocity at the local particles.
  subroutine evaluate_velocity(this, vel_fluid, n)
    class(lpt_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: vel_fluid(3,n)
    real(kind=rp) :: vel_x(n), vel_y(n), vel_z(n)
    logical :: do_interp_on_host

    if (this%particles%n .eq. 0) return

    do_interp_on_host = .false.
    call this%global_interp%evaluate(vel_x, this%u%x, &
         do_interp_on_host)
    call this%global_interp%evaluate(vel_y, this%v%x, &
         do_interp_on_host)
    call this%global_interp%evaluate(vel_z, this%w%x, &
         do_interp_on_host)

    vel_fluid(1,:) = vel_x
    vel_fluid(2,:) = vel_y
    vel_fluid(3,:) = vel_z

  end subroutine evaluate_velocity

  !> Estimate the local particile acceleration
  subroutine evaluate_acceleration(this, vel_fluid, acceleration, n)
    class(lpt_t), intent(inout) :: this
    integer, intent(in) :: n
    real(kind=rp), intent(in) :: vel_fluid(:,:)
    real(kind=rp), intent(out) :: acceleration(3,n)
    real(kind=rp) :: tau_p(n), Re_p(n), f(n), rho_fluid_local(n)
    real(kind=rp) :: mu_fluid_local(n),  nu_fluid_local(n)
    real(kind=rp):: vel_rel(3, n), vel_rel_mag(n)
    real(kind=rp):: wa(n)

    integer :: i
    logical :: do_interp_on_host

    if (this%particles%n .eq. 0) return
    
    do_interp_on_host = .false.
    call this%global_interp%evaluate(mu_fluid_local, this%mu_fluid%x, &
         do_interp_on_host)
    call this%global_interp%evaluate(rho_fluid_local, this%rho_fluid%x, &
         do_interp_on_host)
    call invcol3(nu_fluid_local, mu_fluid_local, rho_fluid_local, n)

    ! compute the time scale
    call cfill(tau_p, 1.0_rp/18.0_rp, n)
    call col2(tau_p, this%particles%rho, n)
    call invcol2(tau_p, mu_fluid_local, n)
    call col2(tau_p, this%particles%d, n)
    call col2(tau_p, this%particles%d, n)

    ! compute the particle Reynolds number
    call sub3(vel_rel, vel_fluid, this%particles%vel, 3*n)
    call vdot3(vel_rel_mag, vel_rel(1,:), vel_rel(2,:), vel_rel(3,:), &
                            vel_rel(1,:), vel_rel(2,:), vel_rel(3,:), n)
    call sqrt_inplace(vel_rel_mag, n)
    call col3(Re_p, vel_rel_mag, this%particles%d, n)
    call invcol2(Re_p, nu_fluid_local, n)
    
    ! compute f
    wa = power(Re_p, 0.687_rp, n)
    call cmult(wa, 0.15_rp, n)
    call cadd2(f, wa, 1.0_rp, n)

    ! assemble compute the acceleration
    do i = 1, 3
         call col3(acceleration(i,:), vel_rel(i,:), f, n)
         call invcol2(acceleration(i,:), tau_p, n)
    end do

  end subroutine evaluate_acceleration

  !> Refresh the particle RHS using the current fluid solution.
  subroutine update_current_rhs(this)
    class(lpt_t), intent(inout) :: this
    real(kind=rp), allocatable :: vel_fluid(:,:)

    call this%redistributor%redistribute_particles(this%global_interp, &
         this%periodic_bc, this%inertia, this%particles%xyz, &
         this%particles%ids, this%particles%vel_lag, this%particles%vel, &
         this%particles%acc, this%particles%d, this%particles%rho, &
         this%particles%acc_lag, this%particles%n, this%particles%n_global)

    allocate(vel_fluid(3, this%particles%n))
    call this%evaluate_velocity(vel_fluid, this%particles%n)

    if (this%inertia) then
       call this%evaluate_acceleration(vel_fluid, this%particles%acc, &
            this%particles%n)
    else
       this%particles%vel = vel_fluid
    end if

    if (allocated(vel_fluid)) deallocate(vel_fluid)
  end subroutine update_current_rhs

  !> Advance particles with local Adams-Bashforth coefficients only.
  subroutine lpt_preprocess(this, time)
    class(lpt_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp), allocatable :: vel_rhs(:,:)
    real(kind=rp), allocatable :: xyz_old(:,:)
    integer :: j

    if (time%t .lt. this%start_time) return
    call this%sync_time_controller(time)
    if (abs(this%lpt_time%dt) .le. epsilon(1.0_rp)) return

    allocate(xyz_old(3, this%particles%n))
    xyz_old = this%particles%xyz
    allocate(vel_rhs(3, this%particles%n))

    ! Advance the particle state from the previously stored RHS.
    if (this%inertia) then
       vel_rhs = this%particles%vel
       call this%ODE_integrate_ab_3c(time, this%particles%vel, &
            this%particles%acc, this%particles%acc_lag, this%particles%n)
    else
       vel_rhs = this%particles%vel
    end if

    ! Advance the coordinates using the velocity history available at step
    ! entry, before the fluid solve refreshes the current RHS.
    call this%ODE_integrate_ab_3c(time, this%particles%xyz, &
         vel_rhs, this%particles%vel_lag, this%particles%n)

    ! Handle the wall collisions with the pre-step RHS.
    call this%handle_elastic_wall_collisions(xyz_old, vel_rhs)

    ! Update lag histories for the next Adams-Bashforth step.
    if (this%lag_len .gt. 0) then
       do j = this%lag_len, 2, -1
          this%particles%vel_lag(:, j, :) = &
               this%particles%vel_lag(:, j - 1, :)
          if (this%inertia) then
             this%particles%acc_lag(:, j, :) = &
                  this%particles%acc_lag(:, j - 1, :)
          end if
       end do
       this%particles%vel_lag(:, 1, :) = vel_rhs
       if (this%inertia) then
          this%particles%acc_lag(:, 1, :) = this%particles%acc
       end if
       this%history_len = min(this%history_len + 1, this%lag_len)
    end if

    if (allocated(vel_rhs)) deallocate(vel_rhs)
    if (allocated(xyz_old)) deallocate(xyz_old)

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
  subroutine ODE_integrate_ab_3c(this, time, solution, rhs, rhslags, n)
    class(lpt_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    real(kind=rp), intent(inout) :: solution(:,:)
    real(kind=rp), intent(in) :: rhs(:,:)
    real(kind=rp), intent(in) :: rhslags(:,:,:)
    integer, intent(in) :: n
    type(ab_time_scheme_t) :: ab_scheme
    real(kind=rp) :: ab_coeffs(4), dt_history(10)
    real(kind=rp) :: dtc
    integer :: i
    integer :: j
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
    do i = 1, 3
       call add2s2(solution(i, :), rhs(i, :), dtc, n)
    end do

    ! contribution from the lagged velocities
    do j = 2, nadv
       dtc = this%lpt_time%dt * ab_coeffs(j)
       do i = 1, 3
          call add2s2(solution(i, :), rhslags(i, j - 1, :), dtc, n)
       end do
    end do

  end subroutine ODE_integrate_ab_3c

  subroutine init_wall_facet_mask(this)
    class(lpt_t), intent(inout) :: this
    integer :: i
    integer :: j
    integer :: facet
    integer :: el
    integer :: zone_id

    if (allocated(this%wall_facet_mask)) deallocate(this%wall_facet_mask)
    allocate(this%wall_facet_mask(2 * this%msh%gdim, this%msh%nelv))
    this%wall_facet_mask = .false.

    do i = 1, size(this%wall_zone_indices)
       zone_id = this%wall_zone_indices(i)
       if (zone_id .lt. 1 .or. zone_id .gt. size(this%msh%labeled_zones)) then
          call neko_error("lpt wall_zone_indices contains an invalid zone id")
       end if

       do j = 1, this%msh%labeled_zones(zone_id)%size
          facet = this%msh%labeled_zones(zone_id)%facet_el(j)%x(1)
          el = this%msh%labeled_zones(zone_id)%facet_el(j)%x(2)
          this%wall_facet_mask(facet, el) = .true.
       end do
    end do
  end subroutine init_wall_facet_mask

  !> Restore particles that crossed a configured wall zone and reflect both
  !! the remaining trajectory and the velocity history for a purely elastic
  !! collision.
  subroutine handle_elastic_wall_collisions(this, xyz_old, vel_rhs)
    class(lpt_t), intent(inout) :: this
    real(kind=rp), intent(inout) :: xyz_old(:, :)
    real(kind=rp), intent(inout) :: vel_rhs(:, :)
    type(matrix_t) :: rst_new
    type(vector_t) :: x_t
    type(vector_t) :: y_t
    type(vector_t) :: z_t
    type(vector_t) :: resx
    type(vector_t) :: resy
    type(vector_t) :: resz
    integer, allocatable :: el_list(:)
    integer :: i
    integer :: j
    integer :: facet
    integer :: el
    integer :: el_mesh
    real(kind=rp) :: normal(3)

    if (.not. this%inertia) return
    if (.not. this%elastic_wall_enabled) return
    if (this%particles%n .eq. 0) return
    allocate(el_list(this%particles%n))
    call rst_new%init(3, this%particles%n)
    call x_t%init(this%particles%n)
    call y_t%init(this%particles%n)
    call z_t%init(this%particles%n)
    call resx%init(this%particles%n)
    call resy%init(this%particles%n)
    call resz%init(this%particles%n)

    do i = 1, this%particles%n
       el_list(i) = this%global_interp%el_owner0_local(i)
       x_t%x(i) = this%particles%xyz(1, i)
       y_t%x(i) = this%particles%xyz(2, i)
       z_t%x(i) = this%particles%xyz(3, i)
    end do

    call this%global_interp%rst_finder%find(rst_new, x_t, y_t, z_t, el_list, &
         this%particles%n, resx, resy, resz)

    do i = 1, this%particles%n
       el = el_list(i)
       if (el .lt. 0) cycle
       el_mesh = el + 1
       if (el_mesh .gt. this%msh%nelv) cycle

       facet = this%identify_wall_facet(this%global_interp%rst_local(:, i), &
            rst_new%x(:, i), el_mesh)
       if (facet .eq. 0) cycle

       call this%wall_facet_normal(el_mesh, facet, normal)
       if (norm2(normal) .le. epsilon(1.0_rp)) cycle

       call this%reflect_position(this%particles%xyz(:, i), &
            this%global_interp%rst_local(:, i), rst_new%x(:, i), el_mesh, &
            facet)

       call reflect_vector(this%particles%vel(:, i), normal)
       call reflect_vector(vel_rhs(:, i), normal)
       do j = 1, this%lag_len
          call reflect_vector(this%particles%vel_lag(:, j, i), normal)
       end do
    end do

    if (allocated(el_list)) deallocate(el_list)
    call rst_new%free()
    call x_t%free()
    call y_t%free()
    call z_t%free()
    call resx%free()
    call resy%free()
    call resz%free()
  end subroutine handle_elastic_wall_collisions

  !> Reflect the post-collision position by mirroring the overshoot in rst and
  !! interpolating the reflected point back to xyz on the same element.
  subroutine reflect_position(this, xyz, rst_old, rst_new, el, facet)
    class(lpt_t), intent(in) :: this
    real(kind=rp), intent(inout) :: xyz(3)
    real(kind=rp), intent(in) :: rst_old(3)
    real(kind=rp), intent(in) :: rst_new(3)
    integer, intent(in) :: el
    integer, intent(in) :: facet
    type(point_interpolator_t) :: interp
    type(point_t) :: rst_point(1)
    type(point_t), allocatable :: xyz_ref(:)
    real(kind=rp) :: collision_rst(3)
    real(kind=rp) :: reflected_rst(3)
    real(kind=rp) :: drst(3)
    real(kind=rp) :: alpha
    real(kind=rp) :: boundary_value
    integer :: dim

    reflected_rst = rst_new
    call facet_rst_info(facet, dim, boundary_value)
    drst = rst_new - rst_old

    if (abs(drst(dim)) .le. epsilon(1.0_rp)) return

    alpha = (boundary_value - rst_old(dim)) / drst(dim)
    alpha = max(0.0_rp, min(1.0_rp, alpha))
    collision_rst = rst_old + alpha * drst
    reflected_rst = collision_rst + (rst_new - collision_rst)
    reflected_rst(dim) = 2.0_rp * boundary_value - rst_new(dim)

    call rst_point(1)%init(real(reflected_rst, kind=kind(rst_point(1)%x)))
    call interp%init(this%dm_Xh%Xh)
    xyz_ref = interp%interpolate(rst_point, this%dm_Xh%x(:, :, :, el), &
         this%dm_Xh%y(:, :, :, el), this%dm_Xh%z(:, :, :, el))
    xyz = real(xyz_ref(1)%x, kind=rp)
    call interp%free()
    if (allocated(xyz_ref)) deallocate(xyz_ref)
  end subroutine reflect_position

  !> Return the wall facet crossed by a particle, or 0 if none matched.
  integer function identify_wall_facet(this, rst_old, rst_new, el) result(facet)
    class(lpt_t), intent(in) :: this
    real(kind=rp), intent(in) :: rst_old(3)
    real(kind=rp), intent(in) :: rst_new(3)
    integer, intent(in) :: el
    real(kind=rp), parameter :: tol = 1.0e-8_rp
    real(kind=rp) :: severity
    real(kind=rp) :: best_severity

    facet = 0
    best_severity = 0.0_rp

    severity = max((-1.0_rp - rst_new(1)), 0.0_rp)
    if (this%wall_facet_mask(1, el) .and. severity .gt. tol .and. &
         rst_new(1) .le. rst_old(1) + tol .and. severity .gt. best_severity) then
       if (this%wall_facet_mask(1, el)) then
          facet = 1
          best_severity = severity
       end if
    end if

    severity = max((rst_new(1) - 1.0_rp), 0.0_rp)
    if (this%wall_facet_mask(2, el) .and. severity .gt. tol .and. &
         rst_new(1) .ge. rst_old(1) - tol .and. severity .gt. best_severity) then
      if (this%wall_facet_mask(2, el)) then
         facet = 2
         best_severity = severity
      end if
    end if

    severity = max((-1.0_rp - rst_new(2)), 0.0_rp)
    if (this%wall_facet_mask(3, el) .and. severity .gt. tol .and. &
         rst_new(2) .le. rst_old(2) + tol .and. severity .gt. best_severity) then
       if (this%wall_facet_mask(3, el)) then
          facet = 3
          best_severity = severity
       end if
    end if

    severity = max((rst_new(2) - 1.0_rp), 0.0_rp)
    if (this%wall_facet_mask(4, el) .and. severity .gt. tol .and. &
         rst_new(2) .ge. rst_old(2) - tol .and. severity .gt. best_severity) then
       if (this%wall_facet_mask(4, el)) then
          facet = 4
          best_severity = severity
       end if
    end if

    if (this%msh%gdim .eq. 3) then
       severity = max((-1.0_rp - rst_new(3)), 0.0_rp)
       if (this%wall_facet_mask(5, el) .and. severity .gt. tol .and. &
            rst_new(3) .le. rst_old(3) + tol .and. severity .gt. best_severity) then
          if (this%wall_facet_mask(5, el)) then
             facet = 5
             best_severity = severity
          end if
       end if

       severity = max((rst_new(3) - 1.0_rp), 0.0_rp)
       if (this%wall_facet_mask(6, el) .and. severity .gt. tol .and. &
            rst_new(3) .ge. rst_old(3) - tol .and. severity .gt. best_severity) then
          if (this%wall_facet_mask(6, el)) then
             facet = 6
          end if
       end if
    end if
  end function identify_wall_facet

  !> Use the face-center SEM normal as the reflection normal.
  subroutine wall_facet_normal(this, el, facet, normal)
    class(lpt_t), intent(in) :: this
    integer, intent(in) :: el
    integer, intent(in) :: facet
    real(kind=rp), intent(out) :: normal(3)
    integer :: ic
    integer :: jc
    integer :: kc

    ic = max(1, (this%coef%Xh%lx + 1) / 2)
    jc = max(1, (this%coef%Xh%ly + 1) / 2)
    kc = max(1, (this%coef%Xh%lz + 1) / 2)

    select case (facet)
    case (1, 2)
       normal = this%coef%get_normal(1, jc, kc, el, facet)
    case (3, 4)
       normal = this%coef%get_normal(ic, 1, kc, el, facet)
    case (5, 6)
       normal = this%coef%get_normal(ic, jc, 1, el, facet)
    case default
       normal = 0.0_rp
    end select
  end subroutine wall_facet_normal

  pure subroutine facet_rst_info(facet, dim, boundary_value)
    integer, intent(in) :: facet
    integer, intent(out) :: dim
    real(kind=rp), intent(out) :: boundary_value

    select case (facet)
    case (1)
       dim = 1
       boundary_value = -1.0_rp
    case (2)
       dim = 1
       boundary_value = 1.0_rp
    case (3)
       dim = 2
       boundary_value = -1.0_rp
    case (4)
       dim = 2
       boundary_value = 1.0_rp
    case (5)
       dim = 3
       boundary_value = -1.0_rp
    case (6)
       dim = 3
       boundary_value = 1.0_rp
    case default
       dim = 1
       boundary_value = 0.0_rp
    end select
  end subroutine facet_rst_info

  pure subroutine reflect_vector(vec, normal)
    real(kind=rp), intent(inout) :: vec(3)
    real(kind=rp), intent(in) :: normal(3)
    real(kind=rp) :: nhat(3)
    real(kind=rp) :: nmag
    real(kind=rp) :: vn

    nmag = norm2(normal)
    if (nmag .le. epsilon(1.0_rp)) return

    nhat = normal / nmag
    vn = dot_product(vec, nhat)
    vec = vec - 2.0_rp * vn * nhat
  end subroutine reflect_vector

  !> Write one trajectory snapshot to CSV by gathering local particle data to
  !! rank 0.
  subroutine write_output(this, time)
    class(lpt_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(matrix_t) :: block
    real(kind=rp), allocatable :: local_data(:,:)
    real(kind=rp), allocatable :: global_data(:,:)
    integer, allocatable :: n_local_particles_per_rank(:)
    integer, allocatable :: recvcounts(:)
    integer, allocatable :: displs(:)
    integer :: n_local
    integer :: total_particles
    integer :: i
    integer :: ierr

    n_local = this%particles%n

    allocate(local_data(9, n_local))
    do i = 1, n_local
       local_data(1,i) = real(time%tstep, rp)
       local_data(2,i) = time%t
       local_data(3,i) = real(this%particles%ids(i), rp)
       local_data(4,i) = this%particles%xyz(1,i)
       local_data(5,i) = this%particles%xyz(2,i)
       local_data(6,i) = this%particles%xyz(3,i)
       local_data(7,i) = this%particles%vel(1,i)
       local_data(8,i) = this%particles%vel(2,i)
       local_data(9,i) = this%particles%vel(3,i)
    end do

    if (pe_rank .eq. 0) then
       allocate(n_local_particles_per_rank(pe_size))
    else
       allocate(n_local_particles_per_rank(0))
    end if
    call MPI_Gather(n_local, 1, MPI_INTEGER, n_local_particles_per_rank, 1, &
         MPI_INTEGER, 0, NEKO_COMM, ierr)

    if (pe_rank .eq. 0) then
       allocate(recvcounts(pe_size))
       allocate(displs(pe_size))
       recvcounts = 0
       displs = 0

       total_particles = 0
       do i = 1, pe_size
          total_particles = total_particles + n_local_particles_per_rank(i)
          recvcounts(i) = 9 * n_local_particles_per_rank(i)
          displs(i) = 9 * (total_particles - n_local_particles_per_rank(i))
       end do

       allocate(global_data(9, total_particles))
    else
       allocate(recvcounts(0))
       allocate(displs(0))
       allocate(global_data(9, 0))
    end if

    call MPI_Gatherv(local_data, 9 * n_local, MPI_REAL_PRECISION, global_data, &
         recvcounts, displs, MPI_REAL_PRECISION, 0, NEKO_COMM, ierr)

    if (pe_rank .eq. 0) then
       call block%init(total_particles, 9)
       call trsp(block%x, total_particles, global_data, 9)
       call this%output_file%write(block)
       call block%free()
    end if

    deallocate(global_data)
    deallocate(n_local_particles_per_rank)
    deallocate(recvcounts)
    deallocate(displs)
    deallocate(local_data)
  end subroutine write_output

  !> Free the component.
  subroutine lpt_free(this)
    class(lpt_t), intent(inout) :: this

    call this%particles%free()
    call this%global_interp%free()
    call this%periodic_bc%free()
    call this%redistributor%free()
    call this%output_file%free()

    this%u => null()
    this%v => null()
    this%w => null()
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
