! Copyright (c) 2025-2026, The Neko Authors
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
module fluid_scheme_compressible_ns
  use comm, only : NEKO_COMM
  use advection, only : advection_t
  use device, only : device_memcpy, HOST_TO_DEVICE
  use operators, only : div, rotate_cyc
  use field_series, only : field_series_t
  use bdf_time_scheme, only : bdf_time_scheme_t
  use time_scheme_controller, only : time_scheme_controller_t
  use field, only : field_t
  use fluid_scheme_compressible, only : fluid_scheme_compressible_t, &
       fluid_scheme_compressible_validate, &
       fluid_scheme_compressible_compute_cfl
  use euler_idp, only : euler_idp_t, EULER_IDP_DIAGNOSTICS_FULL
  use scratch_registry, only : neko_scratch_registry
  use gs_ops, only : GS_OP_MIN, GS_OP_MAX
  use gather_scatter, only : gs_t
  use num_types, only : rp, dp
  use mesh, only : mesh_t
  use checkpoint, only : chkp_t
  use json_module, only : json_file, json_core, json_value
  use json_utils, only : json_get, json_get_or_default, json_extract_item
  use profiler, only : profiler_start_region, profiler_end_region
  use user_intf, only : user_t
  use time_step_controller, only : time_step_controller_t
  use ax_product, only : ax_t, ax_helm_allocator
  use coefs, only : coef_t
  use compressible_residual, only : compressible_rhs_t, compressible_rhs_factory
  use neko_config, only : NEKO_BCKND_DEVICE
  use runge_kutta_time_scheme, only : runge_kutta_time_scheme_t
  use bc_list, only : bc_list_t
  use bc, only : bc_t
  use utils, only : neko_error, neko_type_error
  use logger, only : LOG_SIZE, neko_log, NEKO_LOG_VERBOSE
  use time_state, only : time_state_t
  use compressible_ops_cpu, only : compressible_ops_cpu_update_uvw, &
       compressible_ops_cpu_update_mxyz_p_ruvw, &
       compressible_ops_cpu_update_e
  use compressible_ops_device, only : compressible_ops_device_update_uvw, &
       compressible_ops_device_update_mxyz_p_ruvw, &
       compressible_ops_device_update_e, &
       compressible_ops_device_update_temperature
  use mpi_f08, only : MPI_Allreduce, MPI_INTEGER, MPI_MAX
  use viscous_regularization, only : viscous_regularization_t, &
       viscous_regularization_factory
  use vector_bc_projector, only : vector_bc_projector_t, &
       coupled_vector_bc_projector_t
  implicit none
  private

  type, public, extends(fluid_scheme_compressible_t) &
       :: fluid_scheme_compressible_ns_t
     type(field_t) :: rho_res, m_x_res, m_y_res, m_z_res, m_E_res
     type(field_t) :: drho, dm_x, dm_y, dm_z, dE
     type(euler_idp_t) :: euler_idp_solver
     type(field_t) :: idp_entropy_fraction
     type(field_t), pointer :: idp_wave_speed => null()
     integer :: idp_viscosity_step = -1
     class(advection_t), allocatable :: adv
     class(ax_t), allocatable :: Ax
     class(ax_t), allocatable :: Ax_stress
     class(compressible_rhs_t), allocatable :: compressible_rhs
     type(runge_kutta_time_scheme_t) :: rk_scheme

     class(viscous_regularization_t), allocatable :: viscous_regularization
     logical :: if_viscous_regularization

     !> Boundary conditions projector for velocity constraints.
     type(coupled_vector_bc_projector_t):: bcs_vel_projector
     !> List of boundary conditions for density.
     type(bc_list_t) :: bcs_density

   contains
     procedure, pass(this) :: init => fluid_scheme_compressible_ns_init
     procedure, pass(this) :: free => fluid_scheme_compressible_ns_free
     procedure, pass(this) :: step => fluid_scheme_compressible_ns_step
     procedure, pass(this) :: restart => fluid_scheme_compressible_ns_restart
     procedure, pass(this) :: validate => fluid_scheme_compressible_ns_validate
     procedure, pass(this) :: compute_cfl => &
          fluid_scheme_compressible_ns_compute_cfl
     !> Set up boundary conditions.
     procedure, pass(this) :: setup_bcs &
          => fluid_scheme_compressible_ns_setup_bcs
     procedure, pass(this), private :: setup_viscous_regularization
     procedure, pass(this) :: prepare_idp_graph
     procedure, pass(this), private :: log_euler_idp_diagnostics => &
          fluid_scheme_compressible_ns_log_euler_idp_diagnostics
  end type fluid_scheme_compressible_ns_t

  interface
     !> Boundary condition factory for density.
     !! @details Will mark a mesh zone for the bc and finalize.
     !! @param[inout] object The object to be allocated.
     !! @param[in] scheme The `fluid_scheme_compressible_ns_t` scheme.
     !! @param[inout] json JSON object for initializing the bc.
     !! @param[in] coef SEM coefficients.
     !! @param[in] user The user interface.
     module subroutine density_bc_factory(object, scheme, json, coef, user)
       class(bc_t), pointer, intent(inout) :: object
       type(fluid_scheme_compressible_ns_t), intent(in) :: scheme
       type(json_file), intent(inout) :: json
       type(coef_t), intent(in) :: coef
       type(user_t), intent(in) :: user
     end subroutine density_bc_factory
  end interface

  interface
     !> Boundary condition factory for pressure.
     !! @details Will mark a mesh zone for the bc and finalize.
     !! @param[inout] object The object to be allocated.
     !! @param[in] scheme The `fluid_scheme_compressible_ns_t` scheme.
     !! @param[inout] json JSON object for initializing the bc.
     !! @param[in] coef SEM coefficients.
     !! @param[in] user The user interface.
     module subroutine pressure_bc_factory(object, scheme, json, coef, user)
       class(bc_t), pointer, intent(inout) :: object
       type(fluid_scheme_compressible_ns_t), intent(inout) :: scheme
       type(json_file), intent(inout) :: json
       type(coef_t), intent(in) :: coef
       type(user_t), intent(in) :: user
     end subroutine pressure_bc_factory
  end interface

  interface
     !> Boundary condition factory for velocity
     !! @details Will mark a mesh zone for the bc and finalize.
     !! @param[inout] object The object to be allocated.
     !! @param[in] scheme The `fluid_scheme_compressible_ns_t` scheme.
     !! @param[inout] json JSON object for initializing the bc.
     !! @param[in] coef SEM coefficients.
     !! @param[in] user The user interface.
     module subroutine velocity_bc_factory(object, scheme, json, coef, user)
       class(bc_t), pointer, intent(inout) :: object
       type(fluid_scheme_compressible_ns_t), intent(in) :: scheme
       type(json_file), intent(inout) :: json
       type(coef_t), intent(in) :: coef
       type(user_t), intent(in) :: user
     end subroutine velocity_bc_factory
  end interface

contains
  !> Validate the compressible state and initialise the IDP graph CFL.
  subroutine fluid_scheme_compressible_ns_validate(this)
    class(fluid_scheme_compressible_ns_t), target, intent(inout) :: this

    call fluid_scheme_compressible_validate(this)
    if (this%euler_idp%enabled) then
       call this%euler_idp_solver%invalidate()
       call this%euler_idp_solver%prepare(this%rho, this%m_x, this%m_y, &
            this%m_z, this%E, this%gs_Xh, this%bcs_density, this%bcs_vel, &
            this%bcs_prs, this%gamma)
    end if
  end subroutine fluid_scheme_compressible_ns_validate

  !> Return the graph CFL in Euler-IDP mode and the legacy CFL otherwise.
  function fluid_scheme_compressible_ns_compute_cfl(this, dt) result(cfl)
    class(fluid_scheme_compressible_ns_t), intent(in) :: this
    real(kind=dp), intent(in) :: dt
    real(kind=dp) :: cfl

    if (this%euler_idp%enabled) then
       cfl = this%euler_idp_solver%cfl(real(dt, rp))
    else
       cfl = fluid_scheme_compressible_compute_cfl(this, dt)
    end if
  end function fluid_scheme_compressible_ns_compute_cfl

  !> Initialize the compressible Navier-Stokes fluid scheme
  !! @param this The fluid scheme object
  !! @param msh Mesh data structure
  !! @param lx Polynomial order in x-direction
  !! @param params The fluid parameters, i.e. `case.fluid`
  !! @param numerics_params The numerical parameters, i.e. `case.numerics`
  !! @param user User-defined parameters and functions
  !! @param chkp Checkpoint to write to
  !! @param output_boundary Not supported by this scheme, ignored
  subroutine fluid_scheme_compressible_ns_init(this, msh, lx, params, &
       numerics_params, user, chkp, output_boundary)
    class(fluid_scheme_compressible_ns_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(in) :: lx
    type(json_file), target, intent(inout) :: params
    type(json_file), target, intent(inout) :: numerics_params
    type(user_t), target, intent(in) :: user
    type(chkp_t), target, intent(inout) :: chkp
    logical, intent(in) :: output_boundary
    character(len=12), parameter :: scheme = 'compressible'
    integer :: rk_order

    call this%free()

    ! Initialize base class
    call this%scheme_init(msh, lx, params, numerics_params, scheme, user)

    call compressible_rhs_factory(this%compressible_rhs, this%gamma)

    associate(Xh_lx => this%Xh%lx, Xh_ly => this%Xh%ly, Xh_lz => this%Xh%lz, &
         dm_Xh => this%dm_Xh, nelv => this%msh%nelv)

      call this%drho%init(dm_Xh, 'drho')
      call this%dm_x%init(dm_Xh, 'dm_x')
      call this%dm_y%init(dm_Xh, 'dm_y')
      call this%dm_z%init(dm_Xh, 'dm_z')
      call this%dE%init(dm_Xh, 'dE')

    end associate

    call json_get_or_default(numerics_params, 'time_order', rk_order, 4)

    if (this%euler_idp%enabled) then
       call this%idp_entropy_fraction%init(this%dm_Xh, "idp_entropy_fraction")
       call this%euler_idp_solver%init(this%dm_Xh, this%c_Xh, this%gs_Xh, &
            this%euler_idp, rk_order)
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       associate(p => this%p, rho => this%rho, &
            u => this%u, v => this%v, w => this%w, &
            m_x => this%m_x, m_y => this%m_y, m_z => this%m_z, &
            artificial_visc => this%artificial_visc)
         call device_memcpy(p%x, p%x_d, p%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(rho%x, rho%x_d, rho%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(u%x, u%x_d, u%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(v%x, v%x_d, v%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(w%x, w%x_d, w%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(m_x%x, m_x%x_d, m_x%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(m_y%x, m_y%x_d, m_y%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(m_z%x, m_z%x_d, m_z%dof%size(), &
              HOST_TO_DEVICE, sync = .false.)
         call device_memcpy(artificial_visc%x, artificial_visc%x_d, &
              artificial_visc%dof%size(), HOST_TO_DEVICE, sync = .false.)
       end associate
    end if

    ! Initialize the diffusion operators
    if (this%c_Xh%geo_compression) then
       call ax_helm_allocator(this%Ax, type_name = "standard_compr")
    else
       call ax_helm_allocator(this%Ax, type_name = "standard")
    end if
    call ax_helm_allocator(this%Ax_stress, type_name = "full")

    ! Initialize the velocity BC projector. The coupled projector builds the
    ! local bases required by non-axis aligned mixed velocity conditions.
    call this%bcs_vel_projector%init(this%c_Xh)

    ! Initialize viscous regularization
    call this%setup_viscous_regularization(params)

    ! Initialize Runge-Kutta scheme
    call this%rk_scheme%init(rk_order)

    call neko_log%section("Fluid boundary conditions")
    ! Set up boundary conditions
    call this%setup_bcs(user, params)
    call neko_log%end_section()

  end subroutine fluid_scheme_compressible_ns_init

  !> Free allocated memory and cleanup
  !> @param this The fluid scheme object to destroy
  subroutine fluid_scheme_compressible_ns_free(this)
    class(fluid_scheme_compressible_ns_t), intent(inout) :: this
    class(bc_t), pointer :: bc
    integer :: i

    if (allocated(this%Ax)) then
       call this%Ax%free()
       deallocate(this%Ax)
    end if

    if (allocated(this%Ax_stress)) then
       call this%Ax_stress%free()
       deallocate(this%Ax_stress)
    end if

    call this%scheme_free()

    if (allocated(this%compressible_rhs)) then
       deallocate(this%compressible_rhs)
    end if

    call this%drho%free()
    call this%dm_x%free()
    call this%dm_y%free()
    call this%dm_z%free()
    call this%dE%free()
    call this%idp_entropy_fraction%free()
    nullify(this%idp_wave_speed)
    this%idp_viscosity_step = -1
    call this%euler_idp_solver%free()
    call this%rk_scheme%free()

    if (allocated(this%viscous_regularization)) then
       call this%viscous_regularization%free()
       deallocate(this%viscous_regularization)
    end if

    do i = 1, this%bcs_density%size()
       bc => this%bcs_density%get(i)
       if (associated(bc)) then
          call bc%free()
          deallocate(bc)
       end if
    end do
    call this%bcs_density%free()
    call this%bcs_vel_projector%free()

  end subroutine fluid_scheme_compressible_ns_free

  !> Advance the fluid simulation one timestep
  !> @param this The fluid scheme object
  !> @param time Current simulation time state
  !> @param ext_bdf Time integration controller
  !> @param dt_controller Timestep size controller
  subroutine fluid_scheme_compressible_ns_step(this, time, dt_controller)
    class(fluid_scheme_compressible_ns_t), target, intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_step_controller_t), intent(in) :: dt_controller
    type(field_t), pointer :: temp
    integer :: temp_indices(1)
    ! number of degrees of freedom
    integer :: n
    integer :: i
    class(bc_t), pointer :: b

    if (this%euler_idp%enabled) then
       call fluid_scheme_compressible_ns_step_idp(this, time)
       return
    end if

    n = this%dm_Xh%size()
    call neko_scratch_registry%request_field(temp, temp_indices(1), .false.)
    b => null()

    call profiler_start_region('Fluid compressible', 1)
    associate(u => this%u, v => this%v, w => this%w, p => this%p, &
         m_x=> this%m_x, m_y => this%m_y, m_z => this%m_z, &
         Xh => this%Xh, msh => this%msh, Ax => this%Ax, &
         c_Xh => this%c_Xh, dm_Xh => this%dm_Xh, gs_Xh => this%gs_Xh, &
         E => this%E, rho => this%rho, &
         f_x => this%f_x, f_y => this%f_y, f_z => this%f_z, &
         drho => this%drho, dm_x => this%dm_x, dm_y => this%dm_y, &
         dm_z => this%dm_z, dE => this%dE, &
         compressible_rhs => this%compressible_rhs, &
         t => time%t, tstep => time%tstep, dt => time%dt, &
         rk_scheme => this%rk_scheme)

      !> Update artificial viscosity
      if (allocated(this%viscous_regularization)) then
         call this%viscous_regularization%update(this%artificial_visc)
      end if

      ! Refresh user-specified physical viscosity/conductivity before RHS.
      call this%update_material_properties(time)

      ! Execute RHS step with artificial viscosity field
      call compressible_rhs%step(rho, m_x, m_y, m_z, E, &
           p, u, v, w, this%Ax, &
           this%Ax_stress, c_Xh, gs_Xh, this%artificial_visc, this%mu, &
           this%kappa, this%bcs_vel, time, rk_scheme, real(dt, kind=rp))

      !> Apply density boundary conditions
      call this%bcs_density%apply(rho, time)

      !> Update variables
      ! Update u, v, w
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call compressible_ops_device_update_uvw(u%x_d, v%x_d, w%x_d, &
              m_x%x_d, m_y%x_d, m_z%x_d, rho%x_d, n)
      else
         call compressible_ops_cpu_update_uvw(u%x, v%x, w%x, &
              m_x%x, m_y%x, m_z%x, rho%x, n)
      end if

      !> Apply velocity boundary conditions
      call this%bcs_vel%apply_vector(u%x, v%x, w%x, &
           dm_Xh%size(), time, strong = .true.)

      !> Update m_x, m_y, m_z, p and temp (ruvw)
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call compressible_ops_device_update_mxyz_p_ruvw(m_x%x_d, m_y%x_d, &
              m_z%x_d, p%x_d, temp%x_d, u%x_d, v%x_d, w%x_d, E%x_d, &
              rho%x_d, this%gamma, n)
      else
         call compressible_ops_cpu_update_mxyz_p_ruvw(m_x%x, m_y%x, m_z%x, &
              p%x, temp%x, u%x, v%x, w%x, E%x, rho%x, this%gamma, n)
      end if

      !> Apply pressure boundary conditions
      call this%bcs_prs%apply(p, time)


      !> Update E
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call compressible_ops_device_update_e(E%x_d, p%x_d, &
              temp%x_d, this%gamma, n)
      else
         call compressible_ops_cpu_update_e(E%x, p%x, temp%x, this%gamma, n)
      end if

      !> Update temperature T = p / (rho * (gamma - 1))
      if (NEKO_BCKND_DEVICE .eq. 1) then
         call compressible_ops_device_update_temperature( &
              this%temperature%x_d, p%x_d, rho%x_d, this%gamma, n)
      else
         !OCL NORECURRENCE, NOVREC, NOALIAS
         !DIR$ CONCURRENT
         !DIR$ IVDEP
         !GCC$ ivdep
         !$omp parallel do simd
         do i = 1, n
            this%temperature%x(i,1,1,1) = p%x(i,1,1,1) / &
                 (rho%x(i,1,1,1) * (this%gamma - 1.0_rp))
         end do
         !$omp end parallel do simd
      end if

      !> Update maximum wave speed for CFL computation
      call this%compute_max_wave_speed()

      do i = 1, this%bcs_vel%size()
         b => this%bcs_vel%get(i)
         b%updated = .false.
      end do

      do i = 1, this%bcs_prs%size()
         b => this%bcs_prs%get(i)
         b%updated = .false.
      end do

      do i = 1, this%bcs_density%size()
         b => this%bcs_density%get(i)
         b%updated = .false.
      end do
      nullify(b)

    end associate
    call profiler_end_region('Fluid compressible', 1)

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine fluid_scheme_compressible_ns_step

  !> Advance one invariant-domain-preserving limited Euler step.
  subroutine fluid_scheme_compressible_ns_step_idp(this, time)
    class(fluid_scheme_compressible_ns_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    class(bc_t), pointer :: b
    integer :: i, n

    call profiler_start_region('Fluid Euler IDP', 1)
    if (this%idp_viscosity_step /= time%tstep) then
       call neko_error('Euler IDP requires an entropy_viscosity model ' // &
            'updating the configured regularization field every step')
    end if
    call this%viscous_regularization%update(this%artificial_visc)
    n = this%dm_Xh%size()
    call this%euler_idp_solver%advance(this%rho, this%m_x, this%m_y, &
         this%m_z, this%E, this%c_Xh, this%gs_Xh, this%bcs_density, &
         this%bcs_vel, this%bcs_prs, this%gamma, real(time%dt, rp), time, &
         this%euler_idp_diagnostics, this%idp_entropy_fraction, &
         this%idp_wave_speed)
    call this%log_euler_idp_diagnostics(time)

    call this%euler_idp_solver%copy_primitives(this%u, this%v, this%w, &
         this%p)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call compressible_ops_device_update_temperature( &
            this%temperature%x_d, this%p%x_d, this%rho%x_d, this%gamma, n)
    else
       do concurrent (i = 1:n)
          this%temperature%x(i,1,1,1) = this%p%x(i,1,1,1) / &
               (this%rho%x(i,1,1,1) * (this%gamma - 1.0_rp))
       end do
    end if
    call this%compute_max_wave_speed()
    if (.not. associated(this%idp_wave_speed)) call this%prepare_idp_graph(time)
    do i = 1, this%bcs_vel%size()
       b => this%bcs_vel%get(i)
       b%updated = .false.
    end do
    do i = 1, this%bcs_prs%size()
       b => this%bcs_prs%get(i)
       b%updated = .false.
    end do
    do i = 1, this%bcs_density%size()
       b => this%bcs_density%get(i)
       b%updated = .false.
    end do
    nullify(b)
    call profiler_end_region('Fluid Euler IDP', 1)
  end subroutine fluid_scheme_compressible_ns_step_idp

  !> Prepare the graph CFL using the wave speeds for the updated state.
  subroutine prepare_idp_graph(this, time)
    class(fluid_scheme_compressible_ns_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    call this%euler_idp_solver%prepare(this%rho, this%m_x, this%m_y, &
         this%m_z, this%E, this%gs_Xh, this%bcs_density, this%bcs_vel, &
         this%bcs_prs, this%gamma, time, this%idp_wave_speed)
  end subroutine prepare_idp_graph

  !> Log limiter activity for each stage of one completed IDP timestep.
  subroutine fluid_scheme_compressible_ns_log_euler_idp_diagnostics(this, &
       time)
    class(fluid_scheme_compressible_ns_t), intent(in) :: this
    type(time_state_t), intent(in) :: time
    character(len=2*LOG_SIZE) :: log_buf
    integer :: nstages, stage

    if (this%euler_idp%diagnostics_level .ne. &
         EULER_IDP_DIAGNOSTICS_FULL) return
    if (mod(time%tstep, this%euler_idp%diagnostics_interval) .ne. 0) return
    nstages = 1
    if (this%rk_scheme%order .eq. 3) nstages = 3
    do stage = 1, nstages
       associate(diag => this%euler_idp_solver%stage_diagnostics(stage))
         write(log_buf, '(A,I0,A,I0,A,ES10.3,A,ES10.3)') &
            'IDP limiter: step = ', time%tstep, ' stage = ', stage, &
            ' amin = ', diag%min_limiter, &
            ' frac = ', diag%limited_edge_fraction
         call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)

         write(log_buf, '(A,ES10.3,A,ES10.3,A,ES10.3)') &
            'IDP state: rho_min = ', diag%min_density, &
            ' p_min = ', diag%min_pressure, &
            ' e_min = ', diag%min_internal_energy
         call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)

         write(log_buf, '(A,ES10.3,A,ES10.3,A,ES10.3)') &
            'IDP graph: t = ', diag%stage_time, &
            ' wave = ', diag%max_graph_wave_speed, &
            ' rate = ', diag%max_graph_rate
         call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)

         write(log_buf, '(A,ES10.3,A,ES10.3,A,ES10.3)') &
            'IDP timestep: nodal_wave = ', diag%max_nodal_wave_speed, &
            ' graph_cfl = ', diag%max_graph_cfl, &
            ' dt_floor = ', diag%maximum_floor_timestep
         call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)

         write(log_buf, '(A,3(A,ES10.3))') 'IDP boundary pre:', &
            ' rho = ', diag%before_boundary%min_density, &
            ' p = ', diag%before_boundary%min_pressure, &
            ' e = ', diag%before_boundary%min_internal_energy
         call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)

         write(log_buf, '(A,3(A,ES10.3))') 'IDP boundary pre:', &
            ' nodal = ', diag%before_boundary%max_nodal_wave_speed, &
            ' graph = ', diag%before_boundary%max_graph_wave_speed, &
            ' rate = ', diag%before_boundary%max_graph_rate
         call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)

         write(log_buf, '(A,3(A,ES10.3))') 'IDP boundary post:', &
            ' rho = ', diag%after_boundary%min_density, &
            ' p = ', diag%after_boundary%min_pressure, &
            ' e = ', diag%after_boundary%min_internal_energy
         call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)

         write(log_buf, '(A,3(A,ES10.3))') 'IDP boundary post:', &
            ' nodal = ', diag%after_boundary%max_nodal_wave_speed, &
            ' graph = ', diag%after_boundary%max_graph_wave_speed, &
            ' rate = ', diag%after_boundary%max_graph_rate
         call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)

         write(log_buf, '(A,I0,A,I0,A,I0)') &
            'IDP edges: rho = ', diag%density_limited_edges, &
            ' energy = ', diag%internal_energy_limited_edges, &
            ' entropy = ', diag%entropy_limited_edges
         call neko_log%message(log_buf, lvl = NEKO_LOG_VERBOSE)
       end associate
    end do
  end subroutine fluid_scheme_compressible_ns_log_euler_idp_diagnostics

  !> Set up boundary conditions for the fluid scheme
  !> @param this The fluid scheme object
  !> @param user User-defined boundary conditions
  !> @param params Configuration parameters
  subroutine fluid_scheme_compressible_ns_setup_bcs(this, user, params)
    class(fluid_scheme_compressible_ns_t), target, intent(inout) :: this
    type(user_t), target, intent(in) :: user
    type(json_file), intent(inout) :: params
    integer :: i, n_bcs, zone_index, j, zone_size, global_zone_size, ierr
    class(bc_t), pointer :: bc_i
    type(json_core) :: core
    type(json_value), pointer :: bc_object
    type(json_file) :: bc_subdict
    logical :: found
    integer, allocatable :: zone_indices(:)
    character(len=LOG_SIZE) :: log_buf

    ! Process boundary conditions
    if (params%valid_path('boundary_conditions')) then
       call params%info('boundary_conditions', n_children = n_bcs)
       call params%get_core(core)
       call params%get('boundary_conditions', bc_object, found)

       !
       ! Velocity bcs
       !
       call this%bcs_vel%init(n_bcs)

       do i = 1, n_bcs
          ! Extract BC configuration
          call json_extract_item(core, bc_object, i, bc_subdict)
          call json_get(bc_subdict, "zone_indices", zone_indices)
          ! Validate zones
          do j = 1, size(zone_indices)
             zone_size = this%msh%labeled_zones(zone_indices(j))%size
             call MPI_Allreduce(zone_size, global_zone_size, 1, &
                  MPI_INTEGER, MPI_MAX, NEKO_COMM, ierr)

             if (global_zone_size .eq. 0) then
                write(log_buf, '(A,I0,A)') "Error: Zone ", zone_indices(j), &
                     " has zero size"
                call neko_error(log_buf)
             end if
          end do

          ! Create BC
          bc_i => null()
          call velocity_bc_factory(bc_i, this, bc_subdict, this%c_Xh, user)

          ! Add to appropriate lists
          if (associated(bc_i)) then
             call this%bcs_vel_projector%mark(bc_i)
             call this%bcs_vel%append(bc_i)
          end if
       end do

       !
       ! Pressure bcs
       !
       call this%bcs_prs%init(n_bcs)

       do i = 1, n_bcs
          ! Create a new json containing just the subdict for this bc
          call json_extract_item(core, bc_object, i, bc_subdict)
          bc_i => null()
          call pressure_bc_factory(bc_i, this, bc_subdict, this%c_Xh, user)

          ! Not all bcs require an allocation for pressure in particular,
          ! so we check.
          if (associated(bc_i)) then
             call this%bcs_prs%append(bc_i)
          end if
       end do

       !
       ! Density bcs
       !
       call this%bcs_density%init(n_bcs)

       do i = 1, n_bcs
          ! Create a new json containing just the subdict for this bc
          call json_extract_item(core, bc_object, i, bc_subdict)
          bc_i => null()
          call density_bc_factory(bc_i, this, bc_subdict, this%c_Xh, user)

          ! Not all bcs require an allocation for pressure in particular,
          ! so we check.
          if (associated(bc_i)) then
             call this%bcs_density%append(bc_i)
          end if
       end do
    else
       ! Check that there are no labeled zones, i.e. all are periodic.
       do i = 1, size(this%msh%labeled_zones)
          if (this%msh%labeled_zones(i)%size .gt. 0) then
             call neko_error("No boundary_conditions entry in the case file!")
          end if
       end do

       ! For a pure periodic case, we still need to initilise the bc lists
       ! to a zero size to avoid issues with apply() in step()
       call this%bcs_prs%init()
       call this%bcs_vel%init()
       call this%bcs_density%init()

    end if

    call this%bcs_vel_projector%finalize(rebuild_mask = .true.)
  end subroutine fluid_scheme_compressible_ns_setup_bcs

  !> Restart the simulation from saved state
  !! @param this The fluid scheme object
  !! @param dtlag Previous timestep sizes
  !! @param tlag Previous time values
  subroutine fluid_scheme_compressible_ns_restart(this, chkp)
    class(fluid_scheme_compressible_ns_t), target, intent(inout) :: this
    type(chkp_t), intent(inout) :: chkp
  end subroutine fluid_scheme_compressible_ns_restart

  subroutine setup_viscous_regularization(this, params)
    class(fluid_scheme_compressible_ns_t), target, intent(inout) :: this
    type(json_file), intent(inout) :: params
    type(json_file) :: reg_json
    logical :: found
    character(len=:), allocatable :: viscous_regularization_type

    found = .false.
    if (params%valid_path('viscous_regularization')) then
       call json_get(params, 'viscous_regularization', &
            reg_json)
       found = .true.
    end if

    if (.not. found) then
       ! No viscous_regularization specified, so we skip setup
       this%if_viscous_regularization = .false.
       return
    else
       this%if_viscous_regularization = .true.
    end if

    call json_get(reg_json, 'type', viscous_regularization_type)
    call viscous_regularization_factory(this%viscous_regularization, &
         viscous_regularization_type, &
         reg_json, this%c_Xh, this%dm_Xh)

    call reg_json%destroy()

    if (allocated(viscous_regularization_type)) &
         deallocate(viscous_regularization_type)

  end subroutine setup_viscous_regularization

end module fluid_scheme_compressible_ns
