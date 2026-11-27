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
!> CPU implementation of the Euler IDP solver.
module euler_idp_cpu
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mpi_f08, only : MPI_Allreduce, MPI_INTEGER, MPI_MAX, MPI_MIN, MPI_SUM
  use num_types, only : rp
  use field, only : field_t
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t
  use bc_list, only : bc_list_t
  use time_state, only : time_state_t
  use gs_ops, only : GS_OP_ADD, GS_OP_MIN, GS_OP_MAX
  use bc, only : bc_t
  use compressible_ops_cpu, only : &
       compressible_ops_cpu_conserved_to_primitive, &
       compressible_ops_cpu_update_uvw, &
       compressible_ops_cpu_update_mxyz_p_ruvw, &
       compressible_ops_cpu_update_e, EULER_STATE_OK
  use euler_idp_backend, only : EULER_IDP_NCOMP, euler_idp_backend_t, &
       euler_idp_diagnostics_t, euler_idp_state_observation_t, &
       EULER_IDP_DIAGNOSTICS_OFF, EULER_IDP_DIAGNOSTICS_FULL, &
       euler_idp_maximum_wave_speed, &
       euler_idp_flux_dot_vector, euler_idp_bar_state, &
       euler_idp_internal_energy, euler_idp_internal_energy_timestep, &
       euler_idp_limit_edge, &
       euler_idp_local_entropy_bounds, euler_idp_relax_density_bounds, &
       euler_idp_relax_entropy_bound, &
       euler_idp_specific_entropy, euler_idp_entropy_tolerance, &
       euler_idp_entropy_is_admissible
  use euler_gll_graph, only : euler_gll_graph_t
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION, pe_rank
  use utils, only : neko_error
  use logger, only : LOG_SIZE
  use profiler, only : profiler_start_region, profiler_end_region
  implicit none
  private

  !> Persistent work fields for the CPU Euler IDP implementation.
  type, public, extends(euler_idp_backend_t) :: euler_idp_cpu_t
     logical :: periodic_graph = .false.
     logical :: low_order_only = .false.
     logical :: relax_density_bounds = .false.
     logical :: relax_entropy_bounds = .false.
     logical :: limit_internal_energy = .true.
     logical :: limit_entropy = .true.
     real(kind=rp) :: density_bound_relaxation_factor = 1.0_rp
     real(kind=rp) :: entropy_bound_relaxation_factor = 1.0_rp
     real(kind=rp) :: correction_tolerance = 1.0e-10_rp
     type(euler_gll_graph_t) :: graph
     type(field_t) :: local_residual(EULER_IDP_NCOMP)
     type(field_t) :: low_local_residual(EULER_IDP_NCOMP)
     type(field_t) :: low_candidate(EULER_IDP_NCOMP)
     type(field_t) :: saved_state(EULER_IDP_NCOMP)
     type(field_t) :: viscosity_sum
     type(field_t) :: density_lower_bound
     type(field_t) :: density_upper_bound
     type(field_t) :: entropy_lower_bound
     type(field_t) :: u
     type(field_t) :: v
     type(field_t) :: w
     type(field_t) :: p
     type(field_t) :: sound_speed
     type(field_t) :: internal_energy
     type(field_t) :: flux_x
     type(field_t) :: flux_y
     type(field_t) :: flux_z
     real(kind=rp), allocatable :: element_residual_sum(:,:)
     real(kind=rp), allocatable :: element_boundary_flux(:,:)
     real(kind=rp), allocatable :: low_element_residual_sum(:,:)
     real(kind=rp), allocatable :: low_element_boundary_flux(:,:)
     real(kind=rp), allocatable :: edge_viscosity(:)
     real(kind=rp), allocatable :: correction_flux(:,:)
     integer, allocatable :: state_status(:)
     real(kind=rp) :: domain_volume = 0.0_rp
     real(kind=rp) :: density_relaxation_mass = 0.0_rp
     real(kind=rp) :: maximum_floor_timestep = huge(1.0_rp)
     real(kind=rp) :: limiter_weight_error = 0.0_rp
   contains
     procedure, pass(this) :: init => euler_idp_cpu_init
     procedure, pass(this) :: init_graph => euler_idp_cpu_init_graph
     procedure, pass(this) :: free => euler_idp_cpu_free
     procedure, pass(this) :: update_graph_viscosity => &
          euler_idp_cpu_update_graph_viscosity
     procedure, pass(this) :: prepare_stage => euler_idp_cpu_prepare_stage
     procedure, pass(this) :: evaluate_low_order => &
          euler_idp_cpu_evaluate_low_order
     procedure, pass(this) :: graph_cfl => euler_idp_cpu_graph_cfl
     procedure, pass(this) :: compute_bounds => euler_idp_cpu_compute_bounds
     procedure, pass(this) :: forward_euler => euler_idp_cpu_forward_euler
     procedure, pass(this) :: apply_boundary_conditions => &
          euler_idp_cpu_apply_boundary_conditions
     procedure, pass(this) :: apply_candidate_boundary => &
          euler_idp_cpu_apply_candidate_boundary
     procedure, pass(this) :: observe_candidate => &
          euler_idp_cpu_observe_candidate
     procedure, pass(this) :: validate_candidate => &
          euler_idp_cpu_validate_candidate
     procedure, pass(this) :: save_state => euler_idp_cpu_save_state
     procedure, pass(this) :: combine_stage => euler_idp_cpu_combine_stage
     procedure, pass(this) :: compute_limiter => &
          euler_idp_cpu_compute_limiter
     procedure, pass(this) :: apply_correction => &
          euler_idp_cpu_apply_correction
     procedure, pass(this) :: copy_primitives => &
          euler_idp_cpu_copy_primitives
  end type euler_idp_cpu_t

contains

  !> Initialise persistent fields.
  subroutine euler_idp_cpu_init(this, dof)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(dofmap_t), target, intent(in) :: dof
    character(len=32) :: name
    integer :: component

    call this%free()

    do component = 1, EULER_IDP_NCOMP
       write(name, '(A,I0)') 'euler_idp_low_candidate_', component
       call this%low_candidate(component)%init(dof, trim(name))
    end do

    call this%u%init(dof, 'euler_idp_u')
    call this%v%init(dof, 'euler_idp_v')
    call this%w%init(dof, 'euler_idp_w')
    call this%p%init(dof, 'euler_idp_p')
    call this%sound_speed%init(dof, 'euler_idp_sound_speed')
    call this%internal_energy%init(dof, 'euler_idp_internal_energy')
    call this%flux_x%init(dof, 'euler_idp_flux_x')
    call this%flux_y%init(dof, 'euler_idp_flux_y')
    call this%flux_z%init(dof, 'euler_idp_flux_z')
    call this%viscosity_sum%init(dof, 'euler_idp_viscosity_sum')
    call this%density_lower_bound%init(dof, 'euler_idp_density_lower_bound')
    call this%density_upper_bound%init(dof, 'euler_idp_density_upper_bound')
    allocate(this%state_status(dof%size()))
    this%state_status = EULER_STATE_OK
    this%max_graph_rate = 0.0_rp
    this%max_graph_wave_speed = 0.0_rp
    this%domain_volume = 0.0_rp
    this%maximum_graph_timestep = huge(1.0_rp)
    this%maximum_floor_timestep = huge(1.0_rp)
    this%limiter_weight_error = 0.0_rp
    this%low_order_only = .false.
    this%relax_density_bounds = .false.
    this%relax_entropy_bounds = .false.
    this%limit_internal_energy = .true.
    this%limit_entropy = .true.
    this%density_bound_relaxation_factor = 1.0_rp
    this%entropy_bound_relaxation_factor = 1.0_rp
    this%initialized = .true.
  end subroutine euler_idp_cpu_init

  !> Initialise the sparse GLL graph after coefficients are available.
  subroutine euler_idp_cpu_init_graph(this, coef, gs, relax_density_bounds, &
       relax_entropy_bounds, low_order_only, limit_internal_energy, &
       limit_entropy, density_bound_relaxation_factor, &
       entropy_bound_relaxation_factor, time_order, diagnostics_level, &
       correction_tolerance)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    logical, intent(in) :: relax_density_bounds, relax_entropy_bounds
    logical, intent(in) :: low_order_only
    logical, intent(in) :: limit_internal_energy, limit_entropy
    real(kind=rp), intent(in) :: density_bound_relaxation_factor
    real(kind=rp), intent(in) :: entropy_bound_relaxation_factor
    real(kind=rp), intent(in) :: correction_tolerance
    integer, intent(in) :: time_order, diagnostics_level
    real(kind=rp) :: local_error, global_error
    real(kind=rp) :: local_mass, global_mass
    character(len=32) :: name
    integer :: component, direction, edge, ierr

    if (.not. this%initialized) then
       call neko_error('Euler IDP CPU object is not initialised')
    end if
    this%correction_tolerance = correction_tolerance
    this%relax_density_bounds = .false.
    this%relax_entropy_bounds = .false.
    this%low_order_only = .false.
    this%limit_internal_energy = .true.
    this%limit_entropy = .true.
    this%density_bound_relaxation_factor = 1.0_rp
    this%entropy_bound_relaxation_factor = 1.0_rp
    this%relax_density_bounds = relax_density_bounds
    this%relax_entropy_bounds = relax_entropy_bounds
    this%low_order_only = low_order_only
    this%limit_internal_energy = limit_internal_energy
    this%limit_entropy = limit_entropy
    this%density_bound_relaxation_factor = &
         density_bound_relaxation_factor
    this%entropy_bound_relaxation_factor = &
         entropy_bound_relaxation_factor
    this%diagnostics_level = diagnostics_level
    do component = 1, EULER_IDP_NCOMP
       call this%local_residual(component)%free()
       call this%saved_state(component)%free()
       if (diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL) then
          write(name, '(A,I0)') 'euler_idp_local_', component
          call this%local_residual(component)%init(coef%dof, trim(name))
       end if
       if (time_order .eq. 3) then
          write(name, '(A,I0)') 'euler_idp_saved_state_', component
          call this%saved_state(component)%init(coef%dof, trim(name))
       end if
    end do
    call this%entropy_lower_bound%free()
    if (limit_entropy) then
       call this%entropy_lower_bound%init(coef%dof, &
            'euler_idp_entropy_lower_bound')
    end if
    if (allocated(this%element_residual_sum)) then
       deallocate(this%element_residual_sum)
    end if
    if (allocated(this%element_boundary_flux)) then
       deallocate(this%element_boundary_flux)
    end if
    if (allocated(this%low_element_residual_sum)) then
       deallocate(this%low_element_residual_sum)
    end if
    if (allocated(this%low_element_boundary_flux)) then
       deallocate(this%low_element_boundary_flux)
    end if
    if (diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL) then
       allocate(this%element_residual_sum(EULER_IDP_NCOMP, coef%msh%nelv))
       allocate(this%element_boundary_flux(EULER_IDP_NCOMP, coef%msh%nelv))
       allocate(this%low_element_residual_sum(EULER_IDP_NCOMP, coef%msh%nelv))
       allocate(this%low_element_boundary_flux(EULER_IDP_NCOMP, coef%msh%nelv))
       this%element_residual_sum = 0.0_rp
       this%element_boundary_flux = 0.0_rp
       this%low_element_residual_sum = 0.0_rp
       this%low_element_boundary_flux = 0.0_rp
    end if
    if (allocated(this%edge_viscosity)) deallocate(this%edge_viscosity)
    if (allocated(this%correction_flux)) deallocate(this%correction_flux)
    call this%graph%init(coef, gs)
    this%domain_volume = coef%volume
    if (.not. ieee_is_finite(this%domain_volume) .or. &
         this%domain_volume .le. 0.0_rp) then
       call neko_error('Euler IDP requires positive domain volume')
    end if
    local_mass = huge(1.0_rp)
    if (this%graph%mass%size() .gt. 0) then
       local_mass = minval(this%graph%mass%x)
    end if
    call MPI_Allreduce(local_mass, global_mass, 1, MPI_REAL_PRECISION, &
         MPI_MIN, NEKO_COMM, ierr)
    if (.not. ieee_is_finite(global_mass) .or. global_mass .le. 0.0_rp .or. &
         global_mass .ge. this%domain_volume) then
       call neko_error('Euler IDP density relaxation requires a valid ' // &
            'minimum nodal mass')
    end if
    this%density_relaxation_mass = global_mass
    this%periodic_graph = this%graph%periodic_facets_covered(coef)
    do component = 1, EULER_IDP_NCOMP
       call this%low_local_residual(component)%free()
       if (.not. this%periodic_graph .or. &
            diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL) then
          write(name, '(A,I0)') 'euler_idp_low_local_', component
          call this%low_local_residual(component)%init(coef%dof, trim(name))
       end if
    end do
    allocate(this%edge_viscosity(this%graph%n_edges))
    allocate(this%correction_flux(EULER_IDP_NCOMP, this%graph%n_edges))
    this%edge_viscosity = 0.0_rp
    this%correction_flux = 0.0_rp

    local_error = 0.0_rp
    if (this%graph%mass%size() .gt. 0) then
       local_error = maxval(abs(this%graph%mass%x * coef%Binv - 1.0_rp))
    end if
    call MPI_Allreduce(local_error, global_error, 1, MPI_REAL_PRECISION, &
         MPI_MAX, NEKO_COMM, ierr)
    if (global_error .gt. 128.0_rp * epsilon(1.0_rp)) then
       call neko_error('Euler IDP assembled mass and inverse are inconsistent')
    end if

    this%flux_x%x = 0.0_rp
    do edge = 1, this%graph%n_edges
       associate(a => this%graph%left(:, edge), &
            b => this%graph%right(:, edge))
         direction = this%graph%direction(edge)
         ! Give every tensor direction the same total nodal weight. This
         ! preserves one-dimensional states across replicated element nodes.
         this%flux_x%x(a(1),a(2),a(3),a(4)) = &
              this%flux_x%x(a(1),a(2),a(3),a(4)) + &
              1.0_rp / (real(this%graph%n_directions, rp) * &
              this%graph%directional_degree(direction)%x( &
              a(1),a(2),a(3),a(4)))
         this%flux_x%x(b(1),b(2),b(3),b(4)) = &
              this%flux_x%x(b(1),b(2),b(3),b(4)) + &
              1.0_rp / (real(this%graph%n_directions, rp) * &
              this%graph%directional_degree(direction)%x( &
              b(1),b(2),b(3),b(4)))
       end associate
    end do
    call gs%op(this%flux_x, GS_OP_ADD)
    local_error = 0.0_rp
    if (this%flux_x%size() .gt. 0) then
       local_error = maxval(abs(this%flux_x%x - 1.0_rp))
    end if
    call MPI_Allreduce(local_error, global_error, 1, MPI_REAL_PRECISION, &
         MPI_MAX, NEKO_COMM, ierr)
    this%limiter_weight_error = global_error
    if (global_error .gt. 128.0_rp * epsilon(1.0_rp)) then
       call neko_error('Euler IDP limiter occurrence weights do not sum to one')
    end if
  end subroutine euler_idp_cpu_init_graph

  !> Free persistent fields.
  subroutine euler_idp_cpu_free(this)
    class(euler_idp_cpu_t), intent(inout) :: this
    integer :: component

    call this%graph%free()
    do component = 1, EULER_IDP_NCOMP
       call this%local_residual(component)%free()
       call this%low_local_residual(component)%free()
       call this%low_candidate(component)%free()
       call this%saved_state(component)%free()
    end do
    call this%u%free()
    call this%v%free()
    call this%w%free()
    call this%p%free()
    call this%sound_speed%free()
    call this%internal_energy%free()
    call this%flux_x%free()
    call this%flux_y%free()
    call this%flux_z%free()
    call this%viscosity_sum%free()
    call this%density_lower_bound%free()
    call this%density_upper_bound%free()
    call this%entropy_lower_bound%free()
    if (allocated(this%element_residual_sum)) then
       deallocate(this%element_residual_sum)
    end if
    if (allocated(this%element_boundary_flux)) then
       deallocate(this%element_boundary_flux)
    end if
    if (allocated(this%low_element_residual_sum)) then
       deallocate(this%low_element_residual_sum)
    end if
    if (allocated(this%low_element_boundary_flux)) then
       deallocate(this%low_element_boundary_flux)
    end if
    if (allocated(this%edge_viscosity)) deallocate(this%edge_viscosity)
    if (allocated(this%correction_flux)) deallocate(this%correction_flux)
    if (allocated(this%state_status)) deallocate(this%state_status)
    this%initialized = .false.
    this%periodic_graph = .false.
    this%low_order_only = .false.
    this%relax_density_bounds = .false.
    this%relax_entropy_bounds = .false.
    this%limit_internal_energy = .true.
    this%limit_entropy = .true.
    this%density_bound_relaxation_factor = 1.0_rp
    this%entropy_bound_relaxation_factor = 1.0_rp
    this%correction_tolerance = 1.0e-10_rp
    this%max_graph_rate = 0.0_rp
    this%max_graph_wave_speed = 0.0_rp
    this%domain_volume = 0.0_rp
    this%density_relaxation_mass = 0.0_rp
    this%maximum_graph_timestep = huge(1.0_rp)
    this%maximum_floor_timestep = huge(1.0_rp)
    this%limiter_weight_error = 0.0_rp
  end subroutine euler_idp_cpu_free

  !> Update symmetric graph viscosity and the graph CFL rate.
  subroutine euler_idp_cpu_update_graph_viscosity(this, rho, m_x, m_y, m_z, &
       energy, gs, gamma, graph_wave_speed)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(field_t), intent(in) :: rho, m_x, m_y, m_z, energy
    type(gs_t), intent(inout) :: gs
    real(kind=rp), intent(in) :: gamma
    type(field_t), intent(in), optional :: graph_wave_speed
    real(kind=rp) :: left_state(EULER_IDP_NCOMP)
    real(kind=rp) :: right_state(EULER_IDP_NCOMP)
    real(kind=rp) :: left_flux(EULER_IDP_NCOMP)
    real(kind=rp) :: right_flux(EULER_IDP_NCOMP)
    real(kind=rp) :: bar_state(EULER_IDP_NCOMP)
    real(kind=rp) :: normal(3), coefficient_norm, flux_scale
    real(kind=rp) :: local_maximum(2), global_maximum(2)
    real(kind=rp) :: left_wave_speed, right_wave_speed, edge_wave_speed
    integer :: edge, ierr
    logical :: validate_bar_state

    if (.not. this%graph%initialized) then
       call neko_error('Euler IDP graph is not initialised')
    end if

    call profiler_start_region('Euler IDP graph viscosity')
    this%viscosity_sum%x = 0.0_rp
    local_maximum = 0.0_rp
    validate_bar_state = &
         this%diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL
    do edge = 1, this%graph%n_edges
       associate(a => this%graph%left(:, edge), &
            b => this%graph%right(:, edge), &
            coefficient => this%graph%coefficient(:, edge))
         if (.not. present(graph_wave_speed) .or. validate_bar_state) then
            left_state = [rho%x(a(1),a(2),a(3),a(4)), &
                 m_x%x(a(1),a(2),a(3),a(4)), &
                 m_y%x(a(1),a(2),a(3),a(4)), &
                 m_z%x(a(1),a(2),a(3),a(4)), &
                 energy%x(a(1),a(2),a(3),a(4))]
            right_state = [rho%x(b(1),b(2),b(3),b(4)), &
                 m_x%x(b(1),b(2),b(3),b(4)), &
                 m_y%x(b(1),b(2),b(3),b(4)), &
                 m_z%x(b(1),b(2),b(3),b(4)), &
                 energy%x(b(1),b(2),b(3),b(4))]
         end if
         coefficient_norm = this%graph%coefficient_norm(edge)
         normal = coefficient / coefficient_norm
         if (present(graph_wave_speed)) then
            left_wave_speed = graph_wave_speed%x(a(1),a(2),a(3),a(4))
            right_wave_speed = graph_wave_speed%x(b(1),b(2),b(3),b(4))
            if (.not. ieee_is_finite(left_wave_speed) .or. &
                 .not. ieee_is_finite(right_wave_speed) .or. &
                 left_wave_speed .lt. 0.0_rp .or. &
                 right_wave_speed .lt. 0.0_rp) then
               call neko_error('User graph wave speed must be finite and ' // &
                    'nonnegative')
            end if
            edge_wave_speed = max(left_wave_speed, right_wave_speed)
         else
            edge_wave_speed = euler_idp_maximum_wave_speed( &
                 left_state, right_state, normal, gamma)
         end if
         this%edge_viscosity(edge) = coefficient_norm * &
              edge_wave_speed
         local_maximum(2) = max(local_maximum(2), edge_wave_speed)
         if (validate_bar_state) then
            call euler_idp_flux_dot_vector(left_state, coefficient, gamma, &
                 left_flux)
            call euler_idp_flux_dot_vector(right_state, coefficient, gamma, &
                 right_flux)
            if (this%edge_viscosity(edge) .le. tiny(1.0_rp)) then
               flux_scale = max(1.0_rp, maxval(abs(left_flux)), &
                    maxval(abs(right_flux)))
               if (maxval(abs(right_flux - left_flux)) .gt. &
                    128.0_rp * epsilon(1.0_rp) * flux_scale) then
                  call neko_error('Zero graph viscosity has a nonzero ' // &
                       'projected flux jump')
               end if
            end if
            call euler_idp_bar_state(left_state, right_state, &
                 right_flux - left_flux, this%edge_viscosity(edge), &
                 bar_state)
            if (bar_state(1) .le. 0.0_rp) then
               call neko_error('Euler IDP graph viscosity produced a bar ' // &
                    'state with nonpositive density')
            end if
            if ((this%limit_internal_energy .or. this%limit_entropy) .and. &
                 euler_idp_internal_energy(bar_state) .le. 0.0_rp) then
               call neko_error('Euler IDP graph viscosity produced a bar ' // &
                    'state with nonpositive internal energy')
            end if
         end if
         this%viscosity_sum%x(a(1),a(2),a(3),a(4)) = &
              this%viscosity_sum%x(a(1),a(2),a(3),a(4)) + &
              this%edge_viscosity(edge)
         this%viscosity_sum%x(b(1),b(2),b(3),b(4)) = &
              this%viscosity_sum%x(b(1),b(2),b(3),b(4)) + &
              this%edge_viscosity(edge)
       end associate
    end do
    call profiler_start_region('Euler IDP gather-scatter')
    call gs%op(this%viscosity_sum, GS_OP_ADD)
    call profiler_end_region('Euler IDP gather-scatter')

    if (rho%size() .gt. 0) then
       local_maximum(1) = maxval(2.0_rp * this%viscosity_sum%x / &
            this%graph%mass%x)
    end if
    call profiler_start_region('Euler IDP MPI reduction')
    call MPI_Allreduce(local_maximum, global_maximum, 2, &
         MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM, ierr)
    call profiler_end_region('Euler IDP MPI reduction')
    this%max_graph_rate = global_maximum(1)
    this%max_graph_wave_speed = global_maximum(2)
    if (global_maximum(1) .gt. 0.0_rp) then
       this%maximum_graph_timestep = 1.0_rp / global_maximum(1)
    else
       this%maximum_graph_timestep = huge(1.0_rp)
    end if
    call profiler_end_region('Euler IDP graph viscosity')
  end subroutine euler_idp_cpu_update_graph_viscosity

  !> Prepare primitive and graph data for one immutable stage state.
  subroutine euler_idp_cpu_prepare_stage(this, rho, m_x, m_y, m_z, energy, &
       gs, gamma, internal_energy_floor, primitives_valid, graph_wave_speed)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(field_t), intent(in) :: rho, m_x, m_y, m_z, energy
    type(gs_t), intent(inout) :: gs
    real(kind=rp), intent(in) :: gamma, internal_energy_floor
    logical, intent(in) :: primitives_valid
    type(field_t), intent(in), optional :: graph_wave_speed

    if (.not. primitives_valid) then
       call euler_idp_cpu_primitives(this, rho, m_x, m_y, m_z, energy, &
            gamma, internal_energy_floor, 'prepared stage state')
    end if
    call this%update_graph_viscosity(rho, m_x, m_y, m_z, energy, gs, gamma, &
         graph_wave_speed)
  end subroutine euler_idp_cpu_prepare_stage

  !> Evaluate the conservative low-order graph residual into candidate storage.
  subroutine euler_idp_cpu_evaluate_low_order(this, rho, m_x, m_y, m_z, &
       energy, coef, gs, gamma, internal_energy_floor, &
       diagnostics, directional_error)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(field_t), intent(in) :: rho, m_x, m_y, m_z, energy
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    real(kind=rp), intent(in) :: gamma, internal_energy_floor
    type(euler_idp_diagnostics_t), intent(inout) :: diagnostics
    real(kind=rp), intent(out) :: directional_error(EULER_IDP_NCOMP)
    real(kind=rp) :: difference, left_value, right_value
    real(kind=rp) :: left_flux_value, right_flux_value, flux_difference
    integer :: component, edge
    logical :: need_local_residual

    call profiler_start_region('Euler IDP low-order graph')

    directional_error = 0.0_rp
    need_local_residual = .not. this%periodic_graph .or. &
         this%diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL
    do component = 1, EULER_IDP_NCOMP
       call euler_idp_cpu_flux(this, component, m_x, m_y, m_z, energy)
       if (this%diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL) then
          call profiler_start_region('Euler IDP high-order')
          call this%graph%high_order_flux_divergence( &
               this%local_residual(component)%x, this%flux_x%x, &
               this%flux_y%x, this%flux_z%x, coef)
          call euler_idp_cpu_element_sums(this, component, coef, &
               this%local_residual(component), this%element_residual_sum, &
               this%element_boundary_flux)
          call profiler_end_region('Euler IDP high-order')
       end if
       if (need_local_residual) then
          call this%graph%local_flux_divergence( &
               this%low_local_residual(component)%x, this%flux_x%x, &
               this%flux_y%x, this%flux_z%x)
       end if
       ! The periodic graph cancels the SBP boundary terms, so assemble the
       ! invariant low-order residual directly in storage that becomes U^L.
       ! The nonperiodic path retains the element-local reference residual.
       if (this%periodic_graph) then
          this%low_candidate(component)%x = 0.0_rp
       end if
       do edge = 1, this%graph%n_edges
          associate(a => this%graph%left(:, edge), &
               b => this%graph%right(:, edge), &
               coefficient => this%graph%coefficient(:, edge))
            left_value = euler_idp_cpu_state_component(component, a, rho, &
                 m_x, m_y, m_z, energy)
            right_value = euler_idp_cpu_state_component(component, b, rho, &
                 m_x, m_y, m_z, energy)
            difference = right_value - left_value
            left_flux_value = coefficient(1) * &
                 this%flux_x%x(a(1),a(2),a(3),a(4)) + coefficient(2) * &
                 this%flux_y%x(a(1),a(2),a(3),a(4)) + coefficient(3) * &
                 this%flux_z%x(a(1),a(2),a(3),a(4))
            right_flux_value = coefficient(1) * &
                 this%flux_x%x(b(1),b(2),b(3),b(4)) + coefficient(2) * &
                 this%flux_y%x(b(1),b(2),b(3),b(4)) + coefficient(3) * &
                 this%flux_z%x(b(1),b(2),b(3),b(4))
            flux_difference = right_flux_value - left_flux_value
            if (need_local_residual) then
               this%low_local_residual(component)%x( &
                    a(1),a(2),a(3),a(4)) = &
                    this%low_local_residual(component)%x( &
                    a(1),a(2),a(3),a(4)) - &
                    this%edge_viscosity(edge) * difference
               this%low_local_residual(component)%x( &
                    b(1),b(2),b(3),b(4)) = &
                    this%low_local_residual(component)%x( &
                    b(1),b(2),b(3),b(4)) + &
                    this%edge_viscosity(edge) * difference
            end if
            if (this%periodic_graph) then
               this%low_candidate(component)%x( &
                    a(1),a(2),a(3),a(4)) = &
                    this%low_candidate(component)%x( &
                    a(1),a(2),a(3),a(4)) + &
                    flux_difference - this%edge_viscosity(edge) * difference
               this%low_candidate(component)%x( &
                    b(1),b(2),b(3),b(4)) = &
                    this%low_candidate(component)%x( &
                    b(1),b(2),b(3),b(4)) + &
                    flux_difference + this%edge_viscosity(edge) * difference
            end if
          end associate
       end do
       if (this%diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL) then
          call euler_idp_cpu_element_sums(this, component, coef, &
               this%low_local_residual(component), &
               this%low_element_residual_sum, &
               this%low_element_boundary_flux)
       end if
       call profiler_start_region('Euler IDP reconstruction')
       call this%graph%directional_correction_flux( &
            this%correction_flux(component,:), this%flux_x%x, &
            this%flux_y%x, this%flux_z%x, coef, &
            directional_error(component))
       call profiler_end_region('Euler IDP reconstruction')
       if (.not. this%periodic_graph) then
          this%low_candidate(component)%x = &
               this%low_local_residual(component)%x
       end if
       call profiler_start_region('Euler IDP gather-scatter')
       call gs%op(this%low_candidate(component), GS_OP_ADD)
       call profiler_end_region('Euler IDP gather-scatter')
       this%low_candidate(component)%x = &
            this%low_candidate(component)%x * coef%Binv
    end do

    call profiler_end_region('Euler IDP low-order graph')
  end subroutine euler_idp_cpu_evaluate_low_order

  !> Return the exact maximum nodal graph CFL for a timestep.
  pure real(kind=rp) function euler_idp_cpu_graph_cfl(this, dt) result(cfl)
    class(euler_idp_cpu_t), intent(in) :: this
    real(kind=rp), intent(in) :: dt

    cfl = dt * this%max_graph_rate
  end function euler_idp_cpu_graph_cfl

  !> Build stage-local density and entropy bounds on the assembled graph.
  subroutine euler_idp_cpu_compute_bounds(this, rho, m_x, m_y, m_z, energy, &
       gs, gamma, diagnostics)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(field_t), intent(in) :: rho, m_x, m_y, m_z, energy
    type(gs_t), intent(inout) :: gs
    real(kind=rp), intent(in) :: gamma
    type(euler_idp_diagnostics_t), intent(inout) :: diagnostics
    real(kind=rp) :: left_density, right_density, bar_density
    real(kind=rp) :: left_mass_flux, right_mass_flux
    real(kind=rp) :: entropy
    real(kind=rp) :: difference, pair_average
    real(kind=rp) :: left_weight, right_weight
    real(kind=rp) :: strict_lower, strict_upper
    real(kind=rp) :: relaxed_lower, relaxed_upper
    real(kind=rp) :: local_violation, local_entropy_excess
    real(kind=rp) :: entropy_tolerance
    real(kind=rp) :: scale, tolerance
    real(kind=rp) :: state(EULER_IDP_NCOMP)
    integer :: direction, edge, i

    if (.not. this%graph%initialized) then
       call neko_error('Euler IDP graph is not initialised')
    end if

    call profiler_start_region('Euler IDP bounds')
    this%density_lower_bound%x = rho%x
    this%density_upper_bound%x = rho%x
    if (this%relax_density_bounds) then
       this%flux_x%x = 0.0_rp
    end if
    if (this%limit_entropy) then
       do i = 1, rho%size()
          state = [rho%x(i,1,1,1), m_x%x(i,1,1,1), m_y%x(i,1,1,1), &
               m_z%x(i,1,1,1), energy%x(i,1,1,1)]
          this%flux_z%x(i,1,1,1) = &
               euler_idp_specific_entropy(state, gamma)
       end do
       call euler_idp_local_entropy_bounds(this%flux_z%x, &
            this%graph%left, this%graph%right, this%entropy_lower_bound%x)
    end if
    do edge = 1, this%graph%n_edges
       associate(a => this%graph%left(:, edge), &
            b => this%graph%right(:, edge), &
            coefficient => this%graph%coefficient(:, edge))
         left_density = rho%x(a(1),a(2),a(3),a(4))
         right_density = rho%x(b(1),b(2),b(3),b(4))
         left_mass_flux = coefficient(1) * m_x%x(a(1),a(2),a(3),a(4)) + &
              coefficient(2) * m_y%x(a(1),a(2),a(3),a(4)) + &
              coefficient(3) * m_z%x(a(1),a(2),a(3),a(4))
         right_mass_flux = coefficient(1) * m_x%x(b(1),b(2),b(3),b(4)) + &
              coefficient(2) * m_y%x(b(1),b(2),b(3),b(4)) + &
              coefficient(3) * m_z%x(b(1),b(2),b(3),b(4))
         bar_density = 0.5_rp * (left_density + right_density)
         if (this%edge_viscosity(edge) .gt. tiny(1.0_rp)) then
            bar_density = bar_density - (right_mass_flux - left_mass_flux) / &
                 (2.0_rp * this%edge_viscosity(edge))
         end if
         this%density_lower_bound%x(a(1),a(2),a(3),a(4)) = min( &
              this%density_lower_bound%x(a(1),a(2),a(3),a(4)), &
              right_density, bar_density)
         this%density_lower_bound%x(b(1),b(2),b(3),b(4)) = min( &
              this%density_lower_bound%x(b(1),b(2),b(3),b(4)), &
              left_density, bar_density)
         this%density_upper_bound%x(a(1),a(2),a(3),a(4)) = max( &
              this%density_upper_bound%x(a(1),a(2),a(3),a(4)), &
              right_density, bar_density)
         this%density_upper_bound%x(b(1),b(2),b(3),b(4)) = max( &
              this%density_upper_bound%x(b(1),b(2),b(3),b(4)), &
              left_density, bar_density)
         if (this%relax_density_bounds) then
            ! Geometry-aware divided differences annihilate affine data on
            ! the nonuniform GLL spacing.  The precomputed endpoint weights
            ! also collapse duplicate element occurrences after assembly.
            left_weight = &
                 this%graph%second_difference_weight_left(edge)
            right_weight = &
                 this%graph%second_difference_weight_right(edge)
            difference = left_density - right_density
            this%flux_x%x(a(1),a(2),a(3),a(4)) = &
                 this%flux_x%x(a(1),a(2),a(3),a(4)) + &
                 left_weight * difference
            this%flux_x%x(b(1),b(2),b(3),b(4)) = &
                 this%flux_x%x(b(1),b(2),b(3),b(4)) - &
                 right_weight * difference
         end if
       end associate
    end do
    call profiler_start_region('Euler IDP gather-scatter')
    call gs%op(this%density_lower_bound, GS_OP_MIN)
    call gs%op(this%density_upper_bound, GS_OP_MAX)
    if (this%limit_entropy) then
       call gs%op(this%entropy_lower_bound, GS_OP_MIN)
    end if
    if (this%relax_density_bounds) then
       call gs%op(this%flux_x, GS_OP_ADD)
    end if
    call profiler_end_region('Euler IDP gather-scatter')

    if (this%diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL) then
       local_violation = 0.0_rp
       if (rho%size() .gt. 0) then
          scale = max(1.0_rp, maxval(abs(this%low_candidate(1)%x)), &
               maxval(abs(this%density_lower_bound%x)), &
               maxval(abs(this%density_upper_bound%x)))
          local_violation = max( &
               maxval(this%density_lower_bound%x - &
               this%low_candidate(1)%x), &
               maxval(this%low_candidate(1)%x - &
               this%density_upper_bound%x), 0.0_rp) / scale
       end if
       tolerance = 256.0_rp * epsilon(1.0_rp)
       if (local_violation .gt. tolerance) then
          call neko_error('Euler IDP low-order density is outside its ' // &
               'local bar-state bounds')
       end if

       if (this%limit_entropy) then
          local_entropy_excess = 0.0_rp
          do i = 1, rho%size()
             state = [this%low_candidate(1)%x(i,1,1,1), &
                  this%low_candidate(2)%x(i,1,1,1), &
                  this%low_candidate(3)%x(i,1,1,1), &
                  this%low_candidate(4)%x(i,1,1,1), &
                  this%low_candidate(5)%x(i,1,1,1)]
             entropy = euler_idp_specific_entropy(state, gamma)
             difference = this%entropy_lower_bound%x(i,1,1,1) - entropy
             entropy_tolerance = euler_idp_entropy_tolerance(state, &
                  this%entropy_lower_bound%x(i,1,1,1))
             local_entropy_excess = max(local_entropy_excess, &
                  difference - entropy_tolerance)
          end do
          if (local_entropy_excess .gt. 0.0_rp) then
             call neko_error('Euler IDP low-order state violates its ' // &
                  'local minimum entropy bound')
          end if
       end if
    end if

    if (this%relax_density_bounds) then
       this%flux_y%x = 0.0_rp
       do edge = 1, this%graph%n_edges
          associate(a => this%graph%left(:, edge), &
               b => this%graph%right(:, edge))
            pair_average = 0.5_rp * ( &
                 this%flux_x%x(a(1),a(2),a(3),a(4)) + &
                 this%flux_x%x(b(1),b(2),b(3),b(4)))
            ! Apply the same occurrence normalization to the neighbour average.
            direction = this%graph%direction(edge)
            left_weight = 2.0_rp / &
                 this%graph%directional_degree(direction)%x( &
                 a(1),a(2),a(3),a(4))
            right_weight = 2.0_rp / &
                 this%graph%directional_degree(direction)%x( &
                 b(1),b(2),b(3),b(4))
            this%flux_y%x( &
                 a(1),a(2),a(3),a(4)) = &
                 this%flux_y%x( &
                 a(1),a(2),a(3),a(4)) + left_weight * pair_average
            this%flux_y%x( &
                 b(1),b(2),b(3),b(4)) = &
                 this%flux_y%x( &
                 b(1),b(2),b(3),b(4)) + right_weight * pair_average
          end associate
       end do
       call profiler_start_region('Euler IDP gather-scatter')
       call gs%op(this%flux_y, GS_OP_ADD)
       call profiler_end_region('Euler IDP gather-scatter')
       this%flux_y%x = this%flux_y%x / &
            (2.0_rp * real(2 * this%graph%n_directions + 1, rp))

       do i = 1, rho%size()
          strict_lower = this%density_lower_bound%x(i,1,1,1)
          strict_upper = this%density_upper_bound%x(i,1,1,1)
          call euler_idp_relax_density_bounds(strict_lower, strict_upper, &
               this%flux_y%x(i,1,1,1), &
               this%density_bound_relaxation_factor, &
               this%density_relaxation_mass, this%domain_volume, &
               this%graph%n_directions, relaxed_lower, relaxed_upper)
          this%density_lower_bound%x(i,1,1,1) = relaxed_lower
          this%density_upper_bound%x(i,1,1,1) = relaxed_upper
       end do
    end if

    if (this%limit_entropy .and. this%relax_entropy_bounds) then
       ! The sparse graph is treated as a collection of linear subedges.
       ! Accumulate the maximum entropy of every conservative midpoint, then
       ! relax the strict one-ring lower bound in the positive variable
       ! exp(s) without ever forming that variable explicitly.
       this%flux_y%x = this%entropy_lower_bound%x
       do edge = 1, this%graph%n_edges
          associate(a => this%graph%left(:, edge), &
               b => this%graph%right(:, edge))
            state = 0.5_rp * [ &
                 rho%x(a(1),a(2),a(3),a(4)) + &
                 rho%x(b(1),b(2),b(3),b(4)), &
                 m_x%x(a(1),a(2),a(3),a(4)) + &
                 m_x%x(b(1),b(2),b(3),b(4)), &
                 m_y%x(a(1),a(2),a(3),a(4)) + &
                 m_y%x(b(1),b(2),b(3),b(4)), &
                 m_z%x(a(1),a(2),a(3),a(4)) + &
                 m_z%x(b(1),b(2),b(3),b(4)), &
                 energy%x(a(1),a(2),a(3),a(4)) + &
                 energy%x(b(1),b(2),b(3),b(4))]
            entropy = euler_idp_specific_entropy(state, gamma)
            this%flux_y%x(a(1),a(2),a(3),a(4)) = max( &
                 this%flux_y%x(a(1),a(2),a(3),a(4)), entropy)
            this%flux_y%x(b(1),b(2),b(3),b(4)) = max( &
                 this%flux_y%x(b(1),b(2),b(3),b(4)), entropy)
          end associate
       end do
       call profiler_start_region('Euler IDP gather-scatter')
       call gs%op(this%flux_y, GS_OP_MAX)
       call profiler_end_region('Euler IDP gather-scatter')
       do i = 1, rho%size()
          strict_lower = this%entropy_lower_bound%x(i,1,1,1)
          call euler_idp_relax_entropy_bound(strict_lower, &
               this%flux_y%x(i,1,1,1), &
               this%entropy_bound_relaxation_factor, &
               this%density_relaxation_mass, this%domain_volume, &
               this%graph%n_directions, relaxed_lower)
          this%entropy_lower_bound%x(i,1,1,1) = relaxed_lower
       end do
    end if
    call profiler_end_region('Euler IDP bounds')
  end subroutine euler_idp_cpu_compute_bounds

  !> Form a complete limited Forward Euler update without changing the input.
  subroutine euler_idp_cpu_forward_euler(this, rho, m_x, m_y, m_z, energy, &
       coef, gs, gamma, internal_energy_floor, dt, time, diagnostics, stage, &
       entropy_viscosity_fraction, graph_wave_speed)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(field_t), intent(in) :: rho, m_x, m_y, m_z, energy
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    real(kind=rp), intent(in) :: gamma, internal_energy_floor
    real(kind=rp), intent(in) :: dt
    type(time_state_t), intent(in) :: time
    type(euler_idp_diagnostics_t), intent(inout) :: diagnostics
    integer, intent(in) :: stage
    type(field_t), intent(in), optional :: entropy_viscosity_fraction
    type(field_t), intent(in), optional :: graph_wave_speed
    real(kind=rp) :: state(EULER_IDP_NCOMP), residual(EULER_IDP_NCOMP)
    real(kind=rp) :: local_floor_timestep, global_floor_timestep
    real(kind=rp) :: local_error, local_scale
    real(kind=rp) :: local_entropy_sum, global_entropy_sum
    real(kind=rp) :: local_entropy_max, global_entropy_max
    real(kind=rp) :: high_order_fraction, edge_fraction
    real(kind=rp) :: state_difference(EULER_IDP_NCOMP)
    real(kind=rp) :: directional_error_local(EULER_IDP_NCOMP)
    real(kind=rp) :: local_wave_speed, global_wave_speed
    character(len=2 * LOG_SIZE) :: message
    integer :: component, edge, global_node_count, i, ierr
    logical :: scalar_density_mode

    call profiler_start_region('Euler IDP Forward Euler')
    call diagnostics%reset()
    diagnostics%stage = stage
    diagnostics%stage_time = time%t
    scalar_density_mode = present(graph_wave_speed) .and. &
         .not. this%limit_internal_energy .and. .not. this%limit_entropy
    diagnostics%entropy_viscosity_enabled = &
         present(entropy_viscosity_fraction)
    if (this%diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL) then
       this%flux_y%x = m_x%x
       this%flux_z%x = m_x%x
       call gs%op(this%flux_y, GS_OP_MIN)
       call gs%op(this%flux_z, GS_OP_MAX)
       local_error = maxval(abs(this%flux_z%x - this%flux_y%x))
       local_scale = max(1.0_rp, maxval(abs(this%flux_z%x)), &
            maxval(abs(this%flux_y%x)))
       if (local_error .gt. 512.0_rp * epsilon(1.0_rp) * local_scale) then
          write(message, '(A,I0,A,I0,A,ES13.6)') &
               'Euler IDP stage input momentum-x differs across shared ' // &
               'copies on rank ', pe_rank, ', SSPRK stage ', stage, &
               ', error ', local_error
          call neko_error(trim(message))
       end if
    end if
    call this%evaluate_low_order(rho, m_x, m_y, m_z, energy, coef, gs, &
         gamma, internal_energy_floor, diagnostics, directional_error_local)

    if (this%diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL) then
       local_wave_speed = 0.0_rp
       if (rho%size() .gt. 0) then
          local_wave_speed = maxval(sqrt(this%u%x**2 + this%v%x**2 + &
               this%w%x**2) + this%sound_speed%x)
       end if
       call profiler_start_region('Euler IDP MPI reduction')
       call MPI_Allreduce(local_wave_speed, global_wave_speed, 1, &
            MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM, ierr)
       call profiler_end_region('Euler IDP MPI reduction')
       diagnostics%max_nodal_wave_speed = global_wave_speed
    end if
    diagnostics%max_graph_wave_speed = this%max_graph_wave_speed
    diagnostics%max_graph_rate = this%max_graph_rate

    call profiler_start_region('Euler IDP timestep limits')
    diagnostics%max_graph_cfl = this%graph_cfl(dt)
    diagnostics%min_convex_weight = 1.0_rp - diagnostics%max_graph_cfl
    if (diagnostics%max_graph_cfl .gt. 1.0_rp + 32.0_rp * &
         epsilon(1.0_rp)) then
       if (diagnostics%stage .gt. 0) then
          write(message, '(A,I0,A,ES13.6,A,ES13.6)') &
               'Euler IDP SSPRK stage ', diagnostics%stage, &
               ' graph CFL exceeds one: ', diagnostics%max_graph_cfl, &
               ', maximum dt: ', this%maximum_graph_timestep
       else
          write(message, '(A,ES13.6,A,ES13.6)') &
               'Euler IDP graph CFL exceeds one: ', &
               diagnostics%max_graph_cfl, ', maximum dt: ', &
               this%maximum_graph_timestep
       end if
       call neko_error(trim(message))
    end if

    if (this%limit_internal_energy) then
       local_floor_timestep = this%maximum_graph_timestep
       do i = 1, rho%size()
          state = [rho%x(i,1,1,1), m_x%x(i,1,1,1), m_y%x(i,1,1,1), &
               m_z%x(i,1,1,1), energy%x(i,1,1,1)]
          do component = 1, EULER_IDP_NCOMP
             residual(component) = &
                  this%low_candidate(component)%x(i,1,1,1)
          end do
          local_floor_timestep = min(local_floor_timestep, &
               euler_idp_internal_energy_timestep(state, residual, &
               internal_energy_floor, this%maximum_graph_timestep))
       end do
       call profiler_start_region('Euler IDP MPI reduction')
       call MPI_Allreduce(local_floor_timestep, global_floor_timestep, 1, &
            MPI_REAL_PRECISION, MPI_MIN, NEKO_COMM, ierr)
       call profiler_end_region('Euler IDP MPI reduction')
    else
       global_floor_timestep = this%maximum_graph_timestep
    end if
    this%maximum_floor_timestep = global_floor_timestep
    diagnostics%maximum_floor_timestep = global_floor_timestep
    if (this%limit_internal_energy .and. &
         dt .gt. global_floor_timestep * (1.0_rp + 32.0_rp * &
         epsilon(1.0_rp))) then
       if (diagnostics%stage .gt. 0) then
          write(message, '(A,I0,A,ES13.6,A,ES13.6)') &
               'Euler IDP SSPRK stage ', diagnostics%stage, &
               ' timestep exceeds the floor-aware limit: ', dt, &
               ', maximum dt: ', global_floor_timestep
       else
          write(message, '(A,ES13.6,A,ES13.6)') &
               'Euler IDP timestep exceeds the floor-aware limit: ', dt, &
               ', maximum dt: ', global_floor_timestep
       end if
       call neko_error(trim(message))
    end if
    call profiler_end_region('Euler IDP timestep limits')

    call profiler_start_region('Euler IDP low update')
    this%low_candidate(1)%x = rho%x - dt * this%low_candidate(1)%x
    this%low_candidate(2)%x = m_x%x - dt * this%low_candidate(2)%x
    this%low_candidate(3)%x = m_y%x - dt * this%low_candidate(3)%x
    this%low_candidate(4)%x = m_z%x - dt * this%low_candidate(4)%x
    this%low_candidate(5)%x = energy%x - dt * this%low_candidate(5)%x
    if (scalar_density_mode) then
       ! User-defined scalar problems are embedded in density. The remaining
       ! conservative fields only carry the prescribed scalar flux and are
       ! projected by the user routine after each complete step.
       this%low_candidate(2)%x = m_x%x
       this%low_candidate(3)%x = m_y%x
       this%low_candidate(4)%x = m_z%x
       this%low_candidate(5)%x = energy%x
    end if
    call profiler_end_region('Euler IDP low update')

    call this%compute_bounds(rho, m_x, m_y, m_z, energy, gs, gamma, &
         diagnostics)

    if (present(entropy_viscosity_fraction) .and. &
         this%diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL) then
       local_entropy_max = 0.0_rp
       local_entropy_sum = 0.0_rp
       if (rho%size() .gt. 0) then
          local_entropy_max = maxval(entropy_viscosity_fraction%x)
          local_entropy_sum = sum(entropy_viscosity_fraction%x)
       end if
       call profiler_start_region('Euler IDP MPI reduction')
       call MPI_Allreduce(local_entropy_max, global_entropy_max, 1, &
            MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM, ierr)
       call MPI_Allreduce(local_entropy_sum, global_entropy_sum, 1, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
       call MPI_Allreduce(rho%size(), global_node_count, 1, MPI_INTEGER, &
            MPI_SUM, NEKO_COMM, ierr)
       call profiler_end_region('Euler IDP MPI reduction')
       diagnostics%max_entropy_viscosity = global_entropy_max
       diagnostics%mean_entropy_viscosity = global_entropy_sum / &
            real(max(1, global_node_count), rp)
    end if

    if (this%diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL) then
       do component = 1, EULER_IDP_NCOMP
          local_error = 0.0_rp
          local_scale = 1.0_rp
          if (coef%msh%nelv .gt. 0) then
             local_error = maxval(abs(dt * &
                  (this%low_element_residual_sum(component,:) - &
                  this%element_residual_sum(component,:))))
             local_scale = max(1.0_rp, maxval(abs(dt * &
                  this%low_element_residual_sum(component,:))), &
                  maxval(abs(dt * this%element_residual_sum(component,:))))
          end if
          diagnostics%correction_element_compatibility(component) = &
               local_error / local_scale

          local_error = 0.0_rp
          local_scale = 1.0_rp
          if (coef%msh%nelv .gt. 0) then
             local_error = maxval(abs(dt * &
                  (this%low_element_boundary_flux(component,:) - &
                  this%element_boundary_flux(component,:))))
             local_scale = max(1.0_rp, maxval(abs(dt * &
                  this%low_element_boundary_flux(component,:))), &
                  maxval(abs(dt * this%element_boundary_flux(component,:))))
          end if
          diagnostics%correction_face_mismatch(component) = &
               local_error / local_scale
       end do
       if (maxval(diagnostics%correction_element_compatibility) .gt. &
            this%correction_tolerance .or. &
            maxval(diagnostics%correction_face_mismatch) .gt. &
            this%correction_tolerance) then
          call neko_error('Euler IDP correction is incompatible with the ' // &
               'selected element-local reconstruction')
       end if
    end if

    ! Reconstruct the conservative low-to-high inviscid correction on the
    ! coordinate-neighbour graph.
    call profiler_start_region('Euler IDP reconstruction')
    do edge = 1, this%graph%n_edges
       associate(a => this%graph%left(:, edge), &
            b => this%graph%right(:, edge))
         state_difference(1) = rho%x(a(1),a(2),a(3),a(4)) - &
              rho%x(b(1),b(2),b(3),b(4))
         state_difference(2) = m_x%x(a(1),a(2),a(3),a(4)) - &
              m_x%x(b(1),b(2),b(3),b(4))
         state_difference(3) = m_y%x(a(1),a(2),a(3),a(4)) - &
              m_y%x(b(1),b(2),b(3),b(4))
         state_difference(4) = m_z%x(a(1),a(2),a(3),a(4)) - &
              m_z%x(b(1),b(2),b(3),b(4))
         state_difference(5) = energy%x(a(1),a(2),a(3),a(4)) - &
              energy%x(b(1),b(2),b(3),b(4))
       end associate

       high_order_fraction = 1.0_rp
       if (this%low_order_only) then
          high_order_fraction = 0.0_rp
       else if (present(entropy_viscosity_fraction)) then
          associate(a => this%graph%left(:, edge), &
               b => this%graph%right(:, edge))
            edge_fraction = min(1.0_rp, max(0.0_rp, &
                 entropy_viscosity_fraction%x(a(1),a(2),a(3),a(4)), &
                 entropy_viscosity_fraction%x(b(1),b(2),b(3),b(4))))
          end associate
          high_order_fraction = 1.0_rp - edge_fraction
       end if

       ! Complete the raw low-to-high correction with the low-order graph
       ! viscosity, then apply the entropy-viscosity blend to the complete
       ! conservative edge vector.
       this%correction_flux(:, edge) = dt * high_order_fraction * &
            (this%correction_flux(:, edge) + &
            this%edge_viscosity(edge) * state_difference)
       if (scalar_density_mode) then
          this%correction_flux(2:EULER_IDP_NCOMP, edge) = 0.0_rp
       end if
    end do

    if (this%diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL) then
       if (maxval(directional_error_local) .gt. &
            10.0_rp * this%correction_tolerance) then
          write(message, '(A,5(1X,ES13.6))') &
               'Euler IDP directional correction is not line-compatible:', &
               directional_error_local
          call neko_error(trim(message))
       end if
       diagnostics%reconstruction_residual = directional_error_local
    end if
    call profiler_end_region('Euler IDP reconstruction')

    ! The low-order candidate is the base state for every edge correction.
    ! Validate it after the stage-state fluxes have been reconstructed so its
    ! primitive conversion does not need to be undone.
    if (scalar_density_mode) then
       if (any(.not. ieee_is_finite(this%low_candidate(1)%x)) .or. &
            minval(this%low_candidate(1)%x) .le. 0.0_rp) then
          write(message, '(A,ES13.6)') &
               'Scalar IDP low-order candidate has invalid density: ', &
               minval(this%low_candidate(1)%x)
          call neko_error(trim(message))
       end if
    else
       call euler_idp_cpu_primitives(this, this%low_candidate(1), &
            this%low_candidate(2), this%low_candidate(3), &
            this%low_candidate(4), this%low_candidate(5), gamma, &
            internal_energy_floor, 'low-order candidate')
    end if

    diagnostics%limiter_weight_error = this%limiter_weight_error
    if (this%diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL) then
       call this%validate_candidate(gamma, diagnostics, &
            'low-order candidate before correction')
    end if
    call this%compute_limiter(gamma, internal_energy_floor, diagnostics)
    call this%apply_correction(gs, diagnostics)
    call profiler_end_region('Euler IDP Forward Euler')
  end subroutine euler_idp_cpu_forward_euler

  !> Compute one symmetric vector limiter coefficient per edge occurrence.
  subroutine euler_idp_cpu_compute_limiter(this, gamma, &
       internal_energy_floor, diagnostics)
    class(euler_idp_cpu_t), intent(inout) :: this
    real(kind=rp), intent(in) :: gamma, internal_energy_floor
    type(euler_idp_diagnostics_t), intent(inout) :: diagnostics
    real(kind=rp) :: left_base(EULER_IDP_NCOMP)
    real(kind=rp) :: right_base(EULER_IDP_NCOMP)
    real(kind=rp) :: left_correction(EULER_IDP_NCOMP)
    real(kind=rp) :: right_correction(EULER_IDP_NCOMP)
    real(kind=rp) :: left_trial(EULER_IDP_NCOMP)
    real(kind=rp) :: right_trial(EULER_IDP_NCOMP)
    real(kind=rp) :: left_entropy_bound, right_entropy_bound
    real(kind=rp) :: edge_limit
    real(kind=rp) :: correction_scale, correction_size
    real(kind=rp) :: local_minimum, global_minimum
    real(kind=rp) :: local_maximum, global_maximum
    real(kind=rp) :: local_sum, global_sum
    integer :: local_count(5), global_count(5)
    integer :: direction, edge, component, ierr
    logical :: density_limited, energy_limited, entropy_limited
    logical :: check_base_constraints
    logical :: negligible_correction
    character(len=2 * LOG_SIZE) :: message

    call profiler_start_region('Euler IDP vector limiter')
    local_minimum = 1.0_rp
    local_maximum = 0.0_rp
    local_sum = 0.0_rp
    local_count = 0
    check_base_constraints = &
         this%diagnostics_level .ne. EULER_IDP_DIAGNOSTICS_OFF
    do edge = 1, this%graph%n_edges
       associate(a => this%graph%left(:, edge), &
            b => this%graph%right(:, edge))
         direction = this%graph%direction(edge)
         do component = 1, EULER_IDP_NCOMP
            left_base(component) = this%low_candidate(component)%x( &
                 a(1),a(2),a(3),a(4))
            right_base(component) = this%low_candidate(component)%x( &
                 b(1),b(2),b(3),b(4))
         end do
         ! This is correction_flux / (mass * occurrence_weight), with the
         ! directional occurrence weight constructed during initialization.
         left_correction = this%correction_flux(:, edge) * &
              real(this%graph%n_directions, rp) * &
              this%graph%directional_degree(direction)%x( &
              a(1),a(2),a(3),a(4)) / &
              this%graph%mass%x(a(1),a(2),a(3),a(4))
         right_correction = -this%correction_flux(:, edge) * &
              real(this%graph%n_directions, rp) * &
              this%graph%directional_degree(direction)%x( &
              b(1),b(2),b(3),b(4)) / &
              this%graph%mass%x(b(1),b(2),b(3),b(4))
         left_entropy_bound = -huge(1.0_rp)
         right_entropy_bound = -huge(1.0_rp)
         if (this%limit_entropy) then
            left_entropy_bound = this%entropy_lower_bound%x( &
                 a(1),a(2),a(3),a(4))
            right_entropy_bound = this%entropy_lower_bound%x( &
                 b(1),b(2),b(3),b(4))
         end if
         correction_scale = max(1.0_rp, maxval(abs(left_base)), &
              maxval(abs(right_base)))
         correction_size = max(maxval(abs(left_correction)), &
              maxval(abs(right_correction)))
         negligible_correction = correction_size .le. &
              512.0_rp * epsilon(1.0_rp) * correction_scale
         if (negligible_correction) then
            edge_limit = 1.0_rp
            density_limited = .false.
            energy_limited = .false.
            entropy_limited = .false.
            this%correction_flux(:, edge) = 0.0_rp
            left_correction = 0.0_rp
            right_correction = 0.0_rp
         else
            call euler_idp_limit_edge(left_base, right_base, &
                 left_correction, right_correction, &
                 this%density_lower_bound%x(a(1),a(2),a(3),a(4)), &
                 this%density_upper_bound%x(a(1),a(2),a(3),a(4)), &
                 this%density_lower_bound%x(b(1),b(2),b(3),b(4)), &
                 this%density_upper_bound%x(b(1),b(2),b(3),b(4)), &
                 left_entropy_bound, right_entropy_bound, gamma, &
                 internal_energy_floor, edge_limit, &
                 density_limited, energy_limited, entropy_limited, &
                 this%limit_internal_energy, this%limit_entropy, &
                 check_base_constraints)
         end if
         if (this%diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL .and. &
              this%limit_entropy) then
            left_trial = left_base + edge_limit * left_correction
            right_trial = right_base + edge_limit * right_correction
            if (.not. euler_idp_entropy_is_admissible(left_trial, gamma, &
                 left_entropy_bound, internal_energy_floor)) then
               write(message, '(A,I0,A,I0,A,I0)') &
                    'Euler IDP limiter returned an inadmissible left ' // &
                    'endpoint on rank ', pe_rank, ', SSPRK stage ', &
                    diagnostics%stage, ', edge ', edge
               call neko_error(trim(message))
            end if
            if (.not. euler_idp_entropy_is_admissible(right_trial, gamma, &
                 right_entropy_bound, internal_energy_floor)) then
               write(message, '(A,I0,A,I0,A,I0)') &
                    'Euler IDP limiter returned an inadmissible right ' // &
                    'endpoint on rank ', pe_rank, ', SSPRK stage ', &
                    diagnostics%stage, ', edge ', edge
               call neko_error(trim(message))
            end if
         end if
       end associate
       if (this%diagnostics_level .ne. EULER_IDP_DIAGNOSTICS_OFF) then
          if (.not. ieee_is_finite(edge_limit) .or. &
               edge_limit .lt. 0.0_rp .or. edge_limit .gt. 1.0_rp) then
             call neko_error('Euler IDP produced an invalid edge limiter')
          end if
       end if
       if (this%diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL) then
          local_minimum = min(local_minimum, edge_limit)
          local_maximum = max(local_maximum, edge_limit)
          local_sum = local_sum + edge_limit
          local_count(1) = local_count(1) + 1
          if (edge_limit .lt. &
               1.0_rp - 32.0_rp * epsilon(1.0_rp)) then
             local_count(2) = local_count(2) + 1
          end if
          if (density_limited) local_count(3) = local_count(3) + 1
          if (energy_limited) local_count(4) = local_count(4) + 1
          if (entropy_limited) local_count(5) = local_count(5) + 1
       end if
       this%correction_flux(:, edge) = edge_limit * &
            this%correction_flux(:, edge)
    end do

    if (this%diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL) then
       call profiler_start_region('Euler IDP MPI reduction')
       call MPI_Allreduce(local_minimum, global_minimum, 1, &
            MPI_REAL_PRECISION, MPI_MIN, NEKO_COMM, ierr)
       call MPI_Allreduce(local_maximum, global_maximum, 1, &
            MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM, ierr)
       call MPI_Allreduce(local_sum, global_sum, 1, &
            MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
       call MPI_Allreduce(local_count, global_count, 5, MPI_INTEGER, MPI_SUM, &
            NEKO_COMM, ierr)
       call profiler_end_region('Euler IDP MPI reduction')

       diagnostics%min_limiter = global_minimum
       if (global_count(1) .gt. 0) then
          diagnostics%mean_limiter = global_sum / real(global_count(1), rp)
          diagnostics%max_limiter = global_maximum
          diagnostics%limited_edge_fraction = &
               real(global_count(2), rp) / real(global_count(1), rp)
       else
          diagnostics%min_limiter = 1.0_rp
          diagnostics%mean_limiter = 1.0_rp
          diagnostics%max_limiter = 1.0_rp
          diagnostics%limited_edge_fraction = 0.0_rp
       end if
       diagnostics%density_limited_edges = global_count(3)
       diagnostics%internal_energy_limited_edges = global_count(4)
       diagnostics%entropy_limited_edges = global_count(5)
    end if
    call profiler_end_region('Euler IDP vector limiter')
  end subroutine euler_idp_cpu_compute_limiter

  !> Assemble the limited correction directly through graph incidence.
  subroutine euler_idp_cpu_apply_correction(this, gs, diagnostics)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(gs_t), intent(inout) :: gs
    type(euler_idp_diagnostics_t), intent(inout) :: diagnostics
    real(kind=rp) :: left_weight, right_weight
    real(kind=rp) :: local_error, local_scale
    real(kind=rp) :: local_conservation, global_conservation
    character(len=2 * LOG_SIZE) :: message
    integer :: component, direction, edge, ierr

    if (this%diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL) then
       do component = 1, EULER_IDP_NCOMP
          this%flux_y%x = this%low_candidate(component)%x
          this%flux_z%x = this%low_candidate(component)%x
          call gs%op(this%flux_y, GS_OP_MIN)
          call gs%op(this%flux_z, GS_OP_MAX)
          local_error = maxval(abs(this%flux_z%x - this%flux_y%x))
          local_scale = max(1.0_rp, maxval(abs(this%flux_z%x)), &
               maxval(abs(this%flux_y%x)))
          if (local_error .gt. 512.0_rp * epsilon(1.0_rp) * local_scale) then
             write(message, '(A,I0,A,I0,A,ES13.6)') &
                  'Euler IDP low-order base differs across shared copies ' // &
                  'on rank ', pe_rank, ', component ', component, ', error ', &
                  local_error
             call neko_error(trim(message))
          end if
       end do
       this%flux_y%x = 0.0_rp
       do component = 1, EULER_IDP_NCOMP
          this%local_residual(component)%x = 0.0_rp
       end do
       do edge = 1, this%graph%n_edges
          associate(a => this%graph%left(:, edge), &
               b => this%graph%right(:, edge))
            direction = this%graph%direction(edge)
            left_weight = 1.0_rp / (real(this%graph%n_directions, rp) * &
                 this%graph%directional_degree(direction)%x( &
                 a(1),a(2),a(3),a(4)))
            right_weight = 1.0_rp / (real(this%graph%n_directions, rp) * &
                 this%graph%directional_degree(direction)%x( &
                 b(1),b(2),b(3),b(4)))
            this%flux_y%x(a(1),a(2),a(3),a(4)) = &
                 this%flux_y%x(a(1),a(2),a(3),a(4)) + left_weight
            this%flux_y%x(b(1),b(2),b(3),b(4)) = &
                 this%flux_y%x(b(1),b(2),b(3),b(4)) + right_weight
            do component = 1, EULER_IDP_NCOMP
               this%local_residual(component)%x(a(1),a(2),a(3),a(4)) = &
                    this%local_residual(component)%x(a(1),a(2),a(3),a(4)) + &
                    left_weight * this%low_candidate(component)%x( &
                    a(1),a(2),a(3),a(4)) + &
                    this%correction_flux(component, edge) / &
                    this%graph%mass%x(a(1),a(2),a(3),a(4))
               this%local_residual(component)%x(b(1),b(2),b(3),b(4)) = &
                    this%local_residual(component)%x(b(1),b(2),b(3),b(4)) + &
                    right_weight * this%low_candidate(component)%x( &
                    b(1),b(2),b(3),b(4)) - &
                    this%correction_flux(component, edge) / &
                    this%graph%mass%x(b(1),b(2),b(3),b(4))
            end do
          end associate
       end do
       do component = 1, EULER_IDP_NCOMP
          call gs%op(this%local_residual(component), GS_OP_ADD)
       end do
       call gs%op(this%flux_y, GS_OP_ADD)
       local_error = maxval(abs(this%flux_y%x - 1.0_rp))
       if (local_error .gt. 512.0_rp * epsilon(1.0_rp)) then
          write(message, '(A,I0,A,ES13.6)') &
               'Euler IDP convex weights changed on rank ', pe_rank, &
               ', error ', local_error
          call neko_error(trim(message))
       end if
       this%flux_y%x = this%graph%mass%x
       this%flux_z%x = this%graph%mass%x
       call gs%op(this%flux_y, GS_OP_MIN)
       call gs%op(this%flux_z, GS_OP_MAX)
       local_error = maxval(abs(this%flux_z%x - this%flux_y%x))
       local_scale = max(1.0_rp, maxval(abs(this%flux_z%x)), &
            maxval(abs(this%flux_y%x)))
       if (local_error .gt. 512.0_rp * epsilon(1.0_rp) * local_scale) then
          write(message, '(A,I0,A,ES13.6)') &
               'Euler IDP assembled mass differs across shared copies on ' // &
               'rank ', pe_rank, ', error ', local_error
          call neko_error(trim(message))
       end if
    end if

    call profiler_start_region('Euler IDP correction asm')
    do component = 1, EULER_IDP_NCOMP
       call this%graph%incidence(this%flux_x%x, &
            this%correction_flux(component,:))
       if (this%diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL) then
          local_conservation = sum(this%flux_x%x)
          call MPI_Allreduce(local_conservation, global_conservation, 1, &
               MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
          diagnostics%limited_conservation(component) = &
               abs(global_conservation)
          diagnostics%correction_global_compatibility(component) = &
               abs(global_conservation)
       end if
       call profiler_start_region('Euler IDP gather-scatter')
       call gs%op(this%flux_x, GS_OP_ADD)
       call profiler_end_region('Euler IDP gather-scatter')
       this%low_candidate(component)%x = this%low_candidate(component)%x + &
            this%flux_x%x / this%graph%mass%x
    end do
    if (this%diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_FULL) then
       do component = 1, EULER_IDP_NCOMP
          local_error = maxval(abs(this%low_candidate(component)%x - &
               this%local_residual(component)%x))
          local_scale = max(1.0_rp, &
               maxval(abs(this%low_candidate(component)%x)), &
               maxval(abs(this%local_residual(component)%x)))
          if (local_error .gt. 512.0_rp * epsilon(1.0_rp) * local_scale) then
             write(message, '(A,I0,A,I0,A,ES13.6)') &
                  'Euler IDP direct correction differs from its convex ' // &
                  'endpoint form on rank ', pe_rank, ', component ', &
                  component, ', error ', local_error
             call neko_error(trim(message))
          end if
       end do
    end if
    call profiler_end_region('Euler IDP correction asm')
  end subroutine euler_idp_cpu_apply_correction

  !> Record diagnostics for the current limited candidate.
  subroutine euler_idp_cpu_observe_candidate(this, gs, gamma, &
       internal_energy_floor, time, label, primitives_valid, observation, &
       graph_wave_speed)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(gs_t), intent(inout) :: gs
    real(kind=rp), intent(in) :: gamma, internal_energy_floor
    type(time_state_t), intent(in) :: time
    character(len=*), intent(in) :: label
    logical, intent(in) :: primitives_valid
    type(euler_idp_state_observation_t), intent(out) :: observation
    type(field_t), intent(in), optional :: graph_wave_speed

    call euler_idp_cpu_observe_state(this, this%low_candidate(1), &
         this%low_candidate(2), this%low_candidate(3), &
         this%low_candidate(4), this%low_candidate(5), gs, gamma, &
         internal_energy_floor, time, label, observation, graph_wave_speed, &
         primitives_valid)
  end subroutine euler_idp_cpu_observe_candidate

  !> Apply the strong boundary map to the current limited candidate.
  subroutine euler_idp_cpu_apply_candidate_boundary(this, density_bcs, &
       velocity_bcs, pressure_bcs, gamma, internal_energy_floor, time, label)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(bc_list_t), intent(inout) :: density_bcs, velocity_bcs, pressure_bcs
    real(kind=rp), intent(in) :: gamma, internal_energy_floor
    type(time_state_t), intent(in) :: time
    character(len=*), intent(in) :: label

    call this%apply_boundary_conditions(this%low_candidate(1), &
         this%low_candidate(2), this%low_candidate(3), &
         this%low_candidate(4), this%low_candidate(5), density_bcs, &
         velocity_bcs, pressure_bcs, gamma, internal_energy_floor, time, label)
  end subroutine euler_idp_cpu_apply_candidate_boundary

  !> Validate local bounds and record global minima after the boundary map.
  subroutine euler_idp_cpu_validate_candidate(this, gamma, diagnostics, label)
    class(euler_idp_cpu_t), intent(inout) :: this
    real(kind=rp), intent(in) :: gamma
    type(euler_idp_diagnostics_t), intent(inout) :: diagnostics
    character(len=*), intent(in) :: label
    real(kind=rp) :: state(EULER_IDP_NCOMP)
    real(kind=rp) :: local_bound_violation(2)
    real(kind=rp) :: local_entropy_excess, entropy, entropy_tolerance
    real(kind=rp) :: worst_entropy, worst_entropy_bound
    real(kind=rp) :: worst_entropy_tolerance
    real(kind=rp) :: local_minimum(4), global_minimum(4)
    real(kind=rp) :: local_scale
    character(len=4 * LOG_SIZE) :: message
    integer :: i, ierr, worst_entropy_node

    if (this%diagnostics_level .eq. EULER_IDP_DIAGNOSTICS_OFF) return

    call profiler_start_region('Euler IDP diagnostics')
    local_bound_violation = 0.0_rp
    local_scale = 1.0_rp
    if (this%low_candidate(1)%size() .gt. 0) then
       local_bound_violation(1) = max(0.0_rp, maxval( &
            this%density_lower_bound%x - this%low_candidate(1)%x))
       local_bound_violation(2) = max(0.0_rp, maxval( &
            this%low_candidate(1)%x - this%density_upper_bound%x))
       local_scale = max(1.0_rp, maxval(abs(this%low_candidate(1)%x)), &
            maxval(abs(this%density_lower_bound%x)), &
            maxval(abs(this%density_upper_bound%x)))
    end if
    diagnostics%max_density_lower_violation = local_bound_violation(1)
    diagnostics%max_density_upper_violation = local_bound_violation(2)
    if (maxval(local_bound_violation) .gt. &
         512.0_rp * epsilon(1.0_rp) * local_scale) then
       call neko_error('Euler IDP limited density violates its local bounds')
    end if

    if (this%limit_entropy) then
       local_entropy_excess = 0.0_rp
       worst_entropy_node = 0
       worst_entropy = huge(1.0_rp)
       worst_entropy_bound = -huge(1.0_rp)
       worst_entropy_tolerance = 0.0_rp
       do i = 1, this%low_candidate(1)%size()
          state = [this%low_candidate(1)%x(i,1,1,1), &
               this%low_candidate(2)%x(i,1,1,1), &
               this%low_candidate(3)%x(i,1,1,1), &
               this%low_candidate(4)%x(i,1,1,1), &
               this%low_candidate(5)%x(i,1,1,1)]
          entropy = euler_idp_specific_entropy(state, gamma)
          entropy_tolerance = euler_idp_entropy_tolerance(state, &
               this%entropy_lower_bound%x(i,1,1,1))
          if (this%entropy_lower_bound%x(i,1,1,1) - entropy - &
               entropy_tolerance .gt. local_entropy_excess) then
             local_entropy_excess = &
                  this%entropy_lower_bound%x(i,1,1,1) - entropy - &
                  entropy_tolerance
             worst_entropy_node = i
             worst_entropy = entropy
             worst_entropy_bound = this%entropy_lower_bound%x(i,1,1,1)
             worst_entropy_tolerance = entropy_tolerance
          end if
       end do
       diagnostics%max_entropy_lower_violation = local_entropy_excess
       if (local_entropy_excess .gt. 0.0_rp) then
          write(message, '(A,A,A,I0,A,I0,A,I0,4(A,ES13.6))') &
               'Euler IDP ', trim(label), ' violates its local minimum ' // &
               'entropy bound on rank ', pe_rank, ', SSPRK stage ', &
               diagnostics%stage, ', local node ', worst_entropy_node, &
               ': entropy ', worst_entropy, ', bound ', worst_entropy_bound, &
               ', violation ', worst_entropy_bound - worst_entropy, &
               ', tolerance ', worst_entropy_tolerance
          call neko_error(trim(message))
       end if
    end if

    if (this%diagnostics_level .ne. EULER_IDP_DIAGNOSTICS_OFF) then
       local_minimum = huge(1.0_rp)
       if (this%low_candidate(1)%size() .gt. 0) then
          local_minimum(1:3) = [minval(this%low_candidate(1)%x), &
               minval(this%internal_energy%x), minval(this%p%x)]
          do i = 1, this%low_candidate(1)%size()
             state = [this%low_candidate(1)%x(i,1,1,1), &
                  this%low_candidate(2)%x(i,1,1,1), &
                  this%low_candidate(3)%x(i,1,1,1), &
                  this%low_candidate(4)%x(i,1,1,1), &
                  this%low_candidate(5)%x(i,1,1,1)]
             local_minimum(4) = min(local_minimum(4), &
                  euler_idp_specific_entropy(state, gamma))
          end do
       end if
       call profiler_start_region('Euler IDP MPI reduction')
       call MPI_Allreduce(local_minimum, global_minimum, 4, &
            MPI_REAL_PRECISION, MPI_MIN, NEKO_COMM, ierr)
       call profiler_end_region('Euler IDP MPI reduction')
       diagnostics%min_density = global_minimum(1)
       diagnostics%min_internal_energy = global_minimum(2)
       diagnostics%min_pressure = global_minimum(3)
       diagnostics%min_specific_entropy = global_minimum(4)
    end if
    call profiler_end_region('Euler IDP diagnostics')
  end subroutine euler_idp_cpu_validate_candidate

  !> Preserve the start-of-step state for SSPRK combinations.
  subroutine euler_idp_cpu_save_state(this, rho, m_x, m_y, m_z, energy)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(field_t), intent(in) :: rho, m_x, m_y, m_z, energy

    this%saved_state(1)%x = rho%x
    this%saved_state(2)%x = m_x%x
    this%saved_state(3)%x = m_y%x
    this%saved_state(4)%x = m_z%x
    this%saved_state(5)%x = energy%x
  end subroutine euler_idp_cpu_save_state

  !> Commit or combine a limited candidate with the preserved state.
  subroutine euler_idp_cpu_combine_stage(this, rho, m_x, m_y, m_z, energy, &
       saved_weight, candidate_weight)
    class(euler_idp_cpu_t), intent(in) :: this
    type(field_t), intent(inout) :: rho, m_x, m_y, m_z, energy
    real(kind=rp), intent(in) :: saved_weight, candidate_weight

    call profiler_start_region('Euler IDP SSP combination')
    if (saved_weight .eq. 0.0_rp) then
       rho%x = candidate_weight * this%low_candidate(1)%x
       m_x%x = candidate_weight * this%low_candidate(2)%x
       m_y%x = candidate_weight * this%low_candidate(3)%x
       m_z%x = candidate_weight * this%low_candidate(4)%x
       energy%x = candidate_weight * this%low_candidate(5)%x
    else
       rho%x = saved_weight * this%saved_state(1)%x + &
            candidate_weight * this%low_candidate(1)%x
       m_x%x = saved_weight * this%saved_state(2)%x + &
            candidate_weight * this%low_candidate(2)%x
       m_y%x = saved_weight * this%saved_state(3)%x + &
            candidate_weight * this%low_candidate(3)%x
       m_z%x = saved_weight * this%saved_state(4)%x + &
            candidate_weight * this%low_candidate(4)%x
       energy%x = saved_weight * this%saved_state(5)%x + &
            candidate_weight * this%low_candidate(5)%x
    end if
    call profiler_end_region('Euler IDP SSP combination')
  end subroutine euler_idp_cpu_combine_stage

  !> Copy cached primitive variables to scheme-owned output fields.
  subroutine euler_idp_cpu_copy_primitives(this, u, v, w, p)
    class(euler_idp_cpu_t), intent(in) :: this
    type(field_t), intent(inout) :: u, v, w, p

    u%x = this%u%x
    v%x = this%v%x
    w%x = this%w%x
    p%x = this%p%x
  end subroutine euler_idp_cpu_copy_primitives

  !> Apply the compressible strong boundary conditions to a conserved state.
  subroutine euler_idp_cpu_apply_boundary_conditions(this, rho, m_x, m_y, &
       m_z, energy, density_bcs, velocity_bcs, pressure_bcs, gamma, &
       internal_energy_floor, time, label, refresh)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(field_t), intent(inout) :: rho, m_x, m_y, m_z, energy
    type(bc_list_t), intent(inout) :: density_bcs, velocity_bcs, pressure_bcs
    real(kind=rp), intent(in) :: gamma, internal_energy_floor
    type(time_state_t), intent(in), optional :: time
    character(len=*), intent(in) :: label
    logical, intent(in), optional :: refresh
    class(bc_t), pointer :: boundary
    integer :: i, n
    logical :: refresh_

    call profiler_start_region('Euler IDP boundary')
    if (density_bcs%size() .eq. 0 .and. velocity_bcs%size() .eq. 0 .and. &
         pressure_bcs%size() .eq. 0) then
       call euler_idp_cpu_primitives(this, rho, m_x, m_y, m_z, energy, &
            gamma, internal_energy_floor, label)
       call profiler_end_region('Euler IDP boundary')
       return
    end if
    refresh_ = .false.
    if (present(refresh)) refresh_ = refresh
    if (refresh_) then
       do i = 1, density_bcs%size()
          boundary => density_bcs%get(i)
          boundary%updated = .false.
       end do
       do i = 1, velocity_bcs%size()
          boundary => velocity_bcs%get(i)
          boundary%updated = .false.
       end do
       do i = 1, pressure_bcs%size()
          boundary => pressure_bcs%get(i)
          boundary%updated = .false.
       end do
       nullify(boundary)
    end if
    n = rho%size()
    if (present(time)) then
       call density_bcs%apply(rho, time = time, strong = .true.)
    else
       call density_bcs%apply(rho, strong = .true.)
    end if

    call compressible_ops_cpu_update_uvw(this%u%x, this%v%x, this%w%x, &
         m_x%x, m_y%x, m_z%x, rho%x, n)
    if (present(time)) then
       call velocity_bcs%apply(this%u, this%v, this%w, time = time, &
            strong = .true.)
    else
       call velocity_bcs%apply(this%u, this%v, this%w, &
            strong = .true.)
    end if

    call compressible_ops_cpu_update_mxyz_p_ruvw(m_x%x, m_y%x, m_z%x, &
         this%p%x, this%internal_energy%x, this%u%x, this%v%x, this%w%x, &
         energy%x, rho%x, gamma, n)
    if (present(time)) then
       call pressure_bcs%apply(this%p, time = time, strong = .true.)
    else
       call pressure_bcs%apply(this%p, strong = .true.)
    end if
    call compressible_ops_cpu_update_e(energy%x, this%p%x, &
         this%internal_energy%x, gamma, n, internal_energy_floor)

    call euler_idp_cpu_primitives(this, rho, m_x, m_y, m_z, energy, gamma, &
         internal_energy_floor, label)
    call profiler_end_region('Euler IDP boundary')
  end subroutine euler_idp_cpu_apply_boundary_conditions

  !> Record admissibility and wave-speed data for one conserved state.
  subroutine euler_idp_cpu_observe_state(this, rho, m_x, m_y, m_z, energy, &
       gs, gamma, internal_energy_floor, time, label, observation, &
       graph_wave_speed, primitives_valid)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(field_t), intent(in) :: rho, m_x, m_y, m_z, energy
    type(gs_t), intent(inout) :: gs
    real(kind=rp), intent(in) :: gamma, internal_energy_floor
    type(time_state_t), intent(in) :: time
    character(len=*), intent(in) :: label
    type(euler_idp_state_observation_t), intent(out) :: observation
    type(field_t), intent(in), optional :: graph_wave_speed
    logical, intent(in), optional :: primitives_valid
    real(kind=rp) :: left_state(EULER_IDP_NCOMP)
    real(kind=rp) :: right_state(EULER_IDP_NCOMP)
    real(kind=rp) :: normal(3), coefficient_norm, wave_speed
    real(kind=rp) :: local_minimum(3), global_minimum(3)
    real(kind=rp) :: local_maximum(3), global_maximum(3)
    integer :: edge, ierr
    logical :: primitives_valid_

    call profiler_start_region('Euler IDP state audit')
    primitives_valid_ = .false.
    if (present(primitives_valid)) primitives_valid_ = primitives_valid
    if (.not. primitives_valid_) then
       call euler_idp_cpu_primitives(this, rho, m_x, m_y, m_z, energy, &
            gamma, internal_energy_floor, label)
    end if

    this%viscosity_sum%x = 0.0_rp
    local_maximum = 0.0_rp
    if (rho%size() .gt. 0) then
       local_minimum = [minval(rho%x), &
            minval(this%internal_energy%x), minval(this%p%x)]
       local_maximum(1) = maxval(sqrt(this%u%x**2 + this%v%x**2 + &
            this%w%x**2) + this%sound_speed%x)
    else
       local_minimum = huge(1.0_rp)
    end if

    do edge = 1, this%graph%n_edges
       associate(a => this%graph%left(:, edge), &
            b => this%graph%right(:, edge), &
            coefficient => this%graph%coefficient(:, edge))
         coefficient_norm = sqrt(dot_product(coefficient, coefficient))
         if (coefficient_norm .le. tiny(1.0_rp)) cycle
         if (present(graph_wave_speed)) then
            wave_speed = max(graph_wave_speed%x(a(1),a(2),a(3),a(4)), &
                 graph_wave_speed%x(b(1),b(2),b(3),b(4)))
            if (.not. ieee_is_finite(wave_speed) .or. &
                 wave_speed .lt. 0.0_rp) then
               call neko_error('User graph wave speed must be finite and ' // &
                    'nonnegative')
            end if
         else
            normal = coefficient / coefficient_norm
            left_state = [rho%x(a(1),a(2),a(3),a(4)), &
                 m_x%x(a(1),a(2),a(3),a(4)), &
                 m_y%x(a(1),a(2),a(3),a(4)), &
                 m_z%x(a(1),a(2),a(3),a(4)), &
                 energy%x(a(1),a(2),a(3),a(4))]
            right_state = [rho%x(b(1),b(2),b(3),b(4)), &
                 m_x%x(b(1),b(2),b(3),b(4)), &
                 m_y%x(b(1),b(2),b(3),b(4)), &
                 m_z%x(b(1),b(2),b(3),b(4)), &
                 energy%x(b(1),b(2),b(3),b(4))]
            wave_speed = euler_idp_maximum_wave_speed(left_state, &
                 right_state, normal, gamma)
         end if
         local_maximum(2) = max(local_maximum(2), wave_speed)
         this%viscosity_sum%x(a(1),a(2),a(3),a(4)) = &
              this%viscosity_sum%x(a(1),a(2),a(3),a(4)) + &
              coefficient_norm * wave_speed
         this%viscosity_sum%x(b(1),b(2),b(3),b(4)) = &
              this%viscosity_sum%x(b(1),b(2),b(3),b(4)) + &
              coefficient_norm * wave_speed
       end associate
    end do

    call profiler_start_region('Euler IDP gather-scatter')
    call gs%op(this%viscosity_sum, GS_OP_ADD)
    call profiler_end_region('Euler IDP gather-scatter')
    if (rho%size() .gt. 0) then
       local_maximum(3) = maxval(2.0_rp * this%viscosity_sum%x / &
            this%graph%mass%x)
    end if

    call profiler_start_region('Euler IDP MPI reduction')
    call MPI_Allreduce(local_minimum, global_minimum, 3, &
         MPI_REAL_PRECISION, MPI_MIN, NEKO_COMM, ierr)
    call MPI_Allreduce(local_maximum, global_maximum, 3, &
         MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM, ierr)
    call profiler_end_region('Euler IDP MPI reduction')

    observation%time = time%t
    observation%min_density = global_minimum(1)
    observation%min_internal_energy = global_minimum(2)
    observation%min_pressure = global_minimum(3)
    observation%max_nodal_wave_speed = global_maximum(1)
    observation%max_graph_wave_speed = global_maximum(2)
    observation%max_graph_rate = global_maximum(3)
    call profiler_end_region('Euler IDP state audit')
  end subroutine euler_idp_cpu_observe_state

  !> Reconstruct primitive fields from one conserved state.
  subroutine euler_idp_cpu_primitives(this, rho, m_x, m_y, m_z, energy, &
       gamma, internal_energy_floor, label)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(field_t), intent(in) :: rho, m_x, m_y, m_z, energy
    real(kind=rp), intent(in) :: gamma, internal_energy_floor
    character(len=*), intent(in) :: label
    character(len=LOG_SIZE) :: message
    integer :: first_invalid, n

    call profiler_start_region('Euler IDP primitives')
    n = rho%size()
    call compressible_ops_cpu_conserved_to_primitive(rho%x, m_x%x, m_y%x, &
         m_z%x, energy%x, gamma, internal_energy_floor, &
         this%flux_x%x, this%u%x, this%v%x, this%w%x, this%p%x, &
         this%sound_speed%x, this%internal_energy%x, this%state_status, n)
    if (any(this%state_status .ne. EULER_STATE_OK)) then
       first_invalid = minloc(this%state_status, dim = 1, &
            mask = this%state_status .ne. EULER_STATE_OK)
       write(message, '(A,A,A,I0,A,I0,A,I0)') 'Euler IDP ', trim(label), &
            ' is invalid on rank ', pe_rank, ' at local node ', &
            first_invalid, ', status ', this%state_status(first_invalid)
       call neko_error(trim(message))
    end if
    call profiler_end_region('Euler IDP primitives')
  end subroutine euler_idp_cpu_primitives

  !> Construct one Cartesian Euler flux vector.
  subroutine euler_idp_cpu_flux(this, component, m_x, m_y, m_z, energy)
    class(euler_idp_cpu_t), intent(inout) :: this
    integer, intent(in) :: component
    type(field_t), intent(in) :: m_x, m_y, m_z, energy
    integer :: i, n

    n = m_x%size()
    select case (component)
    case (1)
       this%flux_x%x = m_x%x
       this%flux_y%x = m_y%x
       this%flux_z%x = m_z%x
    case (2)
       do concurrent (i = 1:n)
          this%flux_x%x(i,1,1,1) = m_x%x(i,1,1,1) * &
               this%u%x(i,1,1,1) + this%p%x(i,1,1,1)
          this%flux_y%x(i,1,1,1) = m_x%x(i,1,1,1) * &
               this%v%x(i,1,1,1)
          this%flux_z%x(i,1,1,1) = m_x%x(i,1,1,1) * &
               this%w%x(i,1,1,1)
       end do
    case (3)
       do concurrent (i = 1:n)
          this%flux_x%x(i,1,1,1) = m_y%x(i,1,1,1) * &
               this%u%x(i,1,1,1)
          this%flux_y%x(i,1,1,1) = m_y%x(i,1,1,1) * &
               this%v%x(i,1,1,1) + this%p%x(i,1,1,1)
          this%flux_z%x(i,1,1,1) = m_y%x(i,1,1,1) * &
               this%w%x(i,1,1,1)
       end do
    case (4)
       do concurrent (i = 1:n)
          this%flux_x%x(i,1,1,1) = m_z%x(i,1,1,1) * &
               this%u%x(i,1,1,1)
          this%flux_y%x(i,1,1,1) = m_z%x(i,1,1,1) * &
               this%v%x(i,1,1,1)
          this%flux_z%x(i,1,1,1) = m_z%x(i,1,1,1) * &
               this%w%x(i,1,1,1) + this%p%x(i,1,1,1)
       end do
    case (5)
       do concurrent (i = 1:n)
          this%flux_x%x(i,1,1,1) = &
               (energy%x(i,1,1,1) + this%p%x(i,1,1,1)) * &
               this%u%x(i,1,1,1)
          this%flux_y%x(i,1,1,1) = &
               (energy%x(i,1,1,1) + this%p%x(i,1,1,1)) * &
               this%v%x(i,1,1,1)
          this%flux_z%x(i,1,1,1) = &
               (energy%x(i,1,1,1) + this%p%x(i,1,1,1)) * &
               this%w%x(i,1,1,1)
       end do
    end select
  end subroutine euler_idp_cpu_flux

  !> Read one conserved component at a graph endpoint.
  real(kind=rp) function euler_idp_cpu_state_component(component, index, &
       rho, m_x, m_y, m_z, energy) result(value)
    integer, intent(in) :: component, index(4)
    type(field_t), intent(in) :: rho, m_x, m_y, m_z, energy

    select case (component)
    case (1)
       value = rho%x(index(1), index(2), index(3), index(4))
    case (2)
       value = m_x%x(index(1), index(2), index(3), index(4))
    case (3)
       value = m_y%x(index(1), index(2), index(3), index(4))
    case (4)
       value = m_z%x(index(1), index(2), index(3), index(4))
    case (5)
       value = energy%x(index(1), index(2), index(3), index(4))
    end select
  end function euler_idp_cpu_state_component

  !> Compute element volume and independent surface residual sums.
  subroutine euler_idp_cpu_element_sums(this, component, coef, &
       local_residual, element_residual_sum, element_boundary_flux)
    class(euler_idp_cpu_t), intent(inout) :: this
    integer, intent(in) :: component
    type(coef_t), intent(in) :: coef
    type(field_t), intent(in) :: local_residual
    real(kind=rp), intent(inout) :: element_residual_sum(:,:)
    real(kind=rp), intent(inout) :: element_boundary_flux(:,:)
    real(kind=rp) :: normal(3), area, boundary_flux
    integer :: e, i, j, k, lx, ly, lz

    lx = coef%Xh%lx
    ly = coef%Xh%ly
    lz = coef%Xh%lz
    do e = 1, coef%msh%nelv
       element_residual_sum(component,e) = &
            sum(local_residual%x(:,:,:,e))
       boundary_flux = 0.0_rp

       do k = 1, lz
          do j = 1, ly
             normal = coef%get_normal(1, j, k, e, 1)
             area = coef%get_area(1, j, k, e, 1)
             boundary_flux = boundary_flux + area * &
                  (this%flux_x%x(1,j,k,e) * normal(1) + &
                  this%flux_y%x(1,j,k,e) * normal(2) + &
                  this%flux_z%x(1,j,k,e) * normal(3))
             normal = coef%get_normal(lx, j, k, e, 2)
             area = coef%get_area(lx, j, k, e, 2)
             boundary_flux = boundary_flux + area * &
                  (this%flux_x%x(lx,j,k,e) * normal(1) + &
                  this%flux_y%x(lx,j,k,e) * normal(2) + &
                  this%flux_z%x(lx,j,k,e) * normal(3))
          end do
       end do
       do k = 1, lz
          do i = 1, lx
             normal = coef%get_normal(i, 1, k, e, 3)
             area = coef%get_area(i, 1, k, e, 3)
             boundary_flux = boundary_flux + area * &
                  (this%flux_x%x(i,1,k,e) * normal(1) + &
                  this%flux_y%x(i,1,k,e) * normal(2) + &
                  this%flux_z%x(i,1,k,e) * normal(3))
             normal = coef%get_normal(i, ly, k, e, 4)
             area = coef%get_area(i, ly, k, e, 4)
             boundary_flux = boundary_flux + area * &
                  (this%flux_x%x(i, ly,k,e) * normal(1) + &
                  this%flux_y%x(i, ly,k,e) * normal(2) + &
                  this%flux_z%x(i, ly,k,e) * normal(3))
          end do
       end do
       do j = 1, ly
          do i = 1, lx
             normal = coef%get_normal(i, j, 1, e, 5)
             area = coef%get_area(i, j, 1, e, 5)
             boundary_flux = boundary_flux + area * &
                  (this%flux_x%x(i,j,1,e) * normal(1) + &
                  this%flux_y%x(i,j,1,e) * normal(2) + &
                  this%flux_z%x(i,j,1,e) * normal(3))
             normal = coef%get_normal(i, j, lz, e, 6)
             area = coef%get_area(i, j, lz, e, 6)
             boundary_flux = boundary_flux + area * &
                  (this%flux_x%x(i,j, lz,e) * normal(1) + &
                  this%flux_y%x(i,j, lz,e) * normal(2) + &
                  this%flux_z%x(i,j, lz,e) * normal(3))
          end do
       end do
       element_boundary_flux(component,e) = boundary_flux
    end do
  end subroutine euler_idp_cpu_element_sums

end module euler_idp_cpu
