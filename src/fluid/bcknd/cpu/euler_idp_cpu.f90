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
!> CPU high-order candidate for the Euler IDP solver.
module euler_idp_cpu
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mpi_f08, only : MPI_Allreduce, MPI_INTEGER, MPI_MAX, MPI_MIN, MPI_SUM
  use num_types, only : rp
  use field, only : field_t
  use dofmap, only : dofmap_t
  use coefs, only : coef_t
  use gather_scatter, only : gs_t
  use bc, only : bc_t
  use bc_list, only : bc_list_t
  use gs_ops, only : GS_OP_ADD, GS_OP_MIN, GS_OP_MAX
  use operators, only : div
  use compressible_ops_cpu, only : &
       compressible_ops_cpu_conserved_to_primitive, EULER_STATE_OK
  use euler_idp, only : EULER_IDP_NCOMP, euler_idp_diagnostics_t
  use euler_idp_low_order, only : euler_idp_maximum_wave_speed, &
       euler_idp_flux_dot_vector, euler_idp_bar_state, &
       euler_idp_internal_energy, euler_idp_internal_energy_timestep
  use euler_idp_limiter, only : euler_idp_limit_edge, &
       euler_idp_local_entropy_bounds, euler_idp_relax_density_bounds, &
       euler_idp_specific_entropy
  use euler_gll_graph, only : euler_gll_graph_t
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION, pe_rank
  use utils, only : neko_error
  use logger, only : LOG_SIZE
  implicit none
  private

  !> Persistent work fields for the CPU high-order Euler candidate.
  type, public :: euler_idp_cpu_t
     logical :: initialized = .false.
     logical :: periodic_graph = .false.
     logical :: wall_mask_initialized = .false.
     logical :: low_order_only = .false.
     logical :: relax_density_bounds = .false.
     real(kind=rp) :: correction_tolerance = 1.0e-12_rp
     type(euler_gll_graph_t) :: graph
     type(field_t) :: local_residual(EULER_IDP_NCOMP)
     type(field_t) :: candidate(EULER_IDP_NCOMP)
     type(field_t) :: low_local_residual(EULER_IDP_NCOMP)
     type(field_t) :: low_assembled_residual(EULER_IDP_NCOMP)
     type(field_t) :: low_candidate(EULER_IDP_NCOMP)
     type(field_t) :: limited_candidate(EULER_IDP_NCOMP)
     type(field_t) :: saved_state(EULER_IDP_NCOMP)
     type(field_t) :: viscosity_sum
     type(field_t) :: limiter_weight_sum
     type(field_t) :: wall_mask
     type(field_t) :: density_lower_bound
     type(field_t) :: density_upper_bound
     type(field_t) :: stage_entropy
     type(field_t) :: entropy_lower_bound
     type(field_t) :: density_second_difference
     type(field_t) :: density_second_difference_average
     type(field_t) :: rho
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
     real(kind=rp), allocatable :: edge_entropy_diffusion(:)
     real(kind=rp), allocatable :: edge_wave_speed(:)
     real(kind=rp), allocatable :: bar_state(:,:)
     real(kind=rp), allocatable :: correction_flux(:,:)
     real(kind=rp), allocatable :: edge_limiter(:)
     integer, allocatable :: state_status(:)
     real(kind=rp) :: max_graph_rate = 0.0_rp
     real(kind=rp) :: domain_volume = 0.0_rp
     real(kind=rp) :: density_relaxation_mass = 0.0_rp
     real(kind=rp) :: maximum_graph_timestep = huge(1.0_rp)
     real(kind=rp) :: maximum_floor_timestep = huge(1.0_rp)
     real(kind=rp) :: limiter_weight_error = 0.0_rp
     type(euler_idp_diagnostics_t) :: stage_diagnostics(3)
   contains
     procedure, pass(this) :: init => euler_idp_cpu_init
     procedure, pass(this) :: init_graph => euler_idp_cpu_init_graph
     procedure, pass(this) :: free => euler_idp_cpu_free
     procedure, pass(this) :: evaluate_high_order => &
          euler_idp_cpu_evaluate_high_order
     procedure, pass(this) :: update_graph_viscosity => &
          euler_idp_cpu_update_graph_viscosity
     procedure, pass(this) :: evaluate_low_order => &
          euler_idp_cpu_evaluate_low_order
     procedure, pass(this) :: graph_cfl => euler_idp_cpu_graph_cfl
     procedure, pass(this) :: compute_bounds => euler_idp_cpu_compute_bounds
     procedure, pass(this) :: forward_euler => euler_idp_cpu_forward_euler
     procedure, pass(this) :: advance => euler_idp_cpu_advance
     procedure, pass(this) :: compute_limiter => &
          euler_idp_cpu_compute_limiter
     procedure, pass(this) :: apply_correction => &
          euler_idp_cpu_apply_correction
     procedure, pass(this) :: commit => euler_idp_cpu_commit
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
       write(name, '(A,I0)') 'euler_idp_local_', component
       call this%local_residual(component)%init(dof, trim(name))
       write(name, '(A,I0)') 'euler_idp_candidate_', component
       call this%candidate(component)%init(dof, trim(name))
       write(name, '(A,I0)') 'euler_idp_low_local_', component
       call this%low_local_residual(component)%init(dof, trim(name))
       write(name, '(A,I0)') 'euler_idp_low_assembled_', component
       call this%low_assembled_residual(component)%init(dof, trim(name))
       write(name, '(A,I0)') 'euler_idp_low_candidate_', component
       call this%low_candidate(component)%init(dof, trim(name))
       write(name, '(A,I0)') 'euler_idp_limited_candidate_', component
       call this%limited_candidate(component)%init(dof, trim(name))
       write(name, '(A,I0)') 'euler_idp_saved_state_', component
       call this%saved_state(component)%init(dof, trim(name))
    end do

    call this%rho%init(dof, 'euler_idp_rho')
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
    call this%limiter_weight_sum%init(dof, 'euler_idp_limiter_weight_sum')
    call this%wall_mask%init(dof, 'euler_idp_wall_mask')
    this%wall_mask%x = 0.0_rp
    call this%density_lower_bound%init(dof, 'euler_idp_density_lower_bound')
    call this%density_upper_bound%init(dof, 'euler_idp_density_upper_bound')
    call this%stage_entropy%init(dof, 'euler_idp_stage_entropy')
    call this%entropy_lower_bound%init(dof, 'euler_idp_entropy_lower_bound')
    call this%density_second_difference%init(dof, &
         'euler_idp_density_second_difference')
    call this%density_second_difference_average%init(dof, &
         'euler_idp_density_second_difference_average')

    allocate(this%element_residual_sum(EULER_IDP_NCOMP, dof%msh%nelv))
    allocate(this%element_boundary_flux(EULER_IDP_NCOMP, dof%msh%nelv))
    allocate(this%low_element_residual_sum(EULER_IDP_NCOMP, dof%msh%nelv))
    allocate(this%low_element_boundary_flux(EULER_IDP_NCOMP, dof%msh%nelv))
    allocate(this%state_status(dof%size()))
    this%element_residual_sum = 0.0_rp
    this%element_boundary_flux = 0.0_rp
    this%low_element_residual_sum = 0.0_rp
    this%low_element_boundary_flux = 0.0_rp
    this%state_status = EULER_STATE_OK
    this%max_graph_rate = 0.0_rp
    this%domain_volume = 0.0_rp
    this%maximum_graph_timestep = huge(1.0_rp)
    this%maximum_floor_timestep = huge(1.0_rp)
    this%limiter_weight_error = 0.0_rp
    this%wall_mask_initialized = .false.
    this%low_order_only = .false.
    do component = 1, 3
       call this%stage_diagnostics(component)%reset()
    end do
    this%initialized = .true.
  end subroutine euler_idp_cpu_init

  !> Initialise the sparse GLL graph after coefficients are available.
  subroutine euler_idp_cpu_init_graph(this, coef, gs, relax_density_bounds, &
       low_order_only)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    logical, intent(in), optional :: relax_density_bounds, low_order_only
    real(kind=rp) :: local_error, global_error
    real(kind=rp) :: local_mass, global_mass
    integer :: direction, edge, ierr

    if (.not. this%initialized) then
       call neko_error('Euler IDP CPU object is not initialised')
    end if
    this%correction_tolerance = 1.0e-12_rp
    this%relax_density_bounds = .false.
    this%low_order_only = .false.
    if (present(relax_density_bounds)) then
       this%relax_density_bounds = relax_density_bounds
    end if
    if (present(low_order_only)) this%low_order_only = low_order_only
    if (allocated(this%edge_viscosity)) deallocate(this%edge_viscosity)
    if (allocated(this%edge_entropy_diffusion)) then
       deallocate(this%edge_entropy_diffusion)
    end if
    if (allocated(this%edge_wave_speed)) deallocate(this%edge_wave_speed)
    if (allocated(this%bar_state)) deallocate(this%bar_state)
    if (allocated(this%correction_flux)) deallocate(this%correction_flux)
    if (allocated(this%edge_limiter)) deallocate(this%edge_limiter)
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
    allocate(this%edge_viscosity(this%graph%n_edges))
    allocate(this%edge_entropy_diffusion(this%graph%n_edges))
    allocate(this%edge_wave_speed(this%graph%n_edges))
    allocate(this%bar_state(EULER_IDP_NCOMP, this%graph%n_edges))
    allocate(this%correction_flux(EULER_IDP_NCOMP, this%graph%n_edges))
    allocate(this%edge_limiter(this%graph%n_edges))
    this%edge_viscosity = 0.0_rp
    this%edge_entropy_diffusion = 0.0_rp
    this%edge_wave_speed = 0.0_rp
    this%bar_state = 0.0_rp
    this%correction_flux = 0.0_rp
    this%edge_limiter = 1.0_rp

    local_error = 0.0_rp
    if (this%graph%mass%size() .gt. 0) then
       local_error = maxval(abs(this%graph%mass%x * coef%Binv - 1.0_rp))
    end if
    call MPI_Allreduce(local_error, global_error, 1, MPI_REAL_PRECISION, &
         MPI_MAX, NEKO_COMM, ierr)
    if (global_error .gt. 128.0_rp * epsilon(1.0_rp)) then
       call neko_error('Euler IDP assembled mass and inverse are inconsistent')
    end if

    this%limiter_weight_sum%x = 0.0_rp
    do edge = 1, this%graph%n_edges
       associate(a => this%graph%left(:,edge), &
            b => this%graph%right(:,edge))
         direction = this%graph%direction(edge)
         ! Give every tensor direction the same total nodal weight. This
         ! preserves one-dimensional states across replicated element nodes.
         this%limiter_weight_sum%x(a(1),a(2),a(3),a(4)) = &
              this%limiter_weight_sum%x(a(1),a(2),a(3),a(4)) + &
              1.0_rp / (real(this%graph%n_directions, rp) * &
              this%graph%directional_degree(direction)%x( &
              a(1),a(2),a(3),a(4)))
         this%limiter_weight_sum%x(b(1),b(2),b(3),b(4)) = &
              this%limiter_weight_sum%x(b(1),b(2),b(3),b(4)) + &
              1.0_rp / (real(this%graph%n_directions, rp) * &
              this%graph%directional_degree(direction)%x( &
              b(1),b(2),b(3),b(4)))
       end associate
    end do
    call gs%op(this%limiter_weight_sum, GS_OP_ADD)
    local_error = 0.0_rp
    if (this%limiter_weight_sum%size() .gt. 0) then
       local_error = maxval(abs(this%limiter_weight_sum%x - 1.0_rp))
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

    if (.not. this%initialized) return

    call this%graph%free()
    do component = 1, EULER_IDP_NCOMP
       call this%local_residual(component)%free()
       call this%candidate(component)%free()
       call this%low_local_residual(component)%free()
       call this%low_assembled_residual(component)%free()
       call this%low_candidate(component)%free()
       call this%limited_candidate(component)%free()
       call this%saved_state(component)%free()
    end do
    call this%rho%free()
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
    call this%limiter_weight_sum%free()
    call this%wall_mask%free()
    call this%density_lower_bound%free()
    call this%density_upper_bound%free()
    call this%stage_entropy%free()
    call this%entropy_lower_bound%free()
    call this%density_second_difference%free()
    call this%density_second_difference_average%free()
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
    if (allocated(this%edge_entropy_diffusion)) then
       deallocate(this%edge_entropy_diffusion)
    end if
    if (allocated(this%edge_wave_speed)) deallocate(this%edge_wave_speed)
    if (allocated(this%bar_state)) deallocate(this%bar_state)
    if (allocated(this%correction_flux)) deallocate(this%correction_flux)
    if (allocated(this%edge_limiter)) deallocate(this%edge_limiter)
    if (allocated(this%state_status)) deallocate(this%state_status)
    this%initialized = .false.
    this%periodic_graph = .false.
    this%wall_mask_initialized = .false.
    this%low_order_only = .false.
    this%relax_density_bounds = .false.
    this%correction_tolerance = 1.0e-12_rp
    this%max_graph_rate = 0.0_rp
    this%domain_volume = 0.0_rp
    this%density_relaxation_mass = 0.0_rp
    this%maximum_graph_timestep = huge(1.0_rp)
    this%maximum_floor_timestep = huge(1.0_rp)
    this%limiter_weight_error = 0.0_rp
  end subroutine euler_idp_cpu_free

  !> Evaluate the mass-weighted high-order residual.
  subroutine euler_idp_cpu_evaluate_high_order(this, rho, m_x, m_y, m_z, &
       energy, coef, gs, gamma, internal_energy_floor, &
       diagnostics, entropy_viscosity_fraction)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(field_t), intent(in) :: rho, m_x, m_y, m_z, energy
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    real(kind=rp), intent(in) :: gamma, internal_energy_floor
    type(euler_idp_diagnostics_t), intent(inout) :: diagnostics
    type(field_t), intent(in), optional :: entropy_viscosity_fraction
    real(kind=rp) :: edge_diffusion, edge_fraction
    integer :: component, edge

    if (.not. this%initialized) then
       call neko_error('Euler IDP CPU object is not initialised')
    end if

    call euler_idp_cpu_primitives(this, rho, m_x, m_y, m_z, energy, &
         gamma, internal_energy_floor, 'high-order state')

    this%edge_entropy_diffusion = 0.0_rp
    if (present(entropy_viscosity_fraction)) then
       diagnostics%entropy_viscosity_enabled = .true.
       call this%update_graph_viscosity(rho, m_x, m_y, m_z, energy, gs, &
            gamma)
       do edge = 1, this%graph%n_edges
          associate(a => this%graph%left(:,edge), &
               b => this%graph%right(:,edge))
               ! Entropy viscosity supplies a dimensionless fraction of the
               ! low-order graph viscosity. Its low-order cap therefore maps
               ! exactly to d_ij on every nonuniform GLL subcell.
               edge_fraction = min(1.0_rp, max(0.0_rp, &
                    entropy_viscosity_fraction%x( &
                    a(1),a(2),a(3),a(4)), &
                    entropy_viscosity_fraction%x( &
                    b(1),b(2),b(3),b(4))))
               edge_diffusion = edge_fraction * this%edge_viscosity(edge)
               ! Keep the first subcell layer at a slip wall on the robust
               ! low-order endpoint; no wall-adjacent antidiffusion is added.
               if (this%wall_mask%x(a(1),a(2),a(3),a(4)) .gt. 0.0_rp .or. &
                    this%wall_mask%x(b(1),b(2),b(3),b(4)) .gt. 0.0_rp) then
                  edge_diffusion = this%edge_viscosity(edge)
               end if
               this%edge_entropy_diffusion(edge) = edge_diffusion
          end associate
       end do
    end if

    do component = 1, EULER_IDP_NCOMP
       call euler_idp_cpu_flux(this, component, m_x, m_y, m_z, energy)
       call div(this%local_residual(component)%x, this%flux_x%x, &
            this%flux_y%x, this%flux_z%x, coef)
       this%local_residual(component)%x = &
            this%local_residual(component)%x * coef%B
       call euler_idp_cpu_element_sums(this, component, coef, &
            this%local_residual(component), this%element_residual_sum, &
            this%element_boundary_flux)

    end do

  end subroutine euler_idp_cpu_evaluate_high_order

  !> Update symmetric graph viscosity, bar states, and the graph CFL rate.
  subroutine euler_idp_cpu_update_graph_viscosity(this, rho, m_x, m_y, m_z, &
       energy, gs, gamma)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(field_t), intent(in) :: rho, m_x, m_y, m_z, energy
    type(gs_t), intent(inout) :: gs
    real(kind=rp), intent(in) :: gamma
    real(kind=rp) :: left_state(EULER_IDP_NCOMP)
    real(kind=rp) :: right_state(EULER_IDP_NCOMP)
    real(kind=rp) :: left_flux(EULER_IDP_NCOMP)
    real(kind=rp) :: right_flux(EULER_IDP_NCOMP)
    real(kind=rp) :: normal(3), coefficient_norm
    real(kind=rp) :: local_rate, global_rate
    integer :: edge, ierr

    if (.not. this%graph%initialized) then
       call neko_error('Euler IDP graph is not initialised')
    end if

    this%viscosity_sum%x = 0.0_rp
    do edge = 1, this%graph%n_edges
       associate(a => this%graph%left(:,edge), &
            b => this%graph%right(:,edge), &
            coefficient => this%graph%coefficient(:,edge))
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
         coefficient_norm = sqrt(dot_product(coefficient, coefficient))
         normal = coefficient / coefficient_norm
         this%edge_wave_speed(edge) = euler_idp_maximum_wave_speed( &
              left_state, right_state, normal, gamma)
         this%edge_viscosity(edge) = coefficient_norm * &
              this%edge_wave_speed(edge)
         call euler_idp_flux_dot_vector(left_state, coefficient, gamma, &
              left_flux)
         call euler_idp_flux_dot_vector(right_state, coefficient, gamma, &
              right_flux)
         call euler_idp_bar_state(left_state, right_state, &
              right_flux - left_flux, this%edge_viscosity(edge), &
              this%bar_state(:,edge))
         if (this%bar_state(1,edge) .le. 0.0_rp .or. &
              euler_idp_internal_energy(this%bar_state(:,edge)) .le. &
              0.0_rp) then
            call neko_error('Euler IDP graph viscosity produced an ' // &
                 'inadmissible bar state')
         end if
         this%viscosity_sum%x(a(1),a(2),a(3),a(4)) = &
              this%viscosity_sum%x(a(1),a(2),a(3),a(4)) + &
              this%edge_viscosity(edge)
         this%viscosity_sum%x(b(1),b(2),b(3),b(4)) = &
              this%viscosity_sum%x(b(1),b(2),b(3),b(4)) + &
              this%edge_viscosity(edge)
       end associate
    end do
    call gs%op(this%viscosity_sum, GS_OP_ADD)

    local_rate = 0.0_rp
    if (rho%size() .gt. 0) then
       local_rate = maxval(2.0_rp * this%viscosity_sum%x / &
            this%graph%mass%x)
    end if
    call MPI_Allreduce(local_rate, global_rate, 1, MPI_REAL_PRECISION, &
         MPI_MAX, NEKO_COMM, ierr)
    this%max_graph_rate = global_rate
    if (global_rate .gt. 0.0_rp) then
       this%maximum_graph_timestep = 1.0_rp / global_rate
    else
       this%maximum_graph_timestep = huge(1.0_rp)
    end if
  end subroutine euler_idp_cpu_update_graph_viscosity

  !> Evaluate the conservative low-order graph residual.
  subroutine euler_idp_cpu_evaluate_low_order(this, rho, m_x, m_y, m_z, &
       energy, coef, gs, gamma, internal_energy_floor, &
       diagnostics, graph_viscosity_current)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(field_t), intent(in) :: rho, m_x, m_y, m_z, energy
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    real(kind=rp), intent(in) :: gamma, internal_energy_floor
    type(euler_idp_diagnostics_t), intent(inout) :: diagnostics
    logical, intent(in), optional :: graph_viscosity_current
    real(kind=rp) :: difference, left_value, right_value
    integer :: component, edge

    call euler_idp_cpu_primitives(this, rho, m_x, m_y, m_z, energy, &
         gamma, internal_energy_floor, 'low-order state')
    if (.not. present(graph_viscosity_current)) then
       call this%update_graph_viscosity(rho, m_x, m_y, m_z, energy, gs, &
            gamma)
    else if (.not. graph_viscosity_current) then
       call this%update_graph_viscosity(rho, m_x, m_y, m_z, energy, gs, &
            gamma)
    end if

    do component = 1, EULER_IDP_NCOMP
       call euler_idp_cpu_flux(this, component, m_x, m_y, m_z, energy)
       call this%graph%local_flux_divergence( &
            this%low_local_residual(component)%x, this%flux_x%x, &
            this%flux_y%x, this%flux_z%x)
       ! The supported periodic assembly cancels the SBP boundary terms, so
       ! the assembled residual can be formed directly as
       ! 2*d_ij*(U_i-bar_U_ij), exposing the convex IDP update without an
       ! additional gather-scatter. Keep the physical element residual
       ! separately for coarse-face compatibility diagnostics.
       if (this%periodic_graph) then
          this%low_assembled_residual(component)%x = 0.0_rp
       end if
       do edge = 1, this%graph%n_edges
          associate(a => this%graph%left(:,edge), &
               b => this%graph%right(:,edge))
            left_value = euler_idp_cpu_state_component(component, a, rho, &
                 m_x, m_y, m_z, energy)
            right_value = euler_idp_cpu_state_component(component, b, rho, &
                 m_x, m_y, m_z, energy)
            difference = right_value - left_value
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
            if (this%periodic_graph) then
               this%low_assembled_residual(component)%x( &
                    a(1),a(2),a(3),a(4)) = &
                    this%low_assembled_residual(component)%x( &
                    a(1),a(2),a(3),a(4)) + &
                    2.0_rp * this%edge_viscosity(edge) * &
                    (left_value - this%bar_state(component,edge))
               this%low_assembled_residual(component)%x( &
                    b(1),b(2),b(3),b(4)) = &
                    this%low_assembled_residual(component)%x( &
                    b(1),b(2),b(3),b(4)) + &
                    2.0_rp * this%edge_viscosity(edge) * &
                    (right_value - this%bar_state(component,edge))
            end if
          end associate
       end do
       call euler_idp_cpu_element_sums(this, component, coef, &
            this%low_local_residual(component), &
            this%low_element_residual_sum, &
            this%low_element_boundary_flux)
       if (.not. this%periodic_graph) then
          this%low_assembled_residual(component)%x = &
               this%low_local_residual(component)%x
       end if
       call gs%op(this%low_assembled_residual(component), GS_OP_ADD)
       this%low_assembled_residual(component)%x = &
            this%low_assembled_residual(component)%x * coef%Binv
    end do

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
    real(kind=rp) :: entropy
    real(kind=rp) :: difference, pair_average
    real(kind=rp) :: left_weight, right_weight
    real(kind=rp) :: strict_lower, strict_upper
    real(kind=rp) :: relaxed_lower, relaxed_upper
    real(kind=rp) :: local_violation, local_entropy_violation
    real(kind=rp) :: scale, tolerance
    real(kind=rp) :: state(EULER_IDP_NCOMP)
    integer :: direction, edge, i

    if (.not. this%graph%initialized) then
       call neko_error('Euler IDP graph is not initialised')
    end if

    this%density_lower_bound%x = rho%x
    this%density_upper_bound%x = rho%x
    this%density_second_difference%x = 0.0_rp
    do i = 1, rho%size()
       state = [rho%x(i,1,1,1), m_x%x(i,1,1,1), m_y%x(i,1,1,1), &
            m_z%x(i,1,1,1), energy%x(i,1,1,1)]
       this%stage_entropy%x(i,1,1,1) = &
            euler_idp_specific_entropy(state, gamma)
    end do
    call euler_idp_local_entropy_bounds(this%stage_entropy%x, &
         this%graph%left, this%graph%right, this%entropy_lower_bound%x)
    do edge = 1, this%graph%n_edges
       associate(a => this%graph%left(:,edge), &
            b => this%graph%right(:,edge))
         left_density = rho%x(a(1),a(2),a(3),a(4))
         right_density = rho%x(b(1),b(2),b(3),b(4))
         bar_density = this%bar_state(1,edge)
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
         ! Collapse duplicate element occurrences to the unique-neighbour
         ! stencil used by the density-bound relaxation.
         direction = this%graph%direction(edge)
         left_weight = 2.0_rp / &
              this%graph%directional_degree(direction)%x( &
              a(1),a(2),a(3),a(4))
         right_weight = 2.0_rp / &
              this%graph%directional_degree(direction)%x( &
              b(1),b(2),b(3),b(4))
         difference = left_density - right_density
         this%density_second_difference%x(a(1),a(2),a(3),a(4)) = &
              this%density_second_difference%x(a(1),a(2),a(3),a(4)) + &
              left_weight * difference
         this%density_second_difference%x(b(1),b(2),b(3),b(4)) = &
              this%density_second_difference%x(b(1),b(2),b(3),b(4)) - &
              right_weight * difference
       end associate
    end do
    call gs%op(this%density_lower_bound, GS_OP_MIN)
    call gs%op(this%density_upper_bound, GS_OP_MAX)
    call gs%op(this%entropy_lower_bound, GS_OP_MIN)
    call gs%op(this%density_second_difference, GS_OP_ADD)

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
       call neko_error('Euler IDP low-order density is outside its local ' // &
            'bar-state bounds')
    end if

    local_entropy_violation = 0.0_rp
    do i = 1, rho%size()
       state = [this%low_candidate(1)%x(i,1,1,1), &
            this%low_candidate(2)%x(i,1,1,1), &
            this%low_candidate(3)%x(i,1,1,1), &
            this%low_candidate(4)%x(i,1,1,1), &
            this%low_candidate(5)%x(i,1,1,1)]
       entropy = euler_idp_specific_entropy(state, gamma)
       local_entropy_violation = max(local_entropy_violation, &
            this%entropy_lower_bound%x(i,1,1,1) - entropy)
    end do
    scale = 1.0_rp
    if (rho%size() .gt. 0) then
       scale = max(scale, maxval(abs(this%entropy_lower_bound%x)))
    end if
    if (local_entropy_violation .gt. tolerance * scale) then
       call neko_error('Euler IDP low-order state violates its local ' // &
            'minimum entropy bound')
    end if

    this%density_second_difference_average%x = 0.0_rp
    do edge = 1, this%graph%n_edges
       associate(a => this%graph%left(:,edge), &
            b => this%graph%right(:,edge))
         pair_average = 0.5_rp * ( &
              this%density_second_difference%x(a(1),a(2),a(3),a(4)) + &
              this%density_second_difference%x(b(1),b(2),b(3),b(4)))
         ! Apply the same occurrence normalization to the neighbour average.
         direction = this%graph%direction(edge)
         left_weight = 2.0_rp / &
              this%graph%directional_degree(direction)%x( &
              a(1),a(2),a(3),a(4))
         right_weight = 2.0_rp / &
              this%graph%directional_degree(direction)%x( &
              b(1),b(2),b(3),b(4))
         this%density_second_difference_average%x( &
              a(1),a(2),a(3),a(4)) = &
              this%density_second_difference_average%x( &
              a(1),a(2),a(3),a(4)) + left_weight * pair_average
         this%density_second_difference_average%x( &
              b(1),b(2),b(3),b(4)) = &
              this%density_second_difference_average%x( &
              b(1),b(2),b(3),b(4)) + right_weight * pair_average
       end associate
    end do
    call gs%op(this%density_second_difference_average, GS_OP_ADD)
    this%density_second_difference_average%x = &
         this%density_second_difference_average%x / &
         (2.0_rp * real(2 * this%graph%n_directions + 1, rp))

    if (this%relax_density_bounds) then
       do i = 1, rho%size()
          strict_lower = this%density_lower_bound%x(i,1,1,1)
          strict_upper = this%density_upper_bound%x(i,1,1,1)
          call euler_idp_relax_density_bounds(strict_lower, strict_upper, &
               this%density_second_difference_average%x(i,1,1,1), &
               this%density_relaxation_mass, this%domain_volume, &
               this%graph%n_directions, relaxed_lower, relaxed_upper)
          this%density_lower_bound%x(i,1,1,1) = relaxed_lower
          this%density_upper_bound%x(i,1,1,1) = relaxed_upper
       end do
    end if
  end subroutine euler_idp_cpu_compute_bounds

  !> Form a high-order Forward Euler candidate without changing the input.
  subroutine euler_idp_cpu_forward_euler(this, rho, m_x, m_y, m_z, energy, &
       coef, gs, slip_bcs, gamma, internal_energy_floor, dt, diagnostics, &
       stage, entropy_viscosity_fraction)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(field_t), intent(in) :: rho, m_x, m_y, m_z, energy
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    type(bc_list_t), intent(inout) :: slip_bcs
    real(kind=rp), intent(in) :: gamma, internal_energy_floor
    real(kind=rp), intent(in) :: dt
    type(euler_idp_diagnostics_t), intent(inout) :: diagnostics
    integer, intent(in), optional :: stage
    type(field_t), intent(in), optional :: entropy_viscosity_fraction
    real(kind=rp) :: local_change(EULER_IDP_NCOMP)
    real(kind=rp) :: local_minimum(4), global_minimum(4)
    real(kind=rp) :: state(EULER_IDP_NCOMP), residual(EULER_IDP_NCOMP)
    real(kind=rp) :: local_floor_timestep, global_floor_timestep
    real(kind=rp) :: local_error, local_scale
    real(kind=rp) :: local_bound_violation(2), local_entropy_violation
    real(kind=rp) :: entropy
    real(kind=rp) :: high_order_fraction, state_difference
    real(kind=rp) :: directional_error_local(EULER_IDP_NCOMP)
    character(len=2 * LOG_SIZE) :: message
    integer :: component, edge, i, ierr

    call diagnostics%reset()
    if (present(stage)) diagnostics%stage = stage
    call this%evaluate_high_order(rho, m_x, m_y, m_z, energy, coef, gs, &
         gamma, internal_energy_floor, diagnostics, &
         entropy_viscosity_fraction)
    call this%evaluate_low_order(rho, m_x, m_y, m_z, energy, coef, gs, &
         gamma, internal_energy_floor, diagnostics, &
         present(entropy_viscosity_fraction))

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

    local_floor_timestep = this%maximum_graph_timestep
    do i = 1, rho%size()
       state = [rho%x(i,1,1,1), m_x%x(i,1,1,1), m_y%x(i,1,1,1), &
            m_z%x(i,1,1,1), energy%x(i,1,1,1)]
       do component = 1, EULER_IDP_NCOMP
          residual(component) = &
               this%low_assembled_residual(component)%x(i,1,1,1)
       end do
       local_floor_timestep = min(local_floor_timestep, &
            euler_idp_internal_energy_timestep(state, residual, &
            internal_energy_floor, this%maximum_graph_timestep))
    end do
    call MPI_Allreduce(local_floor_timestep, global_floor_timestep, 1, &
         MPI_REAL_PRECISION, MPI_MIN, NEKO_COMM, ierr)
    this%maximum_floor_timestep = global_floor_timestep
    diagnostics%maximum_floor_timestep = global_floor_timestep
    if (dt .gt. global_floor_timestep * (1.0_rp + 32.0_rp * &
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

    this%low_candidate(1)%x = rho%x - dt * &
         this%low_assembled_residual(1)%x
    this%low_candidate(2)%x = m_x%x - dt * &
         this%low_assembled_residual(2)%x
    this%low_candidate(3)%x = m_y%x - dt * &
         this%low_assembled_residual(3)%x
    this%low_candidate(4)%x = m_z%x - dt * &
         this%low_assembled_residual(4)%x
    this%low_candidate(5)%x = energy%x - dt * &
         this%low_assembled_residual(5)%x

    call euler_idp_cpu_primitives(this, this%low_candidate(1), &
         this%low_candidate(2), this%low_candidate(3), &
         this%low_candidate(4), this%low_candidate(5), gamma, &
         internal_energy_floor, 'low-order candidate')

    call this%compute_bounds(rho, m_x, m_y, m_z, energy, gs, gamma, &
         diagnostics)

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

    ! The low-candidate admissibility check above overwrites the primitive
    ! work fields. Restore primitives of the stage state before constructing
    ! its directional Euler fluxes.
    call euler_idp_cpu_primitives(this, rho, m_x, m_y, m_z, energy, gamma, &
         internal_energy_floor, 'correction reconstruction state')

    ! Reconstruct the low-to-high inviscid correction independently in each
    ! tensor direction. On affine elements the two SBP operators have the same
    ! line boundary flux, so every line correction is compatible and admits a
    ! nearest-neighbour cumulative flux without a graph-Poisson solve.
    directional_error_local = 0.0_rp
    do component = 1, EULER_IDP_NCOMP
       call euler_idp_cpu_flux(this, component, m_x, m_y, m_z, energy)
       call this%graph%directional_correction_flux( &
            this%correction_flux(component,:), this%flux_x%x, &
            this%flux_y%x, this%flux_z%x, coef, &
            directional_error_local(component))
       this%correction_flux(component,:) = dt * &
            this%correction_flux(component,:)
       do edge = 1, this%graph%n_edges
          associate(a => this%graph%left(:,edge), &
               b => this%graph%right(:,edge))
            state_difference = euler_idp_cpu_state_component(component, a, &
                 rho, m_x, m_y, m_z, energy) - &
                 euler_idp_cpu_state_component(component, b, rho, m_x, m_y, &
                 m_z, energy)
          end associate
          ! Complete the raw low-to-high correction with the low-order graph
          ! viscosity. The directional reconstruction above contains only the
          ! difference between the sparse and dense inviscid operators.
          this%correction_flux(component,edge) = &
               this%correction_flux(component,edge) + dt * &
               this%edge_viscosity(edge) * state_difference
          high_order_fraction = 1.0_rp
          if (this%low_order_only) then
             high_order_fraction = 0.0_rp
          else if (this%edge_viscosity(edge) .gt. tiny(1.0_rp)) then
             high_order_fraction = 1.0_rp - min(1.0_rp, &
                  this%edge_entropy_diffusion(edge) / &
                  this%edge_viscosity(edge))
          end if
          ! Entropy viscosity blends the complete low-to-high subcell flux.
          ! The sparse and dense inviscid operators differ, so damping only
          ! the state-jump term would not recover the low-order endpoint.
          this%correction_flux(component,edge) = high_order_fraction * &
               this%correction_flux(component,edge)
       end do
       call this%graph%incidence(this%flux_x%x, &
            this%correction_flux(component,:))
       call gs%op(this%flux_x, GS_OP_ADD)
       this%candidate(component)%x = this%low_candidate(component)%x + &
            this%flux_x%x / this%graph%mass%x
       local_change(component) = 0.0_rp
       if (rho%size() .gt. 0) then
          local_scale = max(1.0_rp, maxval(abs(this%graph%mass%x * &
               (this%candidate(component)%x - &
               this%low_candidate(component)%x))))
          local_change(component) = &
               maxval(abs(this%flux_x%x - this%graph%mass%x * &
               (this%candidate(component)%x - &
               this%low_candidate(component)%x))) / local_scale
       end if
    end do
    if (maxval(directional_error_local) .gt. &
         10.0_rp * this%correction_tolerance) then
       call neko_error('Euler IDP directional correction is not ' // &
            'line-compatible')
    end if
    diagnostics%reconstruction_residual = local_change
    if (maxval(local_change) .gt. &
         10.0_rp * this%correction_tolerance) then
       call neko_error('Euler IDP directional correction reconstruction ' // &
            'failed')
    end if

    diagnostics%limiter_weight_error = this%limiter_weight_error
    call this%compute_limiter(gamma, internal_energy_floor, diagnostics)
    call this%apply_correction(gs)
    call slip_bcs%apply_vector(this%limited_candidate(2)%x, &
         this%limited_candidate(3)%x, this%limited_candidate(4)%x, &
         rho%size(), strong = .true.)

    local_bound_violation = 0.0_rp
    local_scale = 1.0_rp
    if (rho%size() .gt. 0) then
       local_bound_violation(1) = max(0.0_rp, maxval( &
            this%density_lower_bound%x - &
            this%limited_candidate(1)%x))
       local_bound_violation(2) = max(0.0_rp, maxval( &
            this%limited_candidate(1)%x - &
            this%density_upper_bound%x))
       local_scale = max(1.0_rp, maxval(abs(this%limited_candidate(1)%x)), &
            maxval(abs(this%density_lower_bound%x)), &
            maxval(abs(this%density_upper_bound%x)))
    end if
    diagnostics%max_density_lower_violation = local_bound_violation(1)
    diagnostics%max_density_upper_violation = local_bound_violation(2)
    if (maxval(local_bound_violation) .gt. &
         512.0_rp * epsilon(1.0_rp) * local_scale) then
       call neko_error('Euler IDP limited density violates its local bounds')
    end if

    local_entropy_violation = 0.0_rp
    do i = 1, rho%size()
       state = [this%limited_candidate(1)%x(i,1,1,1), &
            this%limited_candidate(2)%x(i,1,1,1), &
            this%limited_candidate(3)%x(i,1,1,1), &
            this%limited_candidate(4)%x(i,1,1,1), &
            this%limited_candidate(5)%x(i,1,1,1)]
       entropy = euler_idp_specific_entropy(state, gamma)
       local_entropy_violation = max(local_entropy_violation, &
            this%entropy_lower_bound%x(i,1,1,1) - entropy)
    end do
    diagnostics%max_entropy_lower_violation = local_entropy_violation
    local_scale = 1.0_rp
    if (rho%size() .gt. 0) then
       local_scale = max(local_scale, &
            maxval(abs(this%entropy_lower_bound%x)))
    end if
    if (local_entropy_violation .gt. &
         512.0_rp * epsilon(1.0_rp) * local_scale) then
       call neko_error('Euler IDP limited state violates its local ' // &
            'minimum entropy bound')
    end if

    call euler_idp_cpu_primitives(this, this%limited_candidate(1), &
         this%limited_candidate(2), this%limited_candidate(3), &
         this%limited_candidate(4), this%limited_candidate(5), gamma, &
         internal_energy_floor, 'limited candidate')

    local_minimum = huge(1.0_rp)
    if (this%limited_candidate(1)%size() .gt. 0) then
       local_minimum(1:3) = [minval(this%limited_candidate(1)%x), &
            minval(this%internal_energy%x), minval(this%p%x)]
       do i = 1, rho%size()
          state = [this%limited_candidate(1)%x(i,1,1,1), &
               this%limited_candidate(2)%x(i,1,1,1), &
               this%limited_candidate(3)%x(i,1,1,1), &
               this%limited_candidate(4)%x(i,1,1,1), &
               this%limited_candidate(5)%x(i,1,1,1)]
          local_minimum(4) = min(local_minimum(4), &
               euler_idp_specific_entropy(state, gamma))
       end do
    end if
    call MPI_Allreduce(local_minimum, global_minimum, 4, &
         MPI_REAL_PRECISION, MPI_MIN, NEKO_COMM, ierr)
    diagnostics%min_density = global_minimum(1)
    diagnostics%min_internal_energy = global_minimum(2)
    diagnostics%min_pressure = global_minimum(3)
    diagnostics%min_specific_entropy = global_minimum(4)
  end subroutine euler_idp_cpu_forward_euler

  !> Advance with Forward Euler or SSPRK3 using the limited Euler map.
  subroutine euler_idp_cpu_advance(this, rho, m_x, m_y, m_z, energy, coef, &
       gs, slip_bcs, gamma, internal_energy_floor, dt, order, &
       diagnostics, entropy_viscosity_fraction)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(field_t), intent(inout) :: rho, m_x, m_y, m_z, energy
    type(coef_t), intent(inout) :: coef
    type(gs_t), intent(inout) :: gs
    type(bc_list_t), intent(inout) :: slip_bcs
    real(kind=rp), intent(in) :: gamma, internal_energy_floor
    real(kind=rp), intent(in) :: dt
    integer, intent(in) :: order
    type(euler_idp_diagnostics_t), intent(out) :: diagnostics
    type(field_t), intent(in), optional :: entropy_viscosity_fraction
    class(bc_t), pointer :: boundary
    integer :: component, i, j

    if (order .ne. 1 .and. order .ne. 3) then
       call neko_error('Euler IDP requires Forward Euler or SSPRK3')
    end if

    if (.not. this%wall_mask_initialized) then
       this%wall_mask%x = 0.0_rp
       do i = 1, slip_bcs%size()
          boundary => slip_bcs%get(i)
          do j = 1, boundary%msk(0)
             this%wall_mask%x(boundary%msk(j),1,1,1) = 1.0_rp
          end do
       end do
       call gs%op(this%wall_mask, GS_OP_ADD)
       where (this%wall_mask%x .gt. 0.0_rp)
          this%wall_mask%x = 1.0_rp
       end where
       this%wall_mask_initialized = .true.
    end if

    call slip_bcs%apply_vector(m_x%x, m_y%x, m_z%x, rho%size(), &
         strong = .true.)

    do component = 1, 3
       call this%stage_diagnostics(component)%reset()
    end do
    if (order .eq. 3) then
       this%saved_state(1)%x = rho%x
       this%saved_state(2)%x = m_x%x
       this%saved_state(3)%x = m_y%x
       this%saved_state(4)%x = m_z%x
       this%saved_state(5)%x = energy%x
    end if
    call this%forward_euler(rho, m_x, m_y, m_z, energy, coef, gs, slip_bcs, &
         gamma, &
         internal_energy_floor, dt, &
         this%stage_diagnostics(1), 1, entropy_viscosity_fraction)
    if (order .eq. 1) then
       call this%commit(rho, m_x, m_y, m_z, energy)
       diagnostics = this%stage_diagnostics(1)
       return
    end if

    ! U(1) = FE(U(n)).
    call this%commit(rho, m_x, m_y, m_z, energy)
    call this%forward_euler(rho, m_x, m_y, m_z, energy, coef, gs, slip_bcs, &
         gamma, &
         internal_energy_floor, dt, &
         this%stage_diagnostics(2), 2, entropy_viscosity_fraction)

    ! U(2) = 3/4 U(n) + 1/4 FE(U(1)).
    rho%x = 0.75_rp * this%saved_state(1)%x + &
         0.25_rp * this%limited_candidate(1)%x
    m_x%x = 0.75_rp * this%saved_state(2)%x + &
         0.25_rp * this%limited_candidate(2)%x
    m_y%x = 0.75_rp * this%saved_state(3)%x + &
         0.25_rp * this%limited_candidate(3)%x
    m_z%x = 0.75_rp * this%saved_state(4)%x + &
         0.25_rp * this%limited_candidate(4)%x
    energy%x = 0.75_rp * this%saved_state(5)%x + &
         0.25_rp * this%limited_candidate(5)%x
    call euler_idp_cpu_primitives(this, rho, m_x, m_y, m_z, energy, gamma, &
         internal_energy_floor, 'SSPRK3 stage 2 state')
    call this%forward_euler(rho, m_x, m_y, m_z, energy, coef, gs, slip_bcs, &
         gamma, &
         internal_energy_floor, dt, &
         this%stage_diagnostics(3), 3, entropy_viscosity_fraction)

    ! U(n+1) = 1/3 U(n) + 2/3 FE(U(2)).
    rho%x = this%saved_state(1)%x / 3.0_rp + &
         2.0_rp * this%limited_candidate(1)%x / 3.0_rp
    m_x%x = this%saved_state(2)%x / 3.0_rp + &
         2.0_rp * this%limited_candidate(2)%x / 3.0_rp
    m_y%x = this%saved_state(3)%x / 3.0_rp + &
         2.0_rp * this%limited_candidate(3)%x / 3.0_rp
    m_z%x = this%saved_state(4)%x / 3.0_rp + &
         2.0_rp * this%limited_candidate(4)%x / 3.0_rp
    energy%x = this%saved_state(5)%x / 3.0_rp + &
         2.0_rp * this%limited_candidate(5)%x / 3.0_rp
    call euler_idp_cpu_primitives(this, rho, m_x, m_y, m_z, energy, gamma, &
         internal_energy_floor, 'SSPRK3 final state')
    diagnostics = this%stage_diagnostics(3)
  end subroutine euler_idp_cpu_advance

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
    real(kind=rp) :: local_minimum, global_minimum
    integer :: local_count(5), global_count(5)
    integer :: direction, edge, component, ierr
    logical :: density_limited, energy_limited, entropy_limited

    local_minimum = 1.0_rp
    local_count = 0
    do edge = 1, this%graph%n_edges
       associate(a => this%graph%left(:,edge), &
            b => this%graph%right(:,edge))
         direction = this%graph%direction(edge)
         do component = 1, EULER_IDP_NCOMP
            left_base(component) = this%low_candidate(component)%x( &
                 a(1),a(2),a(3),a(4))
            right_base(component) = this%low_candidate(component)%x( &
                 b(1),b(2),b(3),b(4))
         end do
         ! This is correction_flux / (mass * occurrence_weight), with the
         ! directional occurrence weight constructed during initialization.
         left_correction = this%correction_flux(:,edge) * &
              real(this%graph%n_directions, rp) * &
              this%graph%directional_degree(direction)%x( &
              a(1),a(2),a(3),a(4)) / &
              this%graph%mass%x(a(1),a(2),a(3),a(4))
         right_correction = -this%correction_flux(:,edge) * &
              real(this%graph%n_directions, rp) * &
              this%graph%directional_degree(direction)%x( &
              b(1),b(2),b(3),b(4)) / &
              this%graph%mass%x(b(1),b(2),b(3),b(4))
         call euler_idp_limit_edge(left_base, right_base, &
              left_correction, right_correction, &
              this%density_lower_bound%x(a(1),a(2),a(3),a(4)), &
              this%density_upper_bound%x(a(1),a(2),a(3),a(4)), &
              this%density_lower_bound%x(b(1),b(2),b(3),b(4)), &
              this%density_upper_bound%x(b(1),b(2),b(3),b(4)), &
              this%entropy_lower_bound%x(a(1),a(2),a(3),a(4)), &
              this%entropy_lower_bound%x(b(1),b(2),b(3),b(4)), gamma, &
              internal_energy_floor, this%edge_limiter(edge), &
              density_limited, energy_limited, entropy_limited)
       end associate
       if (.not. ieee_is_finite(this%edge_limiter(edge)) .or. &
            this%edge_limiter(edge) .lt. 0.0_rp .or. &
            this%edge_limiter(edge) .gt. 1.0_rp) then
          call neko_error('Euler IDP produced an invalid edge limiter')
       end if
       local_minimum = min(local_minimum, this%edge_limiter(edge))
       local_count(1) = local_count(1) + 1
       if (this%edge_limiter(edge) .lt. &
            1.0_rp - 32.0_rp * epsilon(1.0_rp)) then
          local_count(2) = local_count(2) + 1
       end if
       if (density_limited) local_count(3) = local_count(3) + 1
       if (energy_limited) local_count(4) = local_count(4) + 1
       if (entropy_limited) local_count(5) = local_count(5) + 1
    end do

    call MPI_Allreduce(local_minimum, global_minimum, 1, &
         MPI_REAL_PRECISION, MPI_MIN, NEKO_COMM, ierr)
    call MPI_Allreduce(local_count, global_count, 5, MPI_INTEGER, MPI_SUM, &
         NEKO_COMM, ierr)

    diagnostics%min_limiter = global_minimum
    if (global_count(1) .gt. 0) then
       diagnostics%limited_edge_fraction = &
            real(global_count(2), rp) / real(global_count(1), rp)
    else
       diagnostics%min_limiter = 1.0_rp
       diagnostics%limited_edge_fraction = 0.0_rp
    end if
    diagnostics%density_limited_edges = global_count(3)
    diagnostics%internal_energy_limited_edges = global_count(4)
    diagnostics%entropy_limited_edges = global_count(5)
  end subroutine euler_idp_cpu_compute_limiter

  !> Assemble the limited edge correction on the replicated nodal state.
  subroutine euler_idp_cpu_apply_correction(this, gs)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(gs_t), intent(inout) :: gs
    integer :: component

    do component = 1, EULER_IDP_NCOMP
       call this%graph%incidence(this%flux_x%x, &
            this%edge_limiter * this%correction_flux(component,:))
       call gs%op(this%flux_x, GS_OP_ADD)
       this%limited_candidate(component)%x = &
            this%low_candidate(component)%x + &
            this%flux_x%x / this%graph%mass%x
    end do
  end subroutine euler_idp_cpu_apply_correction

  !> Commit the checked limited Forward Euler candidate.
  subroutine euler_idp_cpu_commit(this, rho, m_x, m_y, m_z, energy)
    class(euler_idp_cpu_t), intent(in) :: this
    type(field_t), intent(inout) :: rho, m_x, m_y, m_z, energy

    rho%x = this%limited_candidate(1)%x
    m_x%x = this%limited_candidate(2)%x
    m_y%x = this%limited_candidate(3)%x
    m_z%x = this%limited_candidate(4)%x
    energy%x = this%limited_candidate(5)%x
  end subroutine euler_idp_cpu_commit

  !> Reconstruct primitive fields from one conserved state.
  subroutine euler_idp_cpu_primitives(this, rho, m_x, m_y, m_z, energy, &
       gamma, internal_energy_floor, label)
    class(euler_idp_cpu_t), intent(inout) :: this
    type(field_t), intent(in) :: rho, m_x, m_y, m_z, energy
    real(kind=rp), intent(in) :: gamma, internal_energy_floor
    character(len=*), intent(in) :: label
    character(len=LOG_SIZE) :: message
    integer :: first_invalid, n

    n = rho%size()
    call compressible_ops_cpu_conserved_to_primitive(rho%x, m_x%x, m_y%x, &
         m_z%x, energy%x, gamma, internal_energy_floor, &
         this%rho%x, this%u%x, this%v%x, this%w%x, this%p%x, &
         this%sound_speed%x, this%internal_energy%x, this%state_status, n)
    if (any(this%state_status .ne. EULER_STATE_OK)) then
       first_invalid = minloc(this%state_status, dim = 1, &
            mask = this%state_status .ne. EULER_STATE_OK)
       write(message, '(A,A,A,I0,A,I0,A,I0)') 'Euler IDP ', trim(label), &
            ' is invalid on rank ', pe_rank, ' at local node ', &
            first_invalid, ', status ', this%state_status(first_invalid)
       call neko_error(trim(message))
    end if
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
       value = rho%x(index(1),index(2),index(3),index(4))
    case (2)
       value = m_x%x(index(1),index(2),index(3),index(4))
    case (3)
       value = m_y%x(index(1),index(2),index(3),index(4))
    case (4)
       value = m_z%x(index(1),index(2),index(3),index(4))
    case (5)
       value = energy%x(index(1),index(2),index(3),index(4))
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
                  (this%flux_x%x(i,ly,k,e) * normal(1) + &
                  this%flux_y%x(i,ly,k,e) * normal(2) + &
                  this%flux_z%x(i,ly,k,e) * normal(3))
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
                  (this%flux_x%x(i,j,lz,e) * normal(1) + &
                  this%flux_y%x(i,j,lz,e) * normal(2) + &
                  this%flux_z%x(i,j,lz,e) * normal(3))
          end do
       end do
       element_boundary_flux(component,e) = boundary_flux
    end do
  end subroutine euler_idp_cpu_element_sums

end module euler_idp_cpu
