! Integration-test driver for the periodic Euler IDP cases.
module user
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use neko
  use comm, only : NEKO_COMM, MPI_REAL_PRECISION
  use device, only : HOST_TO_DEVICE, DEVICE_TO_HOST
  use fluid_scheme_compressible_ns, only : fluid_scheme_compressible_ns_t
  use mpi_f08, only : MPI_Allreduce, MPI_MAX, MPI_MIN
  implicit none

  character(len=:), allocatable :: problem
  integer :: polynomial_order = 0
  real(kind=rp) :: gamma_ref = 1.4_rp
  real(kind=rp) :: energy_floor = 1.0e-12_rp
  real(kind=rp) :: initial_integrals(5) = 0.0_rp
  real(kind=rp) :: minimum_stage_density = huge(1.0_rp)
  real(kind=rp) :: minimum_stage_internal_energy = huge(1.0_rp)
  real(kind=rp) :: minimum_limiter = 1.0_rp
  real(kind=rp) :: minimum_mean_limiter = 1.0_rp
  real(kind=rp) :: maximum_limiter = 0.0_rp
  real(kind=rp) :: maximum_limited_fraction = 0.0_rp
  real(kind=rp) :: maximum_entropy_viscosity = 0.0_rp
  real(kind=rp) :: maximum_mean_entropy_viscosity = 0.0_rp
  real(kind=rp) :: maximum_graph_cfl = 0.0_rp
  real(kind=rp) :: maximum_density_lower_violation = 0.0_rp
  real(kind=rp) :: maximum_density_upper_violation = 0.0_rp
  real(kind=rp) :: maximum_entropy_lower_violation = 0.0_rp
  real(kind=rp) :: maximum_stage_conservation(5) = 0.0_rp
  real(kind=rp) :: maximum_reconstruction_residual(5) = 0.0_rp
  real(kind=rp) :: maximum_correction_compatibility(5) = 0.0_rp
  real(kind=rp) :: boundary_metrics(8) = 0.0_rp
  integer :: stage_count = 0
  integer :: maximum_limiter_counts(3) = 0
  logical :: all_stages_finite = .true.

contains

  !> Register the test callbacks.
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    user%startup => startup
    user%initial_conditions => initial_conditions
    user%initialize => initialize
    user%compute => monitor_stages
    user%finalize => finalize
  end subroutine user_setup

  !> Read test-only settings from the generated case file.
  subroutine startup(params)
    type(json_file), intent(inout) :: params

    call json_get(params, 'test.problem', problem)
    call json_get(params, 'case.numerics.polynomial_order', &
         polynomial_order)
    call json_get(params, 'case.fluid.gamma', gamma_ref)
    call json_get(params, &
         'case.numerics.euler_idp.internal_energy_floor', energy_floor)
  end subroutine startup

  !> Set one of the Euler states exercised by pytest.
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(field_t), pointer :: rho, u, v, w, p
    real(kind=rp) :: x, y, outflow_pressure
    integer :: i

    if (trim(scheme_name) .ne. 'fluid') return

    rho => fields%get_by_name('fluid_rho')
    u => fields%get_by_name('u')
    v => fields%get_by_name('v')
    w => fields%get_by_name('w')
    p => fields%get_by_name('p')

    ! Resolve thermal energy above the O(1) kinetic energy in the moving gas.
    outflow_pressure = max(1.0e-12_rp, 32.0_rp * epsilon(1.0_rp), &
         (gamma_ref - 1.0_rp) * energy_floor)
    do i = 1, rho%size()
       x = rho%dof%x%x(i, 1, 1, 1)
       y = rho%dof%y%x(i, 1, 1, 1)

       select case (trim(problem))
       case ('free_stream')
          rho%x(i, 1, 1, 1) = 1.0_rp
          u%x(i, 1, 1, 1) = 0.7_rp
          v%x(i, 1, 1, 1) = -0.2_rp
          w%x(i, 1, 1, 1) = 0.1_rp
          p%x(i, 1, 1, 1) = 1.0_rp
       case ('low_pressure_free_stream')
          rho%x(i,1,1,1) = 1.0_rp
          u%x(i,1,1,1) = 0.0_rp
          v%x(i,1,1,1) = 0.0_rp
          w%x(i,1,1,1) = 0.0_rp
          p%x(i,1,1,1) = 8.0e-13_rp
       case ('smooth_transport')
          rho%x(i, 1, 1, 1) = 1.0_rp + &
               0.2_rp * sin(2.0_rp * pi * (x + y))
          u%x(i, 1, 1, 1) = 2.5_rp
          v%x(i, 1, 1, 1) = -0.5_rp
          w%x(i, 1, 1, 1) = 0.0_rp
          p%x(i, 1, 1, 1) = 1.0_rp
       case ('periodic_discontinuity')
          u%x(i, 1, 1, 1) = 0.0_rp
          v%x(i, 1, 1, 1) = 0.0_rp
          w%x(i, 1, 1, 1) = 0.0_rp
          if (x .lt. 0.25_rp .or. x .ge. 0.75_rp) then
             rho%x(i, 1, 1, 1) = 1.0_rp
             p%x(i, 1, 1, 1) = 1.0_rp
          else
             rho%x(i, 1, 1, 1) = 0.125_rp
             p%x(i, 1, 1, 1) = 0.1_rp
          end if
       case ('near_vacuum')
          u%x(i, 1, 1, 1) = 0.0_rp
          v%x(i, 1, 1, 1) = 0.0_rp
          w%x(i, 1, 1, 1) = 0.0_rp
          if (x .lt. 0.25_rp .or. x .ge. 0.75_rp) then
             rho%x(i, 1, 1, 1) = 1.0_rp
             p%x(i, 1, 1, 1) = &
                  (gamma_ref - 1.0_rp) * 1.0e-1_rp
          else
             rho%x(i, 1, 1, 1) = 1.0e-3_rp
             p%x(i, 1, 1, 1) = &
                  (gamma_ref - 1.0_rp) * 1.0e-10_rp
          end if
       case ('boundary_prescribed', 'boundary_symmetry', 'boundary_slip', &
            'boundary_invalid_density')
          rho%x(i, 1, 1, 1) = 1.25_rp
          u%x(i, 1, 1, 1) = 0.6_rp
          v%x(i, 1, 1, 1) = -0.4_rp
          w%x(i, 1, 1, 1) = 0.2_rp
          p%x(i, 1, 1, 1) = 0.8_rp
       case ('boundary_outflow', 'boundary_normal_outflow')
          rho%x(i, 1, 1, 1) = 1.25_rp
          u%x(i, 1, 1, 1) = 0.6_rp
          v%x(i, 1, 1, 1) = -0.4_rp
          w%x(i, 1, 1, 1) = 0.2_rp
          p%x(i, 1, 1, 1) = outflow_pressure
       case ('boundary_outflow_floor')
          rho%x(i, 1, 1, 1) = 1.25_rp
          u%x(i, 1, 1, 1) = 0.6_rp
          v%x(i, 1, 1, 1) = -0.4_rp
          w%x(i, 1, 1, 1) = 0.2_rp
          p%x(i, 1, 1, 1) = 2.0_rp * outflow_pressure
       case ('boundary_mixed')
          rho%x(i, 1, 1, 1) = 1.25_rp + 0.15_rp * (1.0_rp - x)
          u%x(i, 1, 1, 1) = 0.6_rp + 0.2_rp * (1.0_rp - x)
          v%x(i, 1, 1, 1) = 0.0_rp
          w%x(i, 1, 1, 1) = 0.2_rp - 0.1_rp * (1.0_rp - x)
          p%x(i, 1, 1, 1) = outflow_pressure + &
               (1.0_rp - outflow_pressure) * (1.0_rp - x)
       case default
          call neko_error('Unknown Euler IDP integration-test problem')
       end select
    end do
    call fields%copy_from(HOST_TO_DEVICE, sync = .true.)
  end subroutine initial_conditions

  !> Save the five initial global conserved integrals.
  subroutine initialize(time)
    type(time_state_t), intent(in) :: time
    type(coef_t), pointer :: coef
    type(field_t), pointer :: rho, m_x, m_y, m_z, energy

    coef => neko_user_access%case%fluid%c_Xh
    rho => neko_registry%get_field('fluid_rho')
    m_x => neko_registry%get_field('m_x')
    m_y => neko_registry%get_field('m_y')
    m_z => neko_registry%get_field('m_z')
    energy => neko_registry%get_field('E')

    call rho%copy_from(DEVICE_TO_HOST, sync = .false.)
    call m_x%copy_from(DEVICE_TO_HOST, sync = .false.)
    call m_y%copy_from(DEVICE_TO_HOST, sync = .false.)
    call m_z%copy_from(DEVICE_TO_HOST, sync = .false.)
    call energy%copy_from(DEVICE_TO_HOST, sync = .true.)

    initial_integrals(1) = glsc2(rho%x, coef%B, rho%size())
    initial_integrals(2) = glsc2(m_x%x, coef%B, rho%size())
    initial_integrals(3) = glsc2(m_y%x, coef%B, rho%size())
    initial_integrals(4) = glsc2(m_z%x, coef%B, rho%size())
    initial_integrals(5) = glsc2(energy%x, coef%B, rho%size())

    if (index(trim(problem), 'boundary_') .eq. 1) then
       call probe_boundary_map(time)
    end if
  end subroutine initialize

  !> Exercise the production boundary map on an isolated conserved state.
  subroutine probe_boundary_map(time)
    type(time_state_t), intent(in) :: time
    type(field_t) :: rho, m_x, m_y, m_z, energy
    real(kind=rp), allocatable :: mapped(:,:)
    real(kind=rp) :: density, velocity(3), pressure, total_energy
    real(kind=rp) :: expected_density, expected_velocity(3)
    real(kind=rp) :: expected_pressure, expected_energy
    real(kind=rp) :: local_errors(7), local_minimum
    real(kind=rp) :: internal_energy, kinetic, x, y
    logical :: is_boundary
    integer :: i, ierr, n

    select type (fluid => neko_user_access%case%fluid)
    type is (fluid_scheme_compressible_ns_t)
       call rho%init(fluid%rho%dof, 'boundary_probe_rho')
       call m_x%init(fluid%rho%dof, 'boundary_probe_m_x')
       call m_y%init(fluid%rho%dof, 'boundary_probe_m_y')
       call m_z%init(fluid%rho%dof, 'boundary_probe_m_z')
       call energy%init(fluid%rho%dof, 'boundary_probe_energy')
       n = rho%size()
       allocate(mapped(5, n))

       density = 1.25_rp
       velocity = [0.6_rp, -0.4_rp, 0.2_rp]
       pressure = 0.8_rp
       total_energy = pressure / (gamma_ref - 1.0_rp) + &
            0.5_rp * density * dot_product(velocity, velocity)
       rho%x = density
       m_x%x = density * velocity(1)
       m_y%x = density * velocity(2)
       m_z%x = density * velocity(3)
       energy%x = total_energy
       call rho%copy_from(HOST_TO_DEVICE, sync = .false.)
       call m_x%copy_from(HOST_TO_DEVICE, sync = .false.)
       call m_y%copy_from(HOST_TO_DEVICE, sync = .false.)
       call m_z%copy_from(HOST_TO_DEVICE, sync = .false.)
       call energy%copy_from(HOST_TO_DEVICE, sync = .true.)

       call fluid%euler_idp_solver%backend%apply_boundary_conditions( &
            rho, m_x, m_y, m_z, energy, fluid%bcs_density, &
            fluid%bcs_vel, fluid%bcs_prs, gamma_ref, energy_floor, time, &
            'boundary regression probe', refresh = .true.)
       call rho%copy_from(DEVICE_TO_HOST, sync = .false.)
       call m_x%copy_from(DEVICE_TO_HOST, sync = .false.)
       call m_y%copy_from(DEVICE_TO_HOST, sync = .false.)
       call m_z%copy_from(DEVICE_TO_HOST, sync = .false.)
       call energy%copy_from(DEVICE_TO_HOST, sync = .true.)

       local_errors = 0.0_rp
       local_minimum = huge(1.0_rp)
       do i = 1, n
          x = rho%dof%x%x(i, 1, 1, 1)
          y = rho%dof%y%x(i, 1, 1, 1)
          call expected_boundary_state(x, y, expected_density, &
               expected_velocity, expected_pressure, expected_energy, &
               is_boundary)
          if (is_boundary) then
             density = rho%x(i, 1, 1, 1)
             velocity = [m_x%x(i, 1, 1, 1), m_y%x(i, 1, 1, 1), &
                  m_z%x(i, 1, 1, 1)] / density
             kinetic = 0.5_rp * density * dot_product(velocity, velocity)
             internal_energy = energy%x(i, 1, 1, 1) - kinetic
             pressure = (gamma_ref - 1.0_rp) * internal_energy
             local_errors(1) = max(local_errors(1), &
                  abs(density - expected_density))
             local_errors(2:4) = max(local_errors(2:4), &
                  abs(velocity - expected_velocity))
             local_errors(5) = max(local_errors(5), &
                  abs(pressure - expected_pressure))
             local_errors(6) = max(local_errors(6), &
                  abs(energy%x(i, 1, 1, 1) - expected_energy))
             local_minimum = min(local_minimum, internal_energy)
          end if
          mapped(:,i) = [rho%x(i, 1, 1, 1), m_x%x(i, 1, 1, 1), &
               m_y%x(i, 1, 1, 1), m_z%x(i, 1, 1, 1), &
               energy%x(i, 1, 1, 1)]
       end do

       call fluid%euler_idp_solver%backend%apply_boundary_conditions( &
            rho, m_x, m_y, m_z, energy, fluid%bcs_density, &
            fluid%bcs_vel, fluid%bcs_prs, gamma_ref, energy_floor, time, &
            'repeated boundary regression probe', refresh = .true.)
       call rho%copy_from(DEVICE_TO_HOST, sync = .false.)
       call m_x%copy_from(DEVICE_TO_HOST, sync = .false.)
       call m_y%copy_from(DEVICE_TO_HOST, sync = .false.)
       call m_z%copy_from(DEVICE_TO_HOST, sync = .false.)
       call energy%copy_from(DEVICE_TO_HOST, sync = .true.)
       do i = 1, n
          local_errors(7) = max(local_errors(7), maxval(abs( &
               [rho%x(i, 1, 1, 1), m_x%x(i, 1, 1, 1), &
               m_y%x(i, 1, 1, 1), m_z%x(i, 1, 1, 1), &
               energy%x(i, 1, 1, 1)] - mapped(:,i))))
       end do

       call MPI_Allreduce(local_errors, boundary_metrics(1:7), 7, &
            MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM, ierr)
       call MPI_Allreduce(local_minimum, boundary_metrics(8), 1, &
            MPI_REAL_PRECISION, MPI_MIN, NEKO_COMM, ierr)

       deallocate(mapped)
       call rho%free()
       call m_x%free()
       call m_y%free()
       call m_z%free()
       call energy%free()
    class default
       call neko_error('Euler IDP test requires the compressible scheme')
    end select
  end subroutine probe_boundary_map

  !> Return the expected result of the strong primitive boundary map.
  subroutine expected_boundary_state(x, y, density, velocity, pressure, &
       energy, is_boundary)
    real(kind=rp), intent(in) :: x, y
    real(kind=rp), intent(out) :: density, velocity(3), pressure, energy
    logical, intent(out) :: is_boundary
    real(kind=rp) :: base_energy, kinetic, pressure_floor, tolerance

    tolerance = 128.0_rp * epsilon(1.0_rp)
    density = 1.25_rp
    velocity = [0.6_rp, -0.4_rp, 0.2_rp]
    pressure = 0.8_rp
    base_energy = pressure / (gamma_ref - 1.0_rp) + &
         0.5_rp * density * dot_product(velocity, velocity)
    pressure_floor = max(1.0e-12_rp, &
         (gamma_ref - 1.0_rp) * energy_floor)
    is_boundary = x .le. tolerance .or. x .ge. 1.0_rp - tolerance

    select case (trim(problem))
    case ('boundary_prescribed', 'boundary_invalid_density')
       density = 1.4_rp
       velocity = [0.8_rp, -0.2_rp, 0.1_rp]
    case ('boundary_symmetry', 'boundary_slip')
       velocity(1) = 0.0_rp
    case ('boundary_outflow', 'boundary_normal_outflow', &
         'boundary_outflow_floor')
       continue
    case ('boundary_mixed')
       is_boundary = is_boundary .or. y .le. tolerance .or. &
            y .ge. 1.0_rp - tolerance
       if (x .le. tolerance) then
          density = 1.4_rp
          velocity = [0.8_rp, 0.0_rp, 0.1_rp]
       end if
       if (y .le. tolerance .or. y .ge. 1.0_rp - tolerance) then
          velocity(2) = 0.0_rp
       end if
    case default
       call neko_error('Unknown Euler IDP boundary regression problem')
    end select

    kinetic = 0.5_rp * density * dot_product(velocity, velocity)
    pressure = (gamma_ref - 1.0_rp) * (base_energy - kinetic)
    if (trim(problem) .eq. 'boundary_prescribed' .or. &
         trim(problem) .eq. 'boundary_invalid_density') then
       pressure = 1.0_rp
    else if (trim(problem) .eq. 'boundary_outflow' .or. &
         trim(problem) .eq. 'boundary_normal_outflow' .or. &
         trim(problem) .eq. 'boundary_outflow_floor') then
       pressure = pressure_floor
    else if (trim(problem) .eq. 'boundary_mixed') then
       if (x .le. tolerance) pressure = 1.0_rp
       if (x .ge. 1.0_rp - tolerance) pressure = pressure_floor
    end if
    energy = pressure / (gamma_ref - 1.0_rp) + kinetic
  end subroutine expected_boundary_state

  !> Reduce the diagnostics over every active stage and time step.
  subroutine monitor_stages(time)
    type(time_state_t), intent(in) :: time
    integer :: stage

    select type (fluid => neko_user_access%case%fluid)
    type is (fluid_scheme_compressible_ns_t)
       do stage = 1, fluid%euler_idp_solver%time_order
          associate(diagnostic => &
               fluid%euler_idp_solver%stage_diagnostics(stage))
            stage_count = stage_count + 1
            minimum_stage_density = min(minimum_stage_density, &
                 diagnostic%min_density)
            minimum_stage_internal_energy = &
                 min(minimum_stage_internal_energy, &
                 diagnostic%min_internal_energy)
            minimum_limiter = min(minimum_limiter, &
                 diagnostic%min_limiter)
            minimum_mean_limiter = min(minimum_mean_limiter, &
                 diagnostic%mean_limiter)
            maximum_limiter = max(maximum_limiter, &
                 diagnostic%max_limiter)
            maximum_limited_fraction = max(maximum_limited_fraction, &
                 diagnostic%limited_edge_fraction)
            maximum_entropy_viscosity = max(maximum_entropy_viscosity, &
                 diagnostic%max_entropy_viscosity)
            maximum_mean_entropy_viscosity = &
                 max(maximum_mean_entropy_viscosity, &
                 diagnostic%mean_entropy_viscosity)
            maximum_limiter_counts = max(maximum_limiter_counts, &
                 [diagnostic%density_limited_edges, &
                 diagnostic%internal_energy_limited_edges, &
                 diagnostic%entropy_limited_edges])
            maximum_graph_cfl = max(maximum_graph_cfl, &
                 diagnostic%max_graph_cfl)
            maximum_density_lower_violation = &
                 max(maximum_density_lower_violation, &
                 diagnostic%max_density_lower_violation)
            maximum_density_upper_violation = &
                 max(maximum_density_upper_violation, &
                 diagnostic%max_density_upper_violation)
            maximum_entropy_lower_violation = &
                 max(maximum_entropy_lower_violation, &
                 diagnostic%max_entropy_lower_violation)
            maximum_stage_conservation = max(maximum_stage_conservation, &
                 abs(diagnostic%limited_conservation))
            maximum_reconstruction_residual = &
                 max(maximum_reconstruction_residual, &
                 abs(diagnostic%reconstruction_residual))
            maximum_correction_compatibility = &
                 max(maximum_correction_compatibility, &
                 abs(diagnostic%correction_global_compatibility))
            all_stages_finite = all_stages_finite .and. &
                 diagnostic_is_finite(diagnostic%min_density, &
                 diagnostic%min_internal_energy, &
                 diagnostic%min_limiter, diagnostic%mean_limiter, &
                 diagnostic%max_limiter, &
                 diagnostic%limited_edge_fraction, &
                 diagnostic%max_entropy_viscosity, &
                 diagnostic%mean_entropy_viscosity, &
                 diagnostic%max_graph_cfl, &
                 diagnostic%max_density_lower_violation, &
                 diagnostic%max_density_upper_violation, &
                 diagnostic%max_entropy_lower_violation, &
                 diagnostic%limited_conservation, &
                 diagnostic%reconstruction_residual, &
                 diagnostic%correction_global_compatibility)
          end associate
       end do
    class default
       call neko_error('Euler IDP test requires the compressible scheme')
    end select
  end subroutine monitor_stages

  !> Report exact errors, conservation drift, and aggregated diagnostics.
  subroutine finalize(time)
    type(time_state_t), intent(in) :: time
    type(coef_t), pointer :: coef
    type(field_t), pointer :: rho, m_x, m_y, m_z, energy, temperature
    real(kind=rp), allocatable :: work(:)
    real(kind=rp) :: state(5), exact(5), final_integrals(5)
    real(kind=rp) :: errors(5), drifts(5)
    real(kind=rp) :: minimum_density, minimum_internal_energy, temperature_error
    integer :: component, i, n

    coef => neko_user_access%case%fluid%c_Xh
    rho => neko_registry%get_field('fluid_rho')
    m_x => neko_registry%get_field('m_x')
    m_y => neko_registry%get_field('m_y')
    m_z => neko_registry%get_field('m_z')
    temperature => neko_registry%get_field('temperature')
    call temperature%copy_from(DEVICE_TO_HOST, sync = .false.)
    energy => neko_registry%get_field('E')
    call rho%copy_from(DEVICE_TO_HOST, sync = .false.)
    call m_x%copy_from(DEVICE_TO_HOST, sync = .false.)
    call m_y%copy_from(DEVICE_TO_HOST, sync = .false.)
    call m_z%copy_from(DEVICE_TO_HOST, sync = .false.)
    call energy%copy_from(DEVICE_TO_HOST, sync = .true.)
    n = rho%size()
    allocate(work(n))

    final_integrals(1) = glsc2(rho%x, coef%B, n)
    final_integrals(2) = glsc2(m_x%x, coef%B, n)
    final_integrals(3) = glsc2(m_y%x, coef%B, n)
    final_integrals(4) = glsc2(m_z%x, coef%B, n)
    final_integrals(5) = glsc2(energy%x, coef%B, n)
    drifts = abs(final_integrals - initial_integrals)

    do i = 1, n
       work(i) = rho%x(i, 1, 1, 1)
    end do
    minimum_density = glmin(work, n)
    do i = 1, n
       work(i) = energy%x(i, 1, 1, 1) - 0.5_rp * &
            (m_x%x(i, 1, 1, 1)**2 + m_y%x(i, 1, 1, 1)**2 + &
            m_z%x(i, 1, 1, 1)**2) / rho%x(i, 1, 1, 1)
    end do
    minimum_internal_energy = glmin(work, n)
    do i = 1, n
       work(i) = abs(temperature%x(i,1,1,1) - work(i) / rho%x(i,1,1,1))
    end do
    temperature_error = glmax(work, n)

    errors = 0.0_rp
    if (trim(problem) .eq. 'free_stream' .or. &
         trim(problem) .eq. 'low_pressure_free_stream' .or. &
         trim(problem) .eq. 'smooth_transport') then
       do component = 1, 5
          do i = 1, n
             state = [rho%x(i, 1, 1, 1), m_x%x(i, 1, 1, 1), &
                  m_y%x(i, 1, 1, 1), m_z%x(i, 1, 1, 1), &
                  energy%x(i, 1, 1, 1)]
             call exact_state(rho%dof%x%x(i, 1, 1, 1), &
                  rho%dof%y%x(i, 1, 1, 1), real(time%t, rp), exact)
             work(i) = abs(state(component) - exact(component))
          end do
          errors(component) = glmax(work, n)
       end do
    end if

    all_stages_finite = all_stages_finite .and. &
         fields_are_finite(rho, m_x, m_y, m_z, energy)

    if (pe_rank .eq. 0) then
       write(*, '(A,1X,A,4(1X,I0),1X,L1)') 'EULER_IDP_RESULT', &
            trim(problem), polynomial_order, pe_size, time%tstep, &
            stage_count, all_stages_finite
       write(*, '(A,5(1X,ES25.16E3))') 'EULER_IDP_ERROR', errors
       write(*, '(A,1X,ES25.16E3)') 'EULER_IDP_TEMPERATURE_ERROR', &
            temperature_error
       write(*, '(A,5(1X,ES25.16E3))') 'EULER_IDP_DRIFT', drifts
       write(*, '(A,3(1X,ES25.16E3))') 'EULER_IDP_LIMITER', &
            minimum_limiter, maximum_limited_fraction, maximum_graph_cfl
       write(*, '(A,4(1X,ES25.16E3))') 'EULER_IDP_LIMITER_STATS', &
            minimum_limiter, minimum_mean_limiter, maximum_limiter, &
            maximum_limited_fraction
       write(*, '(A,3(1X,I0))') 'EULER_IDP_LIMITER_COUNTS', &
            maximum_limiter_counts
       write(*, '(A,2(1X,ES25.16E3))') 'EULER_IDP_ENTROPY_BLEND', &
            maximum_entropy_viscosity, maximum_mean_entropy_viscosity
       write(*, '(A,4(1X,ES25.16E3))') 'EULER_IDP_STATE', &
            minimum_stage_density, minimum_stage_internal_energy, &
            minimum_density, minimum_internal_energy
       write(*, '(A,3(1X,ES25.16E3))') 'EULER_IDP_BOUNDS', &
            maximum_density_lower_violation, &
            maximum_density_upper_violation, &
            maximum_entropy_lower_violation
       write(*, '(A,5(1X,ES25.16E3))') 'EULER_IDP_CONSERVATION', &
            maximum_stage_conservation
       write(*, '(A,10(1X,ES25.16E3))') 'EULER_IDP_CORRECTION', &
            maximum_reconstruction_residual, &
            maximum_correction_compatibility
       write(*, '(A,8(1X,ES25.16E3))') 'EULER_IDP_BOUNDARY', &
            boundary_metrics
    end if

    deallocate(work)
  end subroutine finalize

  !> Exact conserved state for the constant and smooth transport cases.
  subroutine exact_state(x, y, time, state)
    real(kind=rp), intent(in) :: x, y, time
    real(kind=rp), intent(out) :: state(5)
    real(kind=rp) :: density, velocity(3), pressure

    pressure = 1.0_rp
    if (trim(problem) .eq. 'free_stream') then
       density = 1.0_rp
       velocity = [0.7_rp, -0.2_rp, 0.1_rp]
    else if (trim(problem) .eq. 'low_pressure_free_stream') then
       density = 1.0_rp
       velocity = 0.0_rp
       pressure = 8.0e-13_rp
    else
       density = 1.0_rp + 0.2_rp * sin(2.0_rp * pi * &
            (x + y - 2.0_rp * time))
       velocity = [2.5_rp, -0.5_rp, 0.0_rp]
    end if

    state(1) = density
    state(2:4) = density * velocity
    state(5) = pressure / (gamma_ref - 1.0_rp) + &
         0.5_rp * density * dot_product(velocity, velocity)
  end subroutine exact_state

  !> Check every scalar carried by one stage diagnostic.
  logical function diagnostic_is_finite(rho_min, energy_min, limiter_min, &
       limiter_mean, limiter_max, limited_fraction, entropy_max, &
       entropy_mean, graph_cfl, &
       density_lower_violation, density_upper_violation, entropy_violation, &
       conservation, reconstruction, correction_compatibility) &
       result(finite)
    real(kind=rp), intent(in) :: rho_min, energy_min, limiter_min
    real(kind=rp), intent(in) :: limiter_mean, limiter_max
    real(kind=rp), intent(in) :: limited_fraction, entropy_max
    real(kind=rp), intent(in) :: entropy_mean, graph_cfl
    real(kind=rp), intent(in) :: density_lower_violation
    real(kind=rp), intent(in) :: density_upper_violation, entropy_violation
    real(kind=rp), intent(in) :: conservation(5), reconstruction(5)
    real(kind=rp), intent(in) :: correction_compatibility(5)

    finite = ieee_is_finite(rho_min) .and. &
         ieee_is_finite(energy_min) .and. &
         ieee_is_finite(limiter_min) .and. &
         ieee_is_finite(limiter_mean) .and. &
         ieee_is_finite(limiter_max) .and. &
         ieee_is_finite(limited_fraction) .and. &
         ieee_is_finite(entropy_max) .and. &
         ieee_is_finite(entropy_mean) .and. &
         ieee_is_finite(graph_cfl) .and. &
         ieee_is_finite(density_lower_violation) .and. &
         ieee_is_finite(density_upper_violation) .and. &
         ieee_is_finite(entropy_violation) .and. &
         all(ieee_is_finite(conservation)) .and. &
         all(ieee_is_finite(reconstruction)) .and. &
         all(ieee_is_finite(correction_compatibility))
  end function diagnostic_is_finite

  !> Check the five final conserved fields for NaN or infinity.
  logical function fields_are_finite(rho, m_x, m_y, m_z, energy) &
       result(finite)
    type(field_t), intent(in) :: rho, m_x, m_y, m_z, energy

    finite = all(ieee_is_finite(rho%x)) .and. &
         all(ieee_is_finite(m_x%x)) .and. &
         all(ieee_is_finite(m_y%x)) .and. &
         all(ieee_is_finite(m_z%x)) .and. &
         all(ieee_is_finite(energy%x))
  end function fields_are_finite

end module user
