! Three-body rigid-rotation density transport for the Euler IDP solver.
!
! This is an Euler proxy for the scalar benchmark in Nazarov (2026). The
! density is one plus the scalar profile. After every Euler step, momentum and
! energy are projected back to the prescribed rigid-rotation velocity and a
! constant pressure. Thus, this case tests density transport and the IDP
! machinery; it is not an unforced solution of the Euler equations.
module user
  use neko
  use fluid_scheme_compressible_ns, only : fluid_scheme_compressible_ns_t
  implicit none

  real(kind=rp), parameter :: gamma_ref = 1.4_rp
  real(kind=rp), parameter :: pressure_ref = 1.0_rp
  real(kind=rp), parameter :: density_background = 1.0_rp
  real(kind=rp), parameter :: body_radius = 0.3_rp
  real(kind=rp) :: initial_mass = 0.0_rp
  real(kind=rp) :: minimum_limiter = 1.0_rp
  real(kind=rp) :: maximum_limited_fraction = 0.0_rp
  real(kind=rp) :: maximum_graph_cfl = 0.0_rp

contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    user%initial_conditions => initial_conditions
    user%initialize => initialize
    user%compute => enforce_rotation
    user%entropy_pair => user_entropy_pair
    user%finalize => finalize
  end subroutine user_setup

  !> Initialize density from the three-body scalar profile.
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(field_t), pointer :: rho, u, v, w, p
    real(kind=rp) :: x, y
    integer :: i

    if (trim(scheme_name) .ne. 'fluid') return

    rho => fields%get_by_name('fluid_rho')
    u => fields%get_by_name('u')
    v => fields%get_by_name('v')
    w => fields%get_by_name('w')
    p => fields%get_by_name('p')

    do i = 1, rho%size()
       x = rho%dof%x(i, 1, 1, 1)
       y = rho%dof%y(i, 1, 1, 1)
       rho%x(i, 1, 1, 1) = density_background + body_profile(x, y)
       u%x(i, 1, 1, 1) = -2.0_rp * pi * y
       v%x(i, 1, 1, 1) = 2.0_rp * pi * x
       w%x(i, 1, 1, 1) = 0.0_rp
       p%x(i, 1, 1, 1) = pressure_ref
    end do
  end subroutine initial_conditions

  !> Store the initial density mass after Neko forms conserved variables.
  subroutine initialize(time)
    type(time_state_t), intent(in) :: time
    type(coef_t), pointer :: coef
    type(field_t), pointer :: rho

    coef => neko_user_access%case%fluid%c_Xh
    rho => neko_registry%get_field('fluid_rho')
    initial_mass = glsc2(rho%x, coef%B, rho%size())
  end subroutine initialize

  !> Project the velocity and pressure after every Euler step.
  subroutine enforce_rotation(time)
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: rho, m_x, m_y, m_z, energy
    real(kind=rp) :: kinetic_energy, u_value, v_value, x, y
    integer :: i

    rho => neko_registry%get_field('fluid_rho')
    m_x => neko_registry%get_field('m_x')
    m_y => neko_registry%get_field('m_y')
    m_z => neko_registry%get_field('m_z')
    energy => neko_registry%get_field('E')

    select type (fluid => neko_user_access%case%fluid)
    type is (fluid_scheme_compressible_ns_t)
       do i = 1, rho%size()
          x = rho%dof%x(i, 1, 1, 1)
          y = rho%dof%y(i, 1, 1, 1)
          u_value = -2.0_rp * pi * y
          v_value = 2.0_rp * pi * x
          kinetic_energy = 0.5_rp * rho%x(i, 1, 1, 1) * &
               (u_value**2 + v_value**2)

          m_x%x(i, 1, 1, 1) = rho%x(i, 1, 1, 1) * u_value
          m_y%x(i, 1, 1, 1) = rho%x(i, 1, 1, 1) * v_value
          m_z%x(i, 1, 1, 1) = 0.0_rp
          energy%x(i, 1, 1, 1) = pressure_ref / (gamma_ref - 1.0_rp) + &
               kinetic_energy
          fluid%u%x(i, 1, 1, 1) = u_value
          fluid%v%x(i, 1, 1, 1) = v_value
          fluid%w%x(i, 1, 1, 1) = 0.0_rp
          fluid%p%x(i, 1, 1, 1) = pressure_ref
          fluid%temperature%x(i, 1, 1, 1) = pressure_ref / &
               (rho%x(i, 1, 1, 1) * (gamma_ref - 1.0_rp))
          fluid%S%x(i, 1, 1, 1) = 0.5_rp * &
               rho%x(i, 1, 1, 1)**2
       end do

       call fluid%compute_max_wave_speed()
       call fluid%euler_idp_cpu%update_graph_viscosity(rho, m_x, m_y, &
            m_z, energy, fluid%gs_Xh, gamma_ref)
       minimum_limiter = min(minimum_limiter, &
            fluid%euler_idp_diagnostics%min_limiter)
       maximum_limited_fraction = max(maximum_limited_fraction, &
            fluid%euler_idp_diagnostics%limited_edge_fraction)
       maximum_graph_cfl = max(maximum_graph_cfl, &
            fluid%euler_idp_diagnostics%max_graph_cfl)
    class default
       call neko_error('Three-body rotation requires compressible Euler')
    end select
  end subroutine enforce_rotation

  !> Evaluate the quadratic entropy pair for rigid-body advection.
  subroutine user_entropy_pair(entropy, flux_x, flux_y, flux_z, &
       wave_speed, time)
    type(field_t), intent(inout) :: entropy
    type(field_t), intent(inout) :: flux_x, flux_y, flux_z
    type(field_t), intent(inout) :: wave_speed
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: rho
    real(kind=rp) :: entropy_value, u_value, v_value, x, y
    integer :: i

    rho => neko_registry%get_field('fluid_rho')
    do i = 1, rho%size()
       x = rho%dof%x(i, 1, 1, 1)
       y = rho%dof%y(i, 1, 1, 1)
       u_value = -2.0_rp * pi * y
       v_value = 2.0_rp * pi * x
       entropy_value = 0.5_rp * rho%x(i, 1, 1, 1)**2
       entropy%x(i, 1, 1, 1) = entropy_value
       flux_x%x(i, 1, 1, 1) = u_value * entropy_value
       flux_y%x(i, 1, 1, 1) = v_value * entropy_value
       flux_z%x(i, 1, 1, 1) = 0.0_rp
       wave_speed%x(i, 1, 1, 1) = &
            sqrt(u_value**2 + v_value**2)
    end do
  end subroutine user_entropy_pair

  !> Report errors after one rotation and write the z=0 density plane.
  subroutine finalize(time)
    type(time_state_t), intent(in) :: time
    type(coef_t), pointer :: coef
    type(field_t), pointer :: rho
    real(kind=rp), allocatable :: error(:), exact_profile(:)
    real(kind=rp) :: angle, cosine_angle, sine_angle
    real(kind=rp) :: current_mass, denominator, l1_error, l2_error
    real(kind=rp) :: minimum_density, maximum_density
    real(kind=rp) :: x, x_initial, y, y_initial
    character(len=80) :: filename
    integer :: e, i, j, k, n, output_unit

    coef => neko_user_access%case%fluid%c_Xh
    rho => neko_registry%get_field('fluid_rho')
    n = rho%size()
    allocate(error(n), exact_profile(n))

    angle = 2.0_rp * pi * time%t
    cosine_angle = cos(angle)
    sine_angle = sin(angle)
    do i = 1, n
       x = rho%dof%x(i, 1, 1, 1)
       y = rho%dof%y(i, 1, 1, 1)
       x_initial = cosine_angle * x + sine_angle * y
       y_initial = -sine_angle * x + cosine_angle * y
       exact_profile(i) = body_profile(x_initial, y_initial)
       error(i) = abs(rho%x(i, 1, 1, 1) - density_background - &
            exact_profile(i))
    end do

    current_mass = glsc2(rho%x, coef%B, n)
    denominator = glsc2(exact_profile, coef%B, n)
    l1_error = glsc2(error, coef%B, n) / denominator
    error = error**2
    l2_error = sqrt(glsc2(error, coef%B, n) / coef%volume)
    minimum_density = glmin(rho%x, n)
    maximum_density = glmax(rho%x, n)

    if (pe_rank .eq. 0) then
       write(*, '(A,I0,9(1X,ES25.16E3))') 'EULER_IDP_THREE_BODY ', &
            time%tstep, time%t, minimum_density, maximum_density, &
            l1_error, l2_error, abs(current_mass - initial_mass), &
            minimum_limiter, maximum_limited_fraction, maximum_graph_cfl
    end if

    write(filename, '(A,I0,A)') 'three_body_density_rank_', pe_rank, '.csv'
    open(newunit = output_unit, file = trim(filename), status = 'replace', &
         action = 'write')
    write(output_unit, '(A)') 'x,y,rho,rho_exact'
    do e = 1, rho%dof%msh%nelv
       do k = 1, rho%dof%Xh%lz
          if (abs(rho%dof%z(1, 1, k, e)) .gt. 1.0e-12_rp) cycle
          do j = 1, rho%dof%Xh%ly
             do i = 1, rho%dof%Xh%lx
                x = rho%dof%x(i, j, k, e)
                y = rho%dof%y(i, j, k, e)
                x_initial = cosine_angle * x + sine_angle * y
                y_initial = -sine_angle * x + cosine_angle * y
                write(output_unit, '(4(ES25.16E3,:,","))') x, y, &
                     rho%x(i, j, k, e), density_background + &
                     body_profile(x_initial, y_initial)
             end do
          end do
       end do
    end do
    close(output_unit)
    deallocate(error, exact_profile)
  end subroutine finalize

  !> Three-body profile from section 6.2 of Nazarov (2026).
  pure real(kind=rp) function body_profile(x, y) result(value)
    real(kind=rp), intent(in) :: x, y
    real(kind=rp) :: r1, r2, r3

    r1 = sqrt(x**2 + (y - 0.5_rp)**2)
    r2 = sqrt((x + 0.5_rp)**2 + y**2)
    r3 = sqrt((x - sqrt(3.0_rp) / 4.0_rp)**2 + &
         (y + 0.25_rp)**2)

    if (r1 .le. body_radius .and. &
         (abs(x) .ge. 0.05_rp .or. y .ge. 0.7_rp)) then
       value = 1.0_rp
    else if (r2 .le. body_radius) then
       value = 0.25_rp * (1.0_rp + cos(pi * r2 / body_radius))
    else if (r3 .le. body_radius) then
       value = 1.0_rp - r3 / body_radius
    else
       value = 0.0_rp
    end if
  end function body_profile

end module user
