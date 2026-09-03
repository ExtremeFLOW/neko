! Two-dimensional Burgers scalar proxy for the Euler IDP solver.
!
! The published signed scalar is shifted by one, rho = u + 1, to keep Euler
! density positive. After every Forward Euler step, momentum is reset to
! f(u) = 0.5 * (u**2, u**2), so the next density residual is exactly the
! Burgers flux. Energy is reset to maintain constant pressure. This is a test
! harness for density transport and limiting, not an unforced Euler flow.
module user
  use neko
  use entropy_viscosity, only : entropy_viscosity_t
  use fluid_scheme_compressible_ns, only : fluid_scheme_compressible_ns_t
  implicit none

  real(kind=rp), parameter :: gamma_ref = 1.4_rp
  real(kind=rp), parameter :: pressure_ref = 1.0_rp
  real(kind=rp), parameter :: density_shift = 1.0_rp
  real(kind=rp), parameter :: scalar_minimum = -0.75_rp
  real(kind=rp), parameter :: scalar_maximum = 1.0_rp
  real(kind=rp) :: initial_mass = 0.0_rp
  real(kind=rp) :: minimum_limiter = 1.0_rp
  real(kind=rp) :: maximum_limited_fraction = 0.0_rp
  real(kind=rp) :: maximum_graph_cfl = 0.0_rp

contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    user%initial_conditions => initial_conditions
    user%initialize => initialize
    user%compute => enforce_burgers_flux
    user%entropy_pair => user_entropy_pair
    user%finalize => finalize
  end subroutine user_setup

  !> Initialize shifted density and the two-dimensional Burgers flux.
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(field_t), pointer :: rho, u, v, w, p
    real(kind=rp) :: flux_value, rho_value, scalar_value, x, y
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
       scalar_value = burgers_initial_scalar(x, y)
       rho_value = scalar_value + density_shift
       flux_value = 0.5_rp * scalar_value**2
       rho%x(i, 1, 1, 1) = rho_value
       u%x(i, 1, 1, 1) = flux_value / rho_value
       v%x(i, 1, 1, 1) = flux_value / rho_value
       w%x(i, 1, 1, 1) = 0.0_rp
       p%x(i, 1, 1, 1) = pressure_ref
    end do
  end subroutine initial_conditions

  !> Store the initial shifted-density mass.
  subroutine initialize(time)
    type(time_state_t), intent(in) :: time
    type(coef_t), pointer :: coef
    type(field_t), pointer :: rho

    coef => neko_user_access%case%fluid%c_Xh
    rho => neko_registry%get_field('fluid_rho')
    initial_mass = glsc2(rho%x, coef%B, rho%size())
  end subroutine initialize

  !> Reset momentum and energy so the next density flux is the Burgers flux.
  subroutine enforce_burgers_flux(time)
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: rho, m_x, m_y, m_z, energy
    real(kind=rp) :: flux_value, kinetic_energy, rho_value, scalar_value
    integer :: i

    rho => neko_registry%get_field('fluid_rho')
    m_x => neko_registry%get_field('m_x')
    m_y => neko_registry%get_field('m_y')
    m_z => neko_registry%get_field('m_z')
    energy => neko_registry%get_field('E')

    select type (fluid => neko_user_access%case%fluid)
    type is (fluid_scheme_compressible_ns_t)
       do i = 1, rho%size()
          rho_value = rho%x(i, 1, 1, 1)
          if (rho_value .le. 0.0_rp) then
             call neko_error('Burgers shifted density must remain positive')
          end if

          scalar_value = rho_value - density_shift
          flux_value = 0.5_rp * scalar_value**2
          kinetic_energy = flux_value**2 / rho_value

          m_x%x(i, 1, 1, 1) = flux_value
          m_y%x(i, 1, 1, 1) = flux_value
          m_z%x(i, 1, 1, 1) = 0.0_rp
          energy%x(i, 1, 1, 1) = pressure_ref / (gamma_ref - 1.0_rp) + &
               kinetic_energy
          fluid%u%x(i, 1, 1, 1) = flux_value / rho_value
          fluid%v%x(i, 1, 1, 1) = flux_value / rho_value
          fluid%w%x(i, 1, 1, 1) = 0.0_rp
          fluid%p%x(i, 1, 1, 1) = pressure_ref
          fluid%temperature%x(i, 1, 1, 1) = pressure_ref / &
               (rho_value * (gamma_ref - 1.0_rp))
          fluid%S%x(i, 1, 1, 1) = scalar_value
       end do

       call fluid%compute_max_wave_speed()
       select type (reg => fluid%regularization)
       type is (entropy_viscosity_t)
          call reg%evaluate_user_entropy_pair(time)
          call fluid%euler_idp_cpu%update_graph_viscosity(rho, m_x, m_y, &
               m_z, energy, fluid%gs_Xh, gamma_ref, &
               reg%entropy_wave_speed)
       class default
          call neko_error('Burgers requires entropy viscosity regularization')
       end select
       minimum_limiter = min(minimum_limiter, &
            fluid%euler_idp_diagnostics%min_limiter)
       maximum_limited_fraction = max(maximum_limited_fraction, &
            fluid%euler_idp_diagnostics%limited_edge_fraction)
       maximum_graph_cfl = max(maximum_graph_cfl, &
            fluid%euler_idp_diagnostics%max_graph_cfl)
    class default
       call neko_error('Burgers proxy requires compressible Euler')
    end select
  end subroutine enforce_burgers_flux

  !> Evaluate the conservation residual and the exact Burgers wave speed.
  subroutine user_entropy_pair(entropy, flux_x, flux_y, flux_z, &
       wave_speed, time)
    type(field_t), intent(inout) :: entropy
    type(field_t), intent(inout) :: flux_x, flux_y, flux_z
    type(field_t), intent(inout) :: wave_speed
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: rho
    real(kind=rp) :: scalar_value
    integer :: i

    rho => neko_registry%get_field('fluid_rho')
    do i = 1, rho%size()
       scalar_value = rho%x(i, 1, 1, 1) - density_shift
       entropy%x(i, 1, 1, 1) = scalar_value
       flux_x%x(i, 1, 1, 1) = 0.5_rp * scalar_value**2
       flux_y%x(i, 1, 1, 1) = 0.5_rp * scalar_value**2
       flux_z%x(i, 1, 1, 1) = 0.0_rp
       wave_speed%x(i, 1, 1, 1) = sqrt(2.0_rp) * abs(scalar_value)
    end do
  end subroutine user_entropy_pair

  !> Report shifted-scalar bounds and mass drift, then write the z=0 plane.
  subroutine finalize(time)
    type(time_state_t), intent(in) :: time
    type(coef_t), pointer :: coef
    type(field_t), pointer :: rho
    real(kind=rp) :: current_mass, lower_violation, upper_violation
    real(kind=rp) :: minimum_scalar, maximum_scalar, x, y
    character(len=80) :: filename
    integer :: e, i, j, k, n, output_unit

    coef => neko_user_access%case%fluid%c_Xh
    rho => neko_registry%get_field('fluid_rho')
    n = rho%size()
    current_mass = glsc2(rho%x, coef%B, n)
    minimum_scalar = glmin(rho%x, n) - density_shift
    maximum_scalar = glmax(rho%x, n) - density_shift
    lower_violation = max(0.0_rp, scalar_minimum - minimum_scalar)
    upper_violation = max(0.0_rp, maximum_scalar - scalar_maximum)

    if (pe_rank .eq. 0) then
       write(*, '(A,I0,9(1X,ES25.16E3))') 'EULER_IDP_BURGERS ', &
            time%tstep, time%t, minimum_scalar, maximum_scalar, &
            abs(current_mass - initial_mass), lower_violation, &
            upper_violation, minimum_limiter, &
            maximum_limited_fraction, maximum_graph_cfl
    end if

    write(filename, '(A,I0,A)') 'burgers_scalar_rank_', pe_rank, '.csv'
    open(newunit = output_unit, file = trim(filename), status = 'replace', &
         action = 'write')
    write(output_unit, '(A)') 'x,y,u,u_initial'
    do e = 1, rho%dof%msh%nelv
       do k = 1, rho%dof%Xh%lz
          if (abs(rho%dof%z(1, 1, k, e)) .gt. 1.0e-12_rp) cycle
          do j = 1, rho%dof%Xh%ly
             do i = 1, rho%dof%Xh%lx
                x = rho%dof%x(i, j, k, e)
                y = rho%dof%y(i, j, k, e)
                write(output_unit, '(4(ES25.16E3,:,","))') x, y, &
                     rho%x(i, j, k, e) - density_shift, &
                     burgers_initial_scalar(x, y)
             end do
          end do
       end do
    end do
    close(output_unit)
  end subroutine finalize

  !> Initial scalar obtained from the exact solution at t = 0.
  pure real(kind=rp) function burgers_initial_scalar(x, y) result(value)
    real(kind=rp), intent(in) :: x, y

    if (abs(x - 1.0_rp) .le. 0.5_rp .and. &
         abs(y - 1.0_rp) .le. 0.5_rp) then
       value = scalar_maximum
    else
       value = scalar_minimum
    end if
  end function burgers_initial_scalar

end module user
