! KPP nonlinear scalar transport proxy for the Euler IDP solver.
!
! Density represents the scalar u. After every Forward Euler step, momentum
! is reset to (sin(rho), cos(rho), 0), so the next density residual has the
! KPP flux. Energy is reset to maintain a constant pressure. This is a test
! harness for the density transport and limiter, not an unforced Euler flow.
module user
  use neko
  use fluid_scheme_compressible_ns, only : fluid_scheme_compressible_ns_t
  implicit none

  real(kind=rp), parameter :: gamma_ref = 1.4_rp
  real(kind=rp), parameter :: pressure_ref = 1.0_rp
  real(kind=rp) :: initial_mass = 0.0_rp
  real(kind=rp) :: minimum_limiter = 1.0_rp
  real(kind=rp) :: maximum_limited_fraction = 0.0_rp
  real(kind=rp) :: maximum_graph_cfl = 0.0_rp

contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    user%initial_conditions => initial_conditions
    user%initialize => initialize
    user%compute => enforce_kpp_flux
    user%entropy_pair => user_entropy_pair
    user%finalize => finalize
  end subroutine user_setup

  !> Initialize the KPP scalar as density and set its nonlinear flux.
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(field_t), pointer :: rho, u, v, w, p
    real(kind=rp) :: rho_value, x, y
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
       rho_value = kpp_initial_density(x, y)
       rho%x(i, 1, 1, 1) = rho_value
       u%x(i, 1, 1, 1) = sin(rho_value) / rho_value
       v%x(i, 1, 1, 1) = cos(rho_value) / rho_value
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

  !> Reset momentum and energy so the next density flux is the KPP flux.
  subroutine enforce_kpp_flux(time)
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: rho, m_x, m_y, m_z, energy
    real(kind=rp) :: kinetic_energy, rho_value
    real(kind=rp) :: momentum_x, momentum_y
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
             call neko_error('KPP density must remain positive')
          end if

          momentum_x = sin(rho_value)
          momentum_y = cos(rho_value)
          kinetic_energy = 0.5_rp * &
               (momentum_x**2 + momentum_y**2) / rho_value

          m_x%x(i, 1, 1, 1) = momentum_x
          m_y%x(i, 1, 1, 1) = momentum_y
          m_z%x(i, 1, 1, 1) = 0.0_rp
          energy%x(i, 1, 1, 1) = pressure_ref / (gamma_ref - 1.0_rp) + &
               kinetic_energy
          fluid%u%x(i, 1, 1, 1) = momentum_x / rho_value
          fluid%v%x(i, 1, 1, 1) = momentum_y / rho_value
          fluid%w%x(i, 1, 1, 1) = 0.0_rp
          fluid%p%x(i, 1, 1, 1) = pressure_ref
          fluid%temperature%x(i, 1, 1, 1) = pressure_ref / &
               (rho_value * (gamma_ref - 1.0_rp))
          fluid%S%x(i, 1, 1, 1) = 0.5_rp * rho_value**2
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
       call neko_error('KPP requires compressible Euler')
    end select
  end subroutine enforce_kpp_flux

  !> Evaluate the quadratic entropy pair for the KPP scalar flux.
  subroutine user_entropy_pair(entropy, flux_x, flux_y, flux_z, &
       wave_speed, time)
    type(field_t), intent(inout) :: entropy
    type(field_t), intent(inout) :: flux_x, flux_y, flux_z
    type(field_t), intent(inout) :: wave_speed
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: rho
    real(kind=rp) :: rho_value
    integer :: i

    rho => neko_registry%get_field('fluid_rho')
    do i = 1, rho%size()
       rho_value = rho%x(i, 1, 1, 1)
       entropy%x(i, 1, 1, 1) = 0.5_rp * rho_value**2
       flux_x%x(i, 1, 1, 1) = &
            rho_value * sin(rho_value) + cos(rho_value)
       flux_y%x(i, 1, 1, 1) = &
            rho_value * cos(rho_value) - sin(rho_value)
       flux_z%x(i, 1, 1, 1) = 0.0_rp
       wave_speed%x(i, 1, 1, 1) = 1.0_rp
    end do
  end subroutine user_entropy_pair

  !> Report scalar bounds and mass drift, then write the z=0 density plane.
  subroutine finalize(time)
    type(time_state_t), intent(in) :: time
    type(coef_t), pointer :: coef
    type(field_t), pointer :: rho
    real(kind=rp) :: current_mass, lower_violation, upper_violation
    real(kind=rp) :: minimum_density, maximum_density, x, y
    real(kind=rp) :: rho_minimum, rho_maximum
    character(len=80) :: filename
    integer :: e, i, j, k, n, output_unit

    coef => neko_user_access%case%fluid%c_Xh
    rho => neko_registry%get_field('fluid_rho')
    n = rho%size()
    rho_minimum = 0.25_rp * pi
    rho_maximum = 3.5_rp * pi
    current_mass = glsc2(rho%x, coef%B, n)
    minimum_density = glmin(rho%x, n)
    maximum_density = glmax(rho%x, n)
    lower_violation = max(0.0_rp, rho_minimum - minimum_density)
    upper_violation = max(0.0_rp, maximum_density - rho_maximum)

    if (pe_rank .eq. 0) then
       write(*, '(A,I0,9(1X,ES25.16E3))') 'EULER_IDP_KPP ', &
            time%tstep, time%t, minimum_density, maximum_density, &
            abs(current_mass - initial_mass), lower_violation, &
            upper_violation, minimum_limiter, &
            maximum_limited_fraction, maximum_graph_cfl
    end if

    write(filename, '(A,I0,A)') 'kpp_density_rank_', pe_rank, '.csv'
    open(newunit = output_unit, file = trim(filename), status = 'replace', &
         action = 'write')
    write(output_unit, '(A)') 'x,y,rho,rho_initial'
    do e = 1, rho%dof%msh%nelv
       do k = 1, rho%dof%Xh%lz
          if (abs(rho%dof%z(1, 1, k, e)) .gt. 1.0e-12_rp) cycle
          do j = 1, rho%dof%Xh%ly
             do i = 1, rho%dof%Xh%lx
                x = rho%dof%x(i, j, k, e)
                y = rho%dof%y(i, j, k, e)
                write(output_unit, '(4(ES25.16E3,:,","))') x, y, &
                     rho%x(i, j, k, e), kpp_initial_density(x, y)
             end do
          end do
       end do
    end do
    close(output_unit)
  end subroutine finalize

  !> Initial density for the KPP rotating-wave benchmark.
  pure real(kind=rp) function kpp_initial_density(x, y) result(value)
    real(kind=rp), intent(in) :: x, y

    if (x**2 + y**2 .lt. 1.0_rp) then
       value = 3.5_rp * pi
    else
       value = 0.25_rp * pi
    end if
  end function kpp_initial_density

end module user
