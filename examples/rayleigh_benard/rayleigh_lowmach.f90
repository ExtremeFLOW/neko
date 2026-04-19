! Rayleigh-Bénard with the low-Mach solver.
!
! This is an adaptation of rayleigh.f90. Two physical changes:
!   * Temperature is shifted to [1, 2] (cold top, hot bottom) so that the
!     EOS rho = P0 / (R * T) used by fluid_lowmach_pnpn stays finite.
!     The driving temperature difference DeltaT = 1 is unchanged.
!   * The Boussinesq buoyancy source is evaluated against the cold-wall
!     reference T_ref = 1, so rhs_w = Ra * Pr * (T - T_ref). This
!     preserves the original zero-buoyancy state at the cold wall.
!
! The low-Mach scheme computes rho and Q_T from this temperature field
! each step in addition to the Boussinesq forcing here.
module user
  use neko
  implicit none

  real(kind=rp) :: Ra
  real(kind=rp) :: Pr
  real(kind=rp) :: ta2 = 0
  real(kind=rp), parameter :: T_ref = 1.0_rp

contains
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initial_conditions => initial_conditions
    user%source_term => forcing_source_term
    user%startup => startup
  end subroutine user_setup

  subroutine startup(params)
    type(json_file), intent(inout) :: params
    real(kind=rp) :: rho, mu, cp, lambda, Re

    call json_get(params, "case.fluid.Ra", Ra)
    call json_get(params, "case.scalar.Pr", Pr)

    Re = 1.0_rp / Pr
    mu = 1.0_rp / Re
    lambda = mu / Pr
    rho = 1.0_rp
    cp = 1.0_rp

    call params%add("case.fluid.mu", mu)
    call params%add("case.fluid.rho", rho)
    call params%add("case.scalar.lambda", lambda)
    call params%add("case.scalar.cp", cp)
  end subroutine startup

  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    integer :: i, e, k, j
    real(kind=rp) :: rand, z
    type(field_t), pointer :: s

    if (scheme_name .ne. 'temperature') return

    s => fields%items(1)%ptr

    ! Linear profile from T=2 at z=0 (bottom) to T=1 at z=1 (top).
    do i = 1, s%dof%size()
       s%x(i,1,1,1) = 2.0_rp - s%dof%z(i,1,1,1)
    end do

    ! Small perturbation on interior nodes to trigger the instability.
    do e = 1, s%msh%nelv
       do k = 2, s%Xh%lx-1
          do j = 2, s%Xh%lx-1
             do i = 2, s%Xh%lx-1
                rand = cos(real(e + s%msh%offset_el, rp) * real(i*j*k, rp))
                z = s%dof%z(i,j,k,e)
                s%x(i,j,k,e) = 2.0_rp - z + 0.0001_rp*rand* &
                     sin(4.0_rp*pi/4.5_rp * s%dof%x(i,j,k,e)) &
                     * sin(4.0_rp*pi/4.5_rp * s%dof%y(i,j,k,e))
             end do
          end do
       end do
    end do

  end subroutine initial_conditions

  subroutine forcing_source_term(scheme_name, rhs, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: rhs
    type(time_state_t), intent(in) :: time

    integer :: i, n
    type(field_t), pointer :: u, v, w, s, rhs_u, rhs_v, rhs_w
    real(kind=rp) :: rapr, ta2pr

    if (scheme_name .eq. 'fluid') then
       u => neko_registry%get_field('u')
       v => neko_registry%get_field('v')
       w => neko_registry%get_field('w')
       s => neko_registry%get_field('temperature')

       rhs_u => rhs%get_by_index(1)
       rhs_v => rhs%get_by_index(2)
       rhs_w => rhs%get_by_index(3)

       rapr = Ra*Pr
       ta2pr = ta2*Pr
       n = rhs_w%dof%size()

       call field_cmult2(rhs_u, v, Ta2Pr)
       call field_cmult2(rhs_v, u, Ta2Pr)

       ! Buoyancy relative to the cold-wall reference temperature.
       do i = 1, n
          rhs_w%x(i,1,1,1) = rapr * (s%x(i,1,1,1) - T_ref)
       end do
    end if
  end subroutine forcing_source_term
end module user
