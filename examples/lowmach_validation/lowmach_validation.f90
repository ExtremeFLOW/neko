! Neko replication of the Nek5000 low-Mach validation case (short_tests/lowMach_test),
! the Tomboulides et al. (JCP 146, 1998) linearized manufactured solution.
!
!   Domain x in [-1,1], y/z periodic (effectively 1-D in x). delta = 0.2.
!   EOS:        rho = P0/(R T) = 1/T            (low_mach, P0=R=cp=1)
!   Properties: mu = 1 (UDIFF), lambda = 1, cp = 1 (constants)
!   Heat source (manufactured): q(x) = sech^2(x/d)/d * (0.5 + tanh(x/d)/d)
!   Analytical steady solution:  u = T = 0.5*(3 + tanh(x/d)),   v = w = 0
!                                QTL = div(u) = 0.5/d * (1 - tanh^2(x/d))
!
!   The IC equals the analytical steady state; the x=+/-1 Dirichlet BCs hold it.
!   The user%compute hook reports MAX errors of u, T and Q_T against the closed
!   forms (mirroring Nek5000 USERCHK). Q_T is the term the new q-in-Q_T fix
!   affects: without the source it is wrong by O(5) near x=0; with it, it should
!   match QTL to solver tolerance, as in Nek5000.
module user
  use neko
  implicit none

  real(kind=rp), parameter :: delta = 0.2_rp

  !> IC selector: .false. = analytic tanh steady state (default, Nek5000-style
  !! "hold" test); .true. = linear profile between the same x = +/-1 boundary
  !! values, so the run must CONVERGE to the steady state (paper-style test).
  !! Set via the optional case-file key "case.linear_ic".
  logical :: linear_ic = .false.

contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%startup => startup
    user%material_properties => material_properties
    user%initial_conditions => initial_conditions
    user%dirichlet_conditions => bc_tanh
    user%source_term => heat_source
    user%compute => errchk
  end subroutine user_setup

  ! --- analytical profiles (Tomboulides manufactured solution) ---
  pure function fprof(x) result(f)        ! u = T = 0.5(3 + tanh(x/d))
    real(kind=rp), intent(in) :: x
    real(kind=rp) :: f
    f = 0.5_rp * (3.0_rp + tanh(x / delta))
  end function fprof

  pure function flin(x) result(f)         ! linear IC with the same end values
    real(kind=rp), intent(in) :: x
    real(kind=rp) :: f
    f = 1.5_rp + 0.5_rp * tanh(1.0_rp / delta) * x
  end function flin

  !> Profile used for the initial condition (BCs always use fprof).
  pure function icprof(x) result(f)
    real(kind=rp), intent(in) :: x
    real(kind=rp) :: f
    if (linear_ic) then
       f = flin(x)
    else
       f = fprof(x)
    end if
  end function icprof

  pure function qtl_ex(x) result(q)       ! QTL = 0.5/d (1 - tanh^2)
    real(kind=rp), intent(in) :: x
    real(kind=rp) :: q, th
    th = tanh(x / delta)
    q = 0.5_rp / delta * (1.0_rp - th * th)
  end function qtl_ex

  pure function qvol(x) result(q)         ! manufactured heat source
    real(kind=rp), intent(in) :: x
    real(kind=rp) :: q, xd, sech2
    xd = x / delta
    sech2 = 1.0_rp / (cosh(xd) * cosh(xd))
    q = sech2 / delta * (0.5_rp + tanh(xd) / delta)
  end function qvol

  subroutine startup(params)
    type(json_file), intent(inout) :: params
    logical :: found
    character(len=256) :: log_buf
    call params%add("case.fluid.mu", 1.0_rp)
    call params%add("case.fluid.rho", 1.0_rp)
    call params%get("case.linear_ic", linear_ic, found)
    if (.not. found) linear_ic = .false.
    write(log_buf, '(A,L1)') '   [validate] linear_ic = ', linear_ic
    call neko_log%message(log_buf)
  end subroutine startup

  !> Constant properties: mu=1 (fluid), lambda=1, cp=1 (temperature). Setting
  !! them here each step guarantees lambda_tot/mu_tot are populated.
  subroutine material_properties(scheme_name, properties, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: properties
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: mu, lambda, cp
    integer :: i, n

    if (scheme_name .eq. 'fluid') then
       mu => properties%get_by_index(2)
       n = mu%dof%size()
       do i = 1, n
          mu%x(i,1,1,1) = 1.0_rp
       end do
    else if (scheme_name .eq. 'temperature') then
       cp => properties%get_by_index(1)
       lambda => properties%get_by_index(2)
       n = cp%dof%size()
       do i = 1, n
          cp%x(i,1,1,1) = 1.0_rp
          lambda%x(i,1,1,1) = 1.0_rp
       end do
    end if
  end subroutine material_properties

  !> IC: u = T = analytical tanh (or linear profile when case.linear_ic).
  !! Also create + fill the registry field "temperature_qdot" with q(x) so
  !! the low-Mach Q_T routine picks it up.
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(field_t), pointer :: u, v, w, s, qd
    integer :: i, n
    real(kind=rp) :: xx

    if (scheme_name .eq. 'fluid') then
       u => fields%items(1)%ptr
       v => fields%items(2)%ptr
       w => fields%items(3)%ptr
       n = u%dof%size()
       do i = 1, n
          xx = u%dof%x(i,1,1,1)
          u%x(i,1,1,1) = icprof(xx)
          v%x(i,1,1,1) = 0.0_rp
          w%x(i,1,1,1) = 0.0_rp
       end do
    else if (scheme_name .eq. 'temperature') then
       s => fields%items(1)%ptr
       n = s%dof%size()
       do i = 1, n
          s%x(i,1,1,1) = icprof(s%dof%x(i,1,1,1))
       end do
       ! Create + fill the heat-source field read by lowmach_update_Q_T.
       if (.not. neko_registry%field_exists('temperature_qdot')) &
            call neko_registry%add_field(s%dof, 'temperature_qdot')
       qd => neko_registry%get_field('temperature_qdot')
       do i = 1, n
          qd%x(i,1,1,1) = qvol(s%dof%x(i,1,1,1))
       end do
    end if
  end subroutine initial_conditions

  !> x = +/-1 Dirichlet: velocity (3-field list) -> u = tanh, v=w=0;
  !! temperature (1-field list) -> T = tanh. Branches on the field count.
  subroutine bc_tanh(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time
    type(dofmap_t), pointer :: dof
    integer :: i, idx
    real(kind=rp) :: xx, val

    dof => fields%dof(1)
    do i = 1, bc%msk(0)
       idx = bc%msk(i)
       xx = dof%x(idx,1,1,1)
       val = fprof(xx)
       if (fields%size() .eq. 3) then       ! velocity vector BC
          fields%items(1)%ptr%x(idx,1,1,1) = val
          fields%items(2)%ptr%x(idx,1,1,1) = 0.0_rp
          fields%items(3)%ptr%x(idx,1,1,1) = 0.0_rp
       else                                  ! scalar (temperature) BC
          fields%items(1)%ptr%x(idx,1,1,1) = val
       end if
    end do
  end subroutine bc_tanh

  !> Energy-equation volumetric heat source q on the temperature scalar.
  subroutine heat_source(scheme_name, rhs, time)
    character(len=*), intent(in) :: scheme_name
    type(time_state_t), intent(in) :: time
    type(field_list_t), intent(inout) :: rhs
    type(field_t), pointer :: f, qd, qT
    integer :: i, n

    if (scheme_name .ne. 'temperature') return
    qd => neko_registry%get_field('temperature_qdot')
    f => rhs%get_by_index(1)
    n = f%dof%size()
    ! Raw volumetric source q (the core scalar_pnpn low_mach path now keeps the
    ! source OUT of the rho*cp(x) weighting, matching Nek5000 makeuq).
    do i = 1, n
       f%x(i,1,1,1) = qd%x(i,1,1,1)
    end do
  end subroutine heat_source

  !> MAX error of u, T, Q_T against the analytical solution (mirrors USERCHK).
  subroutine errchk(time)
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, s, qt, rho, cp
    integer :: i, n
    real(kind=rp) :: eu, et, eq, red(1), xx
    character(len=256) :: log_buf

    u  => neko_registry%get_field('u')
    s  => neko_registry%get_field('temperature')
    qt => neko_registry%get_field('Q_T')
    n = u%dof%size()
    eu = 0.0_rp; et = 0.0_rp; eq = 0.0_rp
    do i = 1, n
       xx = u%dof%x(i,1,1,1)
       eu = max(eu, abs(u%x(i,1,1,1) - fprof(xx)))
       et = max(et, abs(s%x(i,1,1,1) - fprof(xx)))
       eq = max(eq, abs(qt%x(i,1,1,1) - qtl_ex(xx)))
    end do
    red(1) = eu; eu = glmax(red, 1)
    red(1) = et; et = glmax(red, 1)
    red(1) = eq; eq = glmax(red, 1)
    write(log_buf, '(A,I6,A,ES11.3,A,ES11.3,A,ES11.3)') &
         '   [validate] step', time%tstep, '  MAX err  u=', eu, &
         '  T=', et, '  Q_T=', eq
    call neko_log%message(log_buf)

    ! Diagnostic: ranges of the denominator components (rho, cp, T) and Q_T.
    ! Analytical: T in [1,2], rho=1/T in [0.5,1], cp=1, Q_T in [0, 0.5/delta=2.5].
    ! EOS consistency: max|rho*T - 1| measures the one-step rho/T lag (rho built
    ! from T^n vs current T^{n+1}); with heat-first ordering it should sit at
    ! roundoff during the transient instead of O(dt*dT/dt).
    rho => neko_registry%get_field('fluid_rho')
    eu = 0.0_rp
    do i = 1, n
       eu = max(eu, abs(rho%x(i,1,1,1) * s%x(i,1,1,1) - 1.0_rp))
    end do
    red(1) = eu; eu = glmax(red, 1)
    write(log_buf, '(A,ES11.3)') '   [diag] max|rho*T - 1| = ', eu
    call neko_log%message(log_buf)
    write(log_buf, '(A,2ES11.3,A,2ES11.3,A,2ES11.3)') &
         '   [diag] T[', glmin(s%x,n), glmax(s%x,n), '] rho[', &
         glmin(rho%x,n), glmax(rho%x,n), '] Q_T[', glmin(qt%x,n), glmax(qt%x,n)
    call neko_log%message(log_buf)
    write(log_buf, '(A,2ES11.3,A,ES11.3)') &
         '   [diag] u[', glmin(u%x,n), glmax(u%x,n), &
         '] (analytic u in [1,2])'
    call neko_log%message(log_buf)
    if (neko_registry%field_exists('temperature_cp')) then
       cp => neko_registry%get_field('temperature_cp')
       write(log_buf, '(A,2ES11.3)') '   [diag] cp[', &
            glmin(cp%x,n), glmax(cp%x,n)
       call neko_log%message(log_buf)
    end if
  end subroutine errchk

end module user
