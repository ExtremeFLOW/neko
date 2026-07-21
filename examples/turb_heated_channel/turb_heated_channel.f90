! ===========================================================================
!  Turbulent heated channel -- low-Mach, Re_tau ~ 180 (Re_b = 2800), nondim.
! ---------------------------------------------------------------------------
!  Geometry (delta = half-height = 1):
!     x in [0, 2*pi]  streamwise, PERIODIC
!     y in [0, 2]     wall-normal, NO-SLIP walls (wall-clustered mesh, ydist.csv)
!     z in [0, pi]    spanwise,   PERIODIC
!  Heating (this is the "more heating" vs the 2:1 turb_channel reference):
!     hot wall  y=0 (zone 3)  T = 4
!     cold wall y=2 (zone 4)  T = 1
!     ideal-gas EOS rho = P0/(R T) = 1/T  =>  rho in [0.25, 1.0]  (4x density ratio)
!  Transport: constant molecular mu = 1/Re_b, lambda = mu/Pr (Pr = 0.71). Density
!     varies through the EOS; that is the dominant low-Mach effect here. (To make
!     mu(T)/k(T) temperature-dependent, add a material_properties hook -- the
!     laminar mu_tot/lambda_tot refresh already propagates it.)
!  Drive: constant volumetric flow rate (flow_rate_force in the .case), bulk U_b=1.
!  IC: Reichardt log-law mean + divergence-organised + random perturbations to trip
!     turbulence; linear temperature across y. No buoyancy (clean forced test).
! ===========================================================================
module user
  use neko
  implicit none

  real(kind=rp), parameter :: Re_tau = 180.0_rp
  real(kind=rp), parameter :: Re_b   = 2800.0_rp
  real(kind=rp), parameter :: Pr     = 0.71_rp
  real(kind=rp), parameter :: kk     = 0.41_rp     ! von Karman constant
  real(kind=rp), parameter :: Cplus  = 5.17_rp     ! log-law intercept
  real(kind=rp), parameter :: T_hot  = 4.0_rp      ! hot wall  (y = 0, zone 3)
  real(kind=rp), parameter :: T_cold = 1.0_rp      ! cold wall (y = 2, zone 4)
  real(kind=rp), parameter :: llx = 2.0_rp * pi
  real(kind=rp), parameter :: llz = pi

  !> Buoyancy (optional): reduced-gravity magnitude, read from the case key
  !! "case.gravity_nd" (default 0 = off, base case unchanged). Applied per unit
  !! mass as (rho_ref/rho - 1)*g_nd in +y, i.e. gravity points -y with the
  !! hydrostatic reference subtracted: hot fluid (rho < rho_ref) at the y=0 hot
  !! wall is pushed UP -> unstable (Rayleigh-Benard-Poiseuille) orientation.
  !! Needs "source_terms": [{"type": "user"}] in case.fluid to be active.
  real(kind=rp) :: g_nd = 0.0_rp
  real(kind=rp), parameter :: rho_ref = 1.0_rp             ! cold-wall density

contains

  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%startup => startup
    u%initial_conditions => ics
    u%source_term => fluid_source
  end subroutine user_setup

  !> Constant (reference) molecular properties; the fields are filled from these.
  subroutine startup(params)
    type(json_file), intent(inout) :: params
    logical :: found
    character(len=256) :: log_buf
    call params%add("case.fluid.mu", 1.0_rp / Re_b)         ! mu = 1/Re_b
    call params%add("case.fluid.rho", 1.0_rp)               ! reference rho (EOS overrides)
    call params%add("case.scalar.lambda", (1.0_rp / Re_b) / Pr) ! lambda = mu/Pr
    call params%add("case.scalar.cp", 1.0_rp)
    call params%get("case.gravity_nd", g_nd, found)
    if (.not. found) g_nd = 0.0_rp
    write(log_buf, '(A,ES11.3)') '   [channel] gravity_nd = ', g_nd
    call neko_log%message(log_buf)
  end subroutine startup

  !> Buoyancy body force (reduced gravity, -y sense; see g_nd above).
  subroutine fluid_source(scheme_name, rhs, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: rhs
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: ru, rv, rw, rho
    integer :: i, n

    if (scheme_name .ne. 'fluid') return
    ru => rhs%get_by_index(1); rv => rhs%get_by_index(2); rw => rhs%get_by_index(3)
    rho => neko_registry%get_field('fluid_rho')
    n = ru%dof%size()
    do i = 1, n
       ru%x(i,1,1,1) = 0.0_rp
       rv%x(i,1,1,1) = (rho_ref / rho%x(i,1,1,1) - 1.0_rp) * g_nd
       rw%x(i,1,1,1) = 0.0_rp
    end do
  end subroutine fluid_source

  !> Reichardt mean profile + large/small-scale + random perturbations.
  function chan_ic(x, y, z) result(uvw)
    real(kind=rp), intent(in) :: x, y, z
    real(kind=rp) :: uvw(3), d, yp, ux, eps, kx, kz, alpha, beta, ran

    d  = min(y, 2.0_rp - y)                  ! distance to nearest wall (delta = 1)
    yp = d * Re_tau
    ux = (1.0_rp/kk)*log(1.0_rp + kk*yp) + (Cplus - (1.0_rp/kk)*log(kk)) * &
         (1.0_rp - exp(-yp/11.0_rp) - yp/11.0_rp*exp(-yp/3.0_rp))
    ux = ux * Re_tau / Re_b                  ! scale to bulk units (U_b ~ 1)
    uvw(1) = ux; uvw(2) = 0.0_rp; uvw(3) = 0.0_rp

    ! large-scale divergence-organised perturbation
    eps = 0.05_rp; kx = 3.0_rp; kz = 4.0_rp
    alpha = kx*2.0_rp*pi/llx; beta = kz*2.0_rp*pi/llz
    uvw(1) = uvw(1) + eps*beta*sin(alpha*x)*cos(beta*z)
    uvw(2) = uvw(2) + eps*sin(alpha*x)*sin(beta*z)
    uvw(3) = uvw(3) - eps*alpha*cos(alpha*x)*sin(beta*z)
    ! small-scale perturbation
    eps = 0.005_rp; kx = 17.0_rp; kz = 13.0_rp
    alpha = kx*2.0_rp*pi/llx; beta = kz*2.0_rp*pi/llz
    uvw(1) = uvw(1) + eps*beta*sin(alpha*x)*cos(beta*z)
    uvw(2) = uvw(2) + eps*sin(alpha*x)*sin(beta*z)
    uvw(3) = uvw(3) - eps*alpha*cos(alpha*x)*sin(beta*z)
    ! random kick in y
    ran = sin(-20.0_rp*x*z + y**3*tan(x*z**2) + 100.0_rp*z*y &
         - 20.0_rp*sin(x*y*z)**5)
    uvw(2) = uvw(2) + 0.001_rp*ran
  end function chan_ic

  subroutine ics(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(field_t), pointer :: u, v, w, s
    integer :: i, n
    real(kind=rp) :: uvw(3), x, y, z

    if (scheme_name .eq. 'fluid') then
       u => fields%get('u'); v => fields%get('v'); w => fields%get('w')
       n = u%dof%size()
       do i = 1, n
          x = u%dof%x(i,1,1,1); y = u%dof%y(i,1,1,1); z = u%dof%z(i,1,1,1)
          uvw = chan_ic(x, y, z)
          u%x(i,1,1,1) = uvw(1); v%x(i,1,1,1) = uvw(2); w%x(i,1,1,1) = uvw(3)
       end do
    else if (scheme_name .eq. 'temperature') then
       s => fields%items(1)%ptr
       n = s%dof%size()
       do i = 1, n
          y = s%dof%y(i,1,1,1)
          s%x(i,1,1,1) = T_hot - (T_hot - T_cold) * (y / 2.0_rp)   ! linear hot->cold
       end do
    end if
  end subroutine ics

end module user
