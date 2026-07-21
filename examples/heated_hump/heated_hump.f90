! =====================================================================
!  Flow over a wall-mounted hump, low-Mach, with:
!   * STEADY inlet: a fixed PARABOLIC (fully-developed laminar) profile, zero at
!     the walls, bulk = U_in exactly, smooth near-wall layer (no thin-shear
!     corner). NO recycling. LOW Reynolds number (~500) -> laminar & stable;
!     raise Re later, once the case runs cleanly.
!   * HEATED lower wall for x > x_heat (T ramps T_cold -> T_hot=2), cold
!     elsewhere; ideal-gas rho = P0/(R T) = 1/T.
!   * GRAVITY / buoyancy (reduced-gravity body force in -y).
!
!  Incremental switches (do_heat / do_gravity / do_trip) let us isolate.
! =====================================================================
module user
  use neko
  implicit none

  ! ---- geometry / flow parameters (nondimensional) ----
  real(kind=rp), parameter :: U_in   = 1.0_rp
  real(kind=rp), parameter :: Lx     = 25.0_rp
  real(kind=rp), parameter :: Ly     = 3.0_rp
  real(kind=rp), parameter :: Lz     = 4.0_rp
  real(kind=rp), parameter :: Hhump  = 1.0_rp          ! hump height = length scale
  real(kind=rp), parameter :: Re_in  = 500.0_rp     ! LOW Re -> laminar, stable (raise later once running)
  real(kind=rp), parameter :: Pr     = 0.71_rp
  real(kind=rp), parameter :: rho_ref = 1.0_rp         ! cold reference density
  real(kind=rp), parameter :: mu_ref = rho_ref*U_in*Hhump/Re_in
  real(kind=rp), parameter :: k_ref  = mu_ref/Pr

  ! ---- thermal / gravity ----
  real(kind=rp), parameter :: T_cold = 1.0_rp
  real(kind=rp), parameter :: T_hot  = 2.0_rp          ! "heating a bit, ~2"
  real(kind=rp), parameter :: x_heat = 11.0_rp         ! heated wall starts (hump ends x=8)
  real(kind=rp), parameter :: dx_ramp = 1.0_rp         ! smoothstep ramp width
  real(kind=rp), parameter :: g_nd   = 0.15_rp         ! buoyancy magnitude (reduced for stability)

  ! ---- near-inlet trip (seeds transition for the steady inlet) ----
  real(kind=rp), parameter :: trip_x0 = 0.5_rp, trip_x1 = 1.0_rp
  real(kind=rp), parameter :: trip_yh = 0.4_rp         ! within this height of the lower wall
  real(kind=rp), parameter :: trip_A  = 0.2_rp         ! trip amplitude

  ! ---- incremental switches ----
  logical, parameter :: do_heat    = .true.
  logical, parameter :: do_gravity = .false.   ! OFF to isolate: buoyancy (Rayleigh-Benard off the heated wall) is the prime late-blow-up suspect
  logical, parameter :: do_trip    = .false.   ! OFF: the hump trips the flow; the trip was destabilizing the laminar inlet

contains

  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%startup => startup
    u%material_properties => material_properties
    u%initial_conditions => ics
    u%dirichlet_conditions => inflow_and_wall
    u%source_term => fluid_source
  end subroutine user_setup

  subroutine startup(params)
    type(json_file), intent(inout) :: params
    call params%add("case.fluid.mu", mu_ref)
    call params%add("case.fluid.rho", rho_ref)
    call params%add("case.scalar.lambda", k_ref)
    call params%add("case.scalar.cp", 1.0_rp)
  end subroutine startup

  !> Steady PARABOLIC inlet profile (fully-developed laminar): 0 at the walls,
  !! smooth, peak 1.5 at the centre, bulk = U_in EXACTLY (no thin-shear corner,
  !! no bulk-flux over-supply -- both fix the prior blunt-profile blow-up).
  pure function uin_profile(y) result(uu)
    real(kind=rp), intent(in) :: y
    real(kind=rp) :: uu, yn
    yn = min(max(y/Ly, 0.0_rp), 1.0_rp)
    uu = U_in * 1.5_rp * (1.0_rp - (2.0_rp*yn - 1.0_rp)**2)
  end function uin_profile

  !> Constant transport properties (robust first build): mu, lambda, cp.
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
          mu%x(i,1,1,1) = mu_ref
       end do
    else if (scheme_name .eq. 'temperature') then
       cp => properties%get_by_index(1)
       lambda => properties%get_by_index(2)
       n = cp%dof%size()
       do i = 1, n
          cp%x(i,1,1,1) = 1.0_rp
          lambda%x(i,1,1,1) = k_ref
       end do
    end if
  end subroutine material_properties

  !> IC: fluid -> the steady inlet mean profile + GENTLE perturbations (seed
  !! transition without a violent impulsive start); temperature -> uniform cold
  !! (the wall heating builds the hot layer downstream).
  subroutine ics(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(field_t), pointer :: u, v, w, s
    integer :: i, n
    real(kind=rp) :: xx, yy, zz, um, gy, r1, r2, r3

    if (scheme_name .eq. 'fluid') then
       u => fields%items(1)%ptr; v => fields%items(2)%ptr; w => fields%items(3)%ptr
       n = u%dof%size()
       do i = 1, n
          xx = u%dof%x(i,1,1,1); yy = u%dof%y(i,1,1,1); zz = u%dof%z(i,1,1,1)
          um = uin_profile(yy)
          ! near-wall localisation of the (gentle) perturbations
          gy = exp(-((yy - 0.3_rp)/0.3_rp)**2) + exp(-((yy - (Ly-0.3_rp))/0.3_rp)**2)
          call random_number(r1); call random_number(r2); call random_number(r3)
          u%x(i,1,1,1) = um &
               + 0.03_rp*gy*sin(3.0_rp*2.0_rp*pi*xx/Lx)*cos(4.0_rp*2.0_rp*pi*zz/Lz) &
               + 0.01_rp*(r1 - 0.5_rp)
          v%x(i,1,1,1) = 0.02_rp*gy*cos(3.0_rp*2.0_rp*pi*xx/Lx)*sin(4.0_rp*2.0_rp*pi*zz/Lz) &
               + 0.01_rp*(r2 - 0.5_rp)
          w%x(i,1,1,1) = 0.02_rp*gy*sin(3.0_rp*2.0_rp*pi*xx/Lx)*cos(4.0_rp*2.0_rp*pi*zz/Lz) &
               + 0.01_rp*(r3 - 0.5_rp)
       end do
    else if (scheme_name .eq. 'temperature') then
       s => fields%items(1)%ptr
       n = s%dof%size()
       do i = 1, n
          s%x(i,1,1,1) = T_cold
       end do
    end if
  end subroutine ics

  !> Smoothstep heated-wall / inlet temperature: T_cold for x < x_heat, ramping
  !! to T_hot over [x_heat, x_heat+dx_ramp]. Same rule serves the cold inlet
  !! (x~0 -> T_cold) and the lower wall (heated only downstream of x_heat).
  pure function Twall(x) result(T)
    real(kind=rp), intent(in) :: x
    real(kind=rp) :: T, t01
    if (.not. do_heat) then
       T = T_cold; return
    end if
    t01 = min(max((x - x_heat)/dx_ramp, 0.0_rp), 1.0_rp)
    t01 = t01*t01*(3.0_rp - 2.0_rp*t01)
    T = T_cold + (T_hot - T_cold)*t01
  end function Twall

  !> Dirichlet hook: STEADY velocity inlet (fixed blunt profile) AND temperature
  !! (cold inlet + heated lower wall). Branch on field name (velocity list -> 'u').
  subroutine inflow_and_wall(field_bc_list, bc, time)
    type(field_list_t), intent(inout) :: field_bc_list
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: s
    type(dofmap_t), pointer :: dof
    integer :: i, msk_idx

    dof => field_bc_list%dof(1)
    if (field_bc_list%items(1)%ptr%name .eq. "u") then
       ! ---------- STEADY blunt-profile velocity inlet (no recycling) ----------
       associate(u => field_bc_list%items(1)%ptr, &
            v => field_bc_list%items(2)%ptr, w => field_bc_list%items(3)%ptr)
          do i = 1, bc%msk(0)
             msk_idx = bc%msk(i)
             u%x(msk_idx,1,1,1) = uin_profile(dof%y(msk_idx,1,1,1))
             v%x(msk_idx,1,1,1) = 0.0_rp
             w%x(msk_idx,1,1,1) = 0.0_rp
          end do
       end associate
    else
       ! ---------- TEMPERATURE: cold inlet + heated lower wall (x-ramp) ----------
       s => field_bc_list%items(1)%ptr
       do i = 1, bc%msk(0)
          msk_idx = bc%msk(i)
          s%x(msk_idx,1,1,1) = Twall(dof%x(msk_idx,1,1,1))
       end do
    end if
  end subroutine inflow_and_wall

  !> Momentum source: buoyancy (-y reduced gravity) + near-inlet trip.
  subroutine fluid_source(scheme_name, rhs, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: rhs
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: ru, rv, rw, rho
    type(dofmap_t), pointer :: dof
    integer :: i, n
    real(kind=rp) :: xx, yy, zz, fb, ft

    if (scheme_name .ne. 'fluid') return
    ru => rhs%get_by_index(1); rv => rhs%get_by_index(2); rw => rhs%get_by_index(3)
    rho => neko_registry%get_field('fluid_rho')
    dof => ru%dof
    n = ru%dof%size()
    do i = 1, n
       fb = 0.0_rp; ft = 0.0_rp
       ! buoyancy: hot fluid (rho<rho_ref) feels upward (+y) reduced gravity
       if (do_gravity) fb = (rho_ref/rho%x(i,1,1,1) - 1.0_rp)*g_nd
       ! trip: localised, time- and z-varying wall-normal kick near the inlet
       if (do_trip) then
          xx = dof%x(i,1,1,1); yy = dof%y(i,1,1,1); zz = dof%z(i,1,1,1)
          if (xx .ge. trip_x0 .and. xx .le. trip_x1 .and. yy .le. trip_yh) then
             ft = trip_A * sin(2.0_rp*pi*4.0_rp*time%t) &
                  * sin(2.0_rp*pi*6.0_rp*zz/Lz) * (1.0_rp - yy/trip_yh)
          end if
       end if
       ru%x(i,1,1,1) = 0.0_rp
       rv%x(i,1,1,1) = fb + ft
       rw%x(i,1,1,1) = 0.0_rp
    end do
  end subroutine fluid_source

end module user
