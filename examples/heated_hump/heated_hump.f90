! =====================================================================
!  Flow over a wall-mounted hump, low-Mach, with:
!   * SUSTAINED turbulent inflow  -> recycling (global_interpolation): the
!     flow from a plane x_rec upstream of the hump is re-injected at the
!     inlet every step (bulk-velocity rescaled), so genuine turbulent
!     structures feed the inlet indefinitely. Seeded by a turbulent IC and
!     kept lit by a small near-inlet trip.
!   * HEATED lower wall for x > x_heat (T ramps T_cold -> T_hot=2), cold
!     elsewhere; ideal-gas rho = P0/(R T) = 1/T.
!   * GRAVITY / buoyancy (reduced-gravity body force in -y).
!
!  Incremental switches (do_heat / do_gravity / do_trip) let us validate
!  the recycling isothermally first, then turn on heating and gravity.
! =====================================================================
module user
  use neko
  implicit none

  ! ---- recycling state (mirrors examples/recycling/recycling.f90) ----
  type(global_interpolation_t) :: interp
  type(matrix_t) :: xyz
  type(vector_t) :: res, Bw
  real(kind=rp) :: vol
  integer :: n_pts = 0
  logical :: rec_init = .false.

  ! ---- geometry / flow parameters (nondimensional) ----
  real(kind=rp), parameter :: U_in   = 1.0_rp
  real(kind=rp), parameter :: Lx     = 25.0_rp
  real(kind=rp), parameter :: Ly     = 3.0_rp
  real(kind=rp), parameter :: Lz     = 4.0_rp
  real(kind=rp), parameter :: Hhump  = 1.0_rp          ! hump height = length scale
  real(kind=rp), parameter :: Re_in  = 500.0_rp     ! lower Re -> larger mu_ref (more dissipation)
  real(kind=rp), parameter :: Pr     = 0.71_rp
  real(kind=rp), parameter :: rho_ref = 1.0_rp         ! cold reference density
  real(kind=rp), parameter :: mu_ref = rho_ref*U_in*Hhump/Re_in
  real(kind=rp), parameter :: k_ref  = mu_ref/Pr

  ! ---- thermal / gravity ----
  real(kind=rp), parameter :: T_cold = 1.0_rp
  real(kind=rp), parameter :: T_hot  = 2.0_rp          ! "heating a bit, ~2"
  real(kind=rp), parameter :: x_heat = 11.0_rp         ! heated wall starts (hump ends x=8)
  real(kind=rp), parameter :: dx_ramp = 1.0_rp         ! smoothstep ramp width
  real(kind=rp), parameter :: g_nd   = 0.5_rp          ! buoyancy magnitude (reduced for stability)

  ! ---- recycling plane: x = inlet + x_rec, upstream of hump foot (x=4) ----
  real(kind=rp), parameter :: x_rec  = 3.0_rp

  ! ---- near-inlet trip (insurance for sustained transition) ----
  real(kind=rp), parameter :: trip_x0 = 0.5_rp, trip_x1 = 1.0_rp
  real(kind=rp), parameter :: trip_yh = 0.4_rp         ! within this height of the lower wall
  real(kind=rp), parameter :: trip_A  = 0.2_rp        ! reduced trip amplitude

  ! ---- incremental switches ----
  logical, parameter :: do_heat    = .true.
  logical, parameter :: do_gravity = .true.
  logical, parameter :: do_trip    = .true.

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

  !> IC: fluid -> turbulent (blunt channel mean + multi-scale perturbations +
  !! random kick) so the recycle plane is turbulent from step 1; temperature
  !! -> uniform cold (the wall heating builds the hot layer downstream).
  subroutine ics(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(field_t), pointer :: u, v, w, s
    integer :: i, n
    real(kind=rp) :: xx, yy, zz, yn, um, gy, r1, r2, r3

    if (scheme_name .eq. 'fluid') then
       u => fields%items(1)%ptr; v => fields%items(2)%ptr; w => fields%items(3)%ptr
       n = u%dof%size()
       do i = 1, n
          xx = u%dof%x(i,1,1,1); yy = u%dof%y(i,1,1,1); zz = u%dof%z(i,1,1,1)
          yn = min(max(yy/Ly, 0.0_rp), 1.0_rp)
          ! blunt turbulent-like mean (0 at walls, flat core), bulk ~ U_in
          um = U_in * 1.25_rp * (1.0_rp - (2.0_rp*yn - 1.0_rp)**8)
          ! near-wall localisation of the perturbations
          gy = exp(-((yy - 0.3_rp)/0.3_rp)**2) + exp(-((yy - (Ly-0.3_rp))/0.3_rp)**2)
          call random_number(r1); call random_number(r2); call random_number(r3)
          u%x(i,1,1,1) = um &
               + 0.10_rp*gy*sin(3.0_rp*2.0_rp*pi*xx/Lx)*cos(4.0_rp*2.0_rp*pi*zz/Lz) &
               + 0.03_rp*gy*sin(17.0_rp*2.0_rp*pi*xx/Lx)*cos(13.0_rp*2.0_rp*pi*zz/Lz) &
               + 0.02_rp*(r1 - 0.5_rp)
          v%x(i,1,1,1) = 0.08_rp*gy*cos(3.0_rp*2.0_rp*pi*xx/Lx)*sin(4.0_rp*2.0_rp*pi*zz/Lz) &
               + 0.02_rp*(r2 - 0.5_rp)
          w%x(i,1,1,1) = 0.08_rp*gy*sin(3.0_rp*2.0_rp*pi*xx/Lx)*cos(4.0_rp*2.0_rp*pi*zz/Lz) &
               + 0.02_rp*(r3 - 0.5_rp)
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

  !> Dirichlet hook: velocity inlet (recycling) AND temperature (cold inlet +
  !! heated lower wall). Branch on the field name (velocity list starts with 'u').
  subroutine inflow_and_wall(field_bc_list, bc, time)
    type(field_list_t), intent(inout) :: field_bc_list
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time
    type(coef_t), pointer :: coef
    type(field_t), pointer :: field, s
    type(dofmap_t), pointer :: dof
    real(kind=rp) :: scale
    integer :: i, msk_idx

    if (field_bc_list%items(1)%ptr%name .eq. "u") then
       ! ---------- RECYCLING velocity inflow ----------
       coef => neko_user_access%case%fluid%c_Xh
       if (time%tstep .eq. 1 .and. .not. rec_init) then
          n_pts = bc%msk(0)
          call xyz%init(3, n_pts); call Bw%init(n_pts)
          do i = 1, n_pts
             msk_idx = bc%msk(i)
             xyz%x(1,i) = coef%dof%x(msk_idx,1,1,1) + x_rec   ! recycle plane upstream
             xyz%x(2,i) = coef%dof%y(msk_idx,1,1,1)
             xyz%x(3,i) = coef%dof%z(msk_idx,1,1,1)
             Bw%x(i)    = coef%B(msk_idx,1,1,1)
          end do
          vol = glsum(Bw%x, n_pts)
          call xyz%copy_from(HOST_TO_DEVICE, .false.)
          call Bw%copy_from(HOST_TO_DEVICE, .false.)
          call interp%init(coef%dof)
          call interp%find_points(xyz%x, n_pts)
          call res%init(n_pts)
          rec_init = .true.
       end if
       associate(u => field_bc_list%items(1)%ptr, &
            v => field_bc_list%items(2)%ptr, w => field_bc_list%items(3)%ptr)
          field => neko_registry%get_field('u')
          call interp%evaluate(res%x, field%x, .false.)
          if (NEKO_BCKND_DEVICE .eq. 1) then
             scale = 1.0_rp/(device_glsc2(res%x_d, Bw%x_d, n_pts)/vol)
             call device_cmult(res%x_d, scale, n_pts)
          else
             scale = 1.0_rp/(glsc2(res%x, Bw%x, n_pts)/vol)
             call cmult(res%x, scale, n_pts)
          end if
          call field_masked_scatter_copy_0(u, res%x, bc%msk, u%size(), n_pts)
          field => neko_registry%get_field('v')
          call interp%evaluate(res%x, field%x, .false.)
          call field_masked_scatter_copy_0(v, res%x, bc%msk, u%size(), n_pts)
          field => neko_registry%get_field('w')
          call interp%evaluate(res%x, field%x, .false.)
          call field_masked_scatter_copy_0(w, res%x, bc%msk, u%size(), n_pts)
       end associate
    else
       ! ---------- TEMPERATURE: cold inlet + heated lower wall (x-ramp) ----------
       dof => field_bc_list%dof(1)
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
