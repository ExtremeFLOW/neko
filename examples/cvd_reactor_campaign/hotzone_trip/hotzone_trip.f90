! =====================================================================
!  SiC CVD hot zone + tripped inlet extension -- low-Mach, non-dimensional.
!
!  Structure follows the running sigma case
!  (neko_low_mach_solver/full_run/heated_hump_full.f90); the Schlatter-Orlu
!  trip module is dropped because this mesh trips the flow GEOMETRICALLY
!  with a wall-mounted hump 240 mm upstream of the hot zone.
!
!  Non-dimensionalisation: rho_ref = 1, U_in = 1, T_cold = 1.  The MESH IS
!  LEFT IN METRES, so the code length scale is 1 m and Reynolds numbers are
!  quoted against the duct height H below.  Ideal gas: rho = P0/(R T) = 1/T.
!
!  Physical operating point: H2 at 200 mbar / 300 K
!    rho = 0.01617 kg/m3, mu = 8.9e-6 Pa.s -> nu = 5.504e-4 m2/s
!  Re_H = 2000 therefore means U_phys = Re_H*nu/H = 35.7 m/s.  NOTE this is
!  ~25x the nominal CVD throughput (~1.4 m/s): at 200 mbar the real reactor
!  flow is laminar (Re_H ~ 80), and "a bit turbulent" is a deliberate
!  modelling choice, not the design point.  See README.
!
!  T_hot = 6 puts the walls at 1800 K (T_ref = 300 K) -- the real SiC CVD
!  temperature, and the point at which the gas-phase mechanism is actually
!  active. NOTE the thermal-entry-length caveat: at 35.7 m/s the residence in
!  the heated section is only 8.0 ms against a conduction time of order 50 ms,
!  so the CORE GAS DOES NOT REACH 1800 K -- chemistry would be confined to the
!  near-wall layers. Full heating needs Re_H < ~320 (laminar).
! =====================================================================
module user
  use neko
  implicit none

  ! ---- duct geometry at the inlet plane (metres, from the .geo) ----
  real(kind=rp), parameter :: X_HALF = 0.0802_rp        ! half width
  real(kind=rp), parameter :: Y_LO   = -0.014523_rp     ! floor
  real(kind=rp), parameter :: Y_HI   = 0.0163_rp        ! ceiling
  real(kind=rp), parameter :: Hduct  = Y_HI - Y_LO      ! 0.030823 m
  real(kind=rp), parameter :: Wduct  = 2.0_rp*X_HALF    ! 0.1604 m

  ! ---- flow parameters (non-dimensional) ----
  real(kind=rp), parameter :: U_in    = 1.0_rp
  real(kind=rp), parameter :: rho_ref = 1.0_rp
  real(kind=rp), parameter :: Re_H    = 2000.0_rp       ! -> Re_Dh ~ 3355, Re_tau ~ 72
  real(kind=rp), parameter :: Pr      = 0.70_rp         ! H2 at 300 K
  real(kind=rp), parameter :: mu_ref  = rho_ref*U_in*Hduct/Re_H
  real(kind=rp), parameter :: k_ref   = mu_ref/Pr

  ! ---- thermal ----
  real(kind=rp), parameter :: Z_HEAT = 0.2850_rp    ! leading edge of the heated wall
  real(kind=rp), parameter :: L_RAMP = 0.025_rp     ! spatial smoothstep length (~3.7 elements)
  real(kind=rp), parameter :: T_RAMP = 0.030_rp     ! wall heat-up time (= near-wall conduction time)
  real(kind=rp), parameter :: T_cold = 1.0_rp           ! inlet / cold walls
  real(kind=rp), parameter :: T_hot  = 6.0_rp           ! heated section + susceptor (= 1800 K at T_ref = 300 K)

  ! ---- gravity: FULL standard gravity, not a reduced/tuned value ----
  ! The mesh is in metres, so the code length scale is 1 m, the acceleration
  ! scale is U_phys^2 / (1 m), and 9.81 m/s2 maps to g_nd = 9.81/U_phys^2.
  ! That is the whole of Earth's gravity -- nothing is scaled down here.
  real(kind=rp), parameter :: G_EARTH = 9.81_rp          ! m/s2
  real(kind=rp), parameter :: nu_H2  = 5.5040e-4_rp      ! m2/s, H2 at 200 mbar / 300 K
  real(kind=rp), parameter :: U_phys = Re_H*nu_H2/Hduct  ! 35.7 m/s
  real(kind=rp), parameter :: g_nd   = G_EARTH/(U_phys*U_phys)   ! = 7.69e-3
  ! Buoyancy is weak here NOT
  ! because gravity was reduced, but because Re_H = 2000 forces U = 35.7 m/s in
  ! 200 mbar H2.  At the real CVD throughput (~1.4 m/s) the SAME 9.81 m/s2 gives
  ! g_nd ~ 0.15 and buoyancy dominates -- but that flow is laminar (Re_H ~ 80).
  ! Turbulence and buoyancy are not simultaneously reachable at this pressure.

  logical, parameter :: do_gravity = .true.

contains

  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%startup => startup
    u%initial_conditions => ics
    u%dirichlet_conditions => inflow_profile
    u%source_term => fluid_source
  end subroutine user_setup

  subroutine startup(params)
    type(json_file), intent(inout) :: params
    call params%add("case.fluid.mu", mu_ref)
    call params%add("case.fluid.rho", rho_ref)
    call params%add("case.scalar.lambda", k_ref)
    call params%add("case.scalar.cp", 1.0_rp)
  end subroutine startup

  !> Separable parabolic inlet profile: zero at all four walls, bulk EXACTLY
  !! U_in (mean of 6s(1-s) on [0,1] is 1), peak 2.25.  Deliberately laminar --
  !! the hump downstream is what trips it.  A blunt profile blew the sigma case
  !! up; the parabola is the fix that made it stable.
  pure function uin_profile(x, y) result(uu)
    real(kind=rp), intent(in) :: x, y
    real(kind=rp) :: uu, xn, yn
    xn = min(max((x + X_HALF)/Wduct, 0.0_rp), 1.0_rp)
    yn = min(max((y - Y_LO)/Hduct,   0.0_rp), 1.0_rp)
    uu = U_in * (6.0_rp*xn*(1.0_rp - xn)) * (6.0_rp*yn*(1.0_rp - yn))
  end function uin_profile

  !> Heated-wall temperature with a SMOOTHSTEP leading edge.
  !! A step from adiabatic (zone 3) to T_hot (zone 4) at z = Z_HEAT is a genuine
  !! singularity of the continuous problem: heat flux -> infinity and thermal BL
  !! thickness -> 0, so no mesh resolves it and the spectral element sitting on
  !! it rings. Measured in run 4637223: min T fell 0.94 -> -0.51 between t=0.05
  !! and 0.50, MONOTONICALLY (not a decaying transient), all of it at
  !! z = 0.251..0.278 and 0.6-1.2 mm off the floor -- i.e. exactly at the
  !! leading edge. Ramping T_wall over L_RAMP removes the singularity.
  !! Same device as heated_hump_full.f90 on sigma.
  ! NB: the time argument is 'tnow', not 't' -- Fortran is case-insensitive,
  ! so a dummy 't' collides with the result variable 'T'.
  pure function Twall(z, tnow) result(T)
    real(kind=rp), intent(in) :: z, tnow
    real(kind=rp) :: T, sz, st
    ! SPATIAL: flow runs -z, so gas first meets the heated wall at z = Z_HEAT
    ! and travels toward smaller z -- the ramp descends from Z_HEAT. Removes
    ! the leading-edge singularity. Zone 5 (susceptor, z = 0.1095..0.2465) lies
    ! entirely below Z_HEAT-L_RAMP, so sz = 1 there and the same function serves
    ! both zones -- which is why zones 4 AND 5 are both marked "user". Ramping
    ! only zone 4 would put a fresh discontinuity at their junction.
    sz = min(max((Z_HEAT - z)/L_RAMP, 0.0_rp), 1.0_rp)
    sz = sz*sz*(3.0_rp - 2.0_rp*sz)
    ! TEMPORAL: walls start at T_cold, matching the uniform-cold IC exactly, so
    ! there is no shock at t=0. Run 4637228 jumped the walls 1->6 instantly and
    ! dt collapsed 26x (1.29e-5 -> 4.86e-7) inside ~25 steps. T_RAMP = 0.03 is
    ! matched to the near-wall conduction time h^2/alpha = (0.824mm)^2/2.2e-5
    ! = 0.031, i.e. the wall heats no faster than the adjacent gas can follow.
    ! Smoothstep (not linear) because the low-Mach expansion goes as DT/Dt: a
    ! linear ramp would jump dT/dt at both ends, and each jump is a step in the
    ! divergence source. S'(s) = 6s(1-s) vanishes at both ends instead.
    st = min(max(tnow/T_RAMP, 0.0_rp), 1.0_rp)
    st = st*st*(3.0_rp - 2.0_rp*st)
    T = T_cold + (T_hot - T_cold)*sz*st
  end function Twall


  !> IC: inlet mean profile everywhere (flow runs -z) + gentle near-wall
  !! perturbations to seed transition; temperature uniform cold.
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
          um = uin_profile(xx, yy)
          ! localise the perturbations near the floor and ceiling
          gy = exp(-((yy - (Y_LO + 0.2_rp*Hduct))/(0.2_rp*Hduct))**2) &
             + exp(-((yy - (Y_HI - 0.2_rp*Hduct))/(0.2_rp*Hduct))**2)
          call random_number(r1); call random_number(r2); call random_number(r3)
          u%x(i,1,1,1) = 0.02_rp*gy*sin(6.0_rp*pi*zz/Hduct)*cos(4.0_rp*pi*xx/Wduct) &
               + 0.01_rp*(r1 - 0.5_rp)
          v%x(i,1,1,1) = 0.02_rp*gy*cos(6.0_rp*pi*zz/Hduct)*sin(4.0_rp*pi*xx/Wduct) &
               + 0.01_rp*(r2 - 0.5_rp)
          w%x(i,1,1,1) = -um + 0.01_rp*(r3 - 0.5_rp)
       end do
    else if (scheme_name .eq. 'temperature') then
       s => fields%items(1)%ptr
       n = s%dof%size()
       do i = 1, n
          s%x(i,1,1,1) = T_cold
       end do
    end if
  end subroutine ics

  !> Dirichlet hook. Zone 1 = velocity inlet; zone 4 = heated walls, whose
  !! temperature is ramped (see Twall). Everything else is plain
  !! dirichlet/neumann in the case file.
  subroutine inflow_profile(field_bc_list, bc, time)
    type(field_list_t), intent(inout) :: field_bc_list
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time
    type(dofmap_t), pointer :: dof
    integer :: i, msk_idx

    dof => field_bc_list%dof(1)
    if (field_bc_list%items(1)%ptr%name .eq. "u") then
       associate(u => field_bc_list%items(1)%ptr, &
            v => field_bc_list%items(2)%ptr, w => field_bc_list%items(3)%ptr)
          do i = 1, bc%msk(0)
             msk_idx = bc%msk(i)
             u%x(msk_idx,1,1,1) = 0.0_rp
             v%x(msk_idx,1,1,1) = 0.0_rp
             ! flow runs -z
             w%x(msk_idx,1,1,1) = -uin_profile(dof%x(msk_idx,1,1,1), &
                                               dof%y(msk_idx,1,1,1))
          end do
       end associate
    else
       ! Heated walls AND susceptor (zones 4 and 5), both marked "user":
       ! smoothstep in z (leading edge) times smoothstep in t (heat-up).
       associate(s => field_bc_list%items(1)%ptr)
          do i = 1, bc%msk(0)
             msk_idx = bc%msk(i)
             s%x(msk_idx,1,1,1) = Twall(dof%z(msk_idx,1,1,1), time%t)
          end do
       end associate
    end if
  end subroutine inflow_profile

  !> Momentum source: buoyancy only.  y is UP (the .geo calls these surfaces
  !! floor/ceiling and the susceptor sits on the floor), so hot light fluid
  !! rising off the heated substrate is the Rayleigh-Benard-unstable
  !! configuration -- the classic CVD buoyancy concern.
  subroutine fluid_source(scheme_name, rhs, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: rhs
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: ru, rv, rw, rho
    integer :: i, j, k, e
    real(kind=rp) :: fb

    if (scheme_name .ne. 'fluid') return
    ru => rhs%get_by_index(1); rv => rhs%get_by_index(2); rw => rhs%get_by_index(3)
    rho => neko_registry%get_field('fluid_rho')

    do e = 1, ru%msh%nelv
       do k = 1, ru%Xh%lz
          do j = 1, ru%Xh%ly
             do i = 1, ru%Xh%lx
                fb = 0.0_rp
                ! Buoyancy per unit mass under full gravity g_nd: hot fluid
                ! (rho < rho_ref) is pushed up (+y), cold fluid down.  The
                ! (rho_ref/rho - 1) factor is the density contrast, not a
                ! weakening of g.
                if (do_gravity) fb = (rho_ref/rho%x(i,j,k,e) - 1.0_rp)*g_nd
                ru%x(i,j,k,e) = 0.0_rp
                rv%x(i,j,k,e) = fb
                rw%x(i,j,k,e) = 0.0_rp
             end do
          end do
       end do
    end do
  end subroutine fluid_source

end module user
