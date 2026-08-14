! =====================================================================
!  Full CVD train: transitional jet -> diffusor -> hot zone -> extension.
!  Low-Mach, non-dimensional. Mesh train.nmsh (metres), FLOW RUNS -x.
!
!  Frame (differs from hotzone_trip!): x streamwise (inlet x=0, outlet
!  x=-1.2179), y spanwise (wide axis, +-80.2 mm), z VERTICAL (floor
!  -14.5 mm, ceiling +16.5 mm) -> gravity acts on w, not v.
!
!  Non-dimensionalisation: U_ref = PIPE bulk = 1, rho_ref = 1, T_cold=1.
!  Re_pipe = 4000: self-sustaining turbulent pipe flow -- a genuinely
!  turbulent jet (2500 was transitional: puffs, not guaranteed breakdown).
!  U_phys = 72.9 m/s, crest Mach 0.086 (comfortably low-Mach).
!  Duct Re_H ~ 600: the jet's turbulence decays through the chamber and the
!  heated section relaminarizes it fully (nu x21 at T=6).
!  "Jet in, physics decides after."
!
!  Inflow: LAMINAR PARABOLIC profile, bulk exactly 1 -- the hotzone_trip
!  pattern: the inlet is deliberately smooth and the 36% hump at x=-0.1
!  does the tripping. (A 1/7-power turbulent mean profile is not physical
!  at transitional Re; DFSEM synthetic eddies remain the upgrade path.)
! =====================================================================
module user
  use neko
  implicit none

  ! ---- geometry (metres, from build_train.py / ogrid_cvd.py) ----
  real(kind=rp), parameter :: R_PIPE = 0.0151_rp     ! inlet pipe radius
  real(kind=rp), parameter :: X_HEAT = -0.6644_rp    ! heated wall leading edge
  real(kind=rp), parameter :: X_HEND = -0.9509_rp    ! heated wall end
  real(kind=rp), parameter :: L_RAMP = 0.025_rp      ! spatial smoothstep length
  real(kind=rp), parameter :: Z_FLOOR = -0.014523_rp, Z_CEIL = 0.0165_rp

  ! ---- flow ----
  real(kind=rp), parameter :: U_in    = 1.0_rp       ! pipe bulk
  real(kind=rp), parameter :: rho_ref = 1.0_rp
  real(kind=rp), parameter :: Re_p    = 4000.0_rp    ! pipe Reynolds (D=0.0302)
  real(kind=rp), parameter :: Pr      = 0.70_rp
  real(kind=rp), parameter :: mu_ref  = rho_ref*U_in*2.0_rp*R_PIPE/Re_p
  real(kind=rp), parameter :: k_ref   = mu_ref/Pr

  ! ---- thermal ----
  real(kind=rp), parameter :: T_cold = 1.0_rp
  real(kind=rp), parameter :: T_hot  = 6.0_rp        ! 1800 K at T_ref = 300 K
  ! Wall heat-up time, from the near-wall conduction argument (SETUP_NOTES.md):
  ! T_RAMP ~ (2-3 x first duct element 0.67 mm)^2 / alpha, alpha = k_ref =
  ! 1.08e-5 at Re_p 4000 -> 0.17..0.37. 0.30 sits inside the band, 0.30
  ! (~5% of one flow-through ~6.3): the 0.03-was-4x-too-short lesson from
  ! run 4637232 applied with margin.
  real(kind=rp), parameter :: T_RAMP = 0.30_rp

  ! ---- gravity: full 9.81, non-dimensionalised by U_phys = Re_p*nu/D ----
  real(kind=rp), parameter :: nu_H2  = 5.5040e-4_rp  ! m2/s, H2 200 mbar 300 K
  real(kind=rp), parameter :: U_phys = Re_p*nu_H2/(2.0_rp*R_PIPE)  ! 72.9 m/s
  real(kind=rp), parameter :: g_nd   = 9.81_rp/(U_phys*U_phys)     ! 1.85e-3
  logical, parameter :: do_gravity = .true.

contains

  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%startup => startup
    u%initial_conditions => ics
    u%dirichlet_conditions => bc_hook
    u%source_term => fluid_source
  end subroutine user_setup

  subroutine startup(params)
    type(json_file), intent(inout) :: params
    ! Single-scalar case: case.scalar.* is the correct path here (the
    ! multi-scalar shadowing trap only bites "scalars" array cases).
    call params%add("case.fluid.mu", mu_ref)
    call params%add("case.fluid.rho", rho_ref)
    call params%add("case.scalar.lambda", k_ref)
    call params%add("case.scalar.cp", 1.0_rp)
  end subroutine startup

  !> Laminar parabolic pipe profile (Poiseuille), bulk EXACTLY U_in:
  !! u = 2 Ub (1 - (r/R)^2). Deliberately smooth -- the hump trips it.
  pure function ujet(y, z) result(uu)
    real(kind=rp), intent(in) :: y, z
    real(kind=rp) :: uu, r2
    r2 = min((y*y + z*z)/(R_PIPE*R_PIPE), 1.0_rp)
    uu = 2.0_rp*U_in*(1.0_rp - r2)
  end function ujet

  !> Heated-wall temperature: smoothstep in x (leading edge at X_HEAT,
  !! descending -x with the flow) times smoothstep in t (heat-up).
  pure function Twall(x, tnow) result(T)
    real(kind=rp), intent(in) :: x, tnow
    real(kind=rp) :: T, sx, st
    sx = min(max((X_HEAT - x)/L_RAMP, 0.0_rp), 1.0_rp)
    sx = sx*sx*(3.0_rp - 2.0_rp*sx)
    st = min(max(tnow/T_RAMP, 0.0_rp), 1.0_rp)
    st = st*st*(3.0_rp - 2.0_rp*st)
    T = T_cold + (T_hot - T_cold)*sx*st
  end function Twall

  !> IC: jet profile in the pipe decaying to the local duct bulk further in;
  !! simplest robust start: mean profile scaled by local section, plus gentle
  !! near-axis noise. The hump does the real tripping.
  subroutine ics(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(field_t), pointer :: u, v, w, s
    integer :: i, n
    real(kind=rp) :: xx, yy, zz, um, r1, r2, r3

    if (scheme_name .eq. 'fluid') then
       u => fields%items(1)%ptr; v => fields%items(2)%ptr; w => fields%items(3)%ptr
       n = u%dof%size()
       do i = 1, n
          xx = u%dof%x(i,1,1,1); yy = u%dof%y(i,1,1,1); zz = u%dof%z(i,1,1,1)
          if (xx > -0.15_rp) then
             um = ujet(yy, zz)
          else
             um = 0.147_rp                     ! duct bulk (area ratio 6.8)
          end if
          call random_number(r1); call random_number(r2); call random_number(r3)
          u%x(i,1,1,1) = -um + 0.02_rp*(r1 - 0.5_rp)
          v%x(i,1,1,1) = 0.02_rp*(r2 - 0.5_rp)
          w%x(i,1,1,1) = 0.02_rp*(r3 - 0.5_rp)
       end do
    else if (scheme_name .eq. 'temperature') then
       s => fields%items(1)%ptr
       n = s%dof%size()
       do i = 1, n
          s%x(i,1,1,1) = T_cold
       end do
    end if
  end subroutine ics

  !> Dirichlet hook: velocity inlet (zone 1) and heated walls (zones 4+5).
  subroutine bc_hook(field_bc_list, bc, time)
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
             ! flow runs -x
             u%x(msk_idx,1,1,1) = -ujet(dof%y(msk_idx,1,1,1), &
                                        dof%z(msk_idx,1,1,1))
             v%x(msk_idx,1,1,1) = 0.0_rp
             w%x(msk_idx,1,1,1) = 0.0_rp
          end do
       end associate
    else
       ! temperature on zones 4 and 5 (both marked "user"): ramped wall T.
       ! The substrate (x -0.7029..-0.8399) lies past the ramp -> full T_hot.
       associate(s => field_bc_list%items(1)%ptr)
          do i = 1, bc%msk(0)
             msk_idx = bc%msk(i)
             s%x(msk_idx,1,1,1) = Twall(dof%x(msk_idx,1,1,1), time%t)
          end do
       end associate
    end if
  end subroutine bc_hook

  !> Buoyancy: z is UP here (floor at z=-0.0145). Full gravity: at
  !! U_phys = 72.9 m/s, g_nd = 1.85e-3 -- small but not switched off.
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
                if (do_gravity) fb = (rho_ref/rho%x(i,j,k,e) - 1.0_rp)*g_nd
                ru%x(i,j,k,e) = 0.0_rp
                rv%x(i,j,k,e) = 0.0_rp
                rw%x(i,j,k,e) = fb
             end do
          end do
       end do
    end do
  end subroutine fluid_source

end module user
