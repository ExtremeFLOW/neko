! Rider-Kothe vortex-in-a-box: a disk is stretched into a thin filament by a
! non-uniform straining flow, the flow reverses at t = T/2, and the disk should
! return to its initial shape at t = T. Neko's CDI compression term takes its
! interface normal from a transported signed-distance field (`normal = "psi"`).
!
! phi is the phase field (Neko scalar `s`), psi the signed-distance field. Saini
! et al. (2026) use the two symbols the other way round -- see CDI_METHOD.md
! section 1 before reading their equation numbers against this file.
!
! This is the only case in the repo with genuine strain, so it is the one where
! re-distancing could be needed. SVV on the psi transport is required whenever
! `normal = "psi"` (startup check); periodic re-distancing is off in every
! shipped `.case`.
!
! `case.cdi.psi_init` chooses how psi is built at t = 0: "exact" is the analytic
! periodic distance, "redistance" is Saini's Algorithm 1 line 2 -- the same
! pseudo-time solve that later maintains the field. Pairing "exact" with periodic
! re-distancing is the incoherent combination CDI_METHOD.md section 4 describes,
! so a run with re-distancing on wants "redistance" here.
!
! Re-distancing and re-initialisation are different operations here and the names
! are not interchangeable -- but note that Saini's own usage differs: for them
! Eq. (44) is the "re-distancing equation" and the *operation* is always
! reseed-then-relax, called re-initialization, with no in-place mode. `seed`
! selects between the two here. See CDI_METHOD.md section 4.
!
! Method parameters live in `case.cdi`, not under a scalar: `case.scalar`
! (singular) and `case.scalars` (plural) are mutually exclusive paths in
! src/case.f90, so a two-scalar case cannot reuse the singular key.
module user
  use neko
  use elementwise_filter, only : elementwise_filter_t
  use math, only : sqrt_inplace
  use device_math, only : device_sqrt_inplace
  implicit none

  integer, parameter :: NORMAL_PHI = 1, NORMAL_PSI = 2
  integer, parameter :: RD_SEED_PHI = 1, RD_SEED_PSI = 2
  integer, parameter :: RD_NITER_MAX = 50000
  integer, parameter :: PSI_INIT_EXACT = 1, PSI_INIT_REDIST = 2

  !> Floor on |grad psi| before it is divided out to make the unit normal.
  !>
  !> It must sit ABOVE the ~1e-19 round-off gradient of a numerically flat
  !> field. Below it, a flat psi (a re-distanced psi is flat outside its band)
  !> gives a unit normal in a random direction, whose div(n) ~ 1/h turns the
  !> -gamma*phi(1-phi)*div(n) part of the compression term into an exponential
  !> source on phi (CDI_METHOD.md section 4.1). Override with case.cdi.grad_floor.
  real(kind=rp), parameter :: GRAD_FLOOR_DEFAULT = 1.0e-6_rp
  real(kind=rp) :: grad_floor = GRAD_FLOOR_DEFAULT

  !> Plain disk: same centre and radius as ../zalesak_disk's slotted one
  !> (Saini Eq. 77), without the slot.
  real(kind=rp), parameter :: cx = 0.5_rp, cy = 0.75_rp, rad = 0.15_rp

  !> Period of the deformation. The flow reverses at T/2.
  real(kind=rp), parameter :: RK_T = 8.0_rp

  real(kind=rp) :: eps, gamma, u_max, u_max0
  real(kind=rp), parameter :: lambda_bg = 1.0e-16_rp
  integer :: normal_kind, report_every, psi_init_kind = PSI_INIT_EXACT
  logical :: reported = .false., vel_ready = .false., mesh_ready = .false.
  real(kind=rp) :: mesh_hgll, mesh_helem, mesh_hn

  logical :: rd_on = .false., rd_ready = .false.
  real(kind=rp) :: rd_dt_tls, rd_rf, rd_band, rd_cfl, rd_dtau, rd_next
  real(kind=rp) :: rd_dtau_set
  integer :: rd_niter, rd_seed = RD_SEED_PHI, rd_events = 0

  !> In-plane only: the mesh is one element thick in z with w = 0 and a
  !> z-invariant field, so D_z phi = 0 and a third direction would buy nothing
  !> but a tighter explicit stability bound.
  integer, parameter :: SVV_NDIR = 2

  !> Spectral vanishing viscosity, Saini Eqs. (24)-(29), as a derived type
  !> because Saini set c0 and N_svv per equation and there are two here.
  type :: svv_t
     character(len=24) :: tag = ""
     logical :: on = .false., ready = .false., imp = .false.
     real(kind=rp) :: c0 = 0.0_rp, ratio = 2.0_rp, nsvv = 0.0_rp, rho = 0.0_rp
     type(elementwise_filter_t) :: filt
     real(kind=rp), allocatable :: eye(:,:)
     type(c_ptr) :: eye_d = C_NULL_PTR
     type(field_t) :: w1, w2, nu, bass
     type(field_t) :: cg_r, cg_p, cg_q, cg_z
   contains
     procedure, pass(this) :: read_params => svv_read_params
     procedure, pass(this) :: init => svv_init
     procedure, pass(this) :: local => svv_local
     procedure, pass(this) :: op => svv_op
     procedure, pass(this) :: apply_imp => svv_apply_imp
     procedure, pass(this) :: step_imp => svv_step_imp
  end type svv_t

  !> The SVV machinery is byte-identical in all four user files: this repo
  !> keeps each case to one self-contained user file, so a fix to it has to be
  !> made in all four (CLAUDE.md).
  type(svv_t) :: svv_psi    ! the psi transport equation, Saini Eq. (43)
  type(svv_t) :: svv_rd     ! the psi re-distancing equation, Saini Eq. (44)

  !> coef%mult as a field, so weighted inner products go through
  !> field_math like everything else. Built once by ensure_mult_field.
  type(field_t) :: mult_f
  logical :: mult_ready = .false.

  !> The spatial part of the velocity, built once. The time dependence is the
  !> scalar cos(pi t/T), so each step is a rescale rather than a rebuild -- which
  !> on a device build also keeps the field resident instead of copying it up
  !> every step.
  type(field_t) :: u0, v0

contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%startup => startup
    user%initialize => initialize
    user%compute => compute
    user%source_term => source_term
    user%material_properties => material_properties
    user%initial_conditions => initial_conditions
  end subroutine user_setup

  subroutine startup(params)
    type(json_file), intent(inout) :: params
    character(len=:), allocatable :: str
    logical :: var_dt

    call json_get(params, "case.cdi.gamma", gamma)
    call json_get(params, "case.cdi.epsilon", eps)
    u_max = 1.0_rp   ! until compute() measures it

    call json_get_or_default(params, "case.cdi.normal", str, "psi")
    select case (trim(str))
    case ("phi")
      normal_kind = NORMAL_PHI
    case ("psi")
      normal_kind = NORMAL_PSI
    case default
      call neko_error("case.cdi.normal must be 'phi' or 'psi'")
    end select

    ! Default 1000 steps = one output interval at the settled dt.
    call json_get_or_default(params, "case.cdi.band_report_every", &
         report_every, 1000)
    if (report_every .lt. 1) &
         call neko_error("case.cdi.band_report_every must be >= 1")

    call json_get_or_default(params, "case.cdi.grad_floor", grad_floor, &
         GRAD_FLOOR_DEFAULT)
    if (grad_floor .le. 0.0_rp) &
         call neko_error("case.cdi.grad_floor must be > 0")

    call json_get_or_default(params, "case.cdi.psi_init", str, "exact")
    select case (trim(str))
    case ("exact")
      psi_init_kind = PSI_INIT_EXACT
    case ("redistance")
      psi_init_kind = PSI_INIT_REDIST
    case default
      call neko_error("case.cdi.psi_init must be 'exact' or 'redistance'")
    end select

    call svv_psi%read_params(params, "case.cdi.svv_psi", "psi transport", 2.0_rp)
    ! psi is pure advection; its transport always carries SVV, as Saini's
    ! Eq. (43) does. phi's equation has its own, physical, diffusion.
    if (params%valid_path("case.cdi.svv_phi")) call neko_error( &
         "case.cdi.svv_phi is not supported: SVV is a psi-only knob " // &
         "(CDI_METHOD.md section 5)")
    if (normal_kind .eq. NORMAL_PSI .and. .not. svv_psi%on) &
         call neko_error("case.cdi.normal = 'psi' needs case.cdi.svv_psi " // &
         "with c0 > 0")
    call svv_rd%read_params(params, "case.cdi.redistance.svv", &
         "psi re-distancing", 4.0_rp)
    ! Not a tuning choice: at Saini's c0 = 2 the operator's spectral radius
    ! would otherwise set the pseudo timestep rather than the CFL.
    svv_rd%imp = .true.

    call json_get_or_default(params, "case.cdi.redistance.enabled", rd_on, &
         .false.)
    call json_get_or_default(params, "case.cdi.redistance.dt_tls", rd_dt_tls, &
         0.5_rp)
    call json_get_or_default(params, "case.cdi.redistance.r_f", rd_rf, 0.1_rp)
    call json_get_or_default(params, "case.cdi.redistance.band", rd_band, 2.5_rp)
    call json_get_or_default(params, "case.cdi.redistance.cfl", rd_cfl, 0.1_rp)
    ! An explicit pseudo timestep, overriding cfl*h_gll_min when > 0. Saini set
    ! dtau_tls = H/(N+1) directly, which is a pseudo-CFL above 1 from N = 4 up
    ! and so cannot be expressed through the guarded `cfl` key -- and with
    ! N_tls = ceiling(band*H/dtau) it also reproduces their 2.5(N+1) exactly.
    call json_get_or_default(params, "case.cdi.redistance.dtau", rd_dtau_set, &
         0.0_rp)
    if (rd_dt_tls .le. 0.0_rp) &
         call neko_error("case.cdi.redistance.dt_tls must be > 0")
    if (rd_cfl .le. 0.0_rp .or. rd_cfl .gt. 1.0_rp) &
         call neko_error("case.cdi.redistance.cfl must be in (0, 1]")
    ! psi starts from a field that has just been made distance-like either way,
    ! so the first event belongs one interval in rather than at t = 0.
    rd_next = rd_dt_tls

    call json_get_or_default(params, "case.cdi.redistance.seed", str, "phi")
    select case (trim(str))
    case ("phi")
      rd_seed = RD_SEED_PHI
    case ("psi")
      rd_seed = RD_SEED_PSI
    case default
      call neko_error("case.cdi.redistance.seed must be 'phi' or 'psi'")
    end select

    if (rd_on .and. normal_kind .ne. NORMAL_PSI) call neko_error( &
         "re-distancing only does something when case.cdi.normal = 'psi'")
    if (psi_init_kind .eq. PSI_INIT_REDIST .and. normal_kind .ne. NORMAL_PSI) &
         call neko_error("case.cdi.psi_init = 'redistance' only does " // &
         "something when case.cdi.normal = 'psi'")

    call json_get_or_default(params, "case.time.variable_timestep", var_dt, &
         .false.)
    if (var_dt) call neko_error("variable_timestep must be false: Neko sizes " &
         // "dt from the advective CFL and ignores the tighter CDI " &
         // "compression limit gamma*u_max*dt/h_gll_min <= 0.05")
  end subroutine startup

  ! ------------------------------------------------------------------------
  ! Backend helpers: the few operations field_math does not cover
  ! ------------------------------------------------------------------------

  !> a = a * b, where b is one of coef's raw geometric arrays rather than a
  !> field_t. field_math has no overload for that pairing.
  subroutine col2_raw(a, b, b_d, n)
    type(field_t), intent(inout) :: a
    real(kind=rp), intent(in) :: b(n)
    type(c_ptr), intent(in) :: b_d
    integer, intent(in) :: n

    if (NEKO_BCKND_DEVICE .eq. 1) then
      call device_col2(a%x_d, b_d, n)
    else
      call col2(a%x, b, n)
    end if
  end subroutine col2_raw

  !> a = sqrt(a), elementwise.
  subroutine field_sqrt(a, n)
    type(field_t), intent(inout) :: a
    integer, intent(in) :: n

    if (NEKO_BCKND_DEVICE .eq. 1) then
      call device_sqrt_inplace(a%x_d, n)
    else
      call sqrt_inplace(a%x, n)
    end if
  end subroutine field_sqrt

  !> Pull a field back to the host so a host loop can read it.
  subroutine to_host(a, n)
    type(field_t), intent(inout) :: a
    integer, intent(in) :: n

    if (NEKO_BCKND_DEVICE .eq. 1) &
         call device_memcpy(a%x, a%x_d, n, DEVICE_TO_HOST, sync = .true.)
  end subroutine to_host

  !> Push a field the host has just written back to the device.
  subroutine to_device(a, n)
    type(field_t), intent(inout) :: a
    integer, intent(in) :: n

    if (NEKO_BCKND_DEVICE .eq. 1) &
         call device_memcpy(a%x, a%x_d, n, HOST_TO_DEVICE, sync = .true.)
  end subroutine to_device

  ! ------------------------------------------------------------------------
  ! Geometry
  ! ------------------------------------------------------------------------

  !> Signed distance to the plain disk, positive inside.
  pure function disk_distance(x, y) result(d)
    real(kind=rp), intent(in) :: x, y
    real(kind=rp) :: d

    d = rad - sqrt((x - cx)**2 + (y - cy)**2)
  end function disk_distance

  !> The same distance made periodic on [0,1]^2: the nearest of the nine images.
  !>
  !> The plane distance jumps across the y = 0/1 seam. That is invisible in phi
  !> -- both sides saturate to 0 -- and fatal to a transported psi, because the
  !> gather-scatter averages across the jump and plants a large artificial
  !> gradient on the seam. The disk sits only 0.10 from the seam, so this is not
  !> a far-field nicety.
  pure function disk_distance_periodic(x, y) result(d)
    real(kind=rp), intent(in) :: x, y
    real(kind=rp) :: d, cand
    integer :: ix, iy

    d = disk_distance(x, y)
    do ix = -1, 1
      do iy = -1, 1
        cand = disk_distance(x + real(ix, rp), y + real(iy, rp))
        if (abs(cand) .lt. abs(d)) d = cand
      end do
    end do
  end function disk_distance_periodic

  ! ------------------------------------------------------------------------
  ! Per-step hook
  ! ------------------------------------------------------------------------

  !> u = sin^2(pi x) sin(2 pi y) cos(pi t/T), v = -sin(2 pi x) sin^2(pi y)
  !> cos(pi t/T). Prescribed every step; case.fluid.freeze = true is what stops
  !> the flow solver overwriting it.
  subroutine compute(time)
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, v, psifld
    type(coef_t), pointer :: coef
    integer :: i, n, nband, cg_iters
    real(kind=rp) :: x, y, tfac, tmp(1), gmin, gmean, gmax
    character(len=LOG_SIZE) :: mess

    coef => neko_user_access%case%fluid%c_Xh
    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    n = u%size()

    if (.not. vel_ready) then
      call u0%init(coef%dof)
      call v0%init(coef%dof)
      tmp = 0.0_rp
      do i = 1, n
        x = u%dof%x(i,1,1,1)
        y = u%dof%y(i,1,1,1)
        u0%x(i,1,1,1) = sin(pi*x)**2 * sin(2.0_rp*pi*y)
        v0%x(i,1,1,1) = -sin(2.0_rp*pi*x) * sin(pi*y)**2
        tmp(1) = max(tmp(1), sqrt(u0%x(i,1,1,1)**2 + v0%x(i,1,1,1)**2))
      end do
      u_max0 = glmax(tmp, 1)
      call to_device(u0, n)
      call to_device(v0, n)
      vel_ready = .true.
    end if

    ! |cos(pi t/T)| <= 1, so the peak speed is at t = 0 and a resolution report
    ! taken on the first step is not an underestimate for the rest of the run.
    tfac = cos(pi*time%t/RK_T)
    call field_cmult2(u, u0, tfac, n)
    call field_cmult2(v, v0, tfac, n)
    u_max = abs(tfac)*u_max0

    if (.not. reported) then
      call report_resolution(coef%dof, time%dt)
      reported = .true.
    end if

    ! An explicit SVV instance is initialised lazily by its own source-term
    ! hook, so if the case file omits `source_terms` on that scalar the term is
    ! silently never applied while the header still reports it as on. Fail loudly
    ! instead.
    if (time%tstep .gt. 1 .and. svv_psi%on .and. .not. svv_psi%imp &
         .and. .not. svv_psi%ready) call neko_error( &
         "case.cdi.svv_psi is on but the psi source term never fired -- the " &
         // "psi scalar needs a user source_terms entry in the case file")

    ! These hooks run after the scalar step and after slag%update(), so mutating
    ! a field here is a clean Lie split: the next step's BDF history sees the
    ! damped field, and so does the output.
    if (svv_psi%imp .and. svv_psi%ready) then
      psifld => neko_registry%get_field('psi')
      call svv_psi%step_imp(coef, psifld, time%dt, cg_iters)
    end if

    if (rd_on) then
      if (time%t .ge. rd_next) then
        call redistance(coef, time)
        ! This runs after slag%update(), so psi's BDF lags still hold the
        ! pre-event field and BDF3 would settle at psi_old + 11/6 (psi_new -
        ! psi_old). Restart the order at 1, as Saini's ireset_ls does.
        neko_user_access%case%fluid%ext_bdf%nadv = 0
        neko_user_access%case%fluid%ext_bdf%ndiff = 0
        rd_next = rd_next + rd_dt_tls
      end if
    end if

    ! The go/no-go number this case exists to produce: psi is only usable as a
    ! normal source while |grad psi| stays near 1 across the interface band.
    if (mod(time%tstep, report_every) .eq. 0) then
      call band_grad_stats(coef, gmin, gmean, gmax, nband)
      write(mess, '(A,F9.4,A,I0)') "  |grad psi| t=", time%t, &
           "  band nodes ", nband
      call neko_log%message(mess)
      write(mess, '(A,3(E11.4))') "    min/mean/max:", gmin, gmean, gmax
      call neko_log%message(mess)
    end if
  end subroutine compute

  !> Saini's Algorithm 1 line 2: *build* psi with the same pseudo-time solve that
  !> later maintains it, instead of seeding an analytic distance. Pairing an
  !> analytic global psi with band-limited periodic re-distancing is the
  !> incoherent combination of CDI_METHOD.md section 4; this is what makes the
  !> two coherent -- and it is the only path that carries over to a geometry
  !> with no analytic distance to seed from.
  !>
  !> It cannot live in the initial_conditions hook: makeneko binds
  !> neko_user_access only after neko_init returns, and the initial conditions
  !> are set inside it, so coef is unreachable there. Here the case is fully
  !> built and no step has run. psi's BDF history is safe -- step 1 is BDF1,
  !> which reads only the current field (scalar_rhs_maker_bdf's lag loop runs
  !> from 2 to nbd), and slag%update() then captures the field built here.
  subroutine initialize(time)
    type(time_state_t), intent(in) :: time
    type(coef_t), pointer :: coef

    if (psi_init_kind .ne. PSI_INIT_REDIST) return

    coef => neko_user_access%case%fluid%c_Xh
    call measure_mesh(coef%dof)
    call neko_log%section("psi initial condition: Saini Alg. 1 line 2")
    call redistance(coef, time, RD_SEED_PHI, "psi_init build")
    call neko_log%end_section()
  end subroutine initialize

  !> Shortest distance between adjacent GLL nodes, the element edge, and the
  !> nominal node spacing H/N -- measured rather than assumed, once. In-plane
  !> only: w = 0 and the field is z-invariant, so the z spacing cannot limit
  !> anything. Split out of report_resolution because psi_init needs these
  !> before the first step, and that report only happens on it.
  subroutine measure_mesh(dof)
    type(dofmap_t), intent(in) :: dof
    integer :: i, j, k, e
    real(kind=rp) :: hmin, hel, d2, tmp(1)

    if (mesh_ready) return

    hmin = huge(0.0_rp)
    hel = 0.0_rp
    do e = 1, dof%msh%nelv
      do k = 1, dof%Xh%lz
        do j = 1, dof%Xh%ly
          do i = 1, dof%Xh%lx - 1
            d2 = (dof%x(i+1,j,k,e) - dof%x(i,j,k,e))**2 &
               + (dof%y(i+1,j,k,e) - dof%y(i,j,k,e))**2
            hmin = min(hmin, sqrt(d2))
          end do
        end do
      end do
      do k = 1, dof%Xh%lz
        do j = 1, dof%Xh%ly - 1
          do i = 1, dof%Xh%lx
            d2 = (dof%x(i,j+1,k,e) - dof%x(i,j,k,e))**2 &
               + (dof%y(i,j+1,k,e) - dof%y(i,j,k,e))**2
            hmin = min(hmin, sqrt(d2))
          end do
        end do
      end do
      hel = max(hel, abs(dof%x(dof%Xh%lx,1,1,e) - dof%x(1,1,1,e)))
    end do
    tmp(1) = hmin
    mesh_hgll = glmin(tmp, 1)
    tmp(1) = hel
    mesh_helem = glmax(tmp, 1)
    mesh_hn = mesh_helem/real(dof%Xh%lx - 1, rp)
    mesh_ready = .true.
  end subroutine measure_mesh

  !> Print the resolution and CFL bookkeeping once, and refuse to run past the
  !> CDI compression limit. Neko checks the advective CFL but not this one.
  subroutine report_resolution(dof, dt)
    type(dofmap_t), intent(in) :: dof
    real(kind=rp), intent(in) :: dt
    real(kind=rp) :: cc
    character(len=LOG_SIZE) :: mess

    call measure_mesh(dof)
    cc = gamma*u_max*dt/mesh_hgll

    call neko_log%section("Rider-Kothe vortex, CDI with a transported normal")
    write(mess, '(A,I0)') "  polynomial order  : ", dof%Xh%lx - 1
    call neko_log%message(mess)
    call logval("H (element edge)  : ", mesh_helem)
    call logval("H/N               : ", mesh_hn)
    call logval("h_gll_min         : ", mesh_hgll)
    call logval("epsilon           : ", eps)
    call logval("gamma             : ", gamma)
    call logval("u_max (t=0, peak) : ", u_max)
    call logval("xi = eps*N/H      : ", eps/mesh_hn)
    call logval("deformation period: ", RK_T)
    if (normal_kind .eq. NORMAL_PSI) then
      call neko_log%message("  normal            : psi (transported distance)")
    else
      call neko_log%message("  normal            : phi")
    end if
    if (psi_init_kind .eq. PSI_INIT_REDIST) then
      call neko_log%message("  psi at t=0        : built by Eq. (44), Alg. 1 line 2")
    else
      call neko_log%message("  psi at t=0        : exact periodic distance")
    end if
    call svv_report(svv_psi)
    if (rd_on) then
      call svv_report(svv_rd)
      call logval("re-distance dt_tls: ", rd_dt_tls)
      call logval("re-distance band  : ", rd_band*mesh_helem)
      if (rd_seed .eq. RD_SEED_PSI) then
        call neko_log%message("  re-distance seed  : psi (in place, no reinit)")
      else
        call neko_log%message("  re-distance seed  : phi (reinit, Saini Eq. 47)")
      end if
      write(mess, '(A,E15.7)') &
           "  re-distance fires : on the timer, dt_tls = ", rd_dt_tls
      call neko_log%message(mess)
    else
      call neko_log%message("  re-distancing     : off")
    end if
    call logval("C_comp            : ", cc)
    call logval("C_adv             : ", u_max*dt/mesh_hgll)
    call neko_log%end_section()

    if (cc .gt. 0.05_rp) call neko_error( &
         "compression CFL gamma*u_max*dt/h_gll_min > 0.05 -- reduce the timestep")
  end subroutine report_resolution

  subroutine logval(label, v)
    character(len=*), intent(in) :: label
    real(kind=rp), intent(in) :: v
    character(len=LOG_SIZE) :: mess

    write(mess, '(A,A,E15.7)') "  ", label, v
    call neko_log%message(mess)
  end subroutine logval

  ! ------------------------------------------------------------------------
  ! The compression term and its normal
  ! ------------------------------------------------------------------------

  !> n = grad(src)/|grad(src)|, made continuous the way the solver makes any
  !> element-local quantity continuous -- gather-scatter the sum, then divide by
  !> the node multiplicity -- before it is normalised. `w` comes back holding
  !> |grad(src)| after the floor, which the callers reuse.
  subroutine unit_normal(coef, src, g1, g2, g3, w)
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(in) :: src
    type(field_t), intent(inout) :: g1, g2, g3, w
    integer :: n

    n = coef%dof%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
      call grad(g1%x_d, g2%x_d, g3%x_d, src%x_d, coef)
    else
      call grad(g1%x, g2%x, g3%x, src%x, coef)
    end if
    call coef%gs_h%op(g1, GS_OP_ADD)
    call coef%gs_h%op(g2, GS_OP_ADD)
    call coef%gs_h%op(g3, GS_OP_ADD)
    call col2_raw(g1, coef%mult, coef%mult_d, n)
    call col2_raw(g2, coef%mult, coef%mult_d, n)
    call col2_raw(g3, coef%mult, coef%mult_d, n)

    ! |grad| via col3/addcol3 rather than field_vdot3: that routine declares
    ! its result intent(out) on a field_t, which deallocates the field's own
    ! storage on entry. Nothing in Neko calls it, so the defect is unexercised
    ! upstream -- avoid it here rather than rely on it.
    call field_col3(w, g1, g1, n)
    call field_addcol3(w, g2, g2, n)
    call field_addcol3(w, g3, g3, n)
    call field_sqrt(w, n)
    call field_cpwmax2(w, grad_floor, n)
    call field_invcol2(g1, w, n)
    call field_invcol2(g2, w, n)
    call field_invcol2(g3, w, n)
  end subroutine unit_normal

  !> rhs = gamma*u_max * div( -phi*(1-phi) * n ).
  subroutine compression_rhs(coef, s, src, rhs_s)
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(in) :: s, src
    type(field_t), intent(inout) :: rhs_s
    type(field_t), pointer :: g1, g2, g3, w
    integer :: ind(4), n

    n = coef%dof%size()
    call neko_scratch_registry%request_field(g1, ind(1), .false.)
    call neko_scratch_registry%request_field(g2, ind(2), .false.)
    call neko_scratch_registry%request_field(g3, ind(3), .false.)
    call neko_scratch_registry%request_field(w, ind(4), .false.)

    call unit_normal(coef, src, g1, g2, g3, w)

    ! w = phi*(phi - 1) = -phi*(1-phi), the compression flux magnitude. Built
    ! this way, not as col3+sub2, because field_sub2's second argument is
    ! intent(inout) and `s` may be the same field as `src`.
    call field_copy(w, s, n)
    call field_cadd(w, -1.0_rp, n)
    call field_col2(w, s, n)
    call field_col2(g1, w, n)
    call field_col2(g2, w, n)
    call field_col2(g3, w, n)

    if (NEKO_BCKND_DEVICE .eq. 1) then
      call div(rhs_s%x_d, g1%x_d, g2%x_d, g3%x_d, coef)
    else
      call div(rhs_s%x, g1%x, g2%x, g3%x, coef)
    end if
    call field_cmult(rhs_s, gamma*u_max, n)

    call neko_scratch_registry%relinquish_field(ind)
  end subroutine compression_rhs

  !> min / mean / max of |grad psi| over the interface band, and the band size.
  !> The band is the project's convention, phi(1-phi) > 1e-4.
  subroutine band_grad_stats(coef, gmin, gmean, gmax, nband)
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(out) :: gmin, gmean, gmax
    integer, intent(out) :: nband
    type(field_t), pointer :: s, psifld, g1, g2, g3, w
    integer :: ind(4), i, n
    real(kind=rp) :: gsum, loc(1)

    s => neko_registry%get_field('s')
    psifld => neko_registry%get_field('psi')
    n = psifld%size()

    call neko_scratch_registry%request_field(g1, ind(1), .false.)
    call neko_scratch_registry%request_field(g2, ind(2), .false.)
    call neko_scratch_registry%request_field(g3, ind(3), .false.)
    call neko_scratch_registry%request_field(w, ind(4), .false.)

    if (NEKO_BCKND_DEVICE .eq. 1) then
      call grad(g1%x_d, g2%x_d, g3%x_d, psifld%x_d, coef)
    else
      call grad(g1%x, g2%x, g3%x, psifld%x, coef)
    end if
    call coef%gs_h%op(g1, GS_OP_ADD)
    call coef%gs_h%op(g2, GS_OP_ADD)
    call coef%gs_h%op(g3, GS_OP_ADD)
    call col2_raw(g1, coef%mult, coef%mult_d, n)
    call col2_raw(g2, coef%mult, coef%mult_d, n)
    call col2_raw(g3, coef%mult, coef%mult_d, n)
    call field_col3(w, g1, g1, n)      ! see unit_normal on field_vdot3
    call field_addcol3(w, g2, g2, n)
    call field_addcol3(w, g3, g3, n)
    call field_sqrt(w, n)

    call to_host(w, n)
    call to_host(s, n)

    gmin = huge(0.0_rp)
    gmax = 0.0_rp
    gsum = 0.0_rp
    nband = 0
    do i = 1, n
      if (s%x(i,1,1,1)*(1.0_rp - s%x(i,1,1,1)) .gt. 1.0e-4_rp) then
        gmin = min(gmin, w%x(i,1,1,1))
        gmax = max(gmax, w%x(i,1,1,1))
        gsum = gsum + w%x(i,1,1,1)
        nband = nband + 1
      end if
    end do

    loc(1) = gmin
    gmin = glmin(loc, 1)
    loc(1) = gmax
    gmax = glmax(loc, 1)
    loc(1) = gsum
    gsum = glsum(loc, 1)
    loc(1) = real(nband, rp)
    nband = nint(glsum(loc, 1))

    gmean = 0.0_rp
    if (nband .gt. 0) then
      gmean = gsum/real(nband, rp)
    else
      gmin = 0.0_rp
    end if

    call neko_scratch_registry%relinquish_field(ind)
  end subroutine band_grad_stats

  ! ------------------------------------------------------------------------
  ! Re-distancing: Saini et al. (2026) Eqs. (44)-(47)
  !
  ! redistance and band_grad_stats are identical in the three coupled user
  ! files, so a fix here has to be made in all three.
  ! ------------------------------------------------------------------------

  !> sgn(psi) = tanh(psi/(2 eps)).
  !>
  !> The one operation in this file with no device counterpart, so it round-trips
  !> through the host. It runs 3x per pseudo-step inside a re-distancing event
  !> and not at all otherwise, which over a full run is seconds.
  subroutine rd_sgn(psifld, sgnfld, n)
    type(field_t), intent(inout) :: psifld, sgnfld
    integer, intent(in) :: n
    integer :: i

    call to_host(psifld, n)
    do i = 1, n
      sgnfld%x(i,1,1,1) = tanh(psifld%x(i,1,1,1)/(2.0_rp*eps))
    end do
    call to_device(sgnfld, n)
  end subroutine rd_sgn

  !> L(psi) = sgn(psi)*(1 - |grad psi|), the right-hand side of Saini Eq. (44).
  !>
  !> They write it as `-w.grad(psi) + sgn(psi)` with `w = sgn(psi)*n` and
  !> `n = grad(psi)/|grad psi|`, so `w.grad(psi) = sgn(psi)*|grad psi|` and the
  !> two forms are identical -- this one needs no convective operator. Evaluating
  !> sgn at the current psi pins the zero level set: the source vanishes there.
  subroutine rd_rhs(coef, psifld, lval, g1, g2, g3, w)
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: psifld, lval, g1, g2, g3, w
    integer :: n

    n = coef%dof%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
      call grad(g1%x_d, g2%x_d, g3%x_d, psifld%x_d, coef)
    else
      call grad(g1%x, g2%x, g3%x, psifld%x, coef)
    end if
    call coef%gs_h%op(g1, GS_OP_ADD)
    call coef%gs_h%op(g2, GS_OP_ADD)
    call coef%gs_h%op(g3, GS_OP_ADD)
    call col2_raw(g1, coef%mult, coef%mult_d, n)
    call col2_raw(g2, coef%mult, coef%mult_d, n)
    call col2_raw(g3, coef%mult, coef%mult_d, n)

    call field_col3(w, g1, g1, n)      ! see unit_normal on field_vdot3
    call field_addcol3(w, g2, g2, n)
    call field_addcol3(w, g3, g3, n)
    call field_sqrt(w, n)              ! w = |grad psi|
    call rd_sgn(psifld, lval, n)       ! lval = sgn(psi)
    call field_cmult(w, -1.0_rp, n)
    call field_cadd(w, 1.0_rp, n)      ! w = 1 - |grad psi|
    call field_col2(lval, w, n)
  end subroutine rd_rhs

  !> One re-distancing event: optionally re-seed psi from phi, then iterate
  !> Saini Eq. (44) in pseudo time.
  !>
  !> SSP-RK3, not explicit Euler: linearised, Eq. (44) is advection at unit speed
  !> along the normal, and the collocation derivative on a periodic mesh has a
  !> purely imaginary spectrum that Euler amplifies at every wavenumber. RK3
  !> covers the imaginary axis to 1.73.
  subroutine redistance(coef, time, seed, label)
    type(coef_t), intent(inout) :: coef
    type(time_state_t), intent(in) :: time
    integer, intent(in), optional :: seed
    character(len=*), intent(in), optional :: label
    type(field_t), pointer :: s, psifld
    type(field_t), pointer :: g1, g2, g3, w, lval, d0, nx, ny, nz
    integer :: ind(9), i, it, n, cg_iters, nband, use_seed
    real(kind=rp) :: gmin, gmax, gmean, gmin0, gmax0, gmean0, dn, loc(1)
    character(len=LOG_SIZE) :: mess

    s => neko_registry%get_field('s')
    psifld => neko_registry%get_field('psi')
    n = psifld%size()

    if (.not. rd_ready) then
      ! dtau from the pseudo-CFL on h_gll_min, and the iteration count from
      ! Saini's band distance. Their own automated dtau_tls = H/(N+1) is a
      ! pseudo-CFL of 0.9 at N = 3 up to 2.75 at N = 10, and even at 0.9 it is
      ! unstable under repeated events (CDI_METHOD.md 4.2). redistance.dtau
      ! sets it explicitly for a one-shot build.
      if (rd_dtau_set .gt. 0.0_rp) then
        rd_dtau = rd_dtau_set
      else
        rd_dtau = rd_cfl*mesh_hgll
      end if
      ! The epsilon keeps an exactly-divisible band from rounding up a whole
      ! extra pseudo step -- Saini's dtau = H/(N+1) is exactly divisible, and
      ! without it their N_tls = 2.5(N+1) comes out one too many.
      rd_niter = ceiling(rd_band*mesh_helem/rd_dtau - 1.0e-9_rp)
      if (rd_niter .gt. RD_NITER_MAX) call neko_error( &
           "re-distancing needs too many pseudo steps -- raise " // &
           "case.cdi.redistance.cfl")
      if (svv_rd%on .and. .not. svv_rd%ready) call svv_rd%init(coef, rd_dtau)
      write(mess, '(A,E13.6,A,I0,A,F6.3)') "  re-distancing: dtau=", rd_dtau, &
           ", ", rd_niter, " steps, pseudo-CFL ", rd_dtau/mesh_hgll
      call neko_log%message(mess)
      rd_ready = .true.
    end if

    call neko_scratch_registry%request_field(g1, ind(1), .false.)
    call neko_scratch_registry%request_field(g2, ind(2), .false.)
    call neko_scratch_registry%request_field(g3, ind(3), .false.)
    call neko_scratch_registry%request_field(w, ind(4), .false.)
    call neko_scratch_registry%request_field(lval, ind(5), .false.)
    call neko_scratch_registry%request_field(d0, ind(6), .false.)
    call neko_scratch_registry%request_field(nx, ind(7), .false.)
    call neko_scratch_registry%request_field(ny, ind(8), .false.)
    call neko_scratch_registry%request_field(nz, ind(9), .false.)

    use_seed = rd_seed
    if (present(seed)) use_seed = seed
    if (.not. present(label)) then
      rd_events = rd_events + 1
      ! The normal the compression term is reading right now; ||dn|| below is
      ! how far this event moves it.
      call band_grad_stats(coef, gmin0, gmean0, gmax0, nband)
      call unit_normal(coef, psifld, nx, ny, nz, w)
    end if

    ! Saini Eq. (47): discard the transported field and restart from the phase
    ! field, scaled by r_f. phi's errors re-enter psi here, which closes a loop:
    ! psi sets n, n moves phi, phi re-seeds psi. seed = 'psi' iterates in place
    ! from the transported field instead, which keeps phi out of psi's equation.
    ! sgn is evaluated on the current psi either way, so the zero contour is
    ! pinned in both.
    if (use_seed .eq. RD_SEED_PHI) then
      call field_copy(psifld, s, n)
      call field_cadd(psifld, -0.5_rp, n)
      call field_cmult(psifld, rd_rf, n)
    end if

    ! The initial build has no prior psi to compare against and no normal in
    ! use yet, so it reports the raw Eq. (47) seed as its "before" row -- the
    ! field the solve actually starts from -- and skips ||dn|| entirely.
    if (present(label)) call band_grad_stats(coef, gmin0, gmean0, gmax0, nband)

    do it = 1, rd_niter
      call field_copy(d0, psifld, n)

      call rd_rhs(coef, psifld, lval, g1, g2, g3, w)
      call field_add2s2(psifld, lval, rd_dtau, n)                  ! u1

      call rd_rhs(coef, psifld, lval, g1, g2, g3, w)
      call field_add2s2(psifld, lval, rd_dtau, n)
      call field_cmult(psifld, 0.25_rp, n)
      call field_add2s2(psifld, d0, 0.75_rp, n)                    ! u2

      call rd_rhs(coef, psifld, lval, g1, g2, g3, w)
      call field_add2s2(psifld, lval, rd_dtau, n)
      call field_cmult(psifld, 2.0_rp/3.0_rp, n)
      call field_add2s2(psifld, d0, 1.0_rp/3.0_rp, n)              ! u^{n+1}

      if (svv_rd%on) call svv_rd%step_imp(coef, psifld, rd_dtau, cg_iters)
    end do

    ! |grad psi| before and after says whether psi was worth re-distancing;
    ! ||dn|| is how far the event moves the normal the compression term reads.
    call band_grad_stats(coef, gmin, gmean, gmax, nband)
    if (present(label)) then
      write(mess, '(A,A,A,F9.4,A,I0)') "  ", label, " t=", time%t, &
           " band nodes ", nband
      call neko_log%message(mess)
      write(mess, '(A,3(E11.4))') "    |grad psi| seed  :", gmin0, gmean0, gmax0
      call neko_log%message(mess)
      write(mess, '(A,3(E11.4))') "    |grad psi| built :", gmin, gmean, gmax
      call neko_log%message(mess)
      call neko_scratch_registry%relinquish_field(ind)
      return
    end if

    call unit_normal(coef, psifld, g1, g2, g3, w)
    call to_host(g1, n)
    call to_host(g2, n)
    call to_host(g3, n)
    call to_host(nx, n)
    call to_host(ny, n)
    call to_host(nz, n)
    call to_host(s, n)
    dn = 0.0_rp
    do i = 1, n
      if (s%x(i,1,1,1)*(1.0_rp - s%x(i,1,1,1)) .gt. 1.0e-4_rp) then
        dn = dn + (g1%x(i,1,1,1) - nx%x(i,1,1,1))**2 &
                + (g2%x(i,1,1,1) - ny%x(i,1,1,1))**2 &
                + (g3%x(i,1,1,1) - nz%x(i,1,1,1))**2
      end if
    end do
    loc(1) = dn
    dn = sqrt(glsum(loc, 1))

    write(mess, '(A,I0,A,F9.4,A,I0)') "  redist event ", rd_events, &
         " t=", time%t, " band nodes ", nband
    call neko_log%message(mess)
    write(mess, '(A,3(E11.4))') "    |grad psi| before:", gmin0, gmean0, gmax0
    call neko_log%message(mess)
    write(mess, '(A,3(E11.4))') "    |grad psi| after :", gmin, gmean, gmax
    call neko_log%message(mess)
    write(mess, '(A,E11.4)') "    ||dn|| = ", dn
    call neko_log%message(mess)

    call neko_scratch_registry%relinquish_field(ind)
  end subroutine redistance

  ! ------------------------------------------------------------------------
  ! Scalar hooks
  ! ------------------------------------------------------------------------

  !> Neko zeroes the user field before calling this. `psi` -- Saini Eq. (43) --
  !> carries only its SVV term, if that instance is on and explicit.
  subroutine source_term(scheme_name, rhs, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: rhs
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: rhs_s, s, psifld, src
    type(coef_t), pointer :: coef
    integer :: n

    if (scheme_name .ne. 's' .and. scheme_name .ne. 'psi') return

    coef => neko_user_access%case%fluid%c_Xh
    n = coef%dof%size()
    rhs_s => rhs%get_by_index(1)

    if (scheme_name .eq. 'psi') then
      if (.not. svv_psi%on .or. svv_psi%imp) return
      if (.not. svv_psi%ready) call svv_psi%init(coef, time%dt)
      psifld => neko_registry%get_field('psi')
      call svv_psi%op(coef, psifld, rhs_s)
      call field_cmult(rhs_s, -1.0_rp, n)
      return
    end if

    if (gamma .le. 0.0_rp) return

    s => neko_registry%get_field('s')
    if (normal_kind .eq. NORMAL_PSI) then
      src => neko_registry%get_field('psi')
    else
      src => s
    end if
    call compression_rhs(coef, s, src, rhs_s)
  end subroutine source_term

  subroutine material_properties(scheme_name, properties, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: properties
    type(time_state_t), intent(in) :: time

    if (scheme_name .eq. "fluid") then
      call field_cfill(properties%get("fluid_rho"), 1.0_rp)
      call field_cfill(properties%get("fluid_mu"), 1.0_rp)
    else if (scheme_name .eq. "s") then
      call field_cfill(properties%get('s_cp'), 1.0_rp)
      if (gamma .gt. 0.0_rp) then
        ! Physical, not a stabilisation knob: this is the diffusion that
        ! balances compression to hold a tanh profile of half-width eps
        ! (CDI_METHOD.md 5). u_max varies as |cos(pi t/T)| here, so both halves
        ! of that balance shrink together and the equilibrium width is unchanged.
        ! The solve reads s_lambda_tot, which Neko copies from s_lambda only at
        ! init unless a turbulence model is set; without the second fill the
        ! diffusion stays at eps*gamma*u_max(0) and the filament dissolves.
        call field_cfill(properties%get('s_lambda'), eps*gamma*u_max)
        call field_cfill(neko_registry%get_field('s_lambda_tot'), &
             eps*gamma*u_max)
      else
        call field_cfill(properties%get('s_lambda'), lambda_bg)
      end if
    else if (scheme_name .eq. "psi") then
      call field_cfill(properties%get('psi_cp'), 1.0_rp)
      call field_cfill(properties%get('psi_lambda'), lambda_bg)
    end if
  end subroutine material_properties

  !> phi from the plane distance; psi from the *periodic* distance, which the
  !> phase field does not need and a transported field cannot do without.
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(field_t), pointer :: f
    integer :: i

    if (scheme_name .ne. 's' .and. scheme_name .ne. 'psi') return
    f => fields%items(1)%ptr

    if (scheme_name .eq. 'psi') then
      if (psi_init_kind .eq. PSI_INIT_REDIST) then
        ! Saini's Algorithm 1 builds psi from phi alone -- line 2, the Eq. (44)
        ! solve seeded by Eq. (47) -- and that is the whole point of the option:
        ! a real application has no analytic distance to fall back on. So do not
        ! evaluate one here. initialize() builds the field before the first
        ! step; this only has to leave it empty.
        do i = 1, f%dof%size()
          f%x(i,1,1,1) = 0.0_rp
        end do
      else
        do i = 1, f%dof%size()
          f%x(i,1,1,1) = disk_distance_periodic(f%dof%x(i,1,1,1), &
               f%dof%y(i,1,1,1))
        end do
      end if
    else
      do i = 1, f%dof%size()
        f%x(i,1,1,1) = 0.5_rp*(1.0_rp + tanh( &
             disk_distance(f%dof%x(i,1,1,1), f%dof%y(i,1,1,1))/(2.0_rp*eps)))
      end do
    end if

    if ((NEKO_BCKND_DEVICE .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
         .or. (NEKO_BCKND_OPENCL .eq. 1)) then
      call device_memcpy(f%x, f%x_d, f%dof%size(), HOST_TO_DEVICE, sync = .true.)
    end if
  end subroutine initial_conditions

  !> Reads one SVV instance's parameters. `nsvv_ratio` is the divisor d in
  !> Saini's N_svv = N/d; larger d is a broader, stronger filter. c0 <= 0
  !> switches the instance off, which is the default everywhere.
  subroutine svv_read_params(this, params, path, tag, def_ratio)
    class(svv_t), intent(inout) :: this
    type(json_file), intent(inout) :: params
    character(len=*), intent(in) :: path, tag
    real(kind=rp), intent(in) :: def_ratio

    this%tag = tag
    call json_get_or_default(params, path // ".c0", this%c0, 0.0_rp)
    call json_get_or_default(params, path // ".nsvv_ratio", this%ratio, &
         def_ratio)
    if (this%ratio .le. 0.0_rp) &
         call neko_error(path // ".nsvv_ratio must be > 0")
    this%on = this%c0 .gt. 0.0_rp
    call json_get_or_default(params, path // ".implicit", this%imp, .false.)
  end subroutine svv_read_params

  !> a = b, same pairing as col2_raw.
  subroutine copy_raw(a, b, b_d, n)
    type(field_t), intent(inout) :: a
    real(kind=rp), intent(in) :: b(n)
    type(c_ptr), intent(in) :: b_d
    integer, intent(in) :: n

    if (NEKO_BCKND_DEVICE .eq. 1) then
      call device_copy(a%x_d, b_d, n)
    else
      call copy(a%x, b, n)
    end if
  end subroutine copy_raw

  subroutine svv_report(this)
    type(svv_t), intent(in) :: this
    character(len=LOG_SIZE) :: mess

    if (this%on) then
      write(mess, '(A,A,A,F4.1,A,E11.4,A)') "  SVV ", trim(this%tag), &
           ": N_svv = N/", this%ratio, ", c0 = ", this%c0, &
           merge(" (implicit)", " (explicit)", this%imp)
    else
      write(mess, '(A,A,A)') "  SVV ", trim(this%tag), ": off"
    end if
    call neko_log%message(mess)
  end subroutine svv_report

  !> y = S_vv x, element-local (not yet gathered). Saini Eq. (28),
  !> S_vv = D~^T G D~ with D~_l = B Q B^-1 D_l, restricted to the diagonal
  !> geometric factors -- exact on an orthogonal mesh, which init checks.
  !>
  !> The viscosity sits *inside* the bilinear form rather than left-multiplying
  !> the assembled operator as their Eq. (33) has it. Theirs is not
  !> conservative; this is, exactly, because D annihilates constants -- and mass
  !> conservation is one of the things this remedy is judged on.
  !>
  !> tnsr3d's second and third matrices are the transposes of the operators they
  !> apply, so dyt and fht differentiate/filter while dy and fh in those slots
  !> transpose.
  subroutine svv_local(this, coef, x, y)
    class(svv_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: x, y
    integer :: lx, nel, n, dd

    lx = coef%Xh%lx
    nel = coef%msh%nelv
    n = coef%dof%size()
    call field_rzero(y, n)

    do dd = 1, SVV_NDIR
      select case (dd)
      case (1)
        call tnsr3d(this%w1%x, lx, x%x, lx, coef%Xh%dx, this%eye, this%eye, nel)
        call tnsr3d(this%w2%x, lx, this%w1%x, lx, this%filt%fh, this%eye, &
             this%eye, nel)
        call col2_raw(this%w2, coef%G11, coef%G11_d, n)
      case (2)
        call tnsr3d(this%w1%x, lx, x%x, lx, this%eye, coef%Xh%dyt, this%eye, &
             nel)
        call tnsr3d(this%w2%x, lx, this%w1%x, lx, this%eye, this%filt%fht, &
             this%eye, nel)
        call col2_raw(this%w2, coef%G22, coef%G22_d, n)
      end select

      call field_col2(this%w2, this%nu, n)

      select case (dd)
      case (1)
        call tnsr3d(this%w1%x, lx, this%w2%x, lx, this%filt%fht, this%eye, &
             this%eye, nel)
        call tnsr3d(this%w2%x, lx, this%w1%x, lx, coef%Xh%dxt, this%eye, &
             this%eye, nel)
      case (2)
        call tnsr3d(this%w1%x, lx, this%w2%x, lx, this%eye, this%filt%fh, &
             this%eye, nel)
        call tnsr3d(this%w2%x, lx, this%w1%x, lx, this%eye, coef%Xh%dy, &
             this%eye, nel)
      end select

      call field_add2(y, this%w2, n)
    end do
  end subroutine svv_local

  !> y = B^-1 S_vv x, gathered -- the strong form the scalar source term wants.
  subroutine svv_op(this, coef, x, y)
    class(svv_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: x, y
    integer :: n

    n = coef%dof%size()
    call this%local(coef, x, y)
    call coef%gs_h%op(y, GS_OP_ADD)
    call col2_raw(y, coef%Binv, coef%Binv_d, n)
  end subroutine svv_op

  !> y = (B + dt*S_vv) x, for the implicit sub-step below.
  subroutine svv_apply_imp(this, coef, x, y, dt)
    class(svv_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: x, y
    real(kind=rp), intent(in) :: dt
    integer :: n

    n = coef%dof%size()
    call this%local(coef, x, y)
    call coef%gs_h%op(y, GS_OP_ADD)
    call field_cmult(y, dt, n)
    call field_addcol3(y, this%bass, x, n)
  end subroutine svv_apply_imp

  !> One backward-Euler sub-step of the SVV term, Lie-split from the integrator
  !> that carries everything else:
  !>
  !>     (B + dt*S_vv) phi_new = B phi_old
  !>
  !> The SEM mass matrix is diagonal, so B here is just the assembled Jacobian
  !> weight. Being implicit this is unconditionally stable -- it removes the
  !> dt*rho limit an explicit source term imposes -- at the price of treating
  !> the stabilisation to first order in dt. Mass is conserved *exactly*, not
  !> just to solver tolerance: 1^T S_vv = 0, so 1^T B phi_new = 1^T B phi_old
  !> whatever the iteration does.
  !>
  !> Solved by mass-preconditioned CG. The condition number is 1 + dt*rho, so
  !> the iteration count stays small precisely where the explicit treatment
  !> would have failed.
  subroutine svv_step_imp(this, coef, s, dt, iters)
    class(svv_t), intent(inout) :: this
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: s
    real(kind=rp), intent(in) :: dt
    integer, intent(out) :: iters
    real(kind=rp) :: rz, rz_new, alpha, beta, r0
    integer :: n, it
    integer, parameter :: MAXIT = 200
    real(kind=rp), parameter :: RTOL = 1.0e-12_rp

    n = coef%dof%size()

    ! r = B*phi_old - A(phi_old), with phi_old itself as the initial guess
    call this%apply_imp(coef, s, this%cg_r, dt)
    call field_col3(this%cg_q, this%bass, s, n)
    call field_sub2(this%cg_q, this%cg_r, n)
    call field_copy(this%cg_r, this%cg_q, n)
    r0 = sqrt(field_glsc3(this%cg_r, mult_f, this%cg_r, n))

    call field_invcol3(this%cg_z, this%cg_r, this%bass, n)
    call field_copy(this%cg_p, this%cg_z, n)
    rz = field_glsc3(this%cg_r, mult_f, this%cg_z, n)

    iters = 0
    do it = 1, MAXIT
      call this%apply_imp(coef, this%cg_p, this%cg_q, dt)
      alpha = rz/field_glsc3(this%cg_q, mult_f, this%cg_p, n)
      call field_add2s2(s, this%cg_p, alpha, n)
      call field_add2s2(this%cg_r, this%cg_q, -alpha, n)
      iters = it
      if (sqrt(field_glsc3(this%cg_r, mult_f, this%cg_r, n)) &
           .le. RTOL*max(r0, tiny(1.0_rp))) exit
      call field_invcol3(this%cg_z, this%cg_r, this%bass, n)
      rz_new = field_glsc3(this%cg_r, mult_f, this%cg_z, n)
      beta = rz_new/rz
      rz = rz_new
      call field_add2s1(this%cg_p, this%cg_z, beta, n)
    end do
  end subroutine svv_step_imp

  !> Build the operator once, check it, and refuse to run if it is too stiff.
  !> Saini fold SVV into the Helmholtz operator (their Eq. 35); the explicit
  !> path here goes through the source term instead, so BDF3/EXT3 has to carry
  !> it on the real axis and dt*rho(B^-1 S_vv) is the binding constraint rather
  !> than any CFL.
  subroutine svv_init(this, coef, dt)
    class(svv_t), intent(inout) :: this
    type(coef_t), intent(inout), target :: coef
    real(kind=rp), intent(in) :: dt
    type(field_t), pointer :: p, q, sp, bfld
    real(kind=rp), allocatable :: trans(:)
    real(kind=rp) :: tmp(1), h_elem, gdiag, goff, ref, sym1, sym2, nul, nrm
    integer :: ind(4), i, e, lx, nel, n, it
    character(len=LOG_SIZE) :: mess

    lx = coef%Xh%lx
    nel = coef%msh%nelv
    n = coef%dof%size()
    this%nsvv = real(lx - 1, rp)/this%ratio

    call ensure_mult_field(coef)

    ! The power kernel of Saini Eq. (24) as a modal transfer function; Neko's
    ! elementwise filter builds V diag(sigma) V^-1 from it for us.
    allocate(trans(lx))
    do i = 1, lx
      trans(i) = (real(i - 1, rp)/real(lx - 1, rp))**this%nsvv
    end do
    call this%filt%init_from_components(coef, "nonBoyd", trans)
    deallocate(trans)

    allocate(this%eye(lx, lx))
    this%eye = 0.0_rp
    do i = 1, lx
      this%eye(i,i) = 1.0_rp
    end do
    if (NEKO_BCKND_DEVICE .eq. 1) then
      call device_map(this%eye, this%eye_d, lx*lx)
      call device_memcpy(this%eye, this%eye_d, lx*lx, HOST_TO_DEVICE, &
           sync = .true.)
    end if

    call this%w1%init(coef%dof)
    call this%w2%init(coef%dof)
    call this%nu%init(coef%dof)
    call this%bass%init(coef%dof)
    if (this%imp) then
      call this%cg_r%init(coef%dof)
      call this%cg_p%init(coef%dof)
      call this%cg_q%init(coef%dof)
      call this%cg_z%init(coef%dof)
    end if

    tmp(1) = max(maxval(abs(coef%G12)), maxval(abs(coef%G13)), &
                 maxval(abs(coef%G23)))
    goff = glmax(tmp, 1)
    tmp(1) = max(maxval(abs(coef%G11)), maxval(abs(coef%G22)), &
                 maxval(abs(coef%G33)))
    gdiag = glmax(tmp, 1)
    if (goff .gt. 1.0e-10_rp*gdiag) call neko_error( &
         "svv_local drops G12/G13/G23, which this mesh is not entitled to")

    ! nu = c0 |c| H / N, Saini Eq. (29), with |c| the module's u_max. Their
    ! characteristic length is 2*jac^(1/d), which on this one-element-thick slab
    ! would fold in a z extent that means nothing; the element edge is what they
    ! use on a uniform mesh anyway. For re-distancing their |c| is the pointwise
    ! |w| = |sgn psi| <= 1: redistance_circles sets u_max = 1 and applies
    ! |sgn psi| as D_mu in svv_step_eq31, while the coupled files still pass the
    ! flow's u_max to svv_rd.
    h_elem = 0.0_rp
    do e = 1, nel
      h_elem = max(h_elem, abs(coef%dof%x(lx,1,1,e) - coef%dof%x(1,1,1,e)))
    end do
    tmp(1) = h_elem
    h_elem = glmax(tmp, 1)
    call field_cfill(this%nu, this%c0*u_max*h_elem/real(lx - 1, rp), n)

    ! The assembled mass, needed by the implicit split. The SEM mass matrix is
    ! diagonal, so inverting coef%Binv recovers it exactly rather than
    ! approximately.
    call copy_raw(this%bass, coef%Binv, coef%Binv_d, n)
    call field_invcol1(this%bass, n)

    this%ready = .true.

    ! --- self-checks. Two continuous pseudo-random fields test symmetry (i.e.
    ! that the tensor-product transposes above are the right way round), and a
    ! constant tests the nullspace (i.e. conservation). These run on both
    ! backends and are the main evidence that the device path is correct.
    call neko_scratch_registry%request_field(p, ind(1), .false.)
    call neko_scratch_registry%request_field(q, ind(2), .false.)
    call neko_scratch_registry%request_field(sp, ind(3), .false.)
    call neko_scratch_registry%request_field(bfld, ind(4), .false.)
    call copy_raw(bfld, coef%B, coef%B_d, n)

    do i = 1, n
      p%x(i,1,1,1) = sin(real(i, rp))
      q%x(i,1,1,1) = cos(1.7_rp*real(i, rp))
    end do
    call to_device(p, n)
    call to_device(q, n)
    call coef%gs_h%op(p, GS_OP_ADD)
    call col2_raw(p, coef%mult, coef%mult_d, n)
    call coef%gs_h%op(q, GS_OP_ADD)
    call col2_raw(q, coef%mult, coef%mult_d, n)

    call this%local(coef, p, sp)
    sym1 = field_glsc2(q, sp, n)
    ref = max(abs(field_glmax(sp, n)), abs(field_glmin(sp, n)))
    call this%local(coef, q, sp)
    sym2 = field_glsc2(p, sp, n)

    call field_rone(q, n)
    call this%local(coef, q, sp)
    nul = max(abs(field_glmax(sp, n)), abs(field_glmin(sp, n)))

    ! --- spectral radius of B^-1 S_vv by power iteration. The operator is
    ! self-adjoint in the B inner product, so the Rayleigh quotient converges --
    ! from below, which is why the iteration count is generous and the guard
    ! below leaves margin.
    do it = 1, 300
      call this%local(coef, p, sp)
      this%rho = field_glsc2(p, sp, n)/field_glsc3(p, bfld, p, n)
      call field_copy(p, sp, n)
      call coef%gs_h%op(p, GS_OP_ADD)
      call col2_raw(p, coef%Binv, coef%Binv_d, n)
      nrm = sqrt(field_glsc3(p, bfld, p, n))
      if (nrm .gt. 0.0_rp) call field_cmult(p, 1.0_rp/nrm, n)
    end do
    call neko_scratch_registry%relinquish_field(ind)

    call neko_log%section("SVV operator: " // trim(this%tag))
    write(mess, '(A,F6.2,A,F8.4)') "  N_svv = N/", this%ratio, " = ", this%nsvv
    call neko_log%message(mess)
    call logval("c0               : ", this%c0)
    call logval("nu = c0|c|H/N    : ", this%c0*u_max*h_elem/real(lx - 1, rp))
    call logval("|S_vv 1|/|S_vv p|: ", nul/ref)
    call logval("symmetry defect  : ", &
         abs(sym1 - sym2)/max(abs(sym1), abs(sym2)))
    call logval("rho(B^-1 S_vv)   : ", this%rho)
    ! dt*rho means two different things depending on the path: a stability limit
    ! when SVV rides on the explicit integrator, a condition number
    ! (kappa ~ 1 + dt*rho) when it is solved implicitly.
    if (this%imp) then
      call logval("dt*rho (kappa-1) : ", dt*this%rho)
    else
      call logval("dt*rho (limit .4): ", dt*this%rho)
    end if
    call neko_log%end_section()

    if (nul .gt. 1.0e-10_rp*ref) call neko_error( &
         "SVV operator does not annihilate constants -- it would not conserve mass")
    if (abs(sym1 - sym2) .gt. 1.0e-10_rp*max(abs(sym1), abs(sym2))) &
         call neko_error( &
         "SVV operator is not symmetric -- check the tnsr3d transposes")
    ! BDF3/EXT3 carries a real negative eigenvalue explicitly up to
    ! dt*rho = 0.952 (root of its characteristic polynomial); 0.4 leaves the
    ! rest of that budget to the advection and compression terms. The implicit
    ! split is unconditionally stable, so the limit does not apply to it.
    if (.not. this%imp .and. dt*this%rho .gt. 0.4_rp) call neko_error( &
         "SVV explicit stability dt*rho > 0.4 -- reduce the timestep or set " // &
         "the instance's implicit = true")
  end subroutine svv_init

  !> coef%mult as a field_t, so the CG's weighted inner products can go through
  !> field_math like everything else. Built once, on first use.
  subroutine ensure_mult_field(coef)
    type(coef_t), intent(inout) :: coef

    if (mult_ready) return
    call mult_f%init(coef%dof)
    call copy_raw(mult_f, coef%mult, coef%mult_d, coef%dof%size())
    mult_ready = .true.
  end subroutine ensure_mult_field

end module user
