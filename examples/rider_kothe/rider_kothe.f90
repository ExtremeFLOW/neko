! Rider-Kothe vortex-in-a-box (Saini et al. 2026, section 4.5): a disk is
! stretched into a thin filament, the flow reverses at t = T/2, and the disk
! should return to its initial shape at t = T. Neko's CDI compression term takes
! its interface normal from a transported signed distance (`normal = "psi"`).
! The only case with genuine strain. Configuration and results: README.md in
! this directory.
!
! phi is the phase field (Neko scalar `s`), psi the signed distance. Saini's
! names are the other way round (CDI_METHOD.md section 1).
!
! Method parameters live in `case.cdi`, not under a scalar: `case.scalar` and
! `case.scalars` are mutually exclusive paths in src/case.f90.
!
! Seventeen routines (the SVV routines with ensure_mult_field, rd_sgn, rd_rhs,
! unit_normal, logval and the backend helpers) are byte-identical, in-body
! comments included, in the four user files: fix them in all four and check
! with examples/tools/check_shared_routines.py. redistance, band_grad_stats and
! initialize are also identical across the three coupled files.
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

  !> Floor on |grad psi| in unit_normal (case.cdi.grad_floor). It must stay
  !> above the ~1e-19 round-off gradient of a flat field: below it a flat psi (a
  !> built psi is flat outside its band) gives a random unit normal, whose
  !> div(n) ~ 1/h makes the compression term an exponential source on phi
  !> (CDI_METHOD.md section 3).
  real(kind=rp), parameter :: GRAD_FLOOR_DEFAULT = 1.0e-6_rp
  real(kind=rp) :: grad_floor = GRAD_FLOOR_DEFAULT

  !> The disk of Saini section 4.5: centre (0.5, 0.75), radius 0.15.
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

  !> In-plane only: the mesh is one element thick in z and the field is
  !> z-invariant, so a third direction would only tighten the explicit limit.
  integer, parameter :: SVV_NDIR = 2

  !> Spectral vanishing viscosity, Saini Eqs. (24)-(29): one instance per
  !> equation, each with its own c0 and N_svv.
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

  type(svv_t) :: svv_psi    ! the psi transport equation, Saini Eq. (43)
  type(svv_t) :: svv_rd     ! the psi re-distancing equation, Saini Eq. (44)

  type(field_t) :: mult_f
  logical :: mult_ready = .false.

  !> The spatial part of the velocity, built once; each step rescales it by
  !> cos(pi t/T), which also keeps it resident on a device.
  type(field_t) :: u0, v0

contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%startup => startup
    user%initialize => initialize
    user%preprocess => preprocess
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
    u_max = 1.0_rp   ! until preprocess() measures it
    u_max0 = 1.0_rp

    call json_get_or_default(params, "case.cdi.normal", str, "psi")
    select case (trim(str))
    case ("phi")
      normal_kind = NORMAL_PHI
    case ("psi")
      normal_kind = NORMAL_PSI
    case default
      call neko_error("case.cdi.normal must be 'phi' or 'psi'")
    end select

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
    ! psi's transport always carries SVV, as Saini's Eq. (43) does; phi's
    ! equation has its own, physical, diffusion.
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
    ! An explicit pseudo timestep, overriding cfl*h_gll_min when > 0: e.g.
    ! Saini's H/(N+1), a pseudo-CFL above 1 from N = 4 that `cfl` cannot express.
    call json_get_or_default(params, "case.cdi.redistance.dtau", rd_dtau_set, &
         0.0_rp)
    if (rd_dt_tls .le. 0.0_rp) &
         call neko_error("case.cdi.redistance.dt_tls must be > 0")
    if (rd_cfl .le. 0.0_rp .or. rd_cfl .gt. 1.0_rp) &
         call neko_error("case.cdi.redistance.cfl must be in (0, 1]")
    ! The t = 0 field (exact, or built in initialize) stands in for an event at
    ! t = 0, so the first event is one interval in.
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

  subroutine field_sqrt(a, n)
    type(field_t), intent(inout) :: a
    integer, intent(in) :: n

    if (NEKO_BCKND_DEVICE .eq. 1) then
      call device_sqrt_inplace(a%x_d, n)
    else
      call sqrt_inplace(a%x, n)
    end if
  end subroutine field_sqrt

  subroutine to_host(a, n)
    type(field_t), intent(inout) :: a
    integer, intent(in) :: n

    if (NEKO_BCKND_DEVICE .eq. 1) &
         call device_memcpy(a%x, a%x_d, n, DEVICE_TO_HOST, sync = .true.)
  end subroutine to_host

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
  !> The plane distance jumps across the y = 0/1 seam, 0.10 from the disk. phi
  !> saturates there, but a transported psi would get a large artificial
  !> gradient on the seam.
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
  !> cos(pi t/T), Saini Eqs. (85)-(86), prescribed every step at t_n =
  !> time%tlag(1): the scalar step applies advection and compression to s^n and
  !> extrapolates them. Here, not in compute(), so that step 1 sees u(0);
  !> case.fluid.freeze = true stops the flow solver overwriting it.
  subroutine preprocess(time)
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, v
    type(coef_t), pointer :: coef
    integer :: i, n
    real(kind=rp) :: x, y, tfac, tmp(1)

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

    ! |cos(pi t/T)| <= 1: the peak speed is at t = 0, so the step-1 report
    ! covers the whole run.
    tfac = cos(pi*time%tlag(1)/RK_T)
    call field_cmult2(u, u0, tfac, n)
    call field_cmult2(v, v0, tfac, n)
    u_max = abs(tfac)*u_max0

    if (.not. reported) then
      call report_resolution(coef%dof, time%dt)
      reported = .true.
    end if
  end subroutine preprocess

  !> After the scalar step: the split SVV step, the re-distancing events and
  !> the band report.
  subroutine compute(time)
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: psifld
    type(coef_t), pointer :: coef
    integer :: nband, cg_iters
    real(kind=rp) :: gmin, gmean, gmax
    character(len=LOG_SIZE) :: mess

    coef => neko_user_access%case%fluid%c_Xh

    ! An explicit SVV instance is initialised by its own source-term hook, so
    ! without a `source_terms` entry it never fires while the header reports it.
    if (time%tstep .gt. 1 .and. svv_psi%on .and. .not. svv_psi%imp &
         .and. .not. svv_psi%ready) call neko_error( &
         "case.cdi.svv_psi is on but the psi source term never fired -- the " &
         // "psi scalar needs a user source_terms entry in the case file")

    ! After slag%update(), so changing psi here is a clean Lie split: the next
    ! step's history and the output both see it.
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

    ! Band |grad psi| over time, in two lines: one would overrun LOG_SIZE.
    if (mod(time%tstep, report_every) .eq. 0) then
      call band_grad_stats(coef, gmin, gmean, gmax, nband)
      write(mess, '(A,F9.4,A,I0)') "  |grad psi| t=", time%t, &
           "  band nodes ", nband
      call neko_log%message(mess)
      write(mess, '(A,3(E11.4))') "    min/mean/max:", gmin, gmean, gmax
      call neko_log%message(mess)
    end if
  end subroutine compute

  !> psi_init = "redistance": Saini's Algorithm 1 line 2, building psi from phi
  !> by the Eq. (44) solve instead of an analytic distance (REDISTANCING.md
  !> section 4). Not in initial_conditions: makeneko binds neko_user_access only
  !> after neko_init returns, so coef is unreachable there. Step 1 is BDF1, so
  !> psi's time history starts from the field built here.
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

  !> The shortest GLL spacing, the element edge and H/N, measured once and
  !> in-plane (the field is z-invariant). Separate from report_resolution
  !> because initialize needs them before step 1.
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

  !> n = grad(src)/max(|grad(src)|, grad_floor), the gradient made continuous
  !> first (gather-scatter the sum, divide by the node multiplicity). `w` returns
  !> the floored |grad(src)|.
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

    ! Not field_vdot3: its result is intent(out) on a field_t, which
    ! deallocates the field's storage on entry.
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

    ! w = phi*(phi - 1). Not col3 + sub2: `s` is intent(in) here, and
    ! field_sub2's second argument is intent(inout).
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

  !> min/mean/max of |grad psi| over the band phi(1-phi) > 1e-4, and its size.
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
  ! ------------------------------------------------------------------------

  !> sgn(psi) = tanh(psi/(2 eps)), Saini Eq. (46). Neko has no device tanh, so
  !> it round-trips through the host; it runs only inside re-distancing events.
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

  !> L(psi) = sgn(psi)*(1 - |grad psi|), the right-hand side of Saini Eq. (44):
  !> with n = grad(psi)/|grad psi|, w.grad(psi) = sgn(psi)*|grad psi|, so this is
  !> his -w.grad(psi) + sgn(psi) on the GLL points (REDISTANCING.md section 5.3).
  !> sgn vanishes on the zero set, which holds it in place.
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

  !> One re-distancing event: optionally reseed psi from phi (Saini Eq. 47),
  !> then rd_niter pseudo-steps of Eq. (44), each followed by an implicit SVV
  !> step (a Lie split). SSP-RK3, not explicit Euler: linearised, Eq. (44) is
  !> advection along the normal, whose collocation spectrum is imaginary; Euler
  !> amplifies it at every wavenumber, RK3 covers the imaginary axis to 1.73.
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
      ! dtau from the pseudo-CFL on h_gll_min unless redistance.dtau sets it;
      ! the step count covers band*H of pseudo-time. Saini's H/(N+1) is safe only
      ! with his sign width 0.25: at width eps it fails under repeated events
      ! (advecting_slab_1d/README.md section 3.5).
      if (rd_dtau_set .gt. 0.0_rp) then
        rd_dtau = rd_dtau_set
      else
        rd_dtau = rd_cfl*mesh_hgll
      end if
      ! The 1e-9 stops an exactly divisible band (Saini's H/(N+1)) from
      ! rounding up to one step too many.
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
      ! The normal before the event, for ||dn||.
      call band_grad_stats(coef, gmin0, gmean0, gmax0, nband)
      call unit_normal(coef, psifld, nx, ny, nz, w)
    end if

    ! Saini Eq. (47): restart from the phase field, so phi's interface errors
    ! enter psi here. seed = 'psi' relaxes the transported field in place
    ! instead, which cannot re-register psi to phi.
    if (use_seed .eq. RD_SEED_PHI) then
      call field_copy(psifld, s, n)
      call field_cadd(psifld, -0.5_rp, n)
      call field_cmult(psifld, rd_rf, n)
    end if

    ! The t = 0 build reports the Eq. (47) seed as its "before" row; no ||dn||.
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
    real(kind=rp) :: lam

    if (scheme_name .eq. "fluid") then
      call field_cfill(properties%get("fluid_rho"), 1.0_rp)
      call field_cfill(properties%get("fluid_mu"), 1.0_rp)
    else if (scheme_name .eq. "s") then
      call field_cfill(properties%get('s_cp'), 1.0_rp)
      if (gamma .gt. 0.0_rp) then
        ! Physical, not stabilisation (CDI_METHOD.md section 2). Both halves of
        ! the balance follow u_max = |cos(pi t/T)|, so the width is unchanged;
        ! the implicit diffusion reads it at t_{n+1}. The solve reads
        ! s_lambda_tot, which Neko copies from s_lambda only at init unless a
        ! turbulence model is set: fill both, or the diffusion stays at t = 0.
        lam = eps*gamma*u_max0*abs(cos(pi*time%t/RK_T))
        call field_cfill(properties%get('s_lambda'), lam)
        call field_cfill(neko_registry%get_field('s_lambda_tot'), lam)
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
        ! Never evaluate an analytic distance here: a real application has
        ! none. initialize() builds psi by Eq. (44) before the first step.
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
  !> nu sits inside the bilinear form, constant per instance. Saini's Eqs. (31)
  !> and (33) left-multiply the assembled operator by a pointwise D_mu instead;
  !> the two agree only for constant nu (svv_step_eq31 in redistance_circles.f90
  !> is the printed form).
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
  !> that carries everything else: (B + dt*S_vv) s_new = B s_old. Unconditionally
  !> stable and first order in dt; mass is conserved exactly (1^T S_vv = 0).
  !> Mass-preconditioned CG, condition number 1 + dt*rho.
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

    ! r = B*s_old - A(s_old), with s_old itself as the initial guess
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
  !> An explicit instance goes through the source term, so BDF3/EXT3 carries it
  !> on the real axis and dt*rho(B^-1 S_vv) is the binding limit.
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

    ! Saini Eq. (24)'s power kernel as a modal transfer function; Neko's
    ! elementwise filter builds V diag(sigma) V^-1 from it.
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

    ! nu = c0 |c| H / N, Saini Eq. (29), with |c| the module's u_max when the
    ! instance is built. H is the element edge: his 2*jac^(1/d) with d = 3 would
    ! fold in this one-element-thick slab's z extent.
    h_elem = 0.0_rp
    do e = 1, nel
      h_elem = max(h_elem, abs(coef%dof%x(lx,1,1,e) - coef%dof%x(1,1,1,e)))
    end do
    tmp(1) = h_elem
    h_elem = glmax(tmp, 1)
    call field_cfill(this%nu, this%c0*u_max*h_elem/real(lx - 1, rp), n)

    ! The assembled mass for the implicit step: diagonal, so exactly 1/Binv.
    call copy_raw(this%bass, coef%Binv, coef%Binv_d, n)
    call field_invcol1(this%bass, n)

    this%ready = .true.

    ! Self-checks, on both backends: two pseudo-random fields test symmetry
    ! (the tnsr3d transposes), a constant tests the nullspace (conservation).
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

    ! rho(B^-1 S_vv) by power iteration. The operator is self-adjoint in the B
    ! inner product, so the Rayleigh quotient converges from below; hence 300.
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
    ! dt*rho: a stability limit when explicit, kappa - 1 when implicit.
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
    ! dt*rho = 0.952; 0.4 leaves the rest to advection and compression.
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
