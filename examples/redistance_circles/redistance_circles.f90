! Saini et al. (2026) section 4.4, "Re-distancing around intersecting circles":
! Eq. (44) run standalone (no flow, no phase field) from Eq. (83)'s skewed
! initial condition toward a known distance, against his Table 2 and Figs.
! 12-13. The record of the reproduction is README.md in this directory.
!
! Saini calls the distance field phi; here it is psi, the Neko scalar `psi`
! (CDI_METHOD.md section 1). The analytic distance appears only as the error
! reference and, skewed, as the problem's given initial condition.
!
! Saini's BDF2/EXT2 (Eqs. 34-35), with the three choices of his code that the
! paper does not print (README section 2): epsilon = 0.25 in Eq. (46), C(w) psi
! dealiased, and his sign guard.
!
! Seventeen routines (the SVV routines with ensure_mult_field, rd_sgn, rd_rhs,
! unit_normal, logval and the backend helpers) are byte-identical, in-body
! comments included, in the four user files: fix them in all four and check
! with examples/tools/check_shared_routines.py. rd_rhs, svv_op, svv_apply_imp
! and svv_step_imp are unused here; svv_step_eq31 is this case's own.
module user
  use neko
  use adv_dealias, only : adv_dealias_t
  use elementwise_filter, only : elementwise_filter_t
  use math, only : sqrt_inplace
  use device_math, only : device_sqrt_inplace
  implicit none

  integer, parameter :: RD_NITER_MAX = 50000

  !> In-plane only: the mesh is one element thick in z and the field is
  !> z-invariant, so a third direction would only tighten the explicit limit.
  integer, parameter :: SVV_NDIR = 2

  !> Floor on |grad psi| in the shared unit_normal. It must sit above the ~1e-19
  !> round-off gradient of a flat field (CDI_METHOD.md section 3).
  real(kind=rp), parameter :: grad_floor = 1.0e-6_rp

  !> Two circles of radius r centred at (+/-a, 0), Saini Eq. (83).
  real(kind=rp), parameter :: circ_r = 1.0_rp, circ_a = 0.7_rp

  real(kind=rp) :: eps, u_max
  logical :: mesh_ready = .false.
  real(kind=rp) :: mesh_hgll, mesh_helem, mesh_hn

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

  !> The only SVV instance: Eq. (44)'s.
  type(svv_t) :: svv_rd

  type(field_t) :: mult_f
  logical :: mult_ready = .false.

  type(field_t) :: bmass_f
  logical :: bmass_ready = .false.

  !> The exact solution, Eq. (83), and Eq. (84)'s denominator int psi_e, built
  !> once.
  type(field_t) :: psie_f
  logical :: psie_ready = .false.

  !> Saini's dealiased C(w) psi, on floor(3(N+1)/2) Gauss points per direction
  !> (his SIZE: lxd = lx1*3/2), built once.
  type(adv_dealias_t) :: adv_d
  logical :: adv_ready = .false.
  real(kind=rp) :: den_signed = 0.0_rp, den_abs = 0.0_rp, slab_dz = 1.0_rp

  real(kind=rp) :: rd_dtau, rd_tau_end
  !> case.cdi.ic = "exact" starts from psi_e instead of Eq. (83)'s skewed
  !> IC: the control that shows the scheme drifting from the exact solution.
  logical :: ic_exact = .false.
  integer :: rd_niter, report_every

contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%startup => startup
    user%initialize => initialize
    user%material_properties => material_properties
    user%initial_conditions => initial_conditions
  end subroutine user_setup

  subroutine startup(params)
    type(json_file), intent(inout) :: params
    logical :: var_dt

    call json_get(params, "case.cdi.epsilon", eps)
    if (eps .le. 0.0_rp) call neko_error("case.cdi.epsilon must be > 0")

    ! nu = c0*|c|*H/N with |c| = 1 here; Eq. (31)'s pointwise |c| = |sgn(psi)|
    ! is applied by svv_step_eq31 as D_mu.
    u_max = 1.0_rp

    ! N_svv = N/6, c0 = 2: Saini's section 4.4 values (section 4.5 uses N/4).
    call svv_rd%read_params(params, "case.cdi.redistance.svv", &
         "psi re-distancing", 6.0_rp)
    ! svv_step_eq31 is implicit and uses the CG work fields, which svv_init
    ! allocates only for an implicit instance.
    svv_rd%imp = .true.

    ! Saini fixes dtau for this study, so it is required; there is no `cfl` key.
    call json_get(params, "case.cdi.redistance.dtau", rd_dtau)
    if (rd_dtau .le. 0.0_rp) &
         call neko_error("case.cdi.redistance.dtau must be > 0")
    call json_get_or_default(params, "case.cdi.redistance.tau_end", &
         rd_tau_end, 6.0_rp)
    if (rd_tau_end .le. 0.0_rp) &
         call neko_error("case.cdi.redistance.tau_end must be > 0")

    call json_get_or_default(params, "case.cdi.report_every", report_every, 500)
    if (report_every .lt. 1) &
         call neko_error("case.cdi.report_every must be >= 1")

    block
      character(len=:), allocatable :: ic_name
      call json_get_or_default(params, "case.cdi.ic", ic_name, "skewed")
      if (trim(ic_name) .eq. "exact") then
        ic_exact = .true.
      else if (trim(ic_name) .ne. "skewed") then
        call neko_error("case.cdi.ic must be skewed or exact")
      end if
    end block

    call json_get_or_default(params, "case.time.variable_timestep", var_dt, &
         .false.)
    if (var_dt) call neko_error("variable_timestep must be false")
  end subroutine startup

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

  !> Saini Eq. (83): the exact signed distance to the union of two circles of
  !> radius r centred at (+/-a, 0), negative inside. The first branch clamps to
  !> an intersection corner (0, +/-sqrt(r^2 - a^2)) where the radial projection
  !> onto both circles lands inside the other one; that puts the kinks into the
  !> exact solution. (a-x) >= (a/r)*d1 avoids Eq. (83)'s division when d1 = 0.
  pure function circle_distance(x, y) result(d)
    real(kind=rp), intent(in) :: x, y
    real(kind=rp) :: d, d1, d2, hc

    d1 = sqrt((x - circ_a)**2 + y**2)
    d2 = sqrt((x + circ_a)**2 + y**2)
    hc = sqrt(circ_r**2 - circ_a**2)
    if ((circ_a - x) .ge. (circ_a/circ_r)*d1 .and. &
         (circ_a + x) .ge. (circ_a/circ_r)*d2) then
      d = -min(sqrt(x**2 + (y - hc)**2), sqrt(x**2 + (y + hc)**2))
    else
      d = min(d1, d2) - circ_r
    end if
  end function circle_distance

  !> Saini Eq. (83), second line: the exact distance times a factor >= 0.1,
  !> which leaves the zero isocontour where the exact solution has it.
  pure function skewed_ic(x, y) result(d)
    real(kind=rp), intent(in) :: x, y
    real(kind=rp) :: d

    d = circle_distance(x, y)
    d = ((x - 1.0_rp)**2 + (y - 1.0_rp)**2 + 0.1_rp)*d
  end function skewed_ic

  ! ------------------------------------------------------------------------
  ! Mesh, reporting
  ! ------------------------------------------------------------------------

  !> The shortest GLL spacing, the element edge and H/N, measured once and
  !> in-plane (the field is z-invariant).
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

  subroutine logval(label, v)
    character(len=*), intent(in) :: label
    real(kind=rp), intent(in) :: v
    character(len=LOG_SIZE) :: mess

    write(mess, '(A,A,E15.7)') "  ", label, v
    call neko_log%message(mess)
  end subroutine logval

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

  subroutine report_setup(dof)
    type(dofmap_t), intent(in) :: dof
    character(len=LOG_SIZE) :: mess

    call measure_mesh(dof)

    call neko_log%section("Saini 4.4: re-distancing around two circles")
    write(mess, '(A,I0)') "  polynomial order  : ", dof%Xh%lx - 1
    call neko_log%message(mess)
    call logval("H (element edge)  : ", mesh_helem)
    call logval("H/N               : ", mesh_hn)
    call logval("h_gll_min         : ", mesh_hgll)
    call logval("epsilon           : ", eps)
    call logval("xi = eps*N/H      : ", eps/mesh_hn)
    if (ic_exact) then
      call neko_log%message("  initial condition : psi_e (DIAGNOSTIC)")
    else
      call neko_log%message("  initial condition : Eq. (83), skewed")
    end if
    call svv_report(svv_rd)
    if (svv_rd%on) call neko_log%message( &
         "  SVV viscosity     : Eq. (31), D_mu = |sgn psi| left-multiplied")
    call logval("dtau              : ", rd_dtau)
    call logval("tau_end           : ", rd_tau_end)
    call logval("pseudo-CFL        : ", rd_dtau/mesh_hgll)
    write(mess, '(A,I0)') "  pseudo steps      : ", rd_niter
    call neko_log%message(mess)
    call neko_log%end_section()
  end subroutine report_setup

  ! ------------------------------------------------------------------------
  ! The normal
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

  ! ------------------------------------------------------------------------
  ! Re-distancing: Saini et al. (2026) Eq. (44)
  ! ------------------------------------------------------------------------

  !> sgn(psi) = tanh(psi/(2 eps)), Saini Eq. (46), with eps the sign-function
  !> width (0.25 in the committed cases). Neko has no device tanh, so it
  !> round-trips through the host.
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

  !> L(psi) = sgn(psi)*(1 - |grad psi|), the non-dealiased right-hand side of
  !> Saini Eq. (44) that the coupled cases use. Shared, but unused here:
  !> redistance_standalone dealiases C(w) psi instead.
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

  !> min/mean/max of |grad psi| over the whole domain, and the node count:
  !> band_grad_stats without the band (there is no phase field here).
  subroutine grad_stats(coef, gmin, gmean, gmax)
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(out) :: gmin, gmean, gmax
    type(field_t), pointer :: psifld, g1, g2, g3, w
    integer :: ind(4), i, n
    real(kind=rp) :: gsum, loc(1)

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
    gsum = 0.0_rp
    do i = 1, n
      gsum = gsum + w%x(i,1,1,1)
    end do
    loc(1) = gsum
    gsum = glsum(loc, 1)
    loc(1) = real(n, rp)
    gmean = gsum/glsum(loc, 1)
    gmin = field_glmin(w, n)
    gmax = field_glmax(w, n)

    call neko_scratch_registry%relinquish_field(ind)
  end subroutine grad_stats

  ! ------------------------------------------------------------------------
  ! The exact solution and the error norm
  ! ------------------------------------------------------------------------

  !> Build psi_e once, and with it Eq. (84)'s denominator, int psi_e as printed.
  !> Every integral is divided by the slab thickness: coef%B is a 3D mass matrix
  !> and the paper's definitions are areas (it cancels in E_r).
  subroutine ensure_psie(coef)
    type(coef_t), intent(inout) :: coef
    type(field_t), pointer :: w
    integer :: ind(1), i, n
    real(kind=rp) :: loc(1)

    if (psie_ready) return
    n = coef%dof%size()
    call ensure_bmass_field(coef)
    call psie_f%init(coef%dof)

    do i = 1, n
      psie_f%x(i,1,1,1) = circle_distance(coef%dof%x(i,1,1,1), &
           coef%dof%y(i,1,1,1))
    end do
    call to_device(psie_f, n)

    loc(1) = maxval(coef%dof%z)
    slab_dz = glmax(loc, 1)
    loc(1) = minval(coef%dof%z)
    slab_dz = slab_dz - glmin(loc, 1)

    call neko_scratch_registry%request_field(w, ind(1), .false.)
    den_signed = field_glsc2(psie_f, bmass_f, n)/slab_dz
    call field_col3(w, psie_f, psie_f, n)
    call field_sqrt(w, n)
    den_abs = field_glsc2(w, bmass_f, n)/slab_dz
    call neko_scratch_registry%relinquish_field(ind)

    psie_ready = .true.
  end subroutine ensure_psie

  !> Eq. (83) at five hand-derived points, plus the integrals, printed before
  !> the solve as a check on everything downstream.
  subroutine report_exact(coef)
    type(coef_t), intent(inout) :: coef
    character(len=LOG_SIZE) :: mess

    call neko_log%section("Eq. (83) check -- expected values in brackets")
    write(mess, '(A,F10.6,A)') "  psi_e( 0.0, 0.0) = ", &
         circle_distance(0.0_rp, 0.0_rp), "  [-0.714143]"
    call neko_log%message(mess)
    write(mess, '(A,F10.6,A)') "  psi_e( 0.0, 0.5) = ", &
         circle_distance(0.0_rp, 0.5_rp), "  [-0.214143]"
    call neko_log%message(mess)
    write(mess, '(A,F10.6,A)') "  psi_e( 1.0, 0.0) = ", &
         circle_distance(1.0_rp, 0.0_rp), "  [-0.700000]"
    call neko_log%message(mess)
    write(mess, '(A,F10.6,A)') "  psi_e(-1.2, 0.0) = ", &
         circle_distance(-1.2_rp, 0.0_rp), "  [-0.500000]"
    call neko_log%message(mess)
    write(mess, '(A,F10.6,A)') "  psi_e( 0.0, 1.5) = ", &
         circle_distance(0.0_rp, 1.5_rp), "  [ 0.655295]"
    call neko_log%message(mess)
    write(mess, '(A,F10.6,A)') "  int psi_e        = ", den_signed, &
         "  [ 3.419114]"
    call neko_log%message(mess)
    write(mess, '(A,F10.6,A)') "  int |psi_e|      = ", den_abs, &
         "  [ 7.635082]"
    call neko_log%message(mess)
    call neko_log%end_section()
  end subroutine report_exact

  !> Eq. (84), E_r = num / int psi_e, plus the numerator and psi's range. The
  !> /int|psi_e| column is kept for the log format only.
  subroutine error_report(coef, tau)
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(in) :: tau
    type(field_t), pointer :: psifld, d, w
    integer :: ind(2), n
    real(kind=rp) :: num, gmin, gmean, gmax
    character(len=LOG_SIZE) :: mess

    psifld => neko_registry%get_field('psi')
    n = psifld%size()

    call neko_scratch_registry%request_field(d, ind(1), .false.)
    call neko_scratch_registry%request_field(w, ind(2), .false.)
    call field_copy(d, psifld, n)
    call field_sub2(d, psie_f, n)
    call field_col3(w, d, d, n)
    call field_sqrt(w, n)              ! w = |psi - psi_e|
    num = field_glsc2(w, bmass_f, n)/slab_dz
    call neko_scratch_registry%relinquish_field(ind)

    call grad_stats(coef, gmin, gmean, gmax)

    write(mess, '(A,F8.4,A,E12.5)') "  E_r tau=", tau, "  num=", num
    call neko_log%message(mess)
    write(mess, '(A,E12.5,A,E12.5)') "    /int psi_e:", num/den_signed, &
         "   /int|psi_e|:", num/den_abs
    call neko_log%message(mess)
    write(mess, '(A,2(E12.5))') "    psi range :", field_glmin(psifld, n), &
         field_glmax(psifld, n)
    call neko_log%message(mess)
    write(mess, '(A,3(E11.4))') "    |grad psi|:", gmin, gmean, gmax
    call neko_log%message(mess)
  end subroutine error_report

  !> The standalone Eq. (44) relaxation to tau_end, in place, with Saini's
  !> integrator: BDF2/EXT2, Eqs. (34)-(35), the Eq. (31) SVV unsplit in the
  !> Helmholtz operator, BDF1/EXT1 on the first step:
  !>
  !>     F^m     = sgn(psi^m) - C(w^m) psi^m,   w^m = sgn(psi^m) n(psi^m)
  !>     psi_hat = (2 psi^n - psi^{n-1}/2 + dtau (2 F^n - F^{n-1})) / (3/2)
  !>     (B + dtau/(3/2) D_mu(psi^n) S_vv) psi^{n+1} = B psi_hat
  !>
  !> w is not a velocity: it is re-formed from each level's own psi, so the lag
  !> holds F^{n-1}, never w. A w decoupled from the psi it multiplies leaves
  !> growth at up to 1/(2 eps) on the zero set.
  !>
  !> As in Saini's code (README section 2), C is dealiased, so zero-set nodes
  !> can move and change sign, and his sign guard (constrainTLSR) keeps psi^n at
  !> any node whose psi^n disagrees in sign with psi_0; low-N results then
  !> depend on psi_0's round-off on the zero set.
  subroutine redistance_standalone(coef)
    type(coef_t), intent(inout) :: coef
    type(field_t), pointer :: psifld
    type(field_t), pointer :: nx, ny, nz, w, sg, cw, ref, p1, f0, f1
    type(field_t), pointer :: n0x, n0y, n0z, sgn0
    integer :: ind(14), i, it, n, cg_iters
    real(kind=rp) :: dn, loc(1), b0, beta(2), a(2)
    character(len=LOG_SIZE) :: mess

    psifld => neko_registry%get_field('psi')
    n = psifld%size()

    if (svv_rd%on .and. .not. svv_rd%ready) call svv_rd%init(coef, rd_dtau)

    call neko_scratch_registry%request_field(nx, ind(1), .false.)
    call neko_scratch_registry%request_field(ny, ind(2), .false.)
    call neko_scratch_registry%request_field(nz, ind(3), .false.)
    call neko_scratch_registry%request_field(w, ind(4), .false.)
    call neko_scratch_registry%request_field(sg, ind(5), .false.)
    call neko_scratch_registry%request_field(cw, ind(6), .false.)
    call neko_scratch_registry%request_field(ref, ind(7), .false.)
    call neko_scratch_registry%request_field(p1, ind(8), .false.)
    call neko_scratch_registry%request_field(f0, ind(9), .false.)
    call neko_scratch_registry%request_field(f1, ind(10), .false.)
    call neko_scratch_registry%request_field(n0x, ind(11), .false.)
    call neko_scratch_registry%request_field(n0y, ind(12), .false.)
    call neko_scratch_registry%request_field(n0z, ind(13), .false.)
    call neko_scratch_registry%request_field(sgn0, ind(14), .false.)
    call field_rzero(p1, n)
    call field_rzero(f1, n)

    if (.not. adv_ready) then
      call adv_d%init((3*coef%Xh%lx)/2, coef)
      adv_ready = .true.
    end if

    ! The guard's reference: the sign of psi_0, Fortran sign(1, x). Saini's
    ! code reads it off the phase field, whose signs here are psi_0's.
    call to_host(psifld, n)
    do i = 1, n
      sgn0%x(i,1,1,1) = sign(1.0_rp, psifld%x(i,1,1,1))
    end do

    ! The starting normal, for ||dn||.
    call unit_normal(coef, psifld, n0x, n0y, n0z, w)

    do it = 1, rd_niter
      if (it .eq. 1) then
        b0 = 1.0_rp
        beta = [1.0_rp, 0.0_rp]
        a = [1.0_rp, 0.0_rp]
      else
        b0 = 1.5_rp
        beta = [2.0_rp, -0.5_rp]
        a = [2.0_rp, -1.0_rp]
      end if

      ! F^n = sgn(psi^n) - C(w^n) psi^n, w^n from psi^n itself
      call unit_normal(coef, psifld, nx, ny, nz, w)
      call rd_sgn(psifld, sg, n)
      call field_col2(nx, sg, n)
      call field_col2(ny, sg, n)
      call field_col2(nz, sg, n)
      call field_rzero(cw, n)
      call adv_d%compute_scalar(nx, ny, nz, psifld, cw, coef%Xh, coef, n)
      call field_cmult(cw, -1.0_rp, n)
      call coef%gs_h%op(cw, GS_OP_ADD)
      call col2_raw(cw, coef%Binv, coef%Binv_d, n)
      call field_sub3(f0, sg, cw, n)

      ! psi_hat in place; psi^n kept in ref, for D_mu and the lag
      call field_copy(ref, psifld, n)
      call field_cmult(psifld, beta(1), n)
      if (it .gt. 1) call field_add2s2(psifld, p1, beta(2), n)
      call field_add2s2(psifld, f0, rd_dtau*a(1), n)
      if (it .gt. 1) call field_add2s2(psifld, f1, rd_dtau*a(2), n)
      call field_cmult(psifld, 1.0_rp/b0, n)

      if (svv_rd%on) call svv_step_eq31(coef, psifld, ref, rd_dtau/b0, &
           cg_iters)

      ! Saini's sign guard: a node whose psi^n disagrees in sign with psi_0
      ! keeps psi^n.
      call to_host(psifld, n)
      call to_host(ref, n)
      do i = 1, n
        if (sign(1.0_rp, ref%x(i,1,1,1))*sgn0%x(i,1,1,1) .lt. 0.0_rp) &
             psifld%x(i,1,1,1) = ref%x(i,1,1,1)
      end do
      call to_device(psifld, n)

      call field_copy(p1, ref, n)
      call field_copy(f1, f0, n)

      ! A trace, not just an endpoint: drift past tau = 6 shows only while it
      ! happens.
      if (mod(it, report_every) .eq. 0) &
           call error_report(coef, real(it, rp)*rd_dtau)
    end do

    call unit_normal(coef, psifld, nx, ny, nz, w)
    call to_host(nx, n)
    call to_host(ny, n)
    call to_host(nz, n)
    call to_host(n0x, n)
    call to_host(n0y, n)
    call to_host(n0z, n)
    dn = 0.0_rp
    do i = 1, n
      dn = dn + (nx%x(i,1,1,1) - n0x%x(i,1,1,1))**2 &
              + (ny%x(i,1,1,1) - n0y%x(i,1,1,1))**2 &
              + (nz%x(i,1,1,1) - n0z%x(i,1,1,1))**2
    end do
    loc(1) = dn
    dn = sqrt(glsum(loc, 1))
    write(mess, '(A,E11.4)') "  ||dn|| over the solve (global) = ", dn
    call neko_log%message(mess)

    call neko_scratch_registry%relinquish_field(ind)
  end subroutine redistance_standalone

  !> The implicit part of the BDF step: the SVV with Saini's Eq. (31) viscosity
  !> as printed (Eqs. 7, 33, 35), the *assembled* operator left-multiplied at the
  !> GLL points by D_mu = |sgn(psi^n)|, taken from `sref`, the constant c0*H/N
  !> staying in nu:
  !>
  !>     (B + dt*D_mu*S_vv) s_new = B s_old,   dt = dtau/b0
  !>
  !> D_mu vanishes on the zero set, so this SVV cannot move a node on it; the
  !> shared svv_step_imp, with nu inside the bilinear form, does (CDI_METHOD.md
  !> section 5). It does not conserve mass (1^T D_mu S_vv /= 0); a distance does
  !> not need to.
  !>
  !> The operator is not symmetric. s_new = s_old + D_mu z makes it so,
  !>
  !>     (D_mu B + dt D_mu S_vv D_mu) z = -dt D_mu S_vv s_old,
  !>
  !> semi-definite but consistent: rows where D_mu = 0 read 0 = 0, and the
  !> preconditioner is zeroed there so z stays 0. This residual *is* the
  !> unscaled one, B (s_new - s_old) + dt D_mu S_vv s_new.
  subroutine svv_step_eq31(coef, s, sref, dt, iters)
    type(coef_t), intent(inout) :: coef
    type(field_t), intent(inout) :: s, sref
    real(kind=rp), intent(in) :: dt
    integer, intent(out) :: iters
    type(field_t), pointer :: dmu, pinv, t
    real(kind=rp) :: rz, rz_new, alpha, beta, r0
    integer :: ind(3), i, n, it
    integer, parameter :: MAXIT = 200
    real(kind=rp), parameter :: RTOL = 1.0e-12_rp

    n = coef%dof%size()
    call neko_scratch_registry%request_field(dmu, ind(1), .false.)
    call neko_scratch_registry%request_field(pinv, ind(2), .false.)
    call neko_scratch_registry%request_field(t, ind(3), .false.)

    ! D_mu, and the mass preconditioner 1/(D_mu B), zeroed where D_mu = 0
    call to_host(sref, n)
    do i = 1, n
      dmu%x(i,1,1,1) = abs(tanh(sref%x(i,1,1,1)/(2.0_rp*eps)))
      if (dmu%x(i,1,1,1) .ge. tiny(1.0_rp)) then
        pinv%x(i,1,1,1) = 1.0_rp/dmu%x(i,1,1,1)
      else
        pinv%x(i,1,1,1) = 0.0_rp
      end if
    end do
    call to_device(dmu, n)
    call to_device(pinv, n)
    call field_invcol2(pinv, svv_rd%bass, n)

    ! r = -dt D_mu S_vv s_old, the residual at z = 0
    call svv_rd%local(coef, s, svv_rd%cg_r)
    call coef%gs_h%op(svv_rd%cg_r, GS_OP_ADD)
    call field_cmult(svv_rd%cg_r, -dt, n)
    call field_col2(svv_rd%cg_r, dmu, n)
    r0 = sqrt(field_glsc3(svv_rd%cg_r, mult_f, svv_rd%cg_r, n))

    call field_col3(svv_rd%cg_z, svv_rd%cg_r, pinv, n)
    call field_copy(svv_rd%cg_p, svv_rd%cg_z, n)
    rz = field_glsc3(svv_rd%cg_r, mult_f, svv_rd%cg_z, n)

    iters = 0
    if (r0 .gt. 0.0_rp) then
      do it = 1, MAXIT
        ! q = D_mu (B p + dt S_vv t), t = D_mu p, which is also s's increment
        call field_col3(t, dmu, svv_rd%cg_p, n)
        call svv_rd%local(coef, t, svv_rd%cg_q)
        call coef%gs_h%op(svv_rd%cg_q, GS_OP_ADD)
        call field_cmult(svv_rd%cg_q, dt, n)
        call field_addcol3(svv_rd%cg_q, svv_rd%bass, svv_rd%cg_p, n)
        call field_col2(svv_rd%cg_q, dmu, n)

        alpha = rz/field_glsc3(svv_rd%cg_q, mult_f, svv_rd%cg_p, n)
        call field_add2s2(s, t, alpha, n)
        call field_add2s2(svv_rd%cg_r, svv_rd%cg_q, -alpha, n)
        iters = it
        if (sqrt(field_glsc3(svv_rd%cg_r, mult_f, svv_rd%cg_r, n)) &
             .le. RTOL*r0) exit
        call field_col3(svv_rd%cg_z, svv_rd%cg_r, pinv, n)
        rz_new = field_glsc3(svv_rd%cg_r, mult_f, svv_rd%cg_z, n)
        beta = rz_new/rz
        rz = rz_new
        call field_add2s1(svv_rd%cg_p, svv_rd%cg_z, beta, n)
      end do
    end if

    call neko_scratch_registry%relinquish_field(ind)
  end subroutine svv_step_eq31

  ! ------------------------------------------------------------------------
  ! Spectral vanishing viscosity
  ! ------------------------------------------------------------------------

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

  !> Build the operator once, check it (symmetry, nullspace, spectral radius),
  !> and refuse to run if an explicit instance is too stiff. The one instance
  !> here is implicit, used by svv_step_eq31.
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

  !> coef%B as a field_t, so the mass integral is one field_glsc2. Built once,
  !> on first use.
  subroutine ensure_bmass_field(coef)
    type(coef_t), intent(inout) :: coef

    if (bmass_ready) return
    call bmass_f%init(coef%dof)
    call copy_raw(bmass_f, coef%B, coef%B_d, coef%dof%size())
    bmass_ready = .true.
  end subroutine ensure_bmass_field

  ! ------------------------------------------------------------------------
  ! Hooks
  ! ------------------------------------------------------------------------

  !> The whole case runs here, before the first timestep and the first output
  !> (src/simulation.f90 calls user%initialize, then output_controller%execute),
  !> so the one frame on disk is the relaxed field. Not in initial_conditions:
  !> makeneko binds neko_user_access only after neko_init returns. The timestep
  !> that follows is a formality (frozen fluid, psi's diffusivity 1e-16).
  subroutine initialize(time)
    type(time_state_t), intent(in) :: time
    type(coef_t), pointer :: coef

    coef => neko_user_access%case%fluid%c_Xh
    call measure_mesh(coef%dof)
    rd_niter = ceiling(rd_tau_end/rd_dtau - 1.0e-9_rp)
    if (rd_niter .gt. RD_NITER_MAX) call neko_error( &
         "too many pseudo steps -- raise case.cdi.redistance.dtau")

    call ensure_psie(coef)
    call report_setup(coef%dof)
    call report_exact(coef)

    call neko_log%section("Eq. (44) relaxation")
    call error_report(coef, 0.0_rp)
    call redistance_standalone(coef)
    call error_report(coef, rd_tau_end)
    call neko_log%end_section()
  end subroutine initialize

  subroutine material_properties(scheme_name, properties, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: properties
    type(time_state_t), intent(in) :: time

    if (scheme_name .eq. "fluid") then
      call field_cfill(properties%get("fluid_rho"), 1.0_rp)
      call field_cfill(properties%get("fluid_mu"), 1.0_rp)
    else if (scheme_name .eq. "psi") then
      call field_cfill(properties%get('psi_cp'), 1.0_rp)
      call field_cfill(properties%get('psi_lambda'), 1.0e-16_rp)
    end if
  end subroutine material_properties

  !> psi starts as Eq. (83)'s psi_0, the skewed exact distance (or psi_e with
  !> case.cdi.ic = "exact", the drift control): the problem's given initial
  !> condition, which the solve in initialize() relaxes.
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(field_t), pointer :: f
    integer :: i

    if (scheme_name .ne. 'psi') return
    f => fields%items(1)%ptr

    do i = 1, f%dof%size()
      if (ic_exact) then
        f%x(i,1,1,1) = circle_distance(f%dof%x(i,1,1,1), f%dof%y(i,1,1,1))
      else
        f%x(i,1,1,1) = skewed_ic(f%dof%x(i,1,1,1), f%dof%y(i,1,1,1))
      end if
    end do

    if ((NEKO_BCKND_DEVICE .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
         .or. (NEKO_BCKND_OPENCL .eq. 1)) then
      call device_memcpy(f%x, f%x_d, f%dof%size(), HOST_TO_DEVICE, sync = .true.)
    end if
  end subroutine initial_conditions

end module user
