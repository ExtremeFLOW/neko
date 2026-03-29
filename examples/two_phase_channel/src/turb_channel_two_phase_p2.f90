! 3D turbulent two-phase channel flow with CSF surface tension — Phase 2
!
! *** NOTE (2026-03-29): the extra GS on n_hat below is a NO-OP. ***
!
!   The first GS pass (on grad(phi), before normalisation) already makes all
!   shared face nodes identical between neighbouring elements.  After pointwise
!   normalisation, n_hat at face nodes is therefore also identical on both sides.
!   A second GS on n_hat sums these identical values and divides back: unchanged.
!
!   The actual kink is INTRA-ELEMENT: the Lagrange polynomial of n_hat in the
!   face-normal direction connects an averaged face-node value to an
!   element-local first-interior value.  GS cannot reach interior nodes.
!   D[N,N] ~= 14 (Lagrange endpoint derivative, N=7) amplifies this
!   intra-element jump -> kappa_rms >> 2/R regardless of the extra GS pass.
!
!   For sigma=0 runs the CSF force = sigma*kappa*grad(phi) = 0, so even if
!   the GS did smooth n_hat, it would have no effect on those runs.
!
!   Retained here for reference.  A true fix requires acting on interior nodes,
!   e.g. repeated GS on grad(phi) before normalisation (propagating the average
!   inward one ring per pass), or a Helmholtz-type smoother on n_hat.
!
!   Phase 2 case files (mesh 108x18x36, eps=0.09, R=0.4):
!     turb_channel_two_phase_p2_sigma0.case  -- (no-op for sigma=0; same as standard)
!     turb_channel_two_phase_p2_we10.case    -- We=10  (sigma=0.04)
!     turb_channel_two_phase_p2_we1.case     -- We=1   (sigma=0.4)
!
! Background turbulent channel at Re_tau=180, Re_b=2800 with a single drop.
! Drop is initialized with a tanh profile; CDI source term keeps the
! interface sharp; CSF surface tension acts on the fluid.
!
! Domain: 4pi x 2 x 4/3pi  (x x y x z), walls at y = +/-1
! Driving: flow_rate_force in x, value=1 (constant bulk velocity U_b=1)
! Re_tau = 180, Re_b = 2800, matched density/viscosity (rho=mu from Re)
!
! Physical parameters read from case file:
!   case.scalar.epsilon     -- interface thickness
!   case.scalar.gamma       -- CDI mobility coefficient
!   case.scalar.sigma       -- surface tension coefficient
!   case.scalar.drop_radius   -- drop radius
!   case.scalar.drop_center_y -- wall-normal drop centre offset from centreline (default 0)
!   case.fluid.ipostproc      -- diagnostic output frequency (tsteps)
!
module user
  use neko
  implicit none

  ! Global user variables
  type(field_t) :: w1
  type(field_t) :: temp8, temp9
  type(field_t) :: vort1, vort2, vort3

  type(file_t) :: output_file
  type(vector_t) :: vec_out
  integer :: ipostproc

  logical :: turbulent_ic

  real(kind=rp) :: eps
  real(kind=rp) :: gamma, u_max, sigma
  real(kind=rp) :: drop_radius
  real(kind=rp) :: drop_center_y
  real(kind=rp) :: Re_b

contains

  ! Register user-defined functions
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%startup => startup
    user%initialize => initialize
    user%compute => compute
    user%finalize => finalize
    user%source_term => source_term
    user%material_properties => material_properties
    user%initial_conditions => initial_conditions
  end subroutine user_setup

  ! Read simulation parameters from case file
  subroutine startup(params)
    type(json_file), intent(inout) :: params
    character(len=50) :: mess

    call json_get(params, "case.fluid.ipostproc", ipostproc)
    write(mess, *) "postprocessing steps : ", ipostproc
    call neko_log%message(mess)

    call json_get_or_default(params, "case.fluid.turbulent_ic", turbulent_ic, .true.)

    call json_get(params, "case.fluid.Re", Re_b)
    call json_get(params, "case.scalar.epsilon", eps)
    call json_get(params, "case.scalar.gamma", gamma)
    call json_get(params, "case.scalar.sigma", sigma)
    call json_get(params, "case.scalar.drop_radius", drop_radius)
    call json_get_or_default(params, "case.scalar.drop_center_y", drop_center_y, 0.0_rp)
    u_max = 1.0_rp
  end subroutine startup

  ! Initialization: set up diagnostic output and work arrays
  subroutine initialize(time)
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: u

    call output_file%init("ekin.csv", overwrite=.true.)
    call output_file%set_header("# t,Ekin,enst,u_max," // &
         "kappa_max,kappa_min,kappa_rms,Fst_max,phi_min,phi_max")
    call vec_out%init(10)

    u => neko_field_registry%get_field("u")
    call w1%init(u%dof, 'work1')
    call temp8%init(u%dof)
    call temp9%init(u%dof)
    call vort1%init(u%dof)
    call vort2%init(u%dof)
    call vort3%init(u%dof)

    call compute(time)

  end subroutine initialize

  ! Diagnostic output: kinetic energy, enstrophy, CSF diagnostics
  subroutine compute(time)
    type(time_state_t), intent(in) :: time

    integer :: ntot, i
    real(kind=rp) :: ekin, enst
    type(field_t), pointer :: u, v, w
    type(coef_t), pointer :: coef

    real(kind=rp) :: kappa_max, kappa_min, kappa_rms
    real(kind=rp) :: kappa_weight_sum, Fst_max, phi_min, phi_max
    real(kind=rp) :: absgrad
    type(field_t), pointer :: s
    type(field_t), pointer :: t1, t2, t3, t4, t5, t6, t7
    integer :: ind_diag(7)

    if (mod(time%tstep, ipostproc) .ne. 0) return

    coef => neko_user_access%case%fluid%c_Xh
    u => neko_field_registry%get_field("u")
    v => neko_field_registry%get_field("v")
    w => neko_field_registry%get_field("w")

    ntot = u%dof%size()

    call field_col3(w1, u, u, ntot)
    call field_addcol3(w1, v, v, ntot)
    call field_addcol3(w1, w, w, ntot)

    if (NEKO_BCKND_DEVICE .eq. 1) then
      call device_memcpy(w1%x, w1%x_d, w1%size(), &
           DEVICE_TO_HOST, sync=.true.)
    end if
    u_max = sqrt(glmax(w1%x, ntot))

    if (NEKO_BCKND_DEVICE .eq. 1) then
      ekin = 0.5_rp * device_glsc2(w1%x_d, coef%B_d, ntot) / coef%volume
    else
      ekin = 0.5_rp * glsc2(w1%x, coef%B, ntot) / coef%volume
    end if

    call curl(vort1, vort2, vort3, u, v, w, temp8, temp9, coef)
    call field_col3(w1, vort1, vort1, ntot)
    call field_addcol3(w1, vort2, vort2, ntot)
    call field_addcol3(w1, vort3, vort3, ntot)
    if (NEKO_BCKND_DEVICE .eq. 1) then
      enst = 0.5_rp * device_glsc2(w1%x_d, coef%B_d, ntot) / coef%volume
    else
      enst = 0.5_rp * glsc2(w1%x, coef%B, ntot) / coef%volume
    end if

    ! ==================== CSF DIAGNOSTICS ====================
    s => neko_field_registry%get_field('phase')

    call neko_scratch_registry%request_field(t1, ind_diag(1))
    call neko_scratch_registry%request_field(t2, ind_diag(2))
    call neko_scratch_registry%request_field(t3, ind_diag(3))
    call neko_scratch_registry%request_field(t4, ind_diag(4))
    call neko_scratch_registry%request_field(t5, ind_diag(5))
    call neko_scratch_registry%request_field(t6, ind_diag(6))
    call neko_scratch_registry%request_field(t7, ind_diag(7))

    call grad(t1%x, t2%x, t3%x, s%x, coef)
    call coef%gs_h%op(t1, GS_OP_ADD)
    call coef%gs_h%op(t2, GS_OP_ADD)
    call coef%gs_h%op(t3, GS_OP_ADD)
    call col2(t1%x, coef%mult, ntot)
    call col2(t2%x, coef%mult, ntot)
    call col2(t3%x, coef%mult, ntot)

    call copy(t5%x, t1%x, ntot)
    call copy(t6%x, t2%x, ntot)
    call copy(t7%x, t3%x, ntot)

    do i = 1, ntot
      absgrad = sqrt(t1%x(i,1,1,1)**2 + t2%x(i,1,1,1)**2 + t3%x(i,1,1,1)**2)
      if (absgrad < 1.0e-12_rp) then
        t1%x(i,1,1,1) = 0.0_rp
        t2%x(i,1,1,1) = 0.0_rp
        t3%x(i,1,1,1) = 0.0_rp
      else
        t1%x(i,1,1,1) = t1%x(i,1,1,1) / absgrad
        t2%x(i,1,1,1) = t2%x(i,1,1,1) / absgrad
        t3%x(i,1,1,1) = t3%x(i,1,1,1) / absgrad
      end if
    end do

    call div(t4%x, t1%x, t2%x, t3%x, coef)
    call coef%gs_h%op(t4, GS_OP_ADD)
    call col2(t4%x, coef%mult, ntot)
    call cmult(t4%x, -1.0_rp, ntot)

    kappa_max = glmax(t4%x, ntot)
    kappa_min = glmin(t4%x, ntot)

    do i = 1, ntot
      w1%x(i,1,1,1) = sqrt(t5%x(i,1,1,1)**2 + t6%x(i,1,1,1)**2 &
                          + t7%x(i,1,1,1)**2)
    end do
    do i = 1, ntot
      t1%x(i,1,1,1) = t4%x(i,1,1,1)**2 * w1%x(i,1,1,1)
    end do
    kappa_rms = glsc2(t1%x, coef%B, ntot)
    kappa_weight_sum = glsc2(w1%x, coef%B, ntot)
    if (kappa_weight_sum > 0.0_rp) then
      kappa_rms = sqrt(kappa_rms / kappa_weight_sum)
    else
      kappa_rms = 0.0_rp
    end if

    do i = 1, ntot
      w1%x(i,1,1,1) = sigma * abs(t4%x(i,1,1,1)) * w1%x(i,1,1,1)
    end do
    Fst_max = glmax(w1%x, ntot)

    phi_min = glmin(s%x, ntot)
    phi_max = glmax(s%x, ntot)

    call neko_scratch_registry%relinquish_field(ind_diag)
    ! ==================== END CSF DIAGNOSTICS ====================

    call neko_log%message("Writing csv file")
    vec_out%x = (/ekin, enst, u_max, kappa_max, kappa_min, kappa_rms, &
                  Fst_max, phi_min, phi_max, 0.0_rp/)
    call output_file%write(vec_out, time%t)

  end subroutine compute

  subroutine finalize(time)
    type(time_state_t), intent(in) :: time

    call w1%free()
    call file_free(output_file)
    call vec_out%free

  end subroutine finalize

  ! Source terms: CDI for phase field, CSF for fluid momentum
  subroutine source_term(scheme_name, rhs, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: rhs
    type(time_state_t), intent(in) :: time

    integer :: i
    type(field_t), pointer :: rhs_s, s

    type(field_t), pointer :: u, v, w, rhs_u, rhs_v, rhs_w
    type(field_t), pointer :: temp1, temp2, temp3, temp4, temp5, temp6, temp7

    real(kind=rp) :: absgrad
    integer :: ind(7), ind_work(4)
    type(field_t), pointer :: work1, work2, work3, work4
    type(coef_t), pointer :: coef

    ! CDI source term: keeps the interface sharp
    if (scheme_name .eq. 'phase') then

      rhs_s => rhs%get_by_index(1)
      s => neko_field_registry%get_field('phase')
      coef => neko_user_access%case%fluid%c_Xh

      call neko_scratch_registry%request_field(work1, ind_work(1))
      call neko_scratch_registry%request_field(work2, ind_work(2))
      call neko_scratch_registry%request_field(work3, ind_work(3))
      call neko_scratch_registry%request_field(work4, ind_work(4))

      call grad(work1%x, work2%x, work3%x, s%x, coef)

      call coef%gs_h%op(work1, GS_OP_ADD)
      call coef%gs_h%op(work2, GS_OP_ADD)
      call coef%gs_h%op(work3, GS_OP_ADD)
      call col2(work1%x, coef%mult, work4%size())
      call col2(work2%x, coef%mult, work4%size())
      call col2(work3%x, coef%mult, work4%size())

      do i = 1, work4%size()
        absgrad = sqrt(work1%x(i,1,1,1)**2 + work2%x(i,1,1,1)**2 &
                      + work3%x(i,1,1,1)**2)
        if (absgrad < 1.0e-12_rp) then
          work1%x(i,1,1,1) = 0.0_rp
          work2%x(i,1,1,1) = 0.0_rp
          work3%x(i,1,1,1) = 0.0_rp
        else
          work1%x(i,1,1,1) = -s%x(i,1,1,1) * (1.0_rp - s%x(i,1,1,1)) &
               * (work1%x(i,1,1,1) / absgrad)
          work2%x(i,1,1,1) = -s%x(i,1,1,1) * (1.0_rp - s%x(i,1,1,1)) &
               * (work2%x(i,1,1,1) / absgrad)
          work3%x(i,1,1,1) = -s%x(i,1,1,1) * (1.0_rp - s%x(i,1,1,1)) &
               * (work3%x(i,1,1,1) / absgrad)
        end if
      end do

      call dudxyz(work4%x, work1%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
      call copy(rhs_s%x, work4%x, work4%size())
      call dudxyz(work4%x, work2%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
      call add2(rhs_s%x, work4%x, work4%size())
      call dudxyz(work4%x, work3%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
      call add2(rhs_s%x, work4%x, work4%size())

      absgrad = gamma * u_max
      call cmult(rhs_s%x, absgrad, work4%size())

      call neko_scratch_registry%relinquish_field(ind_work)

    ! CSF surface tension force on fluid momentum
    else if (scheme_name .eq. 'fluid') then

      s => neko_field_registry%get_field('phase')
      coef => neko_user_access%case%fluid%c_Xh

      rhs_u => rhs%get_by_index(1)
      rhs_v => rhs%get_by_index(2)
      rhs_w => rhs%get_by_index(3)

      call neko_scratch_registry%request_field(temp1, ind(1))
      call neko_scratch_registry%request_field(temp2, ind(2))
      call neko_scratch_registry%request_field(temp3, ind(3))
      call neko_scratch_registry%request_field(temp4, ind(4))
      call neko_scratch_registry%request_field(temp5, ind(5))
      call neko_scratch_registry%request_field(temp6, ind(6))
      call neko_scratch_registry%request_field(temp7, ind(7))

      ! grad(phi) with C0 continuity
      call grad(temp1%x, temp2%x, temp3%x, s%x, coef)
      call coef%gs_h%op(temp1, GS_OP_ADD)
      call coef%gs_h%op(temp2, GS_OP_ADD)
      call coef%gs_h%op(temp3, GS_OP_ADD)
      call col2(temp1%x, coef%mult, temp4%size())
      call col2(temp2%x, coef%mult, temp4%size())
      call col2(temp3%x, coef%mult, temp4%size())

      ! Store grad(phi)
      call copy(temp5%x, temp1%x, temp4%size())
      call copy(temp6%x, temp2%x, temp4%size())
      call copy(temp7%x, temp3%x, temp4%size())

      ! n = grad(phi) / |grad(phi)|
      do i = 1, temp4%size()
        absgrad = sqrt(temp1%x(i,1,1,1)**2 + temp2%x(i,1,1,1)**2 &
                      + temp3%x(i,1,1,1)**2)
        if (absgrad < 1.0e-12_rp) then
          temp1%x(i,1,1,1) = 0.0_rp
          temp2%x(i,1,1,1) = 0.0_rp
          temp3%x(i,1,1,1) = 0.0_rp
        else
          temp1%x(i,1,1,1) = temp1%x(i,1,1,1) / absgrad
          temp2%x(i,1,1,1) = temp2%x(i,1,1,1) / absgrad
          temp3%x(i,1,1,1) = temp3%x(i,1,1,1) / absgrad
        end if
      end do

      ! *** Phase 2 "fix": second GS pass on n_hat — NO-OP (see file header) ***
      ! Face nodes are already identical after GS(grad(phi))+normalize above.
      ! GS here sums equal values and divides back: no change.
      ! The kink is intra-element (interior nodes), unreachable by GS.
      call coef%gs_h%op(temp1, GS_OP_ADD)
      call col2(temp1%x, coef%mult, temp4%size())
      call coef%gs_h%op(temp2, GS_OP_ADD)
      call col2(temp2%x, coef%mult, temp4%size())
      call coef%gs_h%op(temp3, GS_OP_ADD)
      call col2(temp3%x, coef%mult, temp4%size())
      ! *** end fix ***

      ! kappa = -div(n)  (Brackbill CSF convention)
      call div(temp4%x, temp1%x, temp2%x, temp3%x, coef)
      call coef%gs_h%op(temp4, GS_OP_ADD)
      call col2(temp4%x, coef%mult, temp4%size())
      call cmult(temp4%x, -1.0_rp, temp4%size())

      call copy(temp1%x, temp4%x, temp4%size())

      ! F_ST = (sigma / rho) * kappa * grad(phi)
      ! Matched fluids: rho = 1.0
      do i = 1, temp4%size()
        temp5%x(i,1,1,1) = sigma * temp1%x(i,1,1,1) * temp5%x(i,1,1,1)
        temp6%x(i,1,1,1) = sigma * temp1%x(i,1,1,1) * temp6%x(i,1,1,1)
        temp7%x(i,1,1,1) = sigma * temp1%x(i,1,1,1) * temp7%x(i,1,1,1)
      end do

      call add2(rhs_u%x, temp5%x, temp4%size())
      call add2(rhs_v%x, temp6%x, temp4%size())
      call add2(rhs_w%x, temp7%x, temp4%size())

      call neko_scratch_registry%relinquish_field(ind)
    end if
  end subroutine source_term

  ! Material properties: fluid uses constant Re (set in case file);
  ! only the scalar diffusivity lambda = eps * gamma * u_max is needed.
  subroutine material_properties(scheme_name, properties, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: properties
    type(time_state_t), intent(in) :: time
    real(kind=rp) :: delta

    if (scheme_name .eq. "fluid") then
      call field_cfill(properties%get_by_index(1), 1.0_rp)           ! rho
      call field_cfill(properties%get_by_index(2), 1.0_rp / Re_b)    ! mu
    else if (scheme_name .eq. "phase") then
      delta = u_max * gamma
      call field_cfill(properties%get('phase_cp'), 1.0_rp)
      call field_cfill(properties%get('phase_lambda'), eps * delta)
    end if

  end subroutine material_properties

  ! Initial conditions:
  !   fluid -- turbulent channel (Reichardt + perturbations)
  !   phase -- single drop at channel center
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    type(field_t), pointer :: u, v, w, phase
    integer :: i, j, k, e
    real(kind=rp) :: x, y, z, d
    real(kind=rp) :: uvw(3)
    real(kind=rp) :: x_c, y_c, z_c
    real(kind=rp) :: llx, llz

    llx = 4._rp * pi
    llz = 4._rp / 3._rp * pi

    if (scheme_name .eq. 'fluid') then
      u => fields%get("u")
      v => fields%get("v")
      w => fields%get("w")

      do i = 1, u%size()
        x = u%dof%x(i,1,1,1)
        y = u%dof%y(i,1,1,1)
        z = u%dof%z(i,1,1,1)

        uvw = channel_ic(x, y, z)

        u%x(i,1,1,1) = uvw(1)
        v%x(i,1,1,1) = uvw(2)
        w%x(i,1,1,1) = uvw(3)
      end do

    else if (scheme_name .eq. 'phase') then
      phase => fields%items(1)%ptr

      ! Drop centre: streamwise and spanwise at domain mid-point,
      ! wall-normal offset by drop_center_y (0 = centreline, default).
      x_c = llx / 2._rp
      y_c = drop_center_y
      z_c = llz / 2._rp

      do e = 1, phase%msh%nelv
        do k = 1, phase%Xh%lx
          do j = 1, phase%Xh%lx
            do i = 1, phase%Xh%lx
              x = phase%dof%x(i,j,k,e)
              y = phase%dof%y(i,j,k,e)
              z = phase%dof%z(i,j,k,e)

              d = sqrt((x - x_c)**2 + (y - y_c)**2 + &
                       (z - z_c)**2) - drop_radius

              phase%x(i,j,k,e) = 0.5_rp * &
                   (1.0_rp - tanh(d / (2.0_rp * eps)))

            end do
          end do
        end do
      end do

    end if
  end subroutine initial_conditions

  ! Turbulent channel initial condition: Reichardt profile + perturbations
  function channel_ic(x, y, z) result(uvw)
    real(kind=rp) :: x, y, z
    real(kind=rp) :: uvw(3)
    real(kind=rp) :: ux, eps_ic, Re_tau, yp, Re_b, alpha, beta
    real(kind=rp) :: C, k, kx, kz, eps1, ran
    real(kind=rp) :: llx, llz

    llx = 4._rp * pi
    llz = 4._rp / 3._rp * pi

    Re_tau = 180._rp
    C = 5.17_rp
    k = 0.41_rp
    Re_b = 2800._rp

    yp = (1._rp - y) * Re_tau
    if (y .lt. 0._rp) yp = (1._rp + y) * Re_tau

    ux = 1._rp / k * log(1._rp + k * yp) + (C - (1._rp / k) * log(k)) * &
         (1._rp - exp(-yp / 11._rp) - yp / 11._rp * exp(-yp / 3._rp))
    ux = ux * Re_tau / Re_b

    uvw(1) = ux
    uvw(2) = 0._rp
    uvw(3) = 0._rp

    if (turbulent_ic) then

      ! Large-scale perturbation
      eps_ic = 0.05_rp
      kx = 3._rp
      kz = 4._rp
      alpha = kx * 2._rp * pi / llx
      beta = kz * 2._rp * pi / llz
      uvw(1) = uvw(1) + eps_ic * beta * sin(alpha * x) * cos(beta * z)
      uvw(2) = uvw(2) + eps_ic * sin(alpha * x) * sin(beta * z)
      uvw(3) = uvw(3) - eps_ic * alpha * cos(alpha * x) * sin(beta * z)

      ! Small-scale perturbation
      eps_ic = 0.005_rp
      kx = 17._rp
      kz = 13._rp
      alpha = kx * 2._rp * pi / llx
      beta = kz * 2._rp * pi / llz
      uvw(1) = uvw(1) + eps_ic * beta * sin(alpha * x) * cos(beta * z)
      uvw(2) = uvw(2) + eps_ic * sin(alpha * x) * sin(beta * z)
      uvw(3) = uvw(3) - eps_ic * alpha * cos(alpha * x) * sin(beta * z)

      ! Random perturbation in y
      eps1 = 0.001_rp
      ran = sin(-20._rp * x * z + y**3 * tan(x * z**2) + &
                100._rp * z * y - 20._rp * sin(x * y * z)**5)
      uvw(2) = uvw(2) + eps1 * ran

    end if

  end function channel_ic

end module user
