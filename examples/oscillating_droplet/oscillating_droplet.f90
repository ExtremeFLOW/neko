! Oscillating droplet test case for multiphase flow
!
! An initially elliptical drop (mode-2 perturbation with epsilon=0.1) in a
! 1x1 box with slip boundary conditions. Surface tension drives the drop
! back towards a circle, causing oscillations that decay through viscous
! dissipation. This verifies the CSF surface tension implementation.
!
! Physical parameters (same as spurious currents):
!   rho = 300, mu = 0.1, sigma = 1.0, R0 = 0.2, D = 0.4
!   Laplace number La = sigma*rho*D/mu^2 = 12000
!
! Initial shape: r(theta) = R0 * (1 + epsilon * cos(2*theta))
!   with epsilon = 0.1 (10% perturbation)
!
! The domain is quasi-2D (one element in z). The w velocity is set to
! zero at every step to prevent 3D instabilities.
!
module user
  use neko
  implicit none

  ! Global user variables
  type(field_t) :: w1
  type(field_t) :: temp8,temp9
  type(field_t) :: vort1,vort2,vort3

  type(file_t) output_file ! output file
  type(vector_t) :: vec_out ! will store our output data
  integer :: ipostproc ! frequency of the output

  real(kind=rp) :: eps
  real(kind=rp) :: gamma, u_max, sigma

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%startup => startup
    user%compute => compute
    user%initialize =>initialize
    user%finalize => finalize
    user%source_term => source_term
    user%material_properties => material_properties
    user%initial_conditions => initial_conditions
  end subroutine user_setup

  subroutine startup(params)
    type(json_file), intent(inout) :: params
    character(len=50) :: mess

    call json_get(params, "case.fluid.ipostproc", ipostproc)
    write(mess,*) "postprocessing steps : ",ipostproc
    call neko_log%message(mess)

    call json_get(params, "case.scalar.epsilon",eps)
    call json_get(params, "case.scalar.gamma", gamma)
    call json_get(params, "case.scalar.sigma", sigma)
    u_max = 1.0_rp
  end subroutine startup

  ! User-defined initialization called just before time loop starts
  subroutine initialize(time)
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: u

    ! Initialize the file object and create the output.csv file
    ! in the working directory (overwrite any existing file)
    call output_file%init("ekin.csv", overwrite=.true.)
    call output_file%set_header("# t,Ekin,enst,u_max," // &
         "kappa_max,kappa_min,kappa_rms,Fst_max,phi_min,phi_max")
    call vec_out%init(10)

    ! initialize work arrays for postprocessing
    u => neko_field_registry%get_field("u")
    call w1%init(u%dof, 'work1')
    call temp8%init(u%dof)
    call temp9%init(u%dof)
    call vort1%init(u%dof)
    call vort2%init(u%dof)
    call vort3%init(u%dof)


    ! call usercheck also for tstep=0
    call compute(time)

  end subroutine initialize

  ! User-defined routine called at the end of every time step
  subroutine compute(time)
    type(time_state_t), intent(in) :: time

    integer :: ntot, i
    integer :: ind_diag(7)
    real(kind=rp) :: ekin, enst
    real(kind=rp) :: absgrad, kappa_max, kappa_min, kappa_rms, kappa_weight_sum
    real(kind=rp) :: Fst_max, phi_min, phi_max
    type(field_t), pointer :: u, v, w, s
    type(field_t), pointer :: t1, t2, t3, t4, t5, t6, t7
    type(coef_t), pointer :: coef


    ! do the postprocessing only every 50 time steps (including step 0)
    if (mod(time%tstep, ipostproc) .ne. 0) return

    coef => neko_user_access%case%fluid%c_Xh
    u => neko_field_registry%get_field("u")
    v => neko_field_registry%get_field("v")
    w => neko_field_registry%get_field("w")

    ntot = u%dof%size()

    ! compute the velocity magnitude in w1 (u^2 + v^2 + w^2)
    call field_col3(w1, u, u, ntot)
    call field_addcol3(w1, v, v, ntot)
    call field_addcol3(w1, w, w, ntot)

    ! compute maximum velocity magnitude over the domain
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(w1%x, w1%x_d, w1%size(), &
            DEVICE_TO_HOST, sync=.true.)
    end if
    u_max = sqrt(glmax(w1%x, ntot))

    ! compute the kinetic energy (reuse w1 which already has velocity magnitude squared)
    ! For glsc2 we need to call the correct backend
    if (NEKO_BCKND_DEVICE .eq. 1) then
       ekin = 0.5_rp * device_glsc2(w1%x_d, coef%B_d,ntot) / coef%volume
    else
       ekin = 0.5_rp * glsc2(w1%x, coef%B, ntot) / coef%volume
    end if

    ! compute enstrophy
    call curl(vort1,vort2,vort3, u, v, w, temp8, temp9, coef)
    call field_col3(w1, vort1, vort1, ntot)
    call field_addcol3(w1, vort2, vort2, ntot)
    call field_addcol3(w1, vort3, vort3, ntot)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       enst = 0.5_rp * device_glsc2(w1%x_d, coef%B_d, ntot) / coef%volume
    else
       enst = 0.5_rp * glsc2(w1%x, coef%B, ntot) / coef%volume
    end if


    ! ==================== CSF DIAGNOSTICS ====================
    ! Replicate the CSF curvature/force computation from source_term
    ! to extract diagnostic quantities.

    s => neko_field_registry%get_field('phase')

    call neko_scratch_registry%request_field(t1, ind_diag(1))
    call neko_scratch_registry%request_field(t2, ind_diag(2))
    call neko_scratch_registry%request_field(t3, ind_diag(3))
    call neko_scratch_registry%request_field(t4, ind_diag(4))
    call neko_scratch_registry%request_field(t5, ind_diag(5))
    call neko_scratch_registry%request_field(t6, ind_diag(6))
    call neko_scratch_registry%request_field(t7, ind_diag(7))

    ! Step 1: grad(phi) with GS+mult for C0 continuity
    call grad(t1%x, t2%x, t3%x, s%x, coef)
    call coef%gs_h%op(t1, GS_OP_ADD)
    call coef%gs_h%op(t2, GS_OP_ADD)
    call coef%gs_h%op(t3, GS_OP_ADD)
    call col2(t1%x, coef%mult, ntot)
    call col2(t2%x, coef%mult, ntot)
    call col2(t3%x, coef%mult, ntot)

    ! Save grad(phi) for later
    call copy(t5%x, t1%x, ntot)
    call copy(t6%x, t2%x, ntot)
    call copy(t7%x, t3%x, ntot)

    ! Step 2: normalize to n = grad(phi)/|grad(phi)|
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

    ! Step 3: kappa = -div(n) with GS+mult (Brackbill CSF convention)
    call div(t4%x, t1%x, t2%x, t3%x, coef)
    call coef%gs_h%op(t4, GS_OP_ADD)
    call col2(t4%x, coef%mult, ntot)
    call cmult(t4%x, -1.0_rp, ntot)

    ! Step 4: compute diagnostic scalars
    ! kappa extremes
    kappa_max = glmax(t4%x, ntot)
    kappa_min = glmin(t4%x, ntot)

    ! kappa RMS weighted by |grad(phi)| (focuses on interface)
    ! w1 = |grad(phi)|
    do i = 1, ntot
      w1%x(i,1,1,1) = sqrt(t5%x(i,1,1,1)**2 + t6%x(i,1,1,1)**2 &
                          + t7%x(i,1,1,1)**2)
    end do
    ! t1 = kappa^2 * |grad(phi)|  (reuse t1, normals no longer needed)
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

    ! |F_ST| = sigma * |kappa| * |grad(phi)|
    do i = 1, ntot
      w1%x(i,1,1,1) = sigma * abs(t4%x(i,1,1,1)) * w1%x(i,1,1,1)
    end do
    Fst_max = glmax(w1%x, ntot)

    ! Phase field bounds
    phi_min = glmin(s%x, ntot)
    phi_max = glmax(s%x, ntot)

    call neko_scratch_registry%relinquish_field(ind_diag)
    ! ==================== END CSF DIAGNOSTICS ====================

    ! output all this to file
    call neko_log%message("Writing csv file")
    vec_out%x = (/ekin, enst, u_max, kappa_max, kappa_min, kappa_rms, &
                  Fst_max, phi_min, phi_max, 0.0_rp/)
    call output_file%write(vec_out, time%t)

    ! set the w component to zero to avoid any 3D instability
    ! in this quasi-2D flow
    call field_rzero(w)

  end subroutine compute

  ! User-defined finalization routine called at the end of the simulation
  subroutine finalize(time)
    type(time_state_t), intent(in) :: time

    ! Deallocate the fields
    call w1%free()

    ! Deallocate output file and vector
    call file_free(output_file)
    call vec_out%free

  end subroutine finalize

  ! User-defined source term (forcing function)
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

    ! Only apply forcing to the scalar equation (phase)
    if (scheme_name .eq. 'phase') then

      ! Get the right-hand side field for the scalar
      rhs_s => rhs%get_by_index(1)  ! scalar field

      ! Get the scalar field and coefficients
      s => neko_field_registry%get_field('phase')
      coef => neko_user_access%case%fluid%c_Xh

      ! Request scratch fields
      call neko_scratch_registry%request_field(work1, ind_work(1))
      call neko_scratch_registry%request_field(work2, ind_work(2))
      call neko_scratch_registry%request_field(work3, ind_work(3))
      call neko_scratch_registry%request_field(work4, ind_work(4))

      ! Compute gradient of scalar field
      call grad(work1%x, work2%x, work3%x, s%x, coef)

      ! Apply gather-scatter and multiplicity
      call coef%gs_h%op(work1, GS_OP_ADD)
      call coef%gs_h%op(work2, GS_OP_ADD)
      call coef%gs_h%op(work3, GS_OP_ADD)
      call col2(work1%x, coef%mult, work4%size())
      call col2(work2%x, coef%mult, work4%size())
      call col2(work3%x, coef%mult, work4%size())

      ! Compute normalized gradient and apply phase field forcing
      do i = 1, work4%size()
        absgrad = sqrt(work1%x(i,1,1,1)**2+work2%x(i,1,1,1)**2+work3%x(i,1,1,1)**2)
        if (absgrad == 0.0_rp) then
            print *, 'warning, absgrad==', absgrad
            absgrad = 1e21_rp
        end if

        work1%x(i,1,1,1) = - s%x(i,1,1,1)*(1.0_rp-s%x(i,1,1,1))*(work1%x(i,1,1,1)/absgrad)
        work2%x(i,1,1,1) = - s%x(i,1,1,1)*(1.0_rp-s%x(i,1,1,1))*(work2%x(i,1,1,1)/absgrad)
        work3%x(i,1,1,1) = - s%x(i,1,1,1)*(1.0_rp-s%x(i,1,1,1))*(work3%x(i,1,1,1)/absgrad)
      end do

      ! Compute divergence of the normalized gradient
      call dudxyz(work4%x, work1%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
      call copy(rhs_s%x, work4%x, work4%size())
      call dudxyz(work4%x, work2%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
      call add2(rhs_s%x, work4%x, work4%size())
      call dudxyz(work4%x, work3%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
      call add2(rhs_s%x, work4%x, work4%size())

      ! Scale by gamma * u_max
      absgrad = gamma * u_max
      call cmult(rhs_s%x, absgrad, work4%size())

      ! Release scratch fields
      call neko_scratch_registry%relinquish_field(ind_work)

    ! Only apply forcing to the momentum equation (fluid)
    else if (scheme_name .eq. 'fluid') then

      ! Get fields
      s => neko_field_registry%get_field('phase')
      coef => neko_user_access%case%fluid%c_Xh

      rhs_u => rhs%get_by_index(1)
      rhs_v => rhs%get_by_index(2)
      rhs_w => rhs%get_by_index(3)

      ! Request scratch fields
      call neko_scratch_registry%request_field(temp1, ind(1))
      call neko_scratch_registry%request_field(temp2, ind(2))
      call neko_scratch_registry%request_field(temp3, ind(3))
      call neko_scratch_registry%request_field(temp4, ind(4))
      call neko_scratch_registry%request_field(temp5, ind(5))
      call neko_scratch_registry%request_field(temp6, ind(6))
      call neko_scratch_registry%request_field(temp7, ind(7))

      ! Compute gradient of phase field grad(phi)
      call grad(temp1%x, temp2%x, temp3%x, s%x, coef)

      ! Apply gather-scatter and multiplicity for continuity
      call coef%gs_h%op(temp1, GS_OP_ADD)
      call coef%gs_h%op(temp2, GS_OP_ADD)
      call coef%gs_h%op(temp3, GS_OP_ADD)
      call col2(temp1%x, coef%mult, temp4%size())
      call col2(temp2%x, coef%mult, temp4%size())
      call col2(temp3%x, coef%mult, temp4%size())

      ! Store grad(phi) for later (we need it at the end)
      call copy(temp5%x, temp1%x, temp4%size())  ! temp5 = dphi/dx
      call copy(temp6%x, temp2%x, temp4%size())  ! temp6 = dphi/dy
      call copy(temp7%x, temp3%x, temp4%size())  ! temp7 = dphi/dz

      ! Compute normalized gradient n = grad(phi)/|grad(phi)|
      do i = 1, temp4%size()
        absgrad = sqrt(temp1%x(i,1,1,1)**2 + temp2%x(i,1,1,1)**2 + temp3%x(i,1,1,1)**2)
        if (absgrad < 1.0e-12_rp) then
            ! Avoid division by zero in bulk phases
            temp1%x(i,1,1,1) = 0.0_rp
            temp2%x(i,1,1,1) = 0.0_rp
            temp3%x(i,1,1,1) = 0.0_rp
        else
            temp1%x(i,1,1,1) = temp1%x(i,1,1,1) / absgrad  ! nx
            temp2%x(i,1,1,1) = temp2%x(i,1,1,1) / absgrad  ! ny
            temp3%x(i,1,1,1) = temp3%x(i,1,1,1) / absgrad  ! nz
        end if
      end do

      ! Compute curvature kappa = -div(n) (Brackbill CSF convention)
      call div(temp4%x,temp1%x, temp2%x,temp3%x,coef)

      ! Apply gather-scatter and multiplicity for continuity
      call coef%gs_h%op(temp4, GS_OP_ADD)
      call col2(temp4%x, coef%mult, temp4%size())
      call cmult(temp4%x, -1.0_rp, temp4%size())

      call copy(temp1%x, temp4%x, temp4%size())

      ! Now temp1 contains kappa = -div(n) (curvature, positive for convex drop)

      ! Compute surface tension acceleration: a_ST = (sigma/rho) * kappa * grad(phi)
      do i = 1, temp4%size()
        temp5%x(i,1,1,1) = (sigma / 300.0_rp) * temp1%x(i,1,1,1) * temp5%x(i,1,1,1)
        temp6%x(i,1,1,1) = (sigma / 300.0_rp) * temp1%x(i,1,1,1) * temp6%x(i,1,1,1)
        temp7%x(i,1,1,1) = (sigma / 300.0_rp) * temp1%x(i,1,1,1) * temp7%x(i,1,1,1)
      end do

      ! Add surface tension force to momentum equations
      call add2(rhs_u%x, temp5%x, temp4%size())
      call add2(rhs_v%x, temp6%x, temp4%size())
      call add2(rhs_w%x, temp7%x, temp4%size())

      ! Release scratch fields
      call neko_scratch_registry%relinquish_field(ind)
    end if
  end subroutine source_term

  ! User-defined material properties
  subroutine material_properties(scheme_name, properties, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: properties
    type(time_state_t), intent(in) :: time
    real(kind=rp) :: delta, lambda_val, mu_val, rho_val

    if (scheme_name .eq. "fluid") then
      rho_val = 300.0_rp
      mu_val = 0.1_rp
      call field_cfill(properties%get("fluid_rho"), rho_val)
      call field_cfill(properties%get("fluid_mu"), mu_val)
    else if (scheme_name .eq. "phase") then
      delta = u_max * gamma
      call field_cfill(properties%get('phase_cp'), 1.0_rp)
      call field_cfill(properties%get('phase_lambda'), eps*delta)
    end if

  end subroutine material_properties

  ! User-defined initial conditions: elliptical drop (mode-2 perturbation)
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    type(field_t), pointer :: phase
    integer :: i, j, k, e
    real(kind=rp) :: x, y, d_x, theta, r_theta
    real(kind=rp) :: x_center, y_center, radius, epsilon

    if (scheme_name .ne. 'phase') return

    phase => fields%items(1)%ptr

    x_center = 0.5_rp
    y_center = 0.5_rp
    radius = 0.2_rp
    epsilon = 0.1_rp   ! 10% mode-2 perturbation

    do e = 1, phase%msh%nelv
      do k = 1, phase%Xh%lx
        do j = 1, phase%Xh%lx
          do i = 1, phase%Xh%lx
            x = phase%dof%x(i,j,k,e)
            y = phase%dof%y(i,j,k,e)

            ! Angle from drop center
            theta = atan2(y - y_center, x - x_center)

            ! Perturbed radius: elliptical shape via mode-2 cosine
            r_theta = radius * (1.0_rp + epsilon * cos(2.0_rp * theta))

            ! Signed distance from perturbed interface
            d_x = sqrt((x - x_center)**2 + (y - y_center)**2) - r_theta

            ! Smooth phase field via hyperbolic tangent profile
            phase%x(i,j,k,e) = 0.5_rp * (1.0_rp - tanh(d_x / (2.0_rp * eps)))

          end do
        end do
      end do
    end do

  end subroutine initial_conditions

end module user
