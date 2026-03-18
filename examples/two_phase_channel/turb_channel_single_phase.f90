! Single-phase turbulent channel — fluid spin-up for two-phase restart
!
! Runs the channel to statistical turbulence without any scalar.
! Produces checkpoint files fluid00000.chkp (t=0), fluid00001.chkp (t=5), ...
! at regular intervals (every 5 time units). Use fluid00004.chkp (t=20) as the
! restart_file for turb_channel_two_phase_restart.case.
!
! Writes ekin.csv with columns: t, Ekin, u_max
! Use these to confirm turbulence is established before restarting:
!   - u_max should fluctuate around ~1.15-1.20 (turbulent, not ~1.5 laminar)
!   - Ekin should be statistically stationary (flat mean, no trend)
!   Typically achieved by t ≈ 10-15 convective time units from the Reichardt IC.
!
! Run:
!   makeneko turb_channel_single_phase.f90
!   mpirun -np <N> ./neko turb_channel_single_phase.case
!
module user
  use neko
  implicit none

  type(field_t) :: w1
  type(file_t) :: output_file
  type(vector_t) :: vec_out
  integer :: ipostproc

  real(kind=rp) :: Re_b

contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%startup => startup
    user%initialize => initialize
    user%compute => compute
    user%finalize => finalize
    user%initial_conditions => initial_conditions
  end subroutine user_setup

  subroutine startup(params)
    type(json_file), intent(inout) :: params
    call json_get(params, "case.fluid.Re", Re_b)
    call json_get(params, "case.fluid.ipostproc", ipostproc)
  end subroutine startup

  subroutine initialize(time)
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: u

    call output_file%init("ekin.csv", overwrite=.true.)
    call output_file%set_header("# t,Ekin,u_max")
    call vec_out%init(3)

    u => neko_field_registry%get_field("u")
    call w1%init(u%dof, 'work1')

    call compute(time)

  end subroutine initialize

  ! Diagnostic output: volume-averaged kinetic energy and u_max.
  ! u_max is the key turbulence indicator: turbulent flow at Re_tau=180
  ! gives u_max ~ 1.15-1.20 with fluctuations; laminar Poiseuille gives 1.5.
  ! Ekin reaching a stationary mean (no secular trend) confirms statistical
  ! equilibrium. Expect stationarity by t ~ 10-15 from the Reichardt IC.
  subroutine compute(time)
    type(time_state_t), intent(in) :: time

    integer :: ntot
    real(kind=rp) :: ekin, umax_val
    type(field_t), pointer :: u, v, w
    type(coef_t), pointer :: coef

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
      call device_memcpy(w1%x, w1%x_d, w1%size(), DEVICE_TO_HOST, sync=.true.)
    end if
    umax_val = sqrt(glmax(w1%x, ntot))

    if (NEKO_BCKND_DEVICE .eq. 1) then
      ekin = 0.5_rp * device_glsc2(w1%x_d, coef%B_d, ntot) / coef%volume
    else
      ekin = 0.5_rp * glsc2(w1%x, coef%B, ntot) / coef%volume
    end if

    call neko_log%message("Writing csv file")
    vec_out%x(1) = ekin
    vec_out%x(2) = umax_val
    vec_out%x(3) = 0.0_rp
    call output_file%write(vec_out, time%t)

  end subroutine compute

  subroutine finalize(time)
    type(time_state_t), intent(in) :: time
    call w1%free()
    call file_free(output_file)
    call vec_out%free
  end subroutine finalize

  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    type(field_t), pointer :: u, v, w
    integer :: i
    real(kind=rp) :: x, y, z
    real(kind=rp) :: uvw(3)

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
    end if

  end subroutine initial_conditions

  ! Turbulent channel IC: Reichardt profile + perturbations
  function channel_ic(x, y, z) result(uvw)
    real(kind=rp) :: x, y, z
    real(kind=rp) :: uvw(3)
    real(kind=rp) :: ux, eps_ic, Re_tau, yp, alpha, beta
    real(kind=rp) :: C, k, kx, kz, eps1, ran
    real(kind=rp) :: llx, llz

    llx = 4._rp * pi
    llz = 4._rp / 3._rp * pi

    Re_tau = 180._rp
    C = 5.17_rp
    k = 0.41_rp

    yp = (1._rp - y) * Re_tau
    if (y .lt. 0._rp) yp = (1._rp + y) * Re_tau

    ux = 1._rp / k * log(1._rp + k * yp) + (C - (1._rp / k) * log(k)) * &
         (1._rp - exp(-yp / 11._rp) - yp / 11._rp * exp(-yp / 3._rp))
    ux = ux * Re_tau / Re_b

    uvw(1) = ux
    uvw(2) = 0._rp
    uvw(3) = 0._rp

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

  end function channel_ic

end module user
