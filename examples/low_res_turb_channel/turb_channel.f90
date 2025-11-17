! Martin Karp 13/3-2023
! updated initial condition Philipp Schlatter 09/07/2024
module user
  use neko
  implicit none
  type(file_t) output_file ! output file
  type(vector_t) :: vec_out ! will store our output data

contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initial_conditions => initial_conditions
    user%mesh_setup => user_mesh_scale
    user%compute => compute
    user%initialize =>initialize
  end subroutine user_setup

  ! User-defined initialization called just before time loop starts
  subroutine initialize(time)
    type(time_state_t), intent(in) :: time


    ! Initialize the file object and create the output.csv file
    ! in the working directory
    call output_file%init("divergence.csv")
    call output_file%set_header("# t,integral,abs,L2,weak_abs,weak_L2")
    call vec_out%init(6) ! Initialize our vector with 3 elements

    ! call usercheck also for tstep=0
    call compute(time)

  end subroutine initialize
  ! Rescale mesh, we create a mesh with some refinement close to the wall.
  ! initial mesh: 0..4, -1..1, 0..1.5
  ! mesh size (4*pi,2*delta,4/3*pi)
  ! New mesh can easily be genreated with genmeshbox
  ! OBS refinement is not smooth and the constant values are a bit ad hoc.
  ! Stats converge close to reference DNS
  subroutine user_mesh_scale(msh, time)
    type(mesh_t), intent(inout) :: msh
    type(time_state_t), intent(in) :: time
    integer :: i, p, nvert

    real(kind=rp) :: d, y, viscous_layer, visc_el_h, el_h
    real(kind=rp) :: center_el_h, dist_from_wall
    integer :: el_in_visc_lay, el_in_y
    real(kind=rp) :: llx, llz

    ! target mesh size
    llx = 2.*pi
    llz = pi

    ! rescale mesh
    el_in_y = 18
    el_in_visc_lay = 2
    viscous_layer = 0.0888889
    el_h = 2.0_rp/el_in_y
    visc_el_h = viscous_layer/el_in_visc_lay
    center_el_h = (1.0_rp-viscous_layer)/(el_in_y/2-el_in_visc_lay)

    nvert = size(msh%points)
    do i = 1, nvert
       msh%points(i)%x(1) = llx/4.*msh%points(i)%x(1)
       y = msh%points(i)%x(2)
       if ((1-abs(y)) .le. (el_in_visc_lay*el_h)) then
          dist_from_wall = (1-abs(y))/el_h*visc_el_h
       else
          dist_from_wall = viscous_layer + (1-abs(y)- &
               el_in_visc_lay*el_h)/el_h*center_el_h
       end if
       ! if (y .gt. 0) msh%points(i)%x(2) = 1.0_rp - dist_from_wall
       ! if (y .lt. 0) msh%points(i)%x(2) = -1.0_rp + dist_from_wall
       msh%points(i)%x(3) = 2./3.*llz*msh%points(i)%x(3)
    end do

  end subroutine user_mesh_scale

  ! User defined initial condition
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    real(kind=rp) :: uvw(3), x, y, z
    type (field_t), pointer :: u, v, w
    integer :: i

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

  ! Kind of brute force with rather large initial disturbances
  function channel_ic(x, y, z) result(uvw)
    real(kind=rp) :: x, y, z
    real(kind=rp) :: uvw(3)
    real(kind=rp) :: ux, uy, uz, eps, Re_tau, yp, Re_b, alpha, beta
    real(kind=rp) :: C, k, kx, kz, eps1, ran

    real(kind=rp) :: llx, llz

    llx = 4.*pi
    llz = 4./3.*pi

    Re_tau = 180
    C = 5.17
    k = 0.41
    Re_b = 2800

    yp = (1-y)*Re_tau
    if (y .lt. 0) yp = (1+y)*Re_tau

    ! Reichardt function
    ux = 1/k*log(1.0+k*yp) + (C - (1.0/k)*log(k)) * &
         (1.0 - exp(-yp/11.0) - yp/11*exp(-yp/3.0))
    ux = ux * Re_tau/Re_b

    ! actually, sometimes one may not use the turbulent profile, but
    ! rather the parabolic lamianr one
    ! ux = 1.5*(1-y**2)

    ! add perturbations to trigger turbulence
    ! base flow
    uvw(1) = ux
    uvw(2) = 0
    uvw(3) = 0

    ! first, large scale perturbation
    eps = 0.05
    kx = 3
    kz = 4
    alpha = kx * 2*PI/llx
    beta = kz * 2*PI/llz
    uvw(1) = uvw(1) + eps*beta * sin(alpha*x)*cos(beta*z)
    uvw(2) = uvw(2) + eps * sin(alpha*x)*sin(beta*z)
    uvw(3) = uvw(3) -eps*alpha * cos(alpha*x)*sin(beta*z)

    ! second, small scale perturbation
    eps = 0.005
    kx = 17
    kz = 13
    alpha = kx * 2*PI/llx
    beta = kz * 2*PI/llz
    uvw(1) = uvw(1) + eps*beta * sin(alpha*x)*cos(beta*z)
    uvw(2) = uvw(2) + eps * sin(alpha*x)*sin(beta*z)
    uvw(3) = uvw(3) -eps*alpha * cos(alpha*x)*sin(beta*z)

    ! finally, random perturbations only in y
    eps1 = 0.001
    ran = sin(-20*x*z+y**3*tan(x*z**2)+100*z*y-20*sin(x*y*z)**5)
    uvw(2) = uvw(2) + eps1*ran

  end function channel_ic

! User-defined routine called at the end of every time step
  subroutine compute(time)
    type(time_state_t), intent(in) :: time

    integer :: ntot, i
    real(kind=rp) :: integral, L2, abs, weak_abs, weak_L2
    type(field_t), pointer :: u, v, w, divergence, ta1, ta2
    type(coef_t), pointer :: coef
    integer :: temp_indices(4)


    ! do the postprocessing only every 50 time steps (including step 0)
    if (mod(time%tstep, 20) .ne. 0) return

    coef => neko_user_access%case%fluid%c_Xh
    u => neko_field_registry%get_field("u")
    v => neko_field_registry%get_field("v")
    w => neko_field_registry%get_field("w")
    call neko_scratch_registry%request_field(divergence, temp_indices(1))
    call neko_scratch_registry%request_field(ta1, temp_indices(2))
    call neko_scratch_registry%request_field(ta2, temp_indices(3))

    ntot = u%dof%size()

    call div(divergence%x, u%x, v%x, w%x, coef)

    ! integral
    if (NEKO_BCKND_DEVICE .eq. 1) then
       integral = 0.5_rp * device_glsc2(divergence%x_d, coef%B_d, ntot) / coef%volume
    else
       integral = 0.5_rp * glsc2(divergence%x, coef%B, ntot) / coef%volume
    end if

    ! absolute value
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_absval(divergence%x_d, ntot)
       abs = 0.5_rp * device_glsc2(divergence%x_d, coef%B_d, ntot) / coef%volume
    else
       call absval(divergence%x, ntot)
       abs = 0.5_rp * glsc2(divergence%x, coef%B, ntot) / coef%volume
    end if

    ! L2
    call field_col2(divergence, divergence)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       L2 = 0.5_rp * device_glsc2(divergence%x_d, coef%B_d, ntot) / coef%volume
    else
       L2 = 0.5_rp * glsc2(divergence%x, coef%B, ntot) / coef%volume
    end if

    ! weak divergence (not sure if that's the correct nomenclature)
    ! I mean \nabla v \cdot u
    ! with v the test functions
    call cdtp(divergence%x, u%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call cdtp(ta1%x, u%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call cdtp(ta2%x, u%x, coef%drdz, coef%dsdz, coef%dtdz, coef)

    call field_add2(divergence, ta1)
    call field_add2(divergence, ta2)

    ! weak absolute value
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_absval(divergence%x_d, ntot)
       weak_abs = 0.5_rp * device_glsc2(divergence%x_d, coef%B_d, ntot) / coef%volume
    else
       call absval(divergence%x, ntot)
       weak_abs = 0.5_rp * glsc2(divergence%x, coef%B, ntot) / coef%volume
    end if

    ! weak L2
    call field_col2(divergence, divergence)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       weak_L2 = 0.5_rp * device_glsc2(divergence%x_d, coef%B_d, ntot) / coef%volume
    else
       weak_L2 = 0.5_rp * glsc2(divergence%x, coef%B, ntot) / coef%volume
    end if




    ! output all this to file
    call neko_log%message("Writing csv file")
    vec_out%x = (/integral, abs, L2/)
    call output_file%write(vec_out, time%t)

    call neko_scratch_registry%relinquish_field(temp_indices)

  end subroutine compute
end module user
