module user
  use neko
  implicit none

contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%fluid_user_ic => user_ic
  end subroutine user_setup
  ! User defined initial condition
  subroutine user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    integer :: i
    real(kind=rp) :: uvw(3)

    do i = 1, u%dof%size()
       uvw = pipe_ic(u%dof%x(i,1,1,1), u%dof%y(i,1,1,1), u%dof%z(i,1,1,1))
       u%x(i,1,1,1) = uvw(1)
       v%x(i,1,1,1) = uvw(2)
       w%x(i,1,1,1) = uvw(3)
    end do
  end subroutine user_ic

  function pipe_ic(x, y, z) result(uvw)
    real(kind=rp) :: x, y, z
    real(kind=rp) :: uvw(3)
    real(kind=rp) :: rand, ux, uy, uz, xr, yr, rr, zo
    real(kind=rp) :: amp_z, freq_z, freq_t, amp_tht, amp_clip, blt
    real(kind=rp) :: phase_z, arg_tht, amp_sin, pi, th

    pi = 4_rp * atan(1.0_rp)
    xr = x
    yr = y
    rr = xr*xr + yr*yr
    if (rr .gt. 0.0_rp) rr = sqrt(rr)
    th = atan2(y,x)
    zo = 2.0_rp * pi * z / 25.0_rp

    uz = 6.0_rp * (1.0_rp - rr**6.0_rp) / 5.0_rp

    ! Assign a wiggly shear layer near the wall
    amp_z = 35e-2_rp ! Fraction of 2pi for z-based phase modification
    freq_z = 5.0_rp ! Number of wiggles in axial- (z-) direction
    freq_t = 6.0_rp ! Frequency of wiggles in azimuthal-direction

    amp_tht = 10.0_rp ! Amplification factor for clipped sine function
    amp_clip = 4e-1_rp ! Clipped amplitude

    blt = 3.5e-1_rp ! Fraction of boundary layer with momentum deficit

    phase_z = amp_z * (2.0_rp * pi) * sin(freq_z * zo)

    arg_tht = freq_t * th + phase_z
    amp_sin = 5.0_rp * sin(arg_tht)
    if (amp_sin .gt. amp_clip) amp_sin = amp_clip
    if (amp_sin .lt. -amp_clip) amp_sin = -amp_clip

    if (rr .gt. (1-blt)) uz = uz + amp_sin
    call random_number(rand)

    ux = 5e-2_rp * rand * rand
    uy = 1e-1_rp * rand * rand * rand
    uz = uz + 1e-2_rp * rand

    uvw(1) = ux
    uvw(2) = uy
    uvw(3) = uz

  end function pipe_ic
end module user
