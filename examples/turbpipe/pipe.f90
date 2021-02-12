module user
  use neko
  implicit none
contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%fluid_usr_ic => user_ic
  end subroutine user_setup
  ! User defined initial condition
  subroutine user_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(param_t), intent(inout) :: params
    integer :: i
    real(kind=dp) :: uvw(3)

    do i = 1, u%dof%size()
       uvw = pipe_ic(u%dof%x(i,1,1,1),u%dof%y(i,1,1,1),u%dof%z(i,1,1,1))
       u%x(i,1,1,1) = uvw(1)
       v%x(i,1,1,1) = uvw(2)
       w%x(i,1,1,1) = uvw(3)
    end do
  end subroutine user_ic
  
  function pipe_ic(x, y, z) result(uvw)
    real(kind=dp) :: x, y, z
    real(kind=dp) :: uvw(3)
    real(kind=dp) :: rand, ux, uy, uz, xr, yr, rr, zo
    real(kind=dp) :: amp_z, freq_z, freq_t, amp_tht, amp_clip, blt
    real(kind=dp) :: phase_z, arg_tht, amp_sin, pi, th

    pi = 4d0 * atan(1d0)
    xr = x
    yr = y
    rr = xr*xr + yr*yr
    if (rr.gt.0) rr=sqrt(rr)
    th = atan2(y,x)
    zo = 2*pi*z/25d0

    uz = 6d0*(1d0-rr**6d0)/5d0

    ! Assign a wiggly shear layer near the wall
    amp_z    = 35d-2  ! Fraction of 2pi for z-based phase modification
    freq_z   = 4d0     ! Number of wiggles in axial- (z-) direction
    freq_t   = 9d0     ! Frequency of wiggles in azimuthal-direction

    amp_tht  = 5d0     ! Amplification factor for clipped sine function
    amp_clip = 2d-1   ! Clipped amplitude

    blt      = 7d-2  ! Fraction of boundary layer with momentum deficit

    phase_z = amp_z*(2d0*pi)*sin(freq_z*zo)

    arg_tht = freq_t*th + phase_z
    amp_sin = 5d0*sin(arg_tht)
    if (amp_sin.gt. amp_clip) amp_sin =  amp_clip
    if (amp_sin.lt.-amp_clip) amp_sin = -amp_clip

    if (rr.gt.(1-blt)) uz = uz + amp_sin
    call random_number(rand)

    ux   = 5d-2*rand*rand
    uy   = 1d-1*rand*rand*rand
    uz   = uz + 1d-2*rand
    
    uvw(1) = ux
    uvw(2) = uy
    uvw(3) = uz

  end function pipe_ic


end module user
