module user
  use neko
  implicit none

  ! Data streamer
  type(data_streamer_t) :: dstream
  integer :: ipostproc ! frequency of the streaming

contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%fluid_user_ic => user_ic
    user%user_check => user_check
    user%user_init_modules => user_initialize
    user%user_finalize_modules => user_finalize
  end subroutine user_setup

  ! User-defined initialization called just before time loop starts
  subroutine user_initialize(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u, v, w, p
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params

    integer tstep
    character(len=50) :: mess

    ! read postprocessing interval
    call json_get(params, "case.istream", ipostproc)
    write(mess,*) "streaming steps : ", ipostproc
    call neko_log%message(mess)

    ! Initialize the streamer
    call dstream%init(coef)

    ! Stream the mesh
    call dstream%stream(coef%dof%x)
    call dstream%stream(coef%dof%y)
    call dstream%stream(coef%dof%z)

  end subroutine user_initialize

  ! User-defined routine called at the end of every time step
  subroutine user_check(t, tstep, u, v, w, p, coef, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    type(coef_t), intent(inout) :: coef
    type(json_file), intent(inout) :: params
    type(field_t), intent(inout) :: u, v, w, p
    integer :: ntot, i, n
    real(kind=rp) :: ekin, enst

    if (mod(tstep,ipostproc).ne.0) return

    n = u%dof%size()

    ! Average over interfaces
    call coef%gs_h%op(u, GS_OP_ADD)
    call device_col2(u%x_d, coef%mult_d, n)
    call coef%gs_h%op(v, GS_OP_ADD)
    call device_col2(v%x_d, coef%mult_d, n)
    call coef%gs_h%op(w, GS_OP_ADD)
    call device_col2(w%x_d, coef%mult_d, n)

    ! Sync the GPU-CPU
    call device_memcpy(u%x, u%x_d, n, DEVICE_TO_HOST, sync=.true.)
    call device_memcpy(v%x, v%x_d, n, DEVICE_TO_HOST, sync=.true.)
    call device_memcpy(w%x, w%x_d, n, DEVICE_TO_HOST, sync=.true.)

    ! Stream the data
    call dstream%stream(u%x)
    call dstream%stream(v%x)
    call dstream%stream(w%x)

  end subroutine user_check

  ! User-defined finalization routine called at the end of the simulation
  subroutine user_finalize(t, params)
    real(kind=rp) :: t
    type(json_file), intent(inout) :: params

    ! Finalize the stream
    call dstream%free()

  end subroutine user_finalize

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
