module user
  use neko
  use import_field_utils, only : import_fields
  implicit none

  real(kind=rp), parameter :: u_ref = 4.467_rp
  real(kind=rp), parameter :: z_ref = 51.0_rp
  real(kind=rp), parameter :: z0 = 1.0_rp
  real(kind=rp), parameter :: alpha = 0.16_rp
  real(kind=rp), parameter :: wind_from_deg = 225.0_rp
  real(kind=rp), parameter :: deg_to_rad = pi / 180.0_rp
  real(kind=rp), parameter :: flow_to_deg = modulo(wind_from_deg + 180.0_rp, 360.0_rp)
  real(kind=rp), parameter :: wind_x = sin(flow_to_deg * deg_to_rad)
  real(kind=rp), parameter :: wind_y = cos(flow_to_deg * deg_to_rad)
  real(kind=rp) :: ramp_time = 0.0_rp
  real(kind=rp) :: initial_wind_scale = 1.0_rp
  real(kind=rp) :: inlet_taper_deg = 20.0_rp
  real(kind=rp) :: inlet_arc_width_deg = 180.0_rp
  character(len=256) :: initial_field = ""

contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    ramp_time = read_env_real("SODERMALM_RAMP_TIME", 0.0_rp)
    initial_wind_scale = read_env_real("SODERMALM_INITIAL_WIND_SCALE", 1.0_rp)
    inlet_taper_deg = read_env_real("SODERMALM_INLET_TAPER_DEG", 20.0_rp)
    inlet_arc_width_deg = read_env_real("SODERMALM_INLET_ARC_WIDTH_DEG", 180.0_rp)
    call get_environment_variable("SODERMALM_INITIAL_FIELD", initial_field)

    user%initial_conditions => initial_conditions
    user%dirichlet_conditions => dirichlet_conditions
  end subroutine user_setup

  function read_env_real(name, default_value) result(value)
    character(len=*), intent(in) :: name
    real(kind=rp), intent(in) :: default_value
    real(kind=rp) :: value
    character(len=128) :: raw
    integer :: raw_len, status

    value = default_value
    call get_environment_variable(name, raw, length = raw_len, status = status)
    if (status .eq. 0 .and. raw_len .gt. 0) then
       read(raw(1:raw_len), *, iostat = status) value
       if (status .ne. 0) value = default_value
    end if
  end function read_env_real

  function wind_speed(z) result(speed)
    real(kind=rp), intent(in) :: z
    real(kind=rp) :: speed
    real(kind=rp) :: z_eff

    z_eff = max(z, z0)
    speed = u_ref * (z_eff / z_ref) ** alpha
  end function wind_speed

  function wind_ramp(t) result(scale)
    real(kind=rp), intent(in) :: t
    real(kind=rp) :: scale
    real(kind=rp) :: s

    if (ramp_time <= 0.0_rp) then
       scale = 1.0_rp
    else
       s = min(max(t / ramp_time, 0.0_rp), 1.0_rp)
       scale = s * s * (3.0_rp - 2.0_rp * s)
    end if
  end function wind_ramp

  function circular_distance_deg(a, b) result(distance)
    real(kind=rp), intent(in) :: a, b
    real(kind=rp) :: distance

    distance = abs(modulo(a - b + 180.0_rp, 360.0_rp) - 180.0_rp)
  end function circular_distance_deg

  function smooth01(s) result(value)
    real(kind=rp), intent(in) :: s
    real(kind=rp) :: value
    real(kind=rp) :: x

    x = min(max(s, 0.0_rp), 1.0_rp)
    value = x * x * (3.0_rp - 2.0_rp * x)
  end function smooth01

  function inlet_edge_ramp(x, y) result(scale)
    real(kind=rp), intent(in) :: x, y
    real(kind=rp) :: scale
    real(kind=rp) :: bearing, distance_from_center, edge_distance

    if (inlet_taper_deg <= 0.0_rp) then
       scale = 1.0_rp
       return
    end if

    bearing = modulo(atan2(x, y) / deg_to_rad + 360.0_rp, 360.0_rp)
    distance_from_center = circular_distance_deg(bearing, wind_from_deg)
    edge_distance = 0.5_rp * inlet_arc_width_deg - distance_from_center
    scale = smooth01(edge_distance / inlet_taper_deg)
  end function inlet_edge_ramp

  subroutine set_profile(u, v, w, idx)
    type(field_t), intent(inout) :: u, v, w
    integer, intent(in) :: idx
    real(kind=rp) :: speed

    speed = wind_speed(u%dof%z(idx, 1, 1, 1))
    u%x(idx, 1, 1, 1) = speed * wind_x
    v%x(idx, 1, 1, 1) = speed * wind_y
    w%x(idx, 1, 1, 1) = 0.0_rp
  end subroutine set_profile

  subroutine set_profile_scaled(u, v, w, idx, scale)
    type(field_t), intent(inout) :: u, v, w
    integer, intent(in) :: idx
    real(kind=rp), intent(in) :: scale

    call set_profile(u, v, w, idx)
    u%x(idx, 1, 1, 1) = scale * u%x(idx, 1, 1, 1)
    v%x(idx, 1, 1, 1) = scale * v%x(idx, 1, 1, 1)
    w%x(idx, 1, 1, 1) = scale * w%x(idx, 1, 1, 1)
  end subroutine set_profile_scaled

  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(field_t), pointer :: u, v, w, p
    integer :: i

    if (trim(scheme_name) .ne. "fluid") return

    u => fields%get_by_name("u")
    v => fields%get_by_name("v")
    w => fields%get_by_name("w")
    p => fields%get_by_name("p")

    if (len_trim(initial_field) .gt. 0) then
       call import_fields(trim(initial_field), u = u, v = v, w = w, &
            p = p, interpolate = .false.)
       return
    end if

    do i = 1, u%size()
       call set_profile(u, v, w, i)
       u%x(i, 1, 1, 1) = initial_wind_scale * u%x(i, 1, 1, 1)
       v%x(i, 1, 1, 1) = initial_wind_scale * v%x(i, 1, 1, 1)
       w%x(i, 1, 1, 1) = initial_wind_scale * w%x(i, 1, 1, 1)
    end do
    p%x = 0.0_rp

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call u%copy_from(HOST_TO_DEVICE, sync = .false.)
       call v%copy_from(HOST_TO_DEVICE, sync = .false.)
       call w%copy_from(HOST_TO_DEVICE, sync = .false.)
       call p%copy_from(HOST_TO_DEVICE, sync = .true.)
    end if
  end subroutine initial_conditions

  subroutine dirichlet_conditions(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, v, w
    integer :: i, idx
    real(kind=rp) :: scale

    u => fields%get_by_name("u")
    v => fields%get_by_name("v")
    w => fields%get_by_name("w")

    scale = wind_ramp(time%t)
    do i = 1, bc%msk(0)
       idx = bc%msk(i)
       call set_profile_scaled(u, v, w, idx, scale * &
            inlet_edge_ramp(u%dof%x(idx, 1, 1, 1), u%dof%y(idx, 1, 1, 1)))
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call u%copy_from(HOST_TO_DEVICE, sync = .false.)
       call v%copy_from(HOST_TO_DEVICE, sync = .false.)
       call w%copy_from(HOST_TO_DEVICE, sync = .true.)
    end if
  end subroutine dirichlet_conditions

end module user
