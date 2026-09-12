module user
  use neko
  use fluid_pnpn, only : fluid_pnpn_t
  use import_field_utils, only : import_fields
  use mpi_f08, only : MPI_Allreduce, MPI_Bcast, MPI_INTEGER, MPI_MAX, &
       MPI_MIN, MPI_SUM
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
  real(kind=rp) :: sponge_radius_m = 2536.52_rp
  real(kind=rp) :: sponge_thickness_m = 800.0_rp
  real(kind=rp) :: sponge_rise_m = 500.0_rp
  real(kind=rp) :: diagnostic_interval = 0.0_rp
  real(kind=rp) :: diagnostic_start_time = 0.0_rp
  real(kind=rp) :: next_diagnostic_time = huge(0.0_rp)
  character(len=256) :: initial_field = ""

contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    ramp_time = read_env_real("SODERMALM_RAMP_TIME", 0.0_rp)
    initial_wind_scale = read_env_real("SODERMALM_INITIAL_WIND_SCALE", 1.0_rp)
    inlet_taper_deg = read_env_real("SODERMALM_INLET_TAPER_DEG", 20.0_rp)
    inlet_arc_width_deg = read_env_real("SODERMALM_INLET_ARC_WIDTH_DEG", 180.0_rp)
    sponge_radius_m = read_env_real("SODERMALM_SPONGE_RADIUS_M", 2536.52_rp)
    sponge_thickness_m = read_env_real("SODERMALM_SPONGE_THICKNESS_M", 800.0_rp)
    sponge_rise_m = read_env_real("SODERMALM_SPONGE_RISE_M", 500.0_rp)
    diagnostic_interval = read_env_real( &
         "SODERMALM_DIAGNOSTIC_INTERVAL", 0.0_rp)
    diagnostic_start_time = read_env_real( &
         "SODERMALM_DIAGNOSTIC_START_TIME", 0.0_rp)
    call get_environment_variable("SODERMALM_INITIAL_FIELD", initial_field)

    user%initialize => user_initialize
    user%initial_conditions => initial_conditions
    user%dirichlet_conditions => dirichlet_conditions
    if (diagnostic_interval .gt. 0.0_rp) then
       user%compute => flow_diagnostics
    end if
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
       scale = initial_wind_scale + (1.0_rp - initial_wind_scale) * &
            s * s * (3.0_rp - 2.0_rp * s)
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

  subroutine user_initialize(time)
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, fringe, ubf, vbf, wbf
    integer :: i
    real(kind=rp) :: x, y, z, downstream, dist_from_upstream_edge, rise

    u => neko_registry%get_field("u")
    call neko_registry%add_field(u%dof, "sponge_fringe")
    call neko_registry%add_field(u%dof, "sponge_bf_u")
    call neko_registry%add_field(u%dof, "sponge_bf_v")
    call neko_registry%add_field(u%dof, "sponge_bf_w")
    fringe => neko_registry%get_field("sponge_fringe")
    ubf => neko_registry%get_field("sponge_bf_u")
    vbf => neko_registry%get_field("sponge_bf_v")
    wbf => neko_registry%get_field("sponge_bf_w")

    fringe%x = 0.0_rp
    rise = max(sponge_rise_m, 1.0_rp)
    do i = 1, fringe%size()
       x = fringe%dof%x(i, 1, 1, 1)
       y = fringe%dof%y(i, 1, 1, 1)
       z = fringe%dof%z(i, 1, 1, 1)
       downstream = x * wind_x + y * wind_y
       dist_from_upstream_edge = downstream + sponge_radius_m

       if (dist_from_upstream_edge >= 0.0_rp .and. &
            dist_from_upstream_edge <= sponge_thickness_m) then
          fringe%x(i, 1, 1, 1) = smooth01( &
               (sponge_thickness_m - dist_from_upstream_edge) / rise)
       end if

       ubf%x(i, 1, 1, 1) = wind_speed(z) * wind_x
       vbf%x(i, 1, 1, 1) = wind_speed(z) * wind_y
       wbf%x(i, 1, 1, 1) = 0.0_rp
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call fringe%copy_from(HOST_TO_DEVICE, sync = .false.)
       call ubf%copy_from(HOST_TO_DEVICE, sync = .false.)
       call vbf%copy_from(HOST_TO_DEVICE, sync = .false.)
       call wbf%copy_from(HOST_TO_DEVICE, sync = .true.)
    end if

    if (diagnostic_interval .gt. 0.0_rp) then
       next_diagnostic_time = max(real(time%t, rp), diagnostic_start_time)
    end if

    nullify(u)
    nullify(fringe)
    nullify(ubf)
    nullify(vbf)
    nullify(wbf)
  end subroutine user_initialize

  subroutine flow_diagnostics(time)
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, v, w, p, indicator, permeability
    type(field_t), pointer :: divergence, p_res
    type(coef_t), pointer :: coef
    real(kind=rp) :: local_max2, global_max2, speed2
    real(kind=rp) :: local_drag_power, global_drag_power
    real(kind=rp) :: local_div_max, global_div_max
    real(kind=rp) :: local_div_l2, global_div_l2
    real(kind=rp) :: local_p_max, global_p_max
    real(kind=rp) :: local_pres_max, global_pres_max
    real(kind=rp) :: local_pres_l2, global_pres_l2
    real(kind=rp) :: details(11)
    integer :: i, local_idx, owner_candidate, owner, ierr, divergence_idx
    logical :: has_brinkman

    if (diagnostic_interval .le. 0.0_rp) return
    if (real(time%t, rp) + 0.5_rp * real(time%dt, rp) .lt. &
         next_diagnostic_time) return

    do while (next_diagnostic_time .le. real(time%t, rp) + &
         0.5_rp * real(time%dt, rp))
       next_diagnostic_time = next_diagnostic_time + diagnostic_interval
    end do

    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")
    p => neko_registry%get_field("p")
    has_brinkman = neko_registry%field_exists("brinkman_indicator") .and. &
         neko_registry%field_exists("brinkman_permeability")
    if (has_brinkman) then
       indicator => neko_registry%get_field("brinkman_indicator")
       permeability => neko_registry%get_field("brinkman_permeability")
    end if
    coef => neko_user_access%case%fluid%c_Xh
    nullify(p_res)
    select type (fluid => neko_user_access%case%fluid)
    type is (fluid_pnpn_t)
       p_res => fluid%p_res
    end select

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call u%copy_from(DEVICE_TO_HOST, sync = .false.)
       call v%copy_from(DEVICE_TO_HOST, sync = .false.)
       call w%copy_from(DEVICE_TO_HOST, sync = .false.)
       call p%copy_from(DEVICE_TO_HOST, sync = .false.)
       if (associated(p_res)) then
          call p_res%copy_from(DEVICE_TO_HOST, sync = .false.)
       end if
       if (has_brinkman) then
          call indicator%copy_from(DEVICE_TO_HOST, sync = .false.)
          call permeability%copy_from(DEVICE_TO_HOST, sync = .true.)
       else
          call p%copy_from(DEVICE_TO_HOST, sync = .true.)
       end if
    end if

    call neko_scratch_registry%request_field(divergence, divergence_idx, &
         .false.)
    call div(divergence%x, u%x, v%x, w%x, coef)

    local_max2 = -huge(0.0_rp)
    local_idx = 1
    local_drag_power = 0.0_rp
    local_div_max = 0.0_rp
    local_div_l2 = 0.0_rp
    local_p_max = 0.0_rp
    local_pres_max = 0.0_rp
    local_pres_l2 = 0.0_rp
    do i = 1, u%size()
       speed2 = u%x(i,1,1,1)**2 + v%x(i,1,1,1)**2 + &
            w%x(i,1,1,1)**2
       if (speed2 .gt. local_max2) then
          local_max2 = speed2
          local_idx = i
       end if
       if (has_brinkman) then
          local_drag_power = local_drag_power - permeability%x(i,1,1,1) * &
               speed2 * coef%B(i,1,1,1)
       end if
       local_div_max = max(local_div_max, abs(divergence%x(i,1,1,1)))
       local_div_l2 = local_div_l2 + divergence%x(i,1,1,1)**2 * &
            coef%B(i,1,1,1)
       local_p_max = max(local_p_max, abs(p%x(i,1,1,1)))
       if (associated(p_res)) then
          local_pres_max = max(local_pres_max, abs(p_res%x(i,1,1,1)))
          local_pres_l2 = local_pres_l2 + p_res%x(i,1,1,1)**2 * &
               coef%B(i,1,1,1)
       end if
    end do

    call MPI_Allreduce(local_max2, global_max2, 1, MPI_REAL_PRECISION, &
         MPI_MAX, NEKO_COMM, ierr)
    owner_candidate = huge(0)
    if (local_max2 .eq. global_max2) owner_candidate = pe_rank
    call MPI_Allreduce(owner_candidate, owner, 1, MPI_INTEGER, MPI_MIN, &
         NEKO_COMM, ierr)
    call MPI_Allreduce(local_drag_power, global_drag_power, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    call MPI_Allreduce(local_div_max, global_div_max, 1, &
         MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM, ierr)
    call MPI_Allreduce(local_div_l2, global_div_l2, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)
    call MPI_Allreduce(local_p_max, global_p_max, 1, &
         MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM, ierr)
    call MPI_Allreduce(local_pres_max, global_pres_max, 1, &
         MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM, ierr)
    call MPI_Allreduce(local_pres_l2, global_pres_l2, 1, &
         MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM, ierr)

    details = 0.0_rp
    if (pe_rank .eq. owner) then
       details(1) = u%dof%x(local_idx,1,1,1)
       details(2) = u%dof%y(local_idx,1,1,1)
       details(3) = u%dof%z(local_idx,1,1,1)
       details(4) = u%x(local_idx,1,1,1)
       details(5) = v%x(local_idx,1,1,1)
       details(6) = w%x(local_idx,1,1,1)
       details(7) = p%x(local_idx,1,1,1)
       if (has_brinkman) then
          details(8) = indicator%x(local_idx,1,1,1)
          details(9) = permeability%x(local_idx,1,1,1)
       end if
       details(10) = coef%jac(local_idx,1,1,1)
       details(11) = real(local_idx, rp)
    end if
    call MPI_Bcast(details, size(details), MPI_REAL_PRECISION, owner, &
         NEKO_COMM, ierr)

    if (pe_rank .eq. 0) then
       write(*,'(A,ES15.7,A,I0,A,ES15.7,A,I0)') &
            'BRINKMAN_DIAG t=', time%t, ' step=', time%tstep, &
            ' max_u=', sqrt(global_max2), ' rank=', owner
       write(*,'(A,3(ES15.7,1X),A,3(ES15.7,1X))') &
            'BRINKMAN_DIAG xyz=', details(1:3), ' uvw=', details(4:6)
       write(*,'(A,5(ES15.7,1X),A,I0)') &
            'BRINKMAN_DIAG p_ind_alpha_jac_power=', details(7:10), &
            global_drag_power, ' local_idx=', nint(details(11))
       write(*,'(A,5(ES15.7,1X))') &
            'BRINKMAN_DIAG div_max_div_rms_p_max_pres_max_pres_rms=', &
            global_div_max, sqrt(global_div_l2 / coef%volume), &
            global_p_max, global_pres_max, &
            sqrt(global_pres_l2 / coef%volume)
    end if

    call neko_scratch_registry%relinquish_field(divergence_idx)
    nullify(u, v, w, p, coef, divergence, p_res)
    if (has_brinkman) nullify(indicator, permeability)
  end subroutine flow_diagnostics

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
