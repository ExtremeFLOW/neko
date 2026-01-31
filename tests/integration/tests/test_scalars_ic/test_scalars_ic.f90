module user
  use neko
  use math, only : NEKO_EPS
  use json_module, only : json_file
  implicit none

  integer :: n_scalars = 0

contains
  !> Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initial_conditions => initial_conditions
    user%preprocess => preprocess
    user%startup => startup
  end subroutine user_setup

  !> User initial condition for the scalar
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    type (field_t), pointer :: u, v, w, s


    if (scheme_name .eq. "fluid") then
       u => fields%get("u")
       v => fields%get("v")
       w => fields%get("w")

       call field_cfill(u, 0.0_rp)
       call field_cfill(v, 0.0_rp)
       call field_cfill(w, 0.0_rp)
    else if (scheme_name .eq. "s") then
       s => fields%get("s")
       call field_cfill(s, 1.0_rp)
    else if (scheme_name .eq. "temperature") then
       s => fields%get("temperature")
       call field_cfill(s, 2.0_rp)
    end if

  end subroutine initial_conditions

  !> Startup: determine number of scalars in case file.
  subroutine startup(params)
    type(json_file), intent(inout) :: params

    n_scalars = 0
    if (params%valid_path("case.scalar")) then
       n_scalars = 1
    else if (params%valid_path("case.scalars")) then
       call params%info("case.scalars", n_children = n_scalars)
    end if

    if (n_scalars <= 0) then
       call neko_error("No scalars defined in the case file")
    end if
  end subroutine startup

  !! Checks the set ic values.
  subroutine preprocess(time)
    type(time_state_t), intent(in) :: time
    type (field_t), pointer :: u, v, w, s, t
    integer :: i

    s => neko_registry%get_field("s")
    call s%copy_from(DEVICE_TO_HOST, .true.)

    if (n_scalars .eq. 2) then
       t => neko_registry%get_field("temperature")
       call t%copy_from(DEVICE_TO_HOST, .true.)
    end if

    if (.not. all(abs(s%x - 1.0_rp) < NEKO_EPS)) then
       call neko_error("Incorrect values in s")
    end if

    if (n_scalars .eq. 2) then
       do i = 1, t%dof%size()
          if (.not. abscmp(t%x(i,1,1,1), 2.0_rp)) then
             call neko_error("Incorrect values in temperature")
          end if
       end do
    end if
  end subroutine preprocess
end module user
