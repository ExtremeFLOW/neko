module user
  use neko
  use math, only : NEKO_EPS
  implicit none

contains
  !> Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initial_conditions => initial_conditions
    user%preprocess         => preprocess
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
    end if

  end subroutine initial_conditions

  subroutine preprocess(time)
    type(time_state_t), intent(in) :: time
    type (field_t), pointer :: u, v, w, s
    integer :: i

    s => neko_registry%get_field("s")

    call s%copy_from(DEVICE_TO_HOST, .true.)

    if (.not. all(abs(s%x - 1.0_rp) < NEKO_EPS)) then
       call neko_error("Incorrect values in s")
    end if
  end subroutine preprocess
end module user
