module user
  use neko
  implicit none

contains
  !> Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initial_conditions => initial_conditions
    user%compute => compute
  end subroutine user_setup

  !> User initial condition for the scalar
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    type (field_t), pointer :: u, v, w, s

    if (scheme_name .eq. 'fluid') then
       u => fields%get("u")
       v => fields%get("v")
       w => fields%get("w")

       call field_rzero(u)
       call field_rzero(v)
       call field_rzero(w)
    else
       s => fields%get("s")

       call field_cfill(s, 0.5_rp)
    end if

  end subroutine initial_conditions

  subroutine compute(time)
    type(time_state_t), intent(in) :: time
    type (field_t), pointer :: s
    integer :: i

    s => neko_field_registry%get_field("s")

    do i = 1, s%size()
       call random_number(s%x(i,1,1,1))
    end do
  end subroutine compute
end module user
