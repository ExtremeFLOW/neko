module user
  use neko
  implicit none

contains
  !> Register the user-defined initial condition callback.
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    user%initial_conditions => initial_conditions
  end subroutine user_setup

  !> Initialize linear fields with analytically known plane averages.
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    type(field_t), pointer :: u, v, w, s
    integer :: i, n
    real(kind=rp) :: x, y, z

    if (scheme_name .eq. 'fluid') then
       u => fields%get("u")
       v => fields%get("v")
       w => fields%get("w")

       n = u%size()
       do i = 1, n
          x = u%dof%x(i,1,1,1)
          y = u%dof%y(i,1,1,1)
          z = u%dof%z(i,1,1,1)

          u%x(i,1,1,1) = x + 2.0_rp * y + 3.0_rp * z
          v%x(i,1,1,1) = 0.0_rp
          w%x(i,1,1,1) = 0.0_rp
       end do
    else
       s => fields%get("s")

       n = s%size()
       do i = 1, n
          x = s%dof%x(i,1,1,1)
          y = s%dof%y(i,1,1,1)
          z = s%dof%z(i,1,1,1)

          s%x(i,1,1,1) = 2.0_rp * x - y + 0.5_rp * z
       end do
    end if
  end subroutine initial_conditions

end module user
