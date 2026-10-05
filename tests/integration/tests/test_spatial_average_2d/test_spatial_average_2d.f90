! Overwrites u, v, w and p after every (frozen) fluid step with polynomial
! fields whose averages in x, y or z are known in closed form, to check the
! averaging in one direction (map_2d) of the spatial_average simcomp.
module user
  use neko
  implicit none
contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%compute => set_fields
  end subroutine user_setup

  subroutine set_fields(time)
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, v, w, p
    integer :: i, n
    real(kind=rp) :: x, y, z

    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")
    p => neko_registry%get_field("p")
    n = u%size()
    do i = 1, n
       x = u%dof%x%x(i,1,1,1)
       y = u%dof%y%x(i,1,1,1)
       z = u%dof%z%x(i,1,1,1)
       u%x(i,1,1,1) = x + 2.0_rp * y + 3.0_rp * z
       v%x(i,1,1,1) = x * y * z
       w%x(i,1,1,1) = x * x * y + z * z * x + 1.0_rp
       p%x(i,1,1,1) = x * x + y * y + z * z
    end do
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(u%x, u%x_d, n, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(v%x, v%x_d, n, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(w%x, w%x_d, n, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(p%x, p%x_d, n, HOST_TO_DEVICE, sync = .true.)
    end if
  end subroutine set_fields

end module user
