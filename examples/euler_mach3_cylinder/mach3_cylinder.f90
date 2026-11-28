! Mach-3 Euler flow around a cylinder.
!
! The primitive state matches the Ryujin reference benchmark:
! density = 1.4, velocity = (3, 0, 0), and pressure = 1.
module user
  use neko
  use device, only : HOST_TO_DEVICE
  implicit none

contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    user%initial_conditions => initial_conditions
  end subroutine user_setup

  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    type(field_t), pointer :: rho, u, v, w, p
    integer :: i

    rho => fields%get_by_name("fluid_rho")
    u => fields%get_by_name("u")
    v => fields%get_by_name("v")
    w => fields%get_by_name("w")
    p => fields%get_by_name("p")

    do i = 1, rho%size()
       rho%x(i,1,1,1) = 1.4_rp
       u%x(i,1,1,1) = 3.0_rp
       v%x(i,1,1,1) = 0.0_rp
       w%x(i,1,1,1) = 0.0_rp
       p%x(i,1,1,1) = 1.0_rp
    end do
    call fields%copy_from(HOST_TO_DEVICE, sync = .true.)
  end subroutine initial_conditions

end module user
