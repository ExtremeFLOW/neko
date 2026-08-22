! Mach 800 astrophysical jet
!
module user
  use neko
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
    real(kind=rp), parameter :: gamma = 1.4_rp
    real(kind=rp), parameter :: nozzle_width = 0.05_rp
    real(kind=rp) :: x, y, tolerance
    integer :: i

    rho => fields%get_by_name('fluid_rho')
    u => fields%get_by_name('u')
    v => fields%get_by_name('v')
    w => fields%get_by_name('w')
    p => fields%get_by_name('p')

    tolerance = 100.0_rp * epsilon(1.0_rp)
    do i = 1, rho%size()
       x = rho%dof%x(i,1,1,1)
       y = rho%dof%y(i,1,1,1)
       rho%x(i,1,1,1) = 0.1_rp * gamma
       u%x(i,1,1,1) = 0.0_rp
       v%x(i,1,1,1) = 0.0_rp
       w%x(i,1,1,1) = 0.0_rp
       p%x(i,1,1,1) = 1.0_rp

       if (y .le. tolerance .and. x .le. nozzle_width + tolerance) then
          rho%x(i,1,1,1) = gamma
          v%x(i,1,1,1) = 800.0_rp
       end if
    end do
  end subroutine initial_conditions

end module user
