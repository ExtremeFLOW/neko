! Exact solution for the smooth wave propagation problem
! with periodic boundary conditions
!
module user
  use neko
  implicit none

  ! Global user variables
  type(field_t) :: w1

  type(file_t) output_file ! output file
  type(vector_t) :: vec_out ! will store our output data

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initial_conditions => initial_conditions
  end subroutine user_setup

  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    type (field_t), pointer :: rho, u, v, w, p
    integer :: i
    real(kind=rp) :: x, y, cone_radius, mux, muy, r, theta

    rho => fields%get_by_name("fluid_rho")
    u => fields%get_by_name("u")
    v => fields%get_by_name("v")
    w => fields%get_by_name("w")
    p => fields%get_by_name("p")

    do i = 1, rho%dof%size()
       x = rho%dof%x(i,1,1,1)
       y = rho%dof%y(i,1,1,1)

       u%x(i,1,1,1) = 2.5
       v%x(i,1,1,1) = -0.5
       w%x(i,1,1,1) = 0.0

       rho%x(i,1,1,1) = 1.0 + 0.2 * sin(2 * pi * (x + y))
       p%x(i,1,1,1) = 1.0
    end do
  end subroutine initial_conditions

end module user
