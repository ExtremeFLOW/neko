! Lid-driven cavity
!
! Time-integration of the lid-driven cavity with smoothened
! belt velocity to fulfil continuity equation.
!
module user
  use neko
  implicit none

  ! Global user variables
  type(field_t) :: w1

  type(file_t) output_file ! output file
  type(vector_t) :: vec_out    ! will store our output data

 contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%fluid_compressible_user_ic => user_ic
  end subroutine user_setup
  
  subroutine user_ic(rho, u, v, w, p, params)
    type(field_t), intent(inout) :: rho
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(json_file), intent(inout) :: params
    integer :: i
    real(kind=rp) :: x, y, cone_radius, mux, muy, r, theta

    do i = 1, rho%dof%size()
      x = rho%dof%x(i,1,1,1)
      y = rho%dof%y(i,1,1,1)

      u%x(i,1,1,1) = 2.5
      v%x(i,1,1,1) = -0.5
      w%x(i,1,1,1) = 0.0

      rho%x(i,1,1,1) = 1.0 + 0.2 * sin(2 * pi * (x + y))
      p%x(i,1,1,1) = 1.0
    end do
  end subroutine user_ic

end module user