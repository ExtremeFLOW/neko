! Euler 2D compressible flow over a cylinder
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

    do i = 1, rho%dof%size()
       u%x(i,1,1,1) = 1.1
       v%x(i,1,1,1) = 0.0
       w%x(i,1,1,1) = 0.0
       rho%x(i,1,1,1) = 1.4
       p%x(i,1,1,1) = 1.0
    end do
  end subroutine user_ic

end module user
