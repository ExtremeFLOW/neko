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

    ! rho_L = 1, u_L = (0, 0), p_L = 1
    ! rho_R = 0.125, u_R = (0, 0), p_R = 0.1
    ! end time = 0.2
    theta = 0.01
    mux = 0.5
    do i = 1, rho%dof%size()
      x = rho%dof%x(i,1,1,1)

      u%x(i,1,1,1) = 0.0
      v%x(i,1,1,1) = 0.0
      w%x(i,1,1,1) = 0.0

      if (x < mux) then
        rho%x(i,1,1,1) = 1.0
        p%x(i,1,1,1) = 1.0
      else
        rho%x(i,1,1,1) = 0.125
        p%x(i,1,1,1) = 0.1
      end if

      ! the smoothen initial condition
      ! rho%x(i,1,1,1) = (1-tanh((x-0.5)/theta))/2.0 * (1.0-0.125) + 0.125
      ! p%x(i,1,1,1) = (1-tanh((x-0.5)/theta))/2.0 * (1.0-0.1) + 0.1
    end do
  end subroutine user_ic

end module user