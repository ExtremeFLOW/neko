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
    real(kind=rp) :: x, y, z, rho0, T0, p0, gamma, &
                  Ma, V0

    call json_get(params, "case.fluid.gamma", gamma)

    V0 = 1.0
    !> Mach number
    Ma = 1.25

    rho0 = 1.0
    p0 = rho0 * V0**2 / (gamma * Ma * Ma)

    do i = 1, rho%dof%size()
      x = rho%dof%x(i,1,1,1)
      y = rho%dof%y(i,1,1,1)
      z = rho%dof%z(i,1,1,1)

      u%x(i,1,1,1) = sin(x)*cos(y)*cos(z)
      v%x(i,1,1,1) = -cos(x)*sin(y)*cos(z)
      w%x(i,1,1,1) = 0.0

      rho%x(i,1,1,1) = rho0
      p%x(i,1,1,1) = p0 + rho0 * V0**2.0 / 16.0 * (cos(2.0*x) + cos(2.0*y)) &
                    * (2.0 + cos(2.0 * z))
    end do
  end subroutine user_ic

end module user