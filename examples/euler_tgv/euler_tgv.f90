! Taylor-Green vortex (TGV)
module user
  use neko
  implicit none

  ! Global user variables
  type(field_t) :: w1
  real(kind=rp) :: gamma

  type(file_t) output_file ! output file
  type(vector_t) :: vec_out ! will store our output data

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%startup => startup
    user%initial_conditions => initial_conditions
  end subroutine user_setup

  subroutine startup(params)
    type(json_file), intent(inout) :: params

    call json_get(params, "case.fluid.gamma", gamma)
  end subroutine startup

  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    integer :: i

    real(kind=rp) :: x, y, z, rho0, T0, p0, Ma, V0
    type (field_t), pointer :: rho, u, v, w, p

    V0 = 1.0
    !> Mach number
    Ma = 1.25

    rho0 = 1.0
    p0 = rho0 * V0**2 / (gamma * Ma * Ma)

    rho => fields%get_by_name("fluid_rho")
    u => fields%get_by_name("u")
    v => fields%get_by_name("v")
    w => fields%get_by_name("w")
    p => fields%get_by_name("p")

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
  end subroutine initial_conditions

end module user
