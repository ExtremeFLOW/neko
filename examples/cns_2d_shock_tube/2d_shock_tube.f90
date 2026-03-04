! Reflected shock-boundary layer interaction in a 2D shock tube
!
module user
  use neko
  implicit none

contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initial_conditions => initial_conditions
    user%material_properties => material_properties
  end subroutine user_setup

  ! Riemann problem: diaphragm at x = 0.5
  ! Left  (x < 0.5): rho = 120.0, p = 120.0 / gamma
  ! Right (x > 0.5): rho = 1.2,   p = 1.2   / gamma
  ! Velocity = 0 everywhere
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    type(field_t), pointer :: rho, u, v, w, p
    integer :: i
    real(kind=rp) :: x, gamma

    gamma = 1.4_rp

    rho => fields%get_by_name("fluid_rho")
    u => fields%get_by_name("u")
    v => fields%get_by_name("v")
    w => fields%get_by_name("w")
    p => fields%get_by_name("p")

    do i = 1, rho%dof%size()
       x = rho%dof%x(i, 1, 1, 1)

       u%x(i, 1, 1, 1) = 0.0_rp
       v%x(i, 1, 1, 1) = 0.0_rp
       w%x(i, 1, 1, 1) = 0.0_rp

       if (x < 0.5_rp) then
          rho%x(i, 1, 1, 1) = 120.0_rp
          p%x(i, 1, 1, 1) = 120.0_rp / gamma
       else
          rho%x(i, 1, 1, 1) = 1.2_rp
          p%x(i, 1, 1, 1) = 1.2_rp / gamma
       end if
    end do
  end subroutine initial_conditions

  ! Set dynamic viscosity: mu = 1 / Re
  subroutine material_properties(scheme_name, properties, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: properties
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: mu, kappa
    integer :: i
    real(kind=rp) :: Re, Prandtl_number, gamma

    Re = 1000.0_rp
    Prandtl_number = 0.73_rp
    gamma = 1.4_rp

    if (scheme_name .eq. "fluid") then
       mu => properties%get_by_name("fluid_mu")
       kappa => properties%get_by_name("fluid_kappa")

       do i = 1, mu%dof%size()
          mu%x(i, 1, 1, 1) = 1.0_rp / Re
          kappa%x(i, 1, 1, 1) = mu%x(i, 1, 1, 1) * gamma / Prandtl_number
       end do
    end if
  end subroutine material_properties

end module user
