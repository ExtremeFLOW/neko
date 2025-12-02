module user
  use neko
  implicit none


real(kind=rp), parameter :: U0 = 1.0_rp
real(kind=rp), parameter :: U1 = 1.0_rp
real(kind=rp), parameter :: V0 = 2.0_rp
real(kind=rp), parameter :: V1 = 2.0_rp
real(kind=rp), parameter :: W0 = 3.0_rp
real(kind=rp), parameter :: W1 = 3.0_rp
real(kind=rp), parameter :: S0 = 4.0_rp
real(kind=rp), parameter :: S1 = 4.0_rp
real(kind=rp), parameter :: omega = 1.0_rp
real(kind=rp), parameter :: phi = 3.14192_rp / 4.0_rp
real(kind=rp), parameter :: psi = 3.14192_rp / 2.0_rp
real(kind=rp), parameter :: xi = 3.14192_rp / 3.0_rp

contains
  !> Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initial_conditions => initial_conditions
    user%compute => compute
  end subroutine user_setup

  !> User initial condition for the scalar
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    type (field_t), pointer :: u, v, w, s

    if (scheme_name .eq. 'fluid') then
       u => fields%get("u")
       v => fields%get("v")
       w => fields%get("w")

       call field_cfill(u, U0 + U1 * cos(omega * 0.0_rp))
       call field_cfill(v, V0 + V1 * cos(omega * 0.0_rp + phi))
       call field_cfill(w, W0 + W1 * cos(omega * 0.0_rp + psi))
    else
       call field_cfill(s, S0 + S1 * sin(omega * 0.0_rp + xi))
    end if

  end subroutine initial_conditions

  subroutine compute(time)
    type(time_state_t), intent(in) :: time
    type (field_t), pointer :: u, v, w, s
    integer :: i

    u => neko_field_registry%get_field("u")
    v => neko_field_registry%get_field("v")
    w => neko_field_registry%get_field("w")
    s => neko_field_registry%get_field("s")

    call field_cfill(u, U0 + U1 * cos(omega * time%t))
    call field_cfill(v, V0 + V1 * cos(omega * time%t + phi))
    call field_cfill(w, W0 + W1 * cos(omega * time%t + psi))
    call field_cfill(s, S0 + S1 * sin(omega * time%t + xi))
  end subroutine compute
end module user
