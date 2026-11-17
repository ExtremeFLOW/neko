! User file for test_source_term_time.py. Adds a spatially uniform forcing
! f(t) = cos(omega t) to the x-momentum and scalar equations.
module user
  use neko
  implicit none

  !> Angular frequency of the forcing.
  real(kind=rp), parameter :: omega = 1.0_rp

contains

  !> Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%source_term => forcing
  end subroutine user_setup

  !> Set the first right-hand side component to cos(omega t). For the fluid
  !! this is the x-component, the other two are left at zero.
  subroutine forcing(scheme_name, rhs, time, eval_time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: rhs
    type(time_state_t), intent(in) :: time
    real(kind=dp), intent(in) :: eval_time

    call field_cfill(rhs%items(1)%ptr, real(cos(omega * eval_time), kind=rp))
  end subroutine forcing

end module user
