! User file for test_fluid_traction_rho.py. Logs the volume average of u at
! every time step.
module user
  use neko
  implicit none

contains

  !> Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%compute => compute
  end subroutine user_setup

  !> Log the volume average of u.
  subroutine compute(time)
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u
    type(coef_t), pointer :: coef
    real(kind=rp) :: u_mean
    character(len=LOG_SIZE) :: log_buf
    integer :: n

    coef => neko_user_access%case%fluid%c_Xh
    u => neko_registry%get_field("u")
    n = u%dof%size()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       u_mean = device_glsc2(u%x_d, coef%B_d, n) / coef%volume
    else
       u_mean = glsc2(u%x, coef%B, n) / coef%volume
    end if

    write(log_buf, '(A,ES24.16,1X,ES24.16)') 'Mean u: ', time%t, u_mean
    call neko_log%message(log_buf)
  end subroutine compute

end module user
