! Two-dimensional flow past circular cylinder cavity
!
! Note that the domain is actually 3D with width one element. In order
! to prevent any instability in the z direction, the w velocity is
! set to zero at every step. This is needed for higher Reynolds numbers.
!
module user
  use neko
  implicit none

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%compute => compute
  end subroutine user_setup

  ! User-defined routine called at the end of every time step
  subroutine compute(time)
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: w

    ! set the w component to zero to avoid any 3D instability
    ! in this quasi-2D flow
    w => neko_field_registry%get_field("w")
    call field_rzero(w)

  end subroutine compute

end module user
