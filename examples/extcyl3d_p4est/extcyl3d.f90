! Two-dimensional flow past circular cylinder cavity
!
! Note that the domain is actually 3D with width one element. In order
! to prevent any instability in the z direction, the w velocity is
! set to zero at every step. This is needed for higher Reynolds numbers.
!
module user
  use neko
  use amr_reconstruct, only : amr_reconstruct_t
  implicit none

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%compute => compute
    user%amr_refine_flag => amr_refine_flag
    user%amr_reconstruct => amr_reconstruct
  end subroutine user_setup

  ! User-defined routine called at the end of every time step
  subroutine compute(time)
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: w

    ! set the w component to zero to avoid any 3D instability
    ! in this quasi-2D flow
    w => neko_registry%get_field("w")
    call field_rzero(w)

  end subroutine compute

  subroutine amr_refine_flag(time, ref_mark, ifrefine)
    type(time_state_t), intent(in) :: time
    integer, dimension(:), intent(out) :: ref_mark
    logical, intent(out) :: ifrefine

    ref_mark(:) = 1
    ifrefine = .true.

  end subroutine amr_refine_flag

  subroutine amr_reconstruct(reconstruct, counter)
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter

  end subroutine amr_reconstruct

end module user
