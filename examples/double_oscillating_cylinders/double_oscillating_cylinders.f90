module user
  use neko
  implicit none

contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    user%ale_rigid_kinematics => user_rigid_kinematics
    user%compute => user_check

  end subroutine user_setup

  ! Use this to modify (add/override) the built-in rigid body kinematics.
  subroutine user_rigid_kinematics(body_id, time, vel_trans, vel_ang)
    integer, intent(in) :: body_id
    type(time_state_t), intent(in) :: time
    real(kind=rp), intent(inout) :: vel_trans(3)
    real(kind=rp), intent(inout) :: vel_ang(3)
    real(kind=rp) :: t

    t = time%t

    ! here, body_id maps to the order which are bodies are registered in JSON.
    ! THEY ARE NOT THE ZONE_IDS of the moving BCs.
    select case (body_id)

    case (2) ! second registered body in ALE section in case file.
       vel_trans(1) = vel_trans(1) - &
            0.25_rp * 2 * pi * 0.1_rp * cos(2*pi*0.1_rp*t)
       vel_ang(3) = vel_ang(3) + &
            (3.0_rp * pi / 180.0_rp) * 2 * pi * 0.1_rp * cos(2*pi*0.1_rp*t)
    end select

  end subroutine user_rigid_kinematics

  subroutine user_check(time)
    type(time_state_t), intent(in) :: time
    integer :: ids_to_log(1)

    ! Can be used to log the rotation angle and pivot info of bodies.
    ! Here we only log the second registered body.
    ! If the following optional argument is not present, logging will be for all ALE bodies.

    ids_to_log = [2]

    if (associated(neko_ale)) then
       call neko_ale%log_rot_angles(time, ids_to_log)
       call neko_ale%log_pivot(time)
    end if

  end subroutine user_check
end module user

