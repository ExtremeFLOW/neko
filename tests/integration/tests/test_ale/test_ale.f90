! User file for the ALE integration test.
!
! The kinematics below are added to the built-in oscillation/rotation
! configured in the "ale" section of the case file. The resulting motion is
! asserted against a closed form in test_oscillating_cylinder.py, so the
! amplitudes and frequencies here must stay in sync with the TRANSLATION and
! ROTATION_Z_DEG constants in that file.
!
!   translation x:  0.10 (case) - 0.15 (here) = -0.05 at f = 0.1
!   translation y:  0.30 (case) - 0.35 (here) = -0.05 at f = 0.2
!   rotation   z:   2.0  (case) + 3.0  (here) =  5.0 deg at f = 0.1
!                                            +  2.0 deg at f = 0.3
module user
  use neko
  implicit none

contains

  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    user%ale_rigid_kinematics => user_rigid_kinematics
    user%compute => user_check

  end subroutine user_setup

  !> Modify (add to / override) the built-in rigid body kinematics.
  subroutine user_rigid_kinematics(body_id, time, vel_trans, vel_ang)
    integer, intent(in) :: body_id
    type(time_state_t), intent(in) :: time
    real(kind=rp), intent(inout) :: vel_trans(3)
    real(kind=rp), intent(inout) :: vel_ang(3)
    real(kind=rp) :: t

    t = time%t

    ! body_id is the index of the body in the "bodies" array of the case file,
    ! NOT the zone_id of the moving boundary condition.
    select case (body_id)

    case (1) ! "cylinder", the only body registered in the case file.
       vel_trans(1) = vel_trans(1) - &
            0.15_rp * 2 * pi * 0.1_rp * cos(2*pi*0.1_rp*t)

       vel_trans(2) = vel_trans(2) - &
            0.35_rp * 2 * pi * 0.2_rp * cos(2*pi*0.2_rp*t)

       vel_ang(3) = vel_ang(3) + &
            (3.0_rp * pi / 180.0_rp) * 2 * pi * 0.1_rp * cos(2*pi*0.1_rp*t) + &
            (2.0_rp * pi / 180.0_rp) * 2 * pi * 0.3_rp * cos(2*pi*0.3_rp*t)

    end select

  end subroutine user_rigid_kinematics

  subroutine user_check(time)
    type(time_state_t), intent(in) :: time
    integer :: ids_to_log(1)

    ! Log the rotation angle and pivot state of the first (only) body. The
    ! optional argument restricts logging to the listed bodies; omitting it
    ! logs all of them.
    ids_to_log = [1]

    if (associated(neko_ale)) then
       call neko_ale%log_rot_angles(time, ids_to_log)
       call neko_ale%log_pivot(time, ids_to_log)
    end if

  end subroutine user_check
end module user
