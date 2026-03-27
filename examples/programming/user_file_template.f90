! A template user file containing the user-defined functions
!
module user
  use neko
  implicit none

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user

    user%startup => startup
    user%initialize => initialize
    user%initial_conditions => initial_conditions
    user%mesh_setup => mesh_setup
    user%compute => compute
    user%finalize => finalize
    user%source_term => source_term
    user%dirichlet_conditions => dirichlet_conditions
    user%material_properties => material_properties
    user%ale_mesh_velocity => user_ale_mesh_motion
    user%ale_rigid_kinematics => user_rigid_kinematics
    user%ale_base_shapes => user_ale_base_shapes


  end subroutine user_setup

  subroutine startup(params)
    type(json_file), intent(inout) :: params

  end subroutine startup

  subroutine initialize(time)
    type(time_state_t), intent(in) :: time

  end subroutine initialize

  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

  end subroutine initial_conditions

  subroutine mesh_setup(msh, time)
    type(mesh_t), intent(inout) :: msh
    type(time_state_t), intent(in) :: time

  end subroutine mesh_setup

  subroutine compute(time)
    type(time_state_t), intent(in) :: time
  end subroutine compute

  subroutine finalize(time)
    type(time_state_t), intent(in) :: time

  end subroutine finalize

  subroutine source_term(scheme_name, rhs, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: rhs
    type(time_state_t), intent(in) :: time

  end subroutine source_term

  subroutine dirichlet_conditions(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time

  end subroutine dirichlet_conditions

  subroutine material_properties(scheme_name, properties, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: properties
    type(time_state_t), intent(in) :: time

  end subroutine material_properties

  subroutine user_ale_mesh_motion(wm_x, wm_y, wm_z, coef, &
       x_ref, y_ref, z_ref, base_shapes, time)
    type(coef_t), intent(in) :: coef
    type(field_t), intent(inout) :: wm_x, wm_y, wm_z
    type(field_t), intent(in) :: x_ref, y_ref, z_ref
    type(field_t), intent(in) :: base_shapes(:)
    type(time_state_t), intent(in) :: time

  end subroutine user_ale_mesh_motion

  subroutine user_rigid_kinematics(body_id, time, vel_trans, vel_ang)
    integer, intent(in) :: body_id
    type(time_state_t), intent(in) :: time
    real(kind=rp), intent(inout) :: vel_trans(3)
    real(kind=rp), intent(inout) :: vel_ang(3)
    real(kind=rp) :: t

  end subroutine user_rigid_kinematics

  subroutine user_ale_base_shapes(base_shapes)
    type(field_t), intent(inout) :: base_shapes(:)

  end subroutine user_ale_base_shapes

end module user
