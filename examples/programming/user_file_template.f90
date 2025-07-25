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
    user%user_dirichlet_update => user_dirichlet_update
    user%material_properties => material_properties

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

  subroutine user_dirichlet_update(dirichlet_field_list, dirichlet_bc, coef, t,&
       tstep)
    type(field_list_t), intent(inout) :: dirichlet_field_list
    type(field_dirichlet_t), intent(in) :: dirichlet_bc
    type(coef_t), intent(inout) :: coef
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

  end subroutine user_dirichlet_update

  subroutine material_properties(scheme_name, properties, time)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: properties
    type(time_state_t), intent(in) :: time

  end subroutine material_properties

end module user
