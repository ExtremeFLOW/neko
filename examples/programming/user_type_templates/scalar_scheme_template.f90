!> Template for a user-defined scalar scheme.
module scalar_scheme_template
  use checkpoint, only : chkp_t
  use coefs, only : coef_t
  use field, only : field_t
  use field_series, only : field_series_t
  use gather_scatter, only : gs_t
  use json_module, only : json_file
  use krylov, only : ksp_monitor_t
  use mesh, only : mesh_t
  use scalar_scheme, only : scalar_scheme_t, scalar_scheme_allocate, &
       register_scalar_scheme
  use time_scheme_controller, only : time_scheme_controller_t
  use time_state, only : time_state_t
  use time_step_controller, only : time_step_controller_t
  use user_intf, only : user_t
  use utils, only : neko_error
  implicit none
  private

  type, extends(scalar_scheme_t) :: scalar_scheme_template_t
   contains
     procedure :: init => scalar_scheme_template_init
     procedure :: free => scalar_scheme_template_free
     procedure :: step => scalar_scheme_template_step
     procedure :: restart => scalar_scheme_template_restart
  end type scalar_scheme_template_t

  public :: scalar_scheme_template_register_types

contains

  subroutine scalar_scheme_template_init(this, msh, coef, gs, params, &
       numerics_params, user, chkp, ulag, vlag, wlag, time_scheme, rho)
    class(scalar_scheme_template_t), target, intent(inout) :: this
    type(mesh_t), target, intent(in) :: msh
    type(coef_t), target, intent(in) :: coef
    type(gs_t), target, intent(inout) :: gs
    type(json_file), target, intent(inout) :: params, numerics_params
    type(user_t), target, intent(in) :: user
    type(chkp_t), target, intent(inout) :: chkp
    type(field_series_t), target, intent(in) :: ulag, vlag, wlag
    type(time_scheme_controller_t), target, intent(in) :: time_scheme
    type(field_t), target, intent(in) :: rho

    ! `scheme_init` sets up common fields, material properties, source terms,
    ! and the configured Krylov solver/preconditioner.
    call this%scheme_init(msh, coef, gs, params, 'scalar_scheme_template', &
         user, rho)
    ! TODO: Construct this scheme's operators and use numerics_params,
    !       lagged velocities, time_scheme, and chkp as needed.
  end subroutine scalar_scheme_template_init

  subroutine scalar_scheme_template_free(this)
    class(scalar_scheme_template_t), intent(inout) :: this
    ! TODO: Release operators and state owned by this derived type.
    call this%scheme_free()
  end subroutine scalar_scheme_template_free

  subroutine scalar_scheme_template_step(this, time, ext_bdf, dt_controller, &
       ksp_results)
    class(scalar_scheme_template_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_scheme_controller_t), intent(in) :: ext_bdf
    type(time_step_controller_t), intent(in) :: dt_controller
    type(ksp_monitor_t), intent(inout) :: ksp_results

    ! TODO: Assemble and solve one scalar time step.
    call neko_error('scalar_scheme_template: step is not implemented')
  end subroutine scalar_scheme_template_step

  subroutine scalar_scheme_template_restart(this, chkp)
    class(scalar_scheme_template_t), target, intent(inout) :: this
    type(chkp_t), intent(inout) :: chkp
    ! TODO: Restore scheme-specific checkpoint state.
  end subroutine scalar_scheme_template_restart

  subroutine scalar_scheme_template_register_types()
    procedure(scalar_scheme_allocate), pointer :: allocator
    allocator => scalar_scheme_template_allocate
    call register_scalar_scheme('scalar_scheme_template', allocator)
  end subroutine scalar_scheme_template_register_types

  subroutine scalar_scheme_template_allocate(obj)
    class(scalar_scheme_t), allocatable, intent(inout) :: obj
    allocate(scalar_scheme_template_t :: obj)
  end subroutine scalar_scheme_template_allocate
end module scalar_scheme_template
