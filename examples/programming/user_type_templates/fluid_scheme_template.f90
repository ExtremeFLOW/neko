!> Template for a user-defined fluid scheme.
module fluid_scheme_template
  use checkpoint, only : chkp_t
  use fluid_scheme_base, only : fluid_scheme_base_t, fluid_scheme_allocate, &
       register_fluid_scheme
  use fluid_scheme_incompressible, only : fluid_scheme_incompressible_t
  use json_module, only : json_file
  use mesh, only : mesh_t
  use registry, only : neko_registry
  use time_state, only : time_state_t
  use time_step_controller, only : time_step_controller_t
  use user_intf, only : user_t
  use utils, only : neko_error
  implicit none
  private

  type, extends(fluid_scheme_incompressible_t) :: fluid_scheme_template_t
   contains
     procedure :: init => fluid_scheme_template_init
     procedure :: free => fluid_scheme_template_free
     procedure :: step => fluid_scheme_template_step
     procedure :: restart => fluid_scheme_template_restart
     procedure :: setup_bcs => fluid_scheme_template_setup_bcs
  end type fluid_scheme_template_t

  public :: fluid_scheme_template_register_types

contains

  subroutine fluid_scheme_template_init(this, msh, lx, params, &
       numerics_params, user, chkp, output_boundary)
    class(fluid_scheme_template_t), target, intent(inout) :: this
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(in) :: lx
    type(json_file), target, intent(inout) :: params
    type(json_file), target, intent(inout) :: numerics_params
    type(user_t), target, intent(in) :: user
    type(chkp_t), target, intent(inout) :: chkp
    logical, intent(in) :: output_boundary

    call this%free()

    ! `init_base` sets up the function space, the velocity fields and their
    ! lagged values, material properties, source terms and the velocity
    ! solver configured in `params`.
    call this%init_base(msh, lx, params, numerics_params, &
         'fluid_scheme_template', user, .true., output_boundary)

    call neko_registry%add_field(this%dm_Xh, 'p')
    this%p => neko_registry%get_field('p')

    ! TODO: Set up the pressure solver (this%ksp_prs), the boundary
    !       conditions through setup_bcs, and the time scheme in
    !       this%ext_bdf, which the scalars and LES models also use.
    !       Register the scheme state with chkp for restarts.
  end subroutine fluid_scheme_template_init

  subroutine fluid_scheme_template_free(this)
    class(fluid_scheme_template_t), intent(inout) :: this
    ! TODO: Release operators and state owned by this derived type.
    call this%scheme_free()
  end subroutine fluid_scheme_template_free

  subroutine fluid_scheme_template_step(this, time, dt_controller)
    class(fluid_scheme_template_t), target, intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_step_controller_t), intent(in) :: dt_controller

    ! TODO: Assemble and solve one fluid time step.
    call neko_error('fluid_scheme_template: step is not implemented')
  end subroutine fluid_scheme_template_step

  subroutine fluid_scheme_template_restart(this, chkp)
    class(fluid_scheme_template_t), target, intent(inout) :: this
    type(chkp_t), intent(inout) :: chkp
    ! TODO: Restore scheme-specific checkpoint state.
  end subroutine fluid_scheme_template_restart

  subroutine fluid_scheme_template_setup_bcs(this, user, params)
    class(fluid_scheme_template_t), target, intent(inout) :: this
    type(user_t), target, intent(in) :: user
    type(json_file), intent(inout) :: params
    ! TODO: Build this%bcs_vel and this%bcs_prs from the
    !       boundary_conditions entry in params.
  end subroutine fluid_scheme_template_setup_bcs

  subroutine fluid_scheme_template_register_types()
    procedure(fluid_scheme_allocate), pointer :: allocator
    allocator => fluid_scheme_template_allocate
    call register_fluid_scheme('fluid_scheme_template', allocator)
  end subroutine fluid_scheme_template_register_types

  subroutine fluid_scheme_template_allocate(obj)
    class(fluid_scheme_base_t), allocatable, intent(inout) :: obj
    allocate(fluid_scheme_template_t :: obj)
  end subroutine fluid_scheme_template_allocate
end module fluid_scheme_template
