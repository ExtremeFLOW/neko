!> Template for a user-defined simulation component.
module simulation_component_template
  use case, only : case_t
  use json_module, only : json_file
  use json_utils, only : json_get_or_default
  use simulation_component, only : simulation_component_t, &
       simulation_component_allocate, register_simulation_component
  use time_state, only : time_state_t
  implicit none
  private

  type, extends(simulation_component_t) :: simulation_component_template_t
   contains
     procedure :: init => simulation_component_template_init
     procedure :: free => simulation_component_template_free
     procedure :: compute_ => simulation_component_template_compute
  end type simulation_component_template_t

  public :: simulation_component_template_register_types

contains

  subroutine simulation_component_template_init(this, json, case)
    class(simulation_component_template_t), target, intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), target, intent(inout) :: case

    call this%init_base(json, case)
    call json_get_or_default(json, 'name', this%name, &
         'simulation_component_template')
    ! TODO: Read component-specific parameters from json.
  end subroutine simulation_component_template_init

  subroutine simulation_component_template_free(this)
    class(simulation_component_template_t), intent(inout) :: this
    ! TODO: Release fields owned by this derived type.
    call this%free_base()
  end subroutine simulation_component_template_free

  subroutine simulation_component_template_compute(this, time)
    class(simulation_component_template_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    ! TODO: Implement the component's scheduled action.
  end subroutine simulation_component_template_compute

  subroutine simulation_component_template_register_types()
    procedure(simulation_component_allocate), pointer :: allocator
    allocator => simulation_component_template_allocate
    call register_simulation_component('simulation_component_template', &
         allocator)
  end subroutine simulation_component_template_register_types

  subroutine simulation_component_template_allocate(obj)
    class(simulation_component_t), allocatable, intent(inout) :: obj
    allocate(simulation_component_template_t :: obj)
  end subroutine simulation_component_template_allocate
end module simulation_component_template
