! Copyright (c) 2023, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
!
!> Simulation components are objects that encapsulate functionality that can be
!! fit to a particular compute pattern.
!! @note
!! The canonical way to abbreviate simulation_component is simcomp.
module simulation_component
  use num_types, only : rp
  use json_module, only : json_file
  use case, only : case_t
  use time_based_controller, only : time_based_controller_t
  use json_utils, only : json_get_or_default, json_get
  use time_state, only : time_state_t
  implicit none
  private

  !> Base abstract class for simulation components.
  type, abstract, public :: simulation_component_t
     !> Pointer to the simulation case.
     type(case_t), pointer :: case
     !> Controller for when to run `preprocess`.
     type(time_based_controller_t) :: preprocess_controller
     !> Controller for when to run `compute`.
     type(time_based_controller_t) :: compute_controller
     !> Controller for when to do output.
     type(time_based_controller_t) :: output_controller
     !> The execution order, lowest excutes first.
     integer :: order
   contains
     !> Constructor for the simulation_component_t (base) class.
     procedure, pass(this) :: init_base => simulation_component_init_base
     !> Constructor for the simulation_component_t (base) class from components.
     generic :: init_base_from_components => &
          init_base_from_controllers_properties, &
          init_base_from_controllers
     !> Constructor for the simulation_component_t (base) class from
     !! time_based_controllers, essentially directly from all components (we
     !! reserve the `_from_components` name for the generic binding).
     procedure, pass(this) :: init_base_from_controllers => &
          simulation_component_init_base_from_controllers
     !> Constructor for the simulation_component_t (base) class from
     !! properties of time_based_controllers, so the latter are
     !! constructed instead of assigned.
     procedure, pass(this) :: init_base_from_controllers_properties => &
          simulation_component_init_base_from_controllers_properties
     !> Destructor for the simulation_component_t (base) class.
     procedure, pass(this) :: free_base => simulation_component_free_base
     !> Wrapper for calling `set_counter` for the time based controllers.
     !! Serves as the public interface.
     procedure, pass(this) :: restart => simulation_component_restart_wrapper
     !> Wrapper for calling `preprocess_` based on the `preprocess_controller`.
     !! Serves as the public interface.
     procedure, pass(this) :: preprocess => &
          simulation_component_preprocess_wrapper
     !> Wrapper for calling `compute_` based on the `compute_controller`.
     !! Serves as the public interface.
     procedure, pass(this) :: compute => simulation_component_compute_wrapper
     !> The common constructor using a JSON dictionary.
     procedure(simulation_component_init), pass(this), deferred :: init
     !> Destructor.
     procedure(simulation_component_free), pass(this), deferred :: free
     !> The preprocessing function to be executed during the run.
     procedure, pass(this) :: preprocess_
     !> The main function to be executed during the run.
     procedure, pass(this) :: compute_
     !> The restart function to be called upon restarting simulation
     procedure, pass(this) :: restart_
     !> JSON parameter parser for the time-based controllers
     procedure, pass(this) :: parse_json => simulation_component_parse_json
  end type simulation_component_t

  !> A helper type that is needed to have an array of polymorphic objects
  type, public :: simulation_component_wrapper_t
     class(simulation_component_t), allocatable :: simcomp
  end type simulation_component_wrapper_t


  abstract interface
     !> The common constructor using a JSON dictionary.
     !! @param json The JSON with properties.
     !! @param case The case_t object.
     subroutine simulation_component_init(this, json, case)
       import simulation_component_t, json_file, case_t
       class(simulation_component_t), intent(inout), target :: this
       type(json_file), intent(inout) :: json
       class(case_t), intent(inout), target :: case
     end subroutine simulation_component_init
  end interface

  abstract interface
     !> Destructor.
     subroutine simulation_component_free(this)
       import simulation_component_t
       class(simulation_component_t), intent(inout) :: this
     end subroutine simulation_component_free
  end interface

  interface
     !> Simulation component factory.
     !! Both constructs and initializes the object.
     !! @param object The object to be created and initialized.
     !! @param json JSON object initializing the simulation component.
     !! @param case The simulation case.
     module subroutine simulation_component_factory(object, json, case)
       class(simulation_component_t), allocatable, intent(inout) :: object
       type(json_file), intent(inout) :: json
       class(case_t), intent(inout), target :: case
     end subroutine simulation_component_factory
  end interface

  interface
     !> Simulation component allocator.
     !! @param object The object to be allocated.
     !! @param type_name The name of the simcomp type.
     module subroutine simulation_component_allocator(object, type_name)
       class(simulation_component_t), allocatable, intent(inout) :: object
       character(len=*), intent(in):: type_name
     end subroutine simulation_component_allocator
  end interface

  !
  ! Machinery for injecting user-defined types
  !

  !> Interface for an object allocator.
  !! Implemented in the user modules, should allocate the `obj` to the custom
  !! user type.
  abstract interface
     subroutine simulation_component_allocate(obj)
       import simulation_component_t
       class(simulation_component_t), allocatable, intent(inout) :: obj
     end subroutine simulation_component_allocate
  end interface

  interface
     !> Called in user modules to add an allocator for custom types.
     module subroutine register_simulation_component(type_name, allocator)
       character(len=*), intent(in) :: type_name
       procedure(simulation_component_allocate), pointer, intent(in) :: &
            allocator
     end subroutine register_simulation_component
  end interface

  ! A name-allocator pair for user-defined types. A helper type to define a
  ! registry of custom allocators.
  type allocator_entry
     character(len=20) :: type_name
     procedure(simulation_component_allocate), pointer, nopass :: allocator
  end type allocator_entry

  !> Registry of allocators for user-defined types
  type(allocator_entry), allocatable :: simcomp_registry(:)

  !> The size of the `simulation_component_registry`
  integer :: simcomp_registry_size = 0

  public :: simulation_component_factory, simulation_component_allocator, &
       register_simulation_component, simulation_component_allocate


contains
  !> Constructor for the `simulation_component_t` (base) class.
  subroutine simulation_component_init_base(this, json, case)
    class(simulation_component_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: preprocess_control, compute_control, &
         output_control
    real(kind=rp) :: preprocess_value, compute_value, output_value
    integer :: order

    call this%parse_json(json, case%params, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, &
         output_value)

    call json_get_or_default(json, "order", order, -1)

    call this%init_base_from_components(case, order, &
         preprocess_control, preprocess_value, compute_control, compute_value, &
         output_control, output_value)

  end subroutine simulation_component_init_base

  !> Constructor for the `simulation_component_t` (base) class via the
  !! properties of the `time_based_controller` components.
  !! @param case The simulation case object.
  !! @param order The execution oder priority of the simcomp.
  !! @param preprocess_controller Control mode for preprocessing.
  !! @param preprocess_value Value parameter for preprocessing.
  !! @param compute_controller Control mode for computing.
  !! @param compute_value Value parameter for computing.
  !! @param output_controller Control mode for output.
  !! @param output_value Value parameter for output.
  subroutine simulation_component_init_base_from_controllers_properties(this, &
       case, order, preprocess_control, preprocess_value, compute_control, &
       compute_value, output_control, output_value)
    class(simulation_component_t), intent(inout) :: this
    class(case_t), intent(inout), target :: case
    integer :: order
    character(len=*), intent(in) :: preprocess_control
    real(kind=rp), intent(in) :: preprocess_value
    character(len=*), intent(in) :: compute_control
    real(kind=rp), intent(in) :: compute_value
    character(len=*), intent(in) :: output_control
    real(kind=rp), intent(in) :: output_value

    this%case => case
    this%order = order

    call this%preprocess_controller%init(case%time%start_time, &
         case%time%end_time, preprocess_control, preprocess_value)
    call this%compute_controller%init(case%time%start_time, case%time%end_time,&
         compute_control, compute_value)
    call this%output_controller%init(case%time%start_time, case%time%end_time, &
         output_control, output_value)

  end subroutine simulation_component_init_base_from_controllers_properties

  !> Constructor for the `simulation_component_t` (base) class from components.
  !! @param case The simulation case object.
  !! @param order The execution oder priority of the simcomp.
  !! @param preprocess_controller The controller for running preprocessing.
  !! @param compute_controller The controller for running compute.
  !! @param output_controller The controller for producing output.
  subroutine simulation_component_init_base_from_controllers(this, case, order,&
       preprocess_controller, compute_controller, output_controller)
    class(simulation_component_t), intent(inout) :: this
    class(case_t), intent(inout), target :: case
    integer :: order
    type(time_based_controller_t), intent(in) :: preprocess_controller
    type(time_based_controller_t), intent(in) :: compute_controller
    type(time_based_controller_t), intent(in) :: output_controller

    this%case => case
    this%order = order
    this%preprocess_controller = preprocess_controller
    this%compute_controller = compute_controller
    this%output_controller = output_controller
  end subroutine simulation_component_init_base_from_controllers

  !> Parse JSON to determine the properties of the `time_based_controllers`.
  !! @param json The JSON dictionary of the simcomp.
  !! @param case_params The entire case configuration JSON.
  !! @param preprocess_value Control mode for preprocessing.
  !! @param preprocess_controller Value parameter for preprocessing.
  !! @param compute_controller Control mode for computing.
  !! @param compute_controller Value parameter for computing.
  !! @param output_controller Control mode for output.
  !! @param output_controller Value parameter for output.
  subroutine simulation_component_parse_json(this, json, case_params, &
       preprocess_control, preprocess_value, compute_control, compute_value, &
       output_control, output_value)
    class(simulation_component_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(json_file), intent(inout) :: case_params
    character(len=:), allocatable, intent(inout) :: preprocess_control
    real(kind=rp), intent(out) :: preprocess_value
    character(len=:), allocatable, intent(inout) :: compute_control
    real(kind=rp), intent(out) :: compute_value
    character(len=:), allocatable, intent(inout) :: output_control
    real(kind=rp), intent(out) :: output_value

    ! We default to preprocess every time-step
    call json_get_or_default(json, "preprocess_control", preprocess_control, &
         "tsteps")
    call json_get_or_default(json, "preprocess_value", preprocess_value, 1.0_rp)

    ! We default to compute every time-step
    call json_get_or_default(json, "compute_control", compute_control, &
         "tsteps")
    call json_get_or_default(json, "compute_value", compute_value, 1.0_rp)

    if (compute_control .eq. "fluid_output") then
       call json_get(case_params, 'case.fluid.output_control', &
            compute_control)
       call json_get(case_params, 'case.fluid.output_value', &
            compute_value)
    end if

    ! We default to output whenever we execute
    call json_get_or_default(json, "output_control", output_control, &
         compute_control)
    call json_get_or_default(json, "output_value", output_value, &
         compute_value)

    if (output_control == "global") then
       call json_get(case_params, 'case.fluid.output_control', &
            output_control)
       call json_get(case_params, 'case.fluid.output_value', &
            output_value)
    end if
  end subroutine simulation_component_parse_json

  !> Destructor for the `simulation_component_t` (base) class.
  subroutine simulation_component_free_base(this)
    class(simulation_component_t), intent(inout) :: this

    nullify(this%case)

    call this%preprocess_controller%free()
    call this%compute_controller%free()
    call this%output_controller%free()
  end subroutine simulation_component_free_base

  !> Wrapper for calling `preprocess_` based on the `preprocess_controller`.
  !! Serves as the public interface.
  !! @param time The current time.
  subroutine simulation_component_preprocess_wrapper(this, time)
    class(simulation_component_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    if (this%preprocess_controller%check(time)) then
       call this%preprocess_(time)
       call this%preprocess_controller%register_execution()
    end if
  end subroutine simulation_component_preprocess_wrapper

  !> Wrapper for calling `compute_` based on the `compute_controller`.
  !! Serves as the public interface.
  !! @param time The current time.
  subroutine simulation_component_compute_wrapper(this, time)
    class(simulation_component_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    if (this%compute_controller%check(time)) then
       call this%compute_(time)
       call this%compute_controller%register_execution()
    end if
  end subroutine simulation_component_compute_wrapper

  !> Wrapper for calling `set_counter_` based for the controllers.
  !! @param time The current time.
  subroutine simulation_component_restart_wrapper(this, time)
    class(simulation_component_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    call this%compute_controller%set_counter(time)
    call this%output_controller%set_counter(time)
    call this%restart_(time)

  end subroutine simulation_component_restart_wrapper

  !> Dummy restart function.
  !! @param time The current time.
  subroutine restart_(this, time)
    class(simulation_component_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    ! Do nothing
  end subroutine restart_

  !> Dummy preprocessing function.
  !! @param time The current time.
  subroutine preprocess_(this, time)
    class(simulation_component_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    ! Do nothing
  end subroutine preprocess_

  !> Dummy compute function.
  !! @param time The current time.
  subroutine compute_(this, time)
    class(simulation_component_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    ! Do nothing
  end subroutine compute_
end module simulation_component
