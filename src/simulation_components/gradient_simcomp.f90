! Copyright (c) 2025, The Neko Authors
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
!> Implements the `gradient_t` type.
module gradient_simcomp
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : grad
  use time_state, only : time_state_t
  use case, only : case_t
  use json_utils, only : json_get, json_get_or_default
  use field_writer, only : field_writer_t
  use time_based_controller, only : time_based_controller_t
  implicit none
  private

  !> A simulation component that computes the gradient of a field.
  !! Wraps the `gradient` operator.
  type, public, extends(simulation_component_t) :: gradient_t
     !> The scalar field to compute the weak gradient of.
     type(field_t), pointer :: u
     !> X weak gradient component.
     type(field_t), pointer :: gradient_x
     !> Y weak gradient component.
     type(field_t), pointer :: gradient_y
     !> Z weak gradient component.
     type(field_t), pointer :: gradient_z
     !> Output writer.
     type(field_writer_t) :: writer

   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => gradient_init_from_json
     !> Generic for constructing from components.
     generic :: init_from_components => &
          init_from_controllers, init_from_controllers_properties
     !> Constructor from components, passing time_based_controllers.
     procedure, pass(this) :: init_from_controllers => &
          gradient_init_from_controllers
     !> Constructor from components, passing the properties of
     !! time_based_controllers.
     procedure, pass(this) :: init_from_controllers_properties => &
          gradient_init_from_controllers_properties
     !> Common part of both constructors.
     procedure, private, pass(this) :: init_common => gradient_init_common
     !> Destructor.
     procedure, pass(this) :: free => gradient_free
     !> Compute the gradient field.
     procedure, pass(this) :: compute_ => gradient_compute
  end type gradient_t

contains

  !> Constructor from json.
  subroutine gradient_init_from_json(this, json, case)
    class(gradient_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: field_name
    character(len=20) :: fields(3)
    character(len=:), allocatable :: computed_field

    call json_get(json, "field", field_name)

    ! Add fields keyword to the json so that the field_writer picks it up.
    ! Will also add fields to the registry.
    call json_get_or_default(json, "computed_field", computed_field, &
         "gradient" // trim(field_name))

    fields(1) = computed_field // "_x"
    fields(2) = computed_field // "_y"
    fields(3) = computed_field // "_z"
    write(*,*) fields(1), fields(2), fields(3)

    call json%add("fields", fields)

    call this%init_base(json, case)
    call this%writer%init(json, case)

    call this%init_common(field_name, computed_field)
  end subroutine gradient_init_from_json

  !> Common part of the constructors.
  subroutine gradient_init_common(this, field_name, computed_field)
    class(gradient_t), intent(inout) :: this
    character(len=*) :: field_name
    character(len=*) :: computed_field

    this%u => neko_field_registry%get_field_by_name(trim(field_name))

    this%gradient_x => neko_field_registry%get_field_by_name( &
         computed_field // "_x")
    this%gradient_y => neko_field_registry%get_field_by_name( &
         computed_field // "_y")
    this%gradient_z => neko_field_registry%get_field_by_name( &
         computed_field // "_z")


  end subroutine gradient_init_common

  !> Constructor from components, passing controllers.
  !! @param case The simulation case object.
  !! @param order The execution oder priority of the simcomp.
  !! @param preprocess_controller The controller for running preprocessing.
  !! @param compute_controller The controller for running compute.
  !! @param output_controller The controller for producing output.
  !! @param field_name The name of the field to compute the gradient of.
  !! @param computed_field The base name of the gradient field components.
  !! @param filename The name of the file save the fields to. Optional, if not
  !! @param precision The real precision of the output data. Optional, defaults
  !! to single precision.
  subroutine gradient_init_from_controllers(this, case, order, &
       preprocess_controller, compute_controller, output_controller, &
       field_name, computed_field, filename, precision)
    class(gradient_t), intent(inout) :: this
    class(case_t), intent(inout), target :: case
    integer :: order
    type(time_based_controller_t), intent(in) :: preprocess_controller
    type(time_based_controller_t), intent(in) :: compute_controller
    type(time_based_controller_t), intent(in) :: output_controller
    character(len=*) :: field_name
    character(len=*) :: computed_field
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision

    character(len=20) :: fields(3)

    fields(1) = trim(computed_field) // "_x"
    fields(2) = trim(computed_field) // "_y"
    fields(3) = trim(computed_field) // "_z"

    call this%init_base_from_components(case, order, preprocess_controller, &
         compute_controller, output_controller)
    call this%writer%init_from_components(case, order, preprocess_controller, &
         compute_controller, output_controller, fields, filename, precision)
    call this%init_common(field_name, computed_field)

  end subroutine gradient_init_from_controllers

  !> Constructor from components, passing properties to the
  !! time_based_controller` components in the base type.
  !! @param case The simulation case object.
  !! @param order The execution oder priority of the simcomp.
  !! @param preprocess_controller Control mode for preprocessing.
  !! @param preprocess_value Value parameter for preprocessing.
  !! @param compute_controller Control mode for computing.
  !! @param compute_value Value parameter for computing.
  !! @param output_controller Control mode for output.
  !! @param output_value Value parameter for output.
  !! @param field_name The name of the field to compute the gradient of.
  !! @param computed_field The base name of the gradient field components.
  !! @param filename The name of the file save the fields to. Optional, if not
  !! provided, fields are added to the main output file.
  !! @param precision The real precision of the output data. Optional, defaults
  !! to single precision.
  subroutine gradient_init_from_controllers_properties(this, &
       case, order, preprocess_control, preprocess_value, compute_control, &
       compute_value, output_control, output_value, field_name, &
       computed_field, filename, precision)
    class(gradient_t), intent(inout) :: this
    class(case_t), intent(inout), target :: case
    integer :: order
    character(len=*), intent(in) :: preprocess_control
    real(kind=rp), intent(in) :: preprocess_value
    character(len=*), intent(in) :: compute_control
    real(kind=rp), intent(in) :: compute_value
    character(len=*), intent(in) :: output_control
    real(kind=rp), intent(in) :: output_value
    character(len=*) :: field_name
    character(len=*) :: computed_field
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision

    character(len=20) :: fields(3)

    fields(1) = trim(computed_field) // "_x"
    fields(2) = trim(computed_field) // "_y"
    fields(3) = trim(computed_field) // "_z"

    call this%init_base_from_components(case, order, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, &
         output_value)
    call this%writer%init_from_components(case, order, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, &
         output_value, fields, filename, precision)
    call this%init_common(field_name, computed_field)

  end subroutine gradient_init_from_controllers_properties

  !> Destructor.
  subroutine gradient_free(this)
    class(gradient_t), intent(inout) :: this
    call this%free_base()
    call this%writer%free()
    nullify(this%gradient_x)
    nullify(this%gradient_y)
    nullify(this%gradient_z)
    nullify(this%u)
  end subroutine gradient_free

  !> Compute the gradient field.
  !! @param time The current time value
  subroutine gradient_compute(this, time)
    class(gradient_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    call grad(this%gradient_x%x, this%gradient_y%x, this%gradient_z%x, this%u%x,&
         this%case%fluid%c_Xh)
  end subroutine gradient_compute

end module gradient_simcomp
