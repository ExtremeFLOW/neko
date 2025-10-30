! Copyright (c) 2024-2025, The Neko Authors
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
!> Implements the `derivative_t` type.

module derivative_simcomp
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use time_state, only : time_state_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : dudxyz
  use case, only : case_t
  use json_utils, only : json_get, json_get_or_default
  use field_writer, only : field_writer_t
  use utils, only : neko_error
  use time_based_controller, only : time_based_controller_t
  implicit none
  private

  !> A simulation component that computes a derivative of a field.
  !! Wraps the `duxyz` operator.
  type, public, extends(simulation_component_t) :: derivative_t
     !> The scalar field to compute the weak gradient of.
     type(field_t), pointer :: u
     !> The derivative field
     type(field_t), pointer :: du
     !> Derivatives of r with respect to the direction of derivation.
     real(kind=rp), pointer, contiguous :: dr(:,:,:,:)
     !> Derivatives of s with respect to the direction of derivation.
     real(kind=rp), pointer, contiguous :: ds(:,:,:,:)
     !> Derivatives of t with respect to the direction of derivation.
     real(kind=rp), pointer, contiguous :: dt(:,:,:,:)
     !> Output writer.
     type(field_writer_t) :: writer

   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => derivative_init_from_json
     !> Generic for constructing from components.
     generic :: init_from_components => &
          init_from_controllers, init_from_controllers_properties
     !> Constructor from components, passing time_based_controllers.
     procedure, pass(this) :: init_from_controllers => &
          derivative_init_from_controllers
     !> Constructor from components, passing the properties of
     !! time_based_controllers.
     procedure, pass(this) :: init_from_controllers_properties => &
          derivative_init_from_controllers_properties
     !> Common part of both constructors.
     procedure, private, pass(this) :: init_common => derivative_init_common
     !> Destructor.
     procedure, pass(this) :: free => derivative_free
     !> Compute the derivative field.
     procedure, pass(this) :: compute_ => derivative_compute
  end type derivative_t

contains

  !> Constructor from json.
  subroutine derivative_init_from_json(this, json, case)
    class(derivative_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: field_name
    character(len=:), allocatable :: direction
    character(len=:), allocatable :: computed_field
    character(len=20) :: fields(1)

    ! Add fields keyword to the json so that the field_writer_t picks it up.
    ! Will also add fields to the registry.
    call json_get(json, "field", field_name)
    call json_get(json, "direction", direction)

    call json_get_or_default(json, "computed_field", computed_field, &
         "d" // trim(field_name) // "_d" // direction)

    fields(1) = trim(computed_field)
    call json%add("fields", fields)

    call this%init_base(json, case)
    call this%writer%init(json, case)

    call this%init_common(field_name, computed_field, direction)
  end subroutine derivative_init_from_json

  !> Common part of constructors from components.
  !! @param field_name The name of the field to compute the derivative of.
  !! @param computed_field The base name of the curl field components.
  !! @param direction The direction of the derivative, one of x, y or z.
  subroutine derivative_init_common(this, field_name, computed_field, &
       direction)
    class(derivative_t), intent(inout) :: this
    character(len=*) :: field_name
    character(len=*) :: computed_field
    character(len=*) :: direction

    this%u => neko_field_registry%get_field_by_name(trim(field_name))

    this%du => neko_field_registry%get_field_by_name(&
         "d" // field_name // "_d" // direction)

    if (direction .eq. "x") then
       this%dr => this%case%fluid%c_Xh%drdx
       this%ds => this%case%fluid%c_Xh%dsdx
       this%dt => this%case%fluid%c_Xh%dtdx
    else if (direction .eq. "y") then
       this%dr => this%case%fluid%c_Xh%drdy
       this%ds => this%case%fluid%c_Xh%dsdy
       this%dt => this%case%fluid%c_Xh%dtdy
    else if (direction .eq. "z") then
       this%dr => this%case%fluid%c_Xh%drdz
       this%ds => this%case%fluid%c_Xh%dsdz
       this%dt => this%case%fluid%c_Xh%dtdz
    else
       call neko_error("The direction of the derivative must be x, y or z")
    end if
  end subroutine derivative_init_common

  !> Constructor from components, passing controllers.
  !! @param case The simulation case object.
  !! @param order The execution oder priority of the simcomp.
  !! @param preprocess_controller The controller for running preprocessing.
  !! @param compute_controller The controller for running compute.
  !! @param output_controller The controller for producing output.
  !! @param field_name The name of the field to compute the derivative of.
  !! @param computed_field The base name of the curl field components.
  !! @param direction The direction of the derivative, one of x, y or z.
  !! @param filename The name of the file save the fields to. Optional, if not
  !! @param precision The real precision of the output data. Optional, defaults
  !! to single precision.
  subroutine derivative_init_from_controllers(this, case, order, &
       preprocess_controller, compute_controller, output_controller, &
       field_name, computed_field, direction, filename, precision)
    class(derivative_t), intent(inout) :: this
    class(case_t), intent(inout), target :: case
    integer :: order
    type(time_based_controller_t), intent(in) :: preprocess_controller
    type(time_based_controller_t), intent(in) :: compute_controller
    type(time_based_controller_t), intent(in) :: output_controller
    character(len=*) :: field_name
    character(len=*) :: computed_field
    character(len=*) :: direction
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision

    character(len=20) :: fields(1)

    fields(1) = trim(computed_field)

    call this%init_base_from_components(case, order, preprocess_controller, &
         compute_controller, output_controller)
    call this%writer%init_from_components(case, order, preprocess_controller, &
         compute_controller, output_controller, fields, filename, precision)
    call this%init_common(field_name, computed_field, direction)

  end subroutine derivative_init_from_controllers

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
  !! @param field_name The name of the field to compute the derivative of.
  !! @param computed_field The base name of the curl field components.
  !! @param direction The direction of the derivative, one of x, y or z.
  !! @param filename The name of the file save the fields to. Optional, if not
  !! provided, fields are added to the main output file.
  !! @param precision The real precision of the output data. Optional, defaults
  !! to single precision.
  subroutine derivative_init_from_controllers_properties(this, &
       case, order, preprocess_control, preprocess_value, compute_control, &
       compute_value, output_control, output_value, field_name, &
       computed_field, direction, filename, precision)
    class(derivative_t), intent(inout) :: this
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
    character(len=*) :: direction
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision

    character(len=20) :: fields(1)

    fields(1) = trim(computed_field)

    call this%init_base_from_components(case, order, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, &
         output_value)
    call this%writer%init_from_components(case, order, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, &
         output_value, fields, filename, precision)
    call this%init_common(field_name, computed_field, direction)

  end subroutine derivative_init_from_controllers_properties


  !> Destructor.
  subroutine derivative_free(this)
    class(derivative_t), intent(inout) :: this
    call this%free_base()
    call this%writer%free()
    nullify(this%du)
    nullify(this%u)
    nullify(this%dr)
    nullify(this%ds)
    nullify(this%dt)
  end subroutine derivative_free

  !> Compute the derivative field.
  !! @param time The current time.
  subroutine derivative_compute(this, time)
    class(derivative_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    call dudxyz(this%du%x, this%u%x, this%dr, this%ds, this%dt,&
         this%case%fluid%c_Xh)
  end subroutine derivative_compute

end module derivative_simcomp
