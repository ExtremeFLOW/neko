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
!> Implements the `divergence_t` type.

module divergence_simcomp
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use time_state, only : time_state_t
  use operators, only : div
  use case, only : case_t
  use json_utils, only : json_get, json_get_or_default
  use field_writer, only : field_writer_t
  use time_based_controller, only : time_based_controller_t
  use utils, only : neko_error
  implicit none
  private

  !> A simulation component that computes the divergence of a vector field.
  !! Added to the field registry as `div` by default, bu can be controlled by
  !! the computed_field keyword.
  type, public, extends(simulation_component_t) :: divergence_t
     !> X input vector field component.
     type(field_t), pointer :: u
     !> Y input vector field component.
     type(field_t), pointer :: v
     !> Z input vector field component.
     type(field_t), pointer :: w

     !> The divergenceergence field.
     type(field_t), pointer :: divergence

     !> Output writer.
     type(field_writer_t) :: writer

   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => divergence_init_from_json
     !> Generic for constructing from components.
     generic :: init_from_components => &
          init_from_controllers, init_from_controllers_properties
     !> Constructor from components, passing time_based_controllers.
     procedure, pass(this) :: init_from_controllers => &
          divergence_init_from_controllers
     !> Constructor from components, passing the properties of
     !! time_based_controllers.
     procedure, pass(this) :: init_from_controllers_properties => &
          divergence_init_from_controllers_properties
     !> Common part of both constructors.
     procedure, private, pass(this) :: init_common => divergence_init_common
     !> Destructor.
     procedure, pass(this) :: free => divergence_free
     !> Compute the divergence field.
     procedure, pass(this) :: compute_ => divergence_compute
  end type divergence_t

contains

  !> Constructor from json.
  subroutine divergence_init_from_json(this, json, case)
    class(divergence_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=20) :: fields(1)
    character(len=20), allocatable :: field_names(:)
    character(len=:), allocatable :: computed_field


    call json_get_or_default(json, "computed_field", computed_field, "divergence")
    call json_get(json, "fields", field_names)

    if (size(field_names) .ne. 3) then
       call neko_error("The divergence simcomp requires exactly 3 entries in " // &
            "fieldes.")
    end if

    fields(1) = trim(computed_field)

    ! This is needed for the field writer to pick up the fields.
    call json%add("fields", fields)

    call this%init_base(json, case)
    call this%writer%init(json, case)

    call divergence_init_common(this, field_names, computed_field)
  end subroutine divergence_init_from_json

  !> Actual constructor.
  !! @param field_names The name of the fields to compute the divergence of.
  !! @param computed_field The base name of the divergence field components.
  subroutine divergence_init_common(this, field_names, computed_field)
    class(divergence_t), intent(inout) :: this
    character(len=*) :: field_names(3)
    character(len=*) :: computed_field

    this%u => neko_field_registry%get_field_by_name(field_names(1))
    this%v => neko_field_registry%get_field_by_name(field_names(2))
    this%w => neko_field_registry%get_field_by_name(field_names(3))

    this%divergence => neko_field_registry%get_field_by_name(computed_field)

  end subroutine divergence_init_common

  !> Constructor from components, passing controllers.
  !! @param case The simulation case object.
  !! @param order The execution oder priority of the simcomp.
  !! @param preprocess_controller The controller for running preprocessing.
  !! @param compute_controller The controller for running compute.
  !! @param output_controller The controller for producing output.
  !! @param field_names The name of the fields to compute the divergence of.
  !! @param computed_field The base name of the divergence field components.
  !! @param filename The name of the file save the fields to. Optional, if not
  !! @param precision The real precision of the output data. Optional, defaults
  !! to single precision.
  subroutine divergence_init_from_controllers(this, case, order, &
       preprocess_controller, compute_controller, output_controller, &
       field_names, computed_field, filename, precision)
    class(divergence_t), intent(inout) :: this
    class(case_t), intent(inout), target :: case
    integer :: order
    type(time_based_controller_t), intent(in) :: preprocess_controller
    type(time_based_controller_t), intent(in) :: compute_controller
    type(time_based_controller_t), intent(in) :: output_controller
    character(len=*) :: field_names(3)
    character(len=*) :: computed_field
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision

    character(len=20) :: fields(1)

    fields(1) = trim(computed_field)

    call this%init_base_from_components(case, order, preprocess_controller, &
         compute_controller, output_controller)
    call this%writer%init_from_components(case, order, preprocess_controller, &
         compute_controller, output_controller, fields, filename, precision)
    call this%init_common(field_names, computed_field)

  end subroutine divergence_init_from_controllers

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
  !! @param field_names The name of the field to compute the divergence of.
  !! @param computed_field The base name of the divergence field components.
  !! @param filename The name of the file save the fields to. Optional, if not
  !! provided, fields are added to the main output file.
  !! @param precision The real precision of the output data. Optional, defaults
  !! to single precision.
  subroutine divergence_init_from_controllers_properties(this, &
       case, order, preprocess_control, preprocess_value, compute_control, &
       compute_value, output_control, output_value, field_names, &
       computed_field, filename, precision)
    class(divergence_t), intent(inout) :: this
    class(case_t), intent(inout), target :: case
    integer :: order
    character(len=*), intent(in) :: preprocess_control
    real(kind=rp), intent(in) :: preprocess_value
    character(len=*), intent(in) :: compute_control
    real(kind=rp), intent(in) :: compute_value
    character(len=*), intent(in) :: output_control
    real(kind=rp), intent(in) :: output_value
    character(len=*) :: field_names(3)
    character(len=*) :: computed_field
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision

    character(len=20) :: fields(1)

    fields(1) = trim(computed_field) // "_x"

    call this%init_base_from_components(case, order, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, &
         output_value)
    call this%writer%init_from_components(case, order, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, &
         output_value, fields, filename, precision)
    call this%init_common(field_names, computed_field)

  end subroutine divergence_init_from_controllers_properties

  !> Destructor.
  subroutine divergence_free(this)
    class(divergence_t), intent(inout) :: this
    call this%free_base()
    call this%writer%free()

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%divergence)
  end subroutine divergence_free

  !> Compute the divergence field.
  subroutine divergence_compute(this, time)
    class(divergence_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    call div(this%divergence%x, this%u%x, this%v%x, this%w%x, this%case%fluid%c_Xh)

  end subroutine divergence_compute

end module divergence_simcomp
