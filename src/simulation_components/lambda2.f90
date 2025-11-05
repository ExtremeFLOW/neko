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
!> A simulation component that computes lambda2
!! The values are stored in the field registry under the name 'lambda2'

module lambda2
  use num_types, only : rp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use time_state, only : time_state_t
  use operators, only : lambda2op
  use case, only : case_t
  use field_writer, only : field_writer_t
  use time_based_controller, only : time_based_controller_t
  use device
  implicit none
  private

  type, public, extends(simulation_component_t) :: lambda2_t
     !> X velocity component.
     type(field_t), pointer :: u
     !> Y velocity component.
     type(field_t), pointer :: v
     !> Z velocity component.
     type(field_t), pointer :: w

     !> X lambda2 component.
     type(field_t), pointer :: lambda2

     !> Work arrays.
     type(field_t) :: temp1
     type(field_t) :: temp2

     !> Output writer.
     type(field_writer_t) :: writer

   contains
     !> Constructor from json.
     procedure, pass(this) :: init => lambda2_init_from_json
     !> Generic for constructing from components.
     generic :: init_from_components => &
          init_from_controllers, init_from_controllers_properties
     !> Constructor from components, passing time_based_controllers.
     procedure, pass(this) :: init_from_controllers => &
          lambda2_init_from_controllers
     !> Constructor from components, passing the properties of
     !! time_based_controllers.
     procedure, pass(this) :: init_from_controllers_properties => &
          lambda2_init_from_controllers_properties
     !> Common part of both constructors.
     procedure, private, pass(this) :: init_common => lambda2_init_common
     !> Destructor.
     procedure, pass(this) :: free => lambda2_free
     !> Compute the lambda2 field
     procedure, pass(this) :: compute_ => lambda2_compute
  end type lambda2_t

contains

  !> Constructor from json.
  subroutine lambda2_init_from_json(this, json, case)
    class(lambda2_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target ::case
    character(len=20) :: fields(1)
    type(field_t), pointer :: u, v, w, lambda2

    ! Add fields keyword to the json so that the field_writer picks it up.
    ! Will also add fields to the registry.
    fields(1) = "lambda2"
    call json%add("fields", fields)

    call this%init_base(json, case)
    call this%writer%init(json, case)

    call this%init_common()
  end subroutine lambda2_init_from_json

  !> Common part of constructors.
  subroutine lambda2_init_common(this)
    class(lambda2_t), intent(inout) :: this

    this%u => neko_field_registry%get_field("u")
    this%v => neko_field_registry%get_field("v")
    this%w => neko_field_registry%get_field("w")
    this%lambda2 => neko_field_registry%get_field("lambda2")

  end subroutine lambda2_init_common

  !> Constructor from components, passing controllers.
  !! @param case The simulation case object.
  !! @param order The execution oder priority of the simcomp.
  !! @param preprocess_controller The controller for running preprocessing.
  !! @param compute_controller The controller for running compute.
  !! @param output_controller The controller for producing output.
  !! @param filename The name of the file save the fields to. Optional, if not
  !! @param precision The real precision of the output data. Optional, defaults
  !! to single precision.
  subroutine lambda2_init_from_controllers(this, case, order, &
       preprocess_controller, compute_controller, output_controller, &
       filename, precision)
    class(lambda2_t), intent(inout) :: this
    class(case_t), intent(inout), target :: case
    integer :: order
    type(time_based_controller_t), intent(in) :: preprocess_controller
    type(time_based_controller_t), intent(in) :: compute_controller
    type(time_based_controller_t), intent(in) :: output_controller
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision

    character(len=20) :: fields(1)
    fields(1) = "lambda2"

    call this%init_base_from_components(case, order, preprocess_controller, &
         compute_controller, output_controller)
    call this%writer%init_from_components(case, order, preprocess_controller, &
         compute_controller, output_controller, fields, filename, precision)
    call this%init_common()

  end subroutine lambda2_init_from_controllers

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
  !! @param filename The name of the file save the fields to. Optional, if not
  !! provided, fields are added to the main output file.
  !! @param precision The real precision of the output data. Optional, defaults
  !! to single precision.
  subroutine lambda2_init_from_controllers_properties(this, &
       case, order, preprocess_control, preprocess_value, compute_control, &
       compute_value, output_control, output_value, filename, precision)
    class(lambda2_t), intent(inout) :: this
    class(case_t), intent(inout), target :: case
    integer :: order
    character(len=*), intent(in) :: preprocess_control
    real(kind=rp), intent(in) :: preprocess_value
    character(len=*), intent(in) :: compute_control
    real(kind=rp), intent(in) :: compute_value
    character(len=*), intent(in) :: output_control
    real(kind=rp), intent(in) :: output_value
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision

    character(len=20) :: fields(1)
    fields(1) = "lambda2"

    call this%init_base_from_components(case, order, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, &
         output_value)
    call this%writer%init_from_components(case, order, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, &
         output_value, fields, filename, precision)
    call this%init_common()

  end subroutine lambda2_init_from_controllers_properties

  !> Destructor.
  subroutine lambda2_free(this)
    class(lambda2_t), intent(inout) :: this
    call this%free_base()

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%lambda2)
  end subroutine lambda2_free

  !> Compute the lambda2 field.
  !! @param time The time state.
  subroutine lambda2_compute(this, time)
    class(lambda2_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    call lambda2op(this%lambda2, this%u, this%v, this%w, this%case%fluid%c_Xh)

  end subroutine lambda2_compute

end module lambda2
