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
!> Implements the `user_stats_t` type.

module user_stats
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use case, only : case_t
  use mean_field_output, only : mean_field_output_t
  use json_utils, only : json_get, json_get_or_default
  use mean_field, only : mean_field_t
  use coefs, only : coef_t
  use time_state, only : time_state_t
  use time_based_controller, only : time_based_controller_t
  implicit none
  private

  !> A simulation component that computes the averages of fields in the registry.
  type, public, extends(simulation_component_t) :: user_stats_t

     !> When to start averaging.
     real(kind=rp) :: start_time
     !> Current time. Uses to compute time delta since last run of compute.
     real(kind=rp) :: time
     !> The averaged fields.
     type(mean_field_t), allocatable :: mean_fields(:)
     !> Number of fields to average.
     integer :: n_avg_fields = 0
     !> The names of the fields to average.
     character(len=20), allocatable :: field_names(:)
     !> Output writer.
     type(mean_field_output_t), private :: output

   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => user_stats_init_from_json
     !> Generic for constructing from components.
     generic :: init_from_components => &
          init_from_controllers, init_from_controllers_properties
     !> Constructor from components, passing time_based_controllers.
     procedure, pass(this) :: init_from_controllers => &
          user_stats_init_from_controllers
     !> Constructor from components, passing the properties of
     !! time_based_controllers.
     procedure, pass(this) :: init_from_controllers_properties => &
          user_stats_init_from_controllers_properties
     !> Common part of both constructors.
     procedure, private, pass(this) :: init_common => user_stats_init_common
     !> Destructor.
     procedure, pass(this) :: free => user_stats_free
     !> Compute the means
     procedure, pass(this) :: compute_ => user_stats_compute
     procedure, pass(this) :: restart_ => user_stats_restart
  end type user_stats_t

contains

  !> Constructor from json.
  !! @param json The json paramter dictionary.
  !! @param case The neko case object.
  subroutine user_stats_init_from_json(this, json, case)
    class(user_stats_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: filename
    character(len=:), allocatable :: avg_dir

    call this%init_base(json, case)

    !> Get the number of stat fields and their names
    call json%info('fields', n_children = this%n_avg_fields)
    call json_get(json, 'fields', this%field_names)
    call json_get_or_default(json, 'start_time', this%start_time, 0.0_rp)
    call json_get_or_default(json, 'avg_direction', avg_dir, 'none')
    call json_get_or_default(json, 'output_file', filename, 'user_stats')

    call user_stats_init_common(this, this%start_time, &
         case%fluid%c_Xh, avg_dir, filename = filename)
  end subroutine user_stats_init_from_json

  subroutine user_stats_restart(this, time)
    class(user_stats_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    if (time%t .gt. this%time) this%time = time%t
  end subroutine user_stats_restart

  !> Constructor from components, passing controllers.
  !! @param case The simulation case object.
  !! @param order The execution oder priority of the simcomp.
  !! @param preprocess_controller The controller for running preprocessing.
  !! @param compute_controller The controller for running compute.
  !! @param output_controller The controller for producing output.
  !! @param start_time The start time for gathering samples for the average.
  !! @param coef The SEM coefficients.
  !! @param avg_dir The averaging direction.
  !! @param filename The name of the file save the fields to. Optional, if not
  !! @param precision The real precision of the output data. Optional, defaults
  !! to single precision.
  subroutine user_stats_init_from_controllers(this, case, order, &
       preprocess_controller, compute_controller, output_controller, &
       start_time, coef, avg_dir, filename, precision)
    class(user_stats_t), intent(inout) :: this
    class(case_t), intent(inout), target :: case
    integer :: order
    type(time_based_controller_t), intent(in) :: preprocess_controller
    type(time_based_controller_t), intent(in) :: compute_controller
    type(time_based_controller_t), intent(in) :: output_controller
    real(kind=rp), intent(in) :: start_time
    character(len=*), intent(in) :: avg_dir
    type(coef_t), intent(inout) :: coef
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision

    call this%init_base_from_components(case, order, preprocess_controller, &
         compute_controller, output_controller)
    call this%init_common(start_time, coef, avg_dir, filename, precision)

  end subroutine user_stats_init_from_controllers

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
  !! @param start_time The start time for gathering samples for the average.
  !! @param coef The SEM coefficients.
  !! @param avg_dir The averaging direction.
  !! @param filename The name of the file save the fields to. Optional, if not
  !! provided, fields are added to the main output file.
  !! @param precision The real precision of the output data. Optional, defaults
  !! to single precision.
  subroutine user_stats_init_from_controllers_properties(this, &
       case, order, preprocess_control, preprocess_value, compute_control, &
       compute_value, output_control, output_value, start_time, coef, avg_dir, &
       filename, precision)
    class(user_stats_t), intent(inout) :: this
    class(case_t), intent(inout), target :: case
    integer :: order
    character(len=*), intent(in) :: preprocess_control
    real(kind=rp), intent(in) :: preprocess_value
    character(len=*), intent(in) :: compute_control
    real(kind=rp), intent(in) :: compute_value
    character(len=*), intent(in) :: output_control
    real(kind=rp), intent(in) :: output_value
    real(kind=rp), intent(in) :: start_time
    character(len=*), intent(in) :: avg_dir
    type(coef_t), intent(inout) :: coef
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision

    call this%init_base_from_components(case, order, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, &
         output_value)
    call this%init_common(start_time, coef, avg_dir, filename, precision)

  end subroutine user_stats_init_from_controllers_properties


  !> Common part of constructors
  !! @param start_time The start time for gathering samples for the average.
  !! @param coef The SEM coefficients.
  !! @param avg_dir The averaging direction.
  subroutine user_stats_init_common(this, start_time, coef, avg_dir, &
         filename, precision)
    class(user_stats_t), intent(inout) :: this
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: precision
    real(kind=rp), intent(in) :: start_time
    character(len=*), intent(in) :: avg_dir
    type(coef_t), intent(inout) :: coef
    integer :: i
    type(field_t), pointer :: field_to_avg

    this%start_time = start_time
    this%time = start_time

    !> Allocate and initialize the mean fields
    allocate(this%mean_fields(this%n_avg_fields))
    do i = 1, this%n_avg_fields
       field_to_avg => neko_field_registry%get_field(trim(this%field_names(i)))
       call this%mean_fields(i)%init(field_to_avg)
    end do

    call this%output%init(this%mean_fields, this%n_avg_fields, &
         this%start_time, coef, avg_dir, name=filename)
    call this%case%output_controller%add(this%output, &
         this%output_controller%control_value, &
         this%output_controller%control_mode)
  end subroutine user_stats_init_common

  !> Destructor.
  subroutine user_stats_free(this)
    class(user_stats_t), intent(inout) :: this
    integer :: i

    call this%free_base()

    if (allocated(this%mean_fields)) then
       do i = 1, this%n_avg_fields
          call this%mean_fields(i)%free()
       end do
       deallocate(this%mean_fields)
    end if

    if (allocated(this%field_names)) then
       deallocate(this%field_names)
    end if

  end subroutine user_stats_free

  !> Update the running averages.
  !! @param time The current time state.
  subroutine user_stats_compute(this, time)
    class(user_stats_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer :: i

    !> Update the running average of the fields
    if (time%t .ge. this%start_time) then
       do i = 1, this%n_avg_fields
          call this%mean_fields(i)%update(time%t - this%time)
       end do
       this%time = time%t
    end if

  end subroutine user_stats_compute

end module user_stats
