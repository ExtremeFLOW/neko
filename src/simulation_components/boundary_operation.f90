! Copyright (c) 2026, The Neko Authors
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
!> Implements `boundary_operation_t`.
module boundary_operation
  use num_types, only : rp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use registry, only : neko_registry
  use field, only : field_t
  use case, only : case_t
  use json_utils, only : json_get, json_get_or_default
  use time_state, only : time_state_t
  use coefs, only : coef_t
  use neumann, only : neumann_t
  use logger, only : neko_log, LOG_SIZE
  use utils, only : neko_error
  use space, only : operator(.ne.)
  use file, only : file_t
  use vector, only : vector_t
  use vector_math, only : vector_masked_gather_copy_0, vector_glsc2, &
       vector_glsum, vector_glmin, vector_glmax, &
       vector_face_masked_gather_copy_0
  use time_based_controller, only : time_based_controller_t
  implicit none
  private

  !> A simulation component for boundary reductions on labelled zones.
  type, public, extends(simulation_component_t) :: boundary_operation_t
     !> Field to sample.
     type(field_t), pointer :: field => null()
     !> Surface coefficients.
     type(coef_t), pointer :: coef => null()
     !> Boundary mask built from the chosen zones.
     type(neumann_t) :: bc
     !> Labelled zones to include.
     integer, allocatable :: zone_indices(:)
     !> Name of the input field.
     character(len=:), allocatable :: field_name
     !> Requested operations.
     character(len=16), allocatable :: operations(:)
     !> Whether to compute the integral.
     logical :: compute_integral = .false.
     !> Whether to compute the average.
     logical :: compute_average = .false.
     !> Whether to compute the minimum.
     logical :: compute_min = .false.
     !> Whether to compute the maximum.
     logical :: compute_max = .false.
     !> Whether to write results to the log.
     logical :: log = .true.
     !> Whether CSV output is enabled.
     logical :: csv_output_enabled = .false.
     !> Optional CSV output file.
     type(file_t) :: csv_output
     !> Reusable row buffer for CSV output.
     type(vector_t) :: csv_row
     !> Precomputed boundary quadrature weights.
     type(vector_t) :: areas
     !> Reusable gathered boundary values.
     type(vector_t) :: surface_values
     !> Most recently computed integral.
     real(kind=rp) :: integral = 0.0_rp
     !> Most recently computed average.
     real(kind=rp) :: average = 0.0_rp
     !> Most recently computed minimum.
     real(kind=rp) :: minimum = huge(0.0_rp)
     !> Most recently computed maximum.
     real(kind=rp) :: maximum = -huge(0.0_rp)
   contains
     !> Construct the component from a case-file JSON object.
     procedure, pass(this) :: init => boundary_operation_init_from_json
     !> Generic constructor from components.
     generic :: init_from_components => &
          init_from_controllers, init_from_controllers_properties
     !> Construct from explicit time-based controllers.
     procedure, pass(this) :: init_from_controllers => &
          boundary_operation_init_from_controllers
     !> Construct from time-based controller properties.
     procedure, pass(this) :: init_from_controllers_properties => &
          boundary_operation_init_from_controllers_properties
     !> Common constructor used by all public constructors.
     procedure, private, pass(this) :: init_common => &
          boundary_operation_init_common
     !> Free the component.
     procedure, pass(this) :: free => boundary_operation_free
     !> Compute the requested boundary operations.
     procedure, pass(this) :: compute_ => boundary_operation_compute
  end type boundary_operation_t

contains

  !> Construct from JSON.
  !! @param json JSON object describing the simcomp.
  !! @param case Simulation case.
  subroutine boundary_operation_init_from_json(this, json, case)
    class(boundary_operation_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: name
    integer, allocatable :: zone_indices(:)
    character(len=:), allocatable :: field_name
    character(len=16), allocatable :: operations(:)
    logical :: log
    character(len=:), allocatable :: output_filename

    call this%free()

    call json_get_or_default(json, "name", name, "boundary_operation")
    call json_get_or_default(json, "log", log, .true.)
    call this%init_base(json, case)

    call json_get(json, "zone_indices", zone_indices)
    call json_get(json, "field_name", field_name)
    call json_get(json, "operations", operations)

    if (json%valid_path("output_filename")) then
       call json_get(json, "output_filename", output_filename)
       call this%init_common(name, case%fluid%c_Xh, zone_indices, field_name, &
            operations, log, output_filename)
    else
       call this%init_common(name, case%fluid%c_Xh, zone_indices, field_name, &
            operations, log)
    end if
  end subroutine boundary_operation_init_from_json

  !> Common constructor shared by all public constructors.
  !! @param name Unique simcomp name.
  !! @param coef SEM coefficients.
  !! @param zone_indices Labelled zones to include.
  !! @param field_name Name of the registered field to process.
  !! @param operations Requested operations in output order.
  !! @param log Whether to emit tabular log output.
  !! @param output_filename Optional CSV output filename.
  subroutine boundary_operation_init_common(this, name, coef, zone_indices, &
       field_name, operations, log, output_filename)
    class(boundary_operation_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    type(coef_t), intent(inout), target :: coef
    integer, intent(in) :: zone_indices(:)
    character(len=*), intent(in) :: field_name
    character(len=*), intent(in) :: operations(:)
    logical, intent(in) :: log
    character(len=*), intent(in), optional :: output_filename
    character(len=:), allocatable :: csv_header
    character(len=LOG_SIZE) :: log_buf
    integer :: i
    integer :: n_pts

    this%name = name
    this%log = log
    this%compute_integral = .false.
    this%compute_average = .false.
    this%compute_min = .false.
    this%compute_max = .false.
    this%csv_output_enabled = .false.
    this%integral = 0.0_rp
    this%average = 0.0_rp
    this%minimum = huge(0.0_rp)
    this%maximum = -huge(0.0_rp)

    this%coef => coef
    this%field => neko_registry%get_field_by_name(field_name)

    if (size(zone_indices) .eq. 0) then
       call neko_error("boundary_operation requires at least one zone index")
    end if

    if (size(operations) .eq. 0) then
       call neko_error("boundary_operation requires at least one operation")
    end if

    allocate(this%zone_indices(size(zone_indices)))
    this%zone_indices = zone_indices
    this%field_name = field_name
    allocate(character(len=16) :: this%operations(size(operations)))
    this%operations = operations

    do i = 1, size(this%operations)
       select case (trim(this%operations(i)))
       case ("integral")
          this%compute_integral = .true.
       case ("average")
          this%compute_average = .true.
       case ("min")
          this%compute_min = .true.
       case ("max")
          this%compute_max = .true.
       case default
          call neko_error("boundary_operation supports only operations = " // &
               "[integral, average, min, max]")
       end select
    end do

    call this%bc%init_from_components(this%coef, 0.0_rp)
    this%bc%zone_indices = this%zone_indices
    do i = 1, size(this%zone_indices)
       call this%bc%mark_zone(this%bc%msh%labeled_zones(this%zone_indices(i)))
    end do
    call this%bc%finalize()

    n_pts = this%bc%msk(0)
    if (n_pts .gt. 0) then
       call this%surface_values%init(n_pts)
       call vector_face_masked_gather_copy_0(this%areas, this%coef%area, &
            this%bc%msk, this%bc%facet, this%coef%Xh%lx, this%coef%Xh%ly, &
            this%coef%Xh%lz, n_pts)
    end if

    if (present(output_filename)) then
       csv_header = "tstep,time"
       do i = 1, size(this%operations)
          csv_header = trim(csv_header) // "," // trim(this%operations(i))
       end do
       call this%csv_output%init(trim(output_filename), header = trim(csv_header), &
            overwrite = .true.)
       call this%csv_row%init(2 + size(this%operations))
       this%csv_output_enabled = .true.
    end if

    call neko_log%section("Boundary operation")
    write(log_buf, '(A,A)') "Name: ", trim(this%name)
    call neko_log%message(log_buf)
    write(log_buf, '(A,A)') "Field: ", trim(this%field_name)
    call neko_log%message(log_buf)
    call neko_log%message("Operations:")
    do i = 1, size(this%operations)
       write(log_buf, '(A,A)') "  ", trim(this%operations(i))
       call neko_log%message(log_buf)
    end do
    write(log_buf, '(A,*(I0,:,", "))') "Zone indices: ", this%zone_indices
    call neko_log%message(log_buf)
    write(log_buf, '(A,I0)') "Marked boundary quadrature points: ", &
         this%bc%msk(0)
    call neko_log%message(log_buf)
    call neko_log%end_section()
  end subroutine boundary_operation_init_common

  !> Construct from explicit time-based controllers.
  !! @param name Unique simcomp name.
  !! @param case Simulation case owning the simcomp.
  !! @param order Execution priority.
  !! @param preprocess_controller Controller for preprocessing.
  !! @param compute_controller Controller for computation.
  !! @param output_controller Controller for output.
  !! @param zone_indices Labelled zones to include.
  !! @param field_name Name of the registered field to process.
  !! @param operations Requested operations in output order.
  !! @param log Optional flag controlling log output.
  !! @param output_filename Optional CSV output filename.
  subroutine boundary_operation_init_from_controllers(this, name, case, order, &
       preprocess_controller, compute_controller, output_controller, &
       zone_indices, field_name, operations, log, output_filename)
    class(boundary_operation_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(case_t), intent(inout), target :: case
    integer, intent(in) :: order
    type(time_based_controller_t), intent(in) :: preprocess_controller
    type(time_based_controller_t), intent(in) :: compute_controller
    type(time_based_controller_t), intent(in) :: output_controller
    integer, intent(in) :: zone_indices(:)
    character(len=*), intent(in) :: field_name
    character(len=*), intent(in) :: operations(:)
    logical, intent(in), optional :: log
    character(len=*), intent(in), optional :: output_filename
    logical :: log_enabled

    call this%free()

    log_enabled = .true.
    if (present(log)) log_enabled = log

    call this%init_base_from_components(case, order, preprocess_controller, &
         compute_controller, output_controller)

    if (present(output_filename)) then
       call this%init_common(name, case%fluid%c_Xh, zone_indices, field_name, &
            operations, log_enabled, output_filename)
    else
       call this%init_common(name, case%fluid%c_Xh, zone_indices, field_name, &
            operations, log_enabled)
    end if
  end subroutine boundary_operation_init_from_controllers

  !> Construct from time-based controller properties.
  !! @param name Unique simcomp name.
  !! @param case Simulation case owning the simcomp.
  !! @param order Execution priority.
  !! @param preprocess_control Control mode for preprocessing.
  !! @param preprocess_value Control value for preprocessing.
  !! @param compute_control Control mode for computation.
  !! @param compute_value Control value for computation.
  !! @param output_control Control mode for output.
  !! @param output_value Control value for output.
  !! @param zone_indices Labelled zones to include.
  !! @param field_name Name of the registered field to process.
  !! @param operations Requested operations in output order.
  !! @param log Optional flag controlling log output.
  !! @param output_filename Optional CSV output filename.
  subroutine boundary_operation_init_from_controllers_properties(this, name, &
       case, order, preprocess_control, preprocess_value, compute_control, &
       compute_value, output_control, output_value, zone_indices, field_name, &
       operations, log, output_filename)
    class(boundary_operation_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(case_t), intent(inout), target :: case
    integer, intent(in) :: order
    character(len=*), intent(in) :: preprocess_control
    real(kind=rp), intent(in) :: preprocess_value
    character(len=*), intent(in) :: compute_control
    real(kind=rp), intent(in) :: compute_value
    character(len=*), intent(in) :: output_control
    real(kind=rp), intent(in) :: output_value
    integer, intent(in) :: zone_indices(:)
    character(len=*), intent(in) :: field_name
    character(len=*), intent(in) :: operations(:)
    logical, intent(in), optional :: log
    character(len=*), intent(in), optional :: output_filename
    logical :: log_enabled

    call this%free()

    log_enabled = .true.
    if (present(log)) log_enabled = log

    call this%init_base_from_components(case, order, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, &
         output_value)

    if (present(output_filename)) then
       call this%init_common(name, case%fluid%c_Xh, zone_indices, field_name, &
            operations, log_enabled, output_filename)
    else
       call this%init_common(name, case%fluid%c_Xh, zone_indices, field_name, &
            operations, log_enabled)
    end if
  end subroutine boundary_operation_init_from_controllers_properties

  !> Free all resources owned by the component.
  subroutine boundary_operation_free(this)
    class(boundary_operation_t), intent(inout) :: this

    call this%bc%free()
    call this%csv_output%free()
    call this%csv_row%free()
    call this%areas%free()
    call this%surface_values%free()
    if (allocated(this%zone_indices)) deallocate(this%zone_indices)
    if (allocated(this%field_name)) deallocate(this%field_name)
    if (allocated(this%operations)) deallocate(this%operations)

    nullify(this%field)
    nullify(this%coef)

    call this%free_base()
  end subroutine boundary_operation_free

  !> Compute and optionally output the requested boundary operations.
  !! @param time Current simulation time state.
  subroutine boundary_operation_compute(this, time)
    class(boundary_operation_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer :: n_pts
    integer :: i
    real(kind=rp) :: area
    character(len=18) :: value_buf
    character(len=12) :: step_str
    character(len=:), allocatable :: section_title, header_line, value_line
    integer :: output_col

    n_pts = this%bc%msk(0)

    this%integral = 0.0_rp
    this%average = 0.0_rp
    this%minimum = huge(0.0_rp)
    this%maximum = -huge(0.0_rp)
    area = 0.0_rp

    if (n_pts .gt. 0) then
       call vector_masked_gather_copy_0(this%surface_values, &
            this%field%x, this%bc%msk, this%field%size(), n_pts)

       if (this%compute_integral .or. this%compute_average) then
          this%integral = vector_glsc2(this%surface_values, this%areas, n_pts)
       end if

       if (this%compute_average) then
          area = vector_glsum(this%areas, n_pts)
          if (area .gt. 0.0_rp) then
             this%average = this%integral / area
          else
             this%average = 0.0_rp
          end if
       end if


       if (this%compute_min) then
          this%minimum = vector_glmin(this%surface_values, n_pts)
       end if

       if (this%compute_max) then
          this%maximum = vector_glmax(this%surface_values, n_pts)
       end if
    end if

    if (this%log) then
       section_title = trim(this%name)
       call neko_log%section(section_title)

       header_line = repeat(' ', 12) // ' |'
       write(step_str, '(I12)') time%tstep
       step_str = adjustl(step_str)
       value_line = step_str // ' |'

       do i = 1, size(this%operations)
          select case (trim(this%operations(i)))
          case ('integral')
             header_line = header_line // left_pad('Integral:', 18)
             write(value_buf, '(ES18.9)') this%integral
             value_line = value_line // value_buf
          case ('average')
             header_line = header_line // left_pad('Average:', 18)
             write(value_buf, '(ES18.9)') this%average
             value_line = value_line // value_buf
          case ('min')
             header_line = header_line // left_pad('Min:', 18)
             write(value_buf, '(ES18.9)') this%minimum
             value_line = value_line // value_buf
          case ('max')
             header_line = header_line // left_pad('Max:', 18)
             write(value_buf, '(ES18.9)') this%maximum
             value_line = value_line // value_buf
          end select
       end do

       call neko_log%message(header_line)
       call neko_log%message(value_line)
       call neko_log%end_section()
    end if

    if (this%csv_output_enabled) then
       if (this%output_controller%check(time)) then
          this%csv_row%x(1) = real(time%tstep, rp)
          this%csv_row%x(2) = time%t
          output_col = 3
          do i = 1, size(this%operations)
             select case (trim(this%operations(i)))
             case ('integral')
                this%csv_row%x(output_col) = this%integral
             case ('average')
                this%csv_row%x(output_col) = this%average
             case ('min')
                this%csv_row%x(output_col) = this%minimum
             case ('max')
                this%csv_row%x(output_col) = this%maximum
             end select
             output_col = output_col + 1
          end do
          call this%csv_output%write(this%csv_row)
          call this%output_controller%register_execution()
       end if
    end if
  end subroutine boundary_operation_compute

  !> Left-pad a string to a fixed width.
  !! @param text Input string.
  !! @param width Width of the padded result.
  pure function left_pad(text, width) result(padded)
    character(len=*), intent(in) :: text
    integer, intent(in) :: width
    character(len=:), allocatable :: padded
    integer :: pad_width

    pad_width = max(0, width - len_trim(text))
    padded = repeat(' ', pad_width) // trim(text)
  end function left_pad

end module boundary_operation
