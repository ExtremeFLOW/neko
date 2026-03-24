! Copyright (c) 2024-2026, The Neko Authors
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
!> Implements the `field_writer_t` type.

module field_writer
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use time_state, only : time_state_t
  use registry, only : neko_registry
  use case, only : case_t
  use field_output, only : field_output_t
  use json_utils, only : json_get, json_get_or_default
  use time_based_controller, only : time_based_controller_t
  use logger, only : neko_log, LOG_SIZE, NEKO_LOG_VERBOSE
  use utils, only : neko_error
  use amr_reconstruct, only : amr_reconstruct_t

  implicit none
  private

  !> A simulation component that writes a 3d field to a file.
  type, public, extends(simulation_component_t) :: field_writer_t
     !> Output writer.
     type(field_output_t), private :: output

     ! Default values for optional parameters in the constructor.
     character(len=20), private :: default_name = "field_writer"
     character(len=20), private :: default_precision = "single"
     character(len=20), private :: default_filename = ""
     character(len=20), private :: default_format = "nek5000"
     logical, private :: default_subdivide = .false.

   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => field_writer_init_from_json
     !> Generic for constructing from components.
     generic :: init_from_components => &
          init_from_controllers, init_from_controllers_properties
     !> Constructor from components, passing time_based_controllers.
     procedure, pass(this) :: init_from_controllers => &
          field_writer_init_from_controllers
     !> Constructor from components, passing the properties of
     !! time_based_controllers.
     procedure, pass(this) :: init_from_controllers_properties => &
          field_writer_init_from_controllers_properties
     !> Common part of both constructors.
     procedure, private, pass(this) :: init_common => field_writer_init_common
     !> Destructor.
     procedure, pass(this) :: free => field_writer_free
     !> Here to compy with the interface, does nothing.
     procedure, pass(this) :: compute_ => field_writer_compute
     !> AMR restart
     procedure, pass(this) :: amr_restart => field_writer_amr_restart
  end type field_writer_t

contains

  !> Constructor from json.
  !> @param json JSON object with the parameters.
  !! @param case The case object.
  subroutine field_writer_init_from_json(this, json, case)
    class(field_writer_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: filename
    character(len=:), allocatable :: name
    character(len=:), allocatable :: precision
    character(len=:), allocatable :: format
    character(len=20), allocatable :: fields(:)
    integer :: precision_value
    logical :: subdivide

    call this%init_base(json, case)
    call json_get(json, "fields", fields)
    call json_get_or_default(json, "name", name, this%default_name)
    call json_get_or_default(json, "output_filename", filename, &
         this%default_filename)
    call json_get_or_default(json, "output_precision", precision, &
         this%default_precision)
    call json_get_or_default(json, "output_format", format, this%default_format)

    if (precision .eq. "single") then
       precision_value = sp
    else if (precision .eq. "double") then
       precision_value = dp
    else
       call neko_error("Invalid precision specified for field_writer: " &
            // trim(precision))
    end if

    call json_get_or_default(json, "output_subdivide", subdivide, &
         this%default_subdivide)

    call this%init_common(name, fields, filename, precision_value, format, &
         subdivide)
  end subroutine field_writer_init_from_json

  !> Constructor from components, passing controllers.
  !! @param name The unique name of the simcomp.
  !! @param case The simulation case object.
  !! @param order The execution oder priority of the simcomp.
  !! @param preprocess_controller The controller for running preprocessing.
  !! @param compute_controller The controller for running compute.
  !! @param output_controller The controller for producing output.
  !! @param fields Array of field names to be sampled.
  !! @param filename The name of the file save the fields to. Optional, if not
  !! provided, fields are added to the main output file.
  !! @param precision The real precision of the output data. Optional, defaults
  !! to single precision.
  !! @param format The output format of the data. Optional, defaults to
  !! "nek5000".
  !! @param subdivide Whether to subdivide spectral elements into linear
  !! sub-cells. Optional, defaults to `.false.`.
  subroutine field_writer_init_from_controllers(this, name, case, order, &
       preprocess_controller, compute_controller, output_controller, &
       fields, filename, precision, format, subdivide)
    class(field_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(case_t), intent(inout), target :: case
    integer :: order
    type(time_based_controller_t), intent(in) :: preprocess_controller
    type(time_based_controller_t), intent(in) :: compute_controller
    type(time_based_controller_t), intent(in) :: output_controller
    character(len=20), intent(in) :: fields(:)
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision
    character(len=20), intent(in), optional :: format
    logical, intent(in), optional :: subdivide

    call this%init_base_from_components(case, order, preprocess_controller, &
         compute_controller, output_controller)
    call this%init_common(name, fields, filename, precision, format, subdivide)

  end subroutine field_writer_init_from_controllers

  !> Constructor from components, passing properties to the
  !! time_based_controller` components in the base type.
  !! @param name The unique name of the simcomp.
  !! @param case The simulation case object.
  !! @param order The execution oder priority of the simcomp.
  !! @param preprocess_controller Control mode for preprocessing.
  !! @param preprocess_value Value parameter for preprocessing.
  !! @param compute_controller Control mode for computing.
  !! @param compute_value Value parameter for computing.
  !! @param output_controller Control mode for output.
  !! @param output_value Value parameter for output.
  !! @param fields Array of field names to be sampled.
  !! @param filename The name of the file save the fields to. Optional, if not
  !! provided, fields are added to the main output file.
  !! @param precision The real precision of the output data. Optional, defaults
  !! to single precision.
  !! @param format The output format of the data. Optional, defaults to
  !! "nek5000".
  !! @param subdivide Whether to subdivide spectral elements into linear
  !! sub-cells. Optional, defaults to `.false.`.
  subroutine field_writer_init_from_controllers_properties(this, name, &
       case, order, preprocess_control, preprocess_value, compute_control, &
       compute_value, output_control, output_value, fields, filename, &
       precision, format, subdivide)
    class(field_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(case_t), intent(inout), target :: case
    integer :: order
    character(len=*), intent(in) :: preprocess_control
    real(kind=rp), intent(in) :: preprocess_value
    character(len=*), intent(in) :: compute_control
    real(kind=rp), intent(in) :: compute_value
    character(len=*), intent(in) :: output_control
    real(kind=rp), intent(in) :: output_value
    character(len=20), intent(in) :: fields(:)
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision
    character(len=*), intent(in), optional :: format
    logical, intent(in), optional :: subdivide

    call this%init_base_from_components(case, order, preprocess_control, &
         preprocess_value, compute_control, compute_value, output_control, &
         output_value)
    call this%init_common(name, fields, filename, precision, format, subdivide)

  end subroutine field_writer_init_from_controllers_properties

  !> Common part of both constructors.
  !! @param name The unique name of the simcomp.
  !! @param fields Array of field names to be sampled.
  !! @param filename The name of the file save the fields to. Optional, if not
  !! provided, fields are added to the main output file.
  !! @param precision The real precision of the output data. Optional, defaults
  !! to single precision.
  !! @param format The output format of the data. Optional, defaults to
  !! "nek5000".
  !! @param subdivide Whether to subdivide spectral elements into linear
  !! sub-cells. Optional, defaults to `.false.`.
  subroutine field_writer_init_common(this, name, fields, filename, precision, &
       format, subdivide)
    class(field_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=20), intent(in) :: fields(:)
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision
    character(len=*), intent(in), optional :: format
    logical, intent(in), optional :: subdivide
    character(len=20) :: fieldi
    logical :: filename_provided
    character(len=120) :: message
    integer :: i

    this%name = name
    ! Register fields if they don't exist.
    do i = 1, size(fields)
       fieldi = trim(fields(i))
       call neko_registry%add_field(this%case%fluid%dm_Xh, fieldi, &
            ignore_existing = .true.)
    end do

    filename_provided = .false.
    if (present(filename)) then
       if (len_trim(filename) .ne. 0) then
          filename_provided = .true.
          call this%output%init(trim(filename), size(fields), &
               precision = precision, format = format)

          if (present(subdivide)) then
             call this%output%file_%set_subdivide(subdivide)
          end if

          do i = 1, size(fields)
             fieldi = trim(fields(i))
             call this%output%fields%assign(i, &
                  neko_registry%get_field(fieldi))
          end do

          call this%case%output_controller%add(this%output, &
               this%output_controller%control_value, &
               this%output_controller%control_mode)

       end if
    end if

    if (.not. filename_provided) then
       do i = 1, size(fields)
          fieldi = trim(fields(i))
          call this%case%f_out%fluid%append( &
               neko_registry%get_field(fieldi))
       end do
    end if

  end subroutine field_writer_init_common

  !> Destructor.
  subroutine field_writer_free(this)
    class(field_writer_t), intent(inout) :: this
    call this%free_base()

    call this%output%free()

    call this%free_amr_base()

  end subroutine field_writer_free

  !> Here to comply with the interface, does nothing.
  subroutine field_writer_compute(this, time)
    class(field_writer_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

  end subroutine field_writer_compute

  !> AMR restart
  !! @param[inout]  reconstruct   data reconstruction type
  !! @param[in]     counter       restart counter
  !! @param[in]     tstep         time step
  subroutine field_writer_amr_restart(this, reconstruct, counter, tstep)
    class(field_writer_t), intent(inout) :: this
    type(amr_reconstruct_t), intent(inout) :: reconstruct
    integer, intent(in) :: counter, tstep
    character(len=LOG_SIZE) :: log_buf

    call neko_error('Nothing done for AMR reconstruction')

    ! Was this component already restarted?
    if (this%counter .eq. counter) return

    this%counter = counter

    log_buf = 'Field writer'
    call neko_log%message(log_buf, NEKO_LOG_VERBOSE)
!    call neko_log%section(log_buf, NEKO_LOG_VERBOSE)
!    call neko_log%end_section(lvl = NEKO_LOG_VERBOSE)

  end subroutine field_writer_amr_restart

end module field_writer
