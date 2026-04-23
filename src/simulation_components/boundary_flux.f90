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
!> Implements `boundary_flux_t`.
module boundary_flux
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
  use file, only : file_t
  use vector, only : vector_t
  use vector_math, only : vector_masked_gather_copy_0, vector_glsc2, &
       vector_cmult
  use time_based_controller, only : time_based_controller_t
  use drag_torque, only : setup_normals
  use neko_config, only : NEKO_BCKND_DEVICE
  use device, only : device_memcpy, HOST_TO_DEVICE
  implicit none
  private

  !> A simulation component for total vector flux through labelled zones.
  !! @note Fluxes inside the domain are considered positive.
  type, public, extends(simulation_component_t) :: boundary_flux_t
     !> X component field.
     type(field_t), pointer :: u => null()
     !> Y component field.
     type(field_t), pointer :: v => null()
     !> Z component field.
     type(field_t), pointer :: w => null()
     !> Surface coefficients.
     type(coef_t), pointer :: coef => null()
     !> Boundary mask built from the chosen zones.
     type(neumann_t) :: bc
     !> Labelled zones to include.
     integer, allocatable :: zone_indices(:)
     !> Names of the input vector component fields.
     character(len=80), allocatable :: field_names(:)
     !> Whether to write results to the log.
     logical :: log = .true.
     !> Whether CSV output is enabled.
     logical :: csv_output_enabled = .false.
     !> Optional CSV output file.
     type(file_t) :: csv_output
     !> Reusable row buffer for CSV output.
     type(vector_t) :: csv_row
     !> Area-weighted x normals.
     type(vector_t) :: n1
     !> Area-weighted y normals.
     type(vector_t) :: n2
     !> Area-weighted z normals.
     type(vector_t) :: n3
     !> Reusable gathered x boundary values.
     type(vector_t) :: surface_u
     !> Reusable gathered y boundary values.
     type(vector_t) :: surface_v
     !> Reusable gathered z boundary values.
     type(vector_t) :: surface_w
     !> Most recently computed total flux.
     real(kind=rp) :: flux = 0.0_rp
   contains
     !> Construct the component from a case-file JSON object.
     procedure, pass(this) :: init => boundary_flux_init_from_json
     !> Generic constructor from components.
     generic :: init_from_components => &
          init_from_controllers, init_from_controllers_properties
     !> Construct from explicit time-based controllers.
     procedure, pass(this) :: init_from_controllers => &
          boundary_flux_init_from_controllers
     !> Construct from time-based controller properties.
     procedure, pass(this) :: init_from_controllers_properties => &
          boundary_flux_init_from_controllers_properties
     !> Common constructor used by all public constructors.
     procedure, private, pass(this) :: init_common => &
          boundary_flux_init_common
     !> Free the component.
     procedure, pass(this) :: free => boundary_flux_free
     !> Compute the total boundary flux.
     procedure, pass(this) :: compute_ => boundary_flux_compute
  end type boundary_flux_t

contains

  !> Construct from JSON.
  !! @param json JSON object describing the simcomp.
  !! @param case Simulation case.
  subroutine boundary_flux_init_from_json(this, json, case)
    class(boundary_flux_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: name
    integer, allocatable :: zone_indices(:)
    character(len=80), allocatable :: field_names(:)
    logical :: log
    character(len=:), allocatable :: output_filename

    call this%free()

    call json_get_or_default(json, "name", name, "boundary_flux")
    call json_get_or_default(json, "log", log, .true.)
    call this%init_base(json, case)

    call json_get(json, "zone_indices", zone_indices)
    call json_get(json, "field_names", field_names)

    if (json%valid_path("output_filename")) then
       call json_get(json, "output_filename", output_filename)
       call this%init_common(name, case%fluid%c_Xh, zone_indices, &
            field_names, log, output_filename)
    else
       call this%init_common(name, case%fluid%c_Xh, zone_indices, &
            field_names, log)
    end if
  end subroutine boundary_flux_init_from_json

  !> Common constructor shared by all public constructors.
  !! @param name Unique simcomp name.
  !! @param coef SEM coefficients.
  !! @param zone_indices Labelled zones to include.
  !! @param field_names Names of the registered vector component fields.
  !! @param log Whether to emit tabular log output.
  !! @param output_filename Optional CSV output filename.
  subroutine boundary_flux_init_common(this, name, coef, zone_indices, &
       field_names, log, output_filename)
    class(boundary_flux_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    type(coef_t), intent(inout), target :: coef
    integer, intent(in) :: zone_indices(:)
    character(len=*), intent(in) :: field_names(:)
    logical, intent(in) :: log
    character(len=*), intent(in), optional :: output_filename
    character(len=LOG_SIZE) :: log_buf
    integer :: i
    integer :: n_pts

    this%name = name
    this%log = log
    this%csv_output_enabled = .false.
    this%flux = 0.0_rp

    this%coef => coef

    if (size(zone_indices) .eq. 0) then
       call neko_error("boundary_flux requires at least one zone index")
    end if

    if (size(field_names) .ne. 3) then
       call neko_error("boundary_flux requires exactly three field names")
    end if

    allocate(this%zone_indices(size(zone_indices)))
    this%zone_indices = zone_indices
    allocate(this%field_names(size(field_names)))
    this%field_names = field_names

    this%u => neko_registry%get_field_by_name(trim(this%field_names(1)))
    this%v => neko_registry%get_field_by_name(trim(this%field_names(2)))
    this%w => neko_registry%get_field_by_name(trim(this%field_names(3)))

    call this%bc%init_from_components(this%coef, 0.0_rp)
    this%bc%zone_indices = this%zone_indices
    do i = 1, size(this%zone_indices)
       call this%bc%mark_zone(this%bc%msh%labeled_zones(this%zone_indices(i)))
    end do
    call this%bc%finalize()

    n_pts = this%bc%facet_msk(0)
    if (n_pts .gt. 0) then
       call this%n1%init(n_pts)
       call this%n2%init(n_pts)
       call this%n3%init(n_pts)
       call this%surface_u%init(n_pts)
       call this%surface_v%init(n_pts)
       call this%surface_w%init(n_pts)
       call setup_normals(this%coef, this%bc%facet_msk, this%bc%facet, &
            this%n1%x, this%n2%x, this%n3%x, n_pts)
       call vector_cmult(this%n1, -1.0_rp, n_pts)
       call vector_cmult(this%n2, -1.0_rp, n_pts)
       call vector_cmult(this%n3, -1.0_rp, n_pts)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(this%n1%x, this%n1%x_d, n_pts, &
               HOST_TO_DEVICE, .false.)
          call device_memcpy(this%n2%x, this%n2%x_d, n_pts, &
               HOST_TO_DEVICE, .false.)
          call device_memcpy(this%n3%x, this%n3%x_d, n_pts, &
               HOST_TO_DEVICE, .true.)
       end if
    end if

    if (present(output_filename)) then
       call this%csv_output%init(trim(output_filename), &
            header = "tstep,time,flux", overwrite = .true.)
       call this%csv_row%init(3)
       this%csv_output_enabled = .true.
    end if

    call neko_log%section("Boundary flux")
    write(log_buf, '(A,A)') "Name: ", trim(this%name)
    call neko_log%message(log_buf)
    write(log_buf, '(A,A,", ",A,", ",A)') "Fields: ", &
         trim(this%field_names(1)), trim(this%field_names(2)), &
         trim(this%field_names(3))
    call neko_log%message(log_buf)
    write(log_buf, '(A,*(I0,:,", "))') "Zone indices: ", this%zone_indices
    call neko_log%message(log_buf)
    write(log_buf, '(A,I0)') "Marked boundary quadrature points: ", &
         this%bc%facet_msk(0)
    call neko_log%message(log_buf)
    call neko_log%end_section()
  end subroutine boundary_flux_init_common

  !> Construct from explicit time-based controllers.
  !! @param name Unique simcomp name.
  !! @param case Simulation case owning the simcomp.
  !! @param order Execution priority.
  !! @param preprocess_controller Controller for preprocessing.
  !! @param compute_controller Controller for computation.
  !! @param output_controller Controller for output.
  !! @param zone_indices Labelled zones to include.
  !! @param field_names Names of the registered vector component fields.
  !! @param log Optional flag controlling log output.
  !! @param output_filename Optional CSV output filename.
  subroutine boundary_flux_init_from_controllers(this, name, case, order, &
       preprocess_controller, compute_controller, output_controller, &
       zone_indices, field_names, log, output_filename)
    class(boundary_flux_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(case_t), intent(inout), target :: case
    integer, intent(in) :: order
    type(time_based_controller_t), intent(in) :: preprocess_controller
    type(time_based_controller_t), intent(in) :: compute_controller
    type(time_based_controller_t), intent(in) :: output_controller
    integer, intent(in) :: zone_indices(:)
    character(len=*), intent(in) :: field_names(:)
    logical, intent(in), optional :: log
    character(len=*), intent(in), optional :: output_filename
    logical :: log_enabled

    call this%free()

    log_enabled = .true.
    if (present(log)) log_enabled = log

    call this%init_base_from_components(case, order, preprocess_controller, &
         compute_controller, output_controller)

    if (present(output_filename)) then
       call this%init_common(name, case%fluid%c_Xh, zone_indices, &
            field_names, log_enabled, output_filename)
    else
       call this%init_common(name, case%fluid%c_Xh, zone_indices, &
            field_names, log_enabled)
    end if
  end subroutine boundary_flux_init_from_controllers

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
  !! @param field_names Names of the registered vector component fields.
  !! @param log Optional flag controlling log output.
  !! @param output_filename Optional CSV output filename.
  subroutine boundary_flux_init_from_controllers_properties(this, name, case, &
       order, preprocess_control, preprocess_value, compute_control, &
       compute_value, output_control, output_value, zone_indices, field_names, &
       log, output_filename)
    class(boundary_flux_t), intent(inout) :: this
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
    character(len=*), intent(in) :: field_names(:)
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
       call this%init_common(name, case%fluid%c_Xh, zone_indices, &
            field_names, log_enabled, output_filename)
    else
       call this%init_common(name, case%fluid%c_Xh, zone_indices, &
            field_names, log_enabled)
    end if
  end subroutine boundary_flux_init_from_controllers_properties

  !> Free all resources owned by the component.
  subroutine boundary_flux_free(this)
    class(boundary_flux_t), intent(inout) :: this

    call this%bc%free()
    call this%csv_output%free()
    call this%csv_row%free()
    call this%n1%free()
    call this%n2%free()
    call this%n3%free()
    call this%surface_u%free()
    call this%surface_v%free()
    call this%surface_w%free()
    if (allocated(this%zone_indices)) deallocate(this%zone_indices)
    if (allocated(this%field_names)) deallocate(this%field_names)

    nullify(this%u)
    nullify(this%v)
    nullify(this%w)
    nullify(this%coef)

    call this%free_base()
  end subroutine boundary_flux_free

  !> Compute and optionally output the total boundary flux.
  !! @param time Current simulation time state.
  subroutine boundary_flux_compute(this, time)
    class(boundary_flux_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer :: n_pts
    character(len=18) :: value_buf
    character(len=12) :: step_str
    character(len=:), allocatable :: header_line, value_line

    n_pts = this%bc%facet_msk(0)
    this%flux = 0.0_rp

    call vector_masked_gather_copy_0(this%surface_u, this%u%x, &
         this%bc%facet_msk, this%u%size(), n_pts)
    call vector_masked_gather_copy_0(this%surface_v, this%v%x, &
         this%bc%facet_msk, this%v%size(), n_pts)
    call vector_masked_gather_copy_0(this%surface_w, this%w%x, &
         this%bc%facet_msk, this%w%size(), n_pts)

    this%flux = vector_glsc2(this%surface_u, this%n1, n_pts) + &
         vector_glsc2(this%surface_v, this%n2, n_pts) + &
         vector_glsc2(this%surface_w, this%n3, n_pts)

    if (this%log) then
       call neko_log%section(trim(this%name))

       header_line = repeat(' ', 12) // ' |' // left_pad('Flux:', 18)
       write(step_str, '(I12)') time%tstep
       step_str = adjustl(step_str)
       write(value_buf, '(ES18.9)') this%flux
       value_line = step_str // ' |' // value_buf

       call neko_log%message(header_line)
       call neko_log%message(value_line)
       call neko_log%end_section()
    end if

    if (this%csv_output_enabled) then
       if (this%output_controller%check(time)) then
          this%csv_row%x(1) = real(time%tstep, rp)
          this%csv_row%x(2) = time%t
          this%csv_row%x(3) = this%flux
          call this%csv_output%write(this%csv_row)
          call this%output_controller%register_execution()
       end if
    end if
  end subroutine boundary_flux_compute

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

end module boundary_flux
