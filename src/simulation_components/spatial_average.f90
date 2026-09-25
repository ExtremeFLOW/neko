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
!> Implements the `spatial_average_t` type.
module spatial_average
  use num_types, only : rp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use case, only : case_t
  use time_state, only : time_state_t
  use json_utils, only : json_get, json_get_or_default
  use coefs, only : coef_t
  use utils, only : NEKO_FNAME_LEN, NEKO_VARNAME_LEN
  use spatial_average_output, only : spatial_average_output_t
  implicit none
  private

  !> A simulation component that writes spatial averages of registry fields.
  type, public, extends(simulation_component_t) :: spatial_average_t
     !> Output writer.
     type(spatial_average_output_t), private :: output
   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => spatial_average_init_from_json
     !> Constructor from components.
     procedure, pass(this) :: init_from_components => &
          spatial_average_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => spatial_average_free
     !> Here to comply with the interface, does nothing.
     procedure, pass(this) :: compute_ => spatial_average_compute
  end type spatial_average_t

contains

  !> Constructor from json.
  !! @param json The json parameter dictionary.
  !! @param case The neko case object.
  subroutine spatial_average_init_from_json(this, json, case)
    class(spatial_average_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: name
    character(len=:), allocatable :: avg_direction
    character(len=:), allocatable :: output_filename
    character(len=NEKO_VARNAME_LEN), allocatable :: fields(:)

    call this%init_base(json, case)

    call json_get_or_default(json, 'name', name, 'spatial_average')
    call json_get(json, 'fields', fields)
    call json_get(json, 'avg_direction', avg_direction)
    call json_get_or_default(json, 'output_filename', output_filename, &
         'spatial_average')

    call this%init_from_components(name, fields, case%fluid%c_Xh, &
         avg_direction, output_filename)

    if (allocated(fields)) deallocate(fields)
    if (allocated(output_filename)) deallocate(output_filename)
    if (allocated(avg_direction)) deallocate(avg_direction)
    if (allocated(name)) deallocate(name)
  end subroutine spatial_average_init_from_json

  !> Constructor from components.
  !! @param name The unique name of the simcomp.
  !! @param fields Names of the fields to average.
  !! @param coef The SEM coefficients.
  !! @param avg_direction Direction(s) to average in.
  !! @param output_filename Base name of the output file.
  subroutine spatial_average_init_from_components(this, name, fields, coef, &
       avg_direction, output_filename)
    class(spatial_average_t), target, intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: fields(:)
    type(coef_t), target, intent(inout) :: coef
    character(len=*), intent(in) :: avg_direction
    character(len=*), intent(in) :: output_filename
    character(len=NEKO_FNAME_LEN) :: resolved_output_filename

    this%name = name

    if (allocated(this%case%output_directory) .and. &
         len_trim(this%case%output_directory) .gt. 0) then
       resolved_output_filename = trim(this%case%output_directory) // &
            trim(output_filename)
    else
       resolved_output_filename = trim(output_filename)
    end if

    call this%output%init(fields, coef, avg_direction, &
         resolved_output_filename)
    call this%case%output_controller%add(this%output, &
         this%output_controller%control_value, &
         this%output_controller%control_mode)
  end subroutine spatial_average_init_from_components

  !> Destructor.
  subroutine spatial_average_free(this)
    class(spatial_average_t), intent(inout) :: this

    call this%free_base()
    call this%output%free()
  end subroutine spatial_average_free

  !> Here to comply with the interface, does nothing.
  !! @param time The current time state.
  subroutine spatial_average_compute(this, time)
    class(spatial_average_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
  end subroutine spatial_average_compute

end module spatial_average
