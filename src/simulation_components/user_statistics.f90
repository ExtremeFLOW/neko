! Copyright (c) 2024, The Neko Authors
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
!> Implements the `vorticity_t` type.

module user_statistics
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use case, only : case_t
  use fld_file_output, only : fld_file_output_t
  use json_utils, only : json_get, json_get_or_default
  use mean_field, only : mean_field_t


  implicit none
  private

  !> A simulation component that computes the averages of fields in the registry.
  type, public, extends(simulation_component_t) :: user_statistics_t

     !> Variables
     real(kind=rp) :: average_starting_t
     real(kind=rp) :: average_last_t
     type(mean_field_t), allocatable :: mean_fields(:)
     integer :: n_fields_in_registry = 0
     character(len=20), allocatable  :: field_in_registry_name(:)

     !> Output writer.
     type(fld_file_output_t), private :: output

   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => user_statistics_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
        user_statistics_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => user_statistics_free
     !> Compute the vorticity field.
     procedure, pass(this) :: compute_ => user_statistics_compute
  end type vorticity_t

contains

  !> Constructor from json.
  subroutine user_statistics_init_from_json(this, json, case)
    class(vorticity_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: filename
    character(len=:), allocatable :: precision

    call this%init_base(json, case)

    !> Get the number of stat fields and their names
    call json%info('fields', n_children=this%n_fields_in_registry)
    call json_get(json, 'fields', this%field_in_registry_name)

    if (json%valid_path("output_filename")) then
       call json_get(json, "output_filename", filename)
       if (json%valid_path("output_precision")) then
           call json_get(json, "output_precision", precision)
           if (precision == "double") then
              call user_statistics_init_from_attributes(this, filename, dp)
           else
              call user_statistics_init_from_attributes(this, filename, sp)
           end if
       else
           call user_statistics_init_from_attributes(this, filename)
       end if
    else
       call user_statistics_init_from_attributes(this)
    end if
  end subroutine user_statistics_init_from_json

  !> Actual constructor.
  subroutine user_statistics_init_from_attributes(this, filename, precision)
    class(user_statistics_t), intent(inout) :: this
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision
    type(fld_file_output_t) :: output
    integer :: i


    !> Allocate and initialize the mean fields
    allocate(this%mean_fields(this%n_fields_in_registry))
    do i = 1, this%n_fields_in_registry
       call this%mean_fields(i)%init(neko_field_registry%get_field(&
                                     trim(this%field_in_registry_name(i))))
    end do

    if (present(filename)) then
       if (present(precision)) then
          call this%output%init(precision, filename, this%n_fields_in_registry)
       else
          call this%output%init(sp, filename, this%n_fields_in_registry)
       end if
        
       do i = 1, this%n_fields_in_registry
          this%output%fields%fields(i)%f => this%mean_fields(i)%mf
       end if

       call this%case%s%add(this%output, this%output_controller%control_value, &
                            this%output_controller%control_mode)
    else

       do i = 1, this%n_fields_in_registry
          call this%case%f_out%fluid%append(this%mean_fields(i)%mf)
       end if

    end if

  end subroutine user_statistics_init_from_attributes

  !> Destructor.
  subroutine user_statistics_free(this)
    class(user_statistics_t), intent(inout) :: this
    integer :: i
    call this%free_base()
    do i = 1, this%n_fields_in_registry
       call this%mean_fields(i)%free
    end if
  end subroutine user_statistics_free

  !> Compute the vorticity field.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine user_statistics_compute(this, t, tstep)
    class(vorticity_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
 
    !> Update the running average of the fields
    do i = 1, this%n_fields_in_registry
       call this%mean_fields(i)%update(t-this%average_last_t) 
    end do
    this%average_last_t = t

  end subroutine user_statistics_compute

end module user_statistics
