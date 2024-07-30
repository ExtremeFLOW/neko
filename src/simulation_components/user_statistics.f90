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


  implicit none
  private

  !> A simulation component that computes the averages of fields in the registry.
  type, public, extends(simulation_component_t) :: user_stats_t

     !> Variables
     real(kind=rp) :: start_time !<When to start average
     real(kind=rp) :: time !< current time
     type(mean_field_t), allocatable :: mean_fields(:) !<mean fields
     integer :: n_avg_fields = 0 !< NUmber of fields to average
     character(len=20), allocatable  :: field_names(:) !< field names

     !> Output writer.
     type(mean_field_output_t), private :: output

   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => user_stats_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
        user_stats_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => user_stats_free
     !> Compute the means
     procedure, pass(this) :: compute_ => user_stats_compute
     procedure, pass(this) :: restart_ => user_stats_restart
  end type user_stats_t

contains

  !> Constructor from json.
  subroutine user_stats_init_from_json(this, json, case)
    class(user_stats_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: filename
    character(len=:), allocatable :: precision

    call this%init_base(json, case)

    !> Get the number of stat fields and their names
    call json%info('fields', n_children=this%n_avg_fields)
    call json_get(json, 'fields', this%field_names)
    call json_get_or_default(json, 'start_time', &
         this%start_time, 0.0_rp)
    this%time = this%start_time
    call user_stats_init_from_attributes(this,this%start_time,filename="user_stats")
   ! if (json%valid_path("output_filename")) then
   !    call json_get_or_default(json, "output_filename", filename,'user_stats')
   !    if (json%valid_path("output_precision")) then
   !        call json_get_or_default(json, "output_precision", precision,"single")
   !        if (precision == "double") then
   !           call user_stats_init_from_attributes(this, filename, dp)
   !        else
   !           call user_stats_init_from_attributes(this, filename, sp)
   !        end if
   !    else
   !        call user_stats_init_from_attributes(this, filename)
   !    end if
   ! else
   !    call user_stats_init_from_attributes(this)
   ! end if
  end subroutine user_stats_init_from_json

  subroutine user_stats_restart(this, t)
    class(user_stats_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    if(t .gt. this%time) this%time = t
  end subroutine user_stats_restart


  !> Actual constructor.
  subroutine user_stats_init_from_attributes(this,start_time, filename, precision)
    class(user_stats_t), intent(inout) :: this
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: precision
    real(kind=rp), intent(in) :: start_time
    integer :: i
    type(field_t), pointer :: field_to_avg

    this%start_time = start_time
    this%time = start_time
    !> Allocate and initialize the mean fields
    allocate(this%mean_fields(this%n_avg_fields))
    do i = 1, this%n_avg_fields
       field_to_avg => neko_field_registry%get_field(&
                                     trim(this%field_names(i)))
       call this%mean_fields(i)%init(field_to_avg)
    end do

    call this%output%init(this%mean_fields, this%n_avg_fields,this%start_time, sp, name=filename)
    call this%case%s%add(this%output, this%output_controller%control_value, &
                         this%output_controller%control_mode)


  end subroutine user_stats_init_from_attributes

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

  !> Compute the vorticity field.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine user_stats_compute(this, t, tstep)
    class(user_stats_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: i
 
    !> Update the running average of the fields
    if (t .ge. this%start_time) then
       do i = 1, this%n_avg_fields
          call this%mean_fields(i)%update(t-this%time) 
       end do
       this%time= t
    end if

  end subroutine user_stats_compute

end module user_stats
