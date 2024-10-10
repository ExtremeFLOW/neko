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
!> Implements the `field_writer_t` type.

module field_writer
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : curl
  use case, only : case_t
  use fld_file_output, only : fld_file_output_t
  use json_utils, only : json_get, json_get_or_default
  implicit none
  private

  !> A simulation component that writes a 3d field to a file.
  type, public, extends(simulation_component_t) :: field_writer_t
     !> Output writer.
     type(fld_file_output_t), private :: output

   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => field_writer_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
        field_writer_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => field_writer_free
     !> Here to compy with the interface, does nothing.
     procedure, pass(this) :: compute_ => field_writer_compute
  end type field_writer_t

contains

  !> Constructor from json.
  !> @param json JSON object with the parameters.
  !! @param case The case object.
  subroutine field_writer_init_from_json(this, json, case)
    class(field_writer_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: filename
    character(len=:), allocatable :: precision
    character(len=20), allocatable :: fields(:)

    call this%init_base(json, case)
    call json_get(json, "fields", fields)

    if (json%valid_path("output_filename")) then
       call json_get(json, "output_filename", filename)
       if (json%valid_path("output_precision")) then
           call json_get(json, "output_precision", precision)
           if (precision == "double") then
              call field_writer_init_from_attributes(this, fields, filename, dp)
           else
              call field_writer_init_from_attributes(this, fields, filename, sp)
           end if
       else
           call field_writer_init_from_attributes(this, fields, filename)
       end if
    else
       call field_writer_init_from_attributes(this, fields)
    end if
  end subroutine field_writer_init_from_json

  !> Actual constructor.
  !! @param fields Array of field names to be sampled.
  !! @param filename The name of the file save the fields to. Optional, if not
  !! provided, fields are added to the main output file.
  !! @param precision The real precision of the output data. Optional, defaults
  !! to single precision.
  subroutine field_writer_init_from_attributes(this, fields, filename, precision)
    class(field_writer_t), intent(inout) :: this
    character(len=20), allocatable, intent(in) :: fields(:)
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision
    character(len=20) :: fieldi
    integer :: i

    ! Regsiter fields if they don't exist.
    do i=1, size(fields)
      fieldi = trim(fields(i))
      call neko_field_registry%add_field(this%case%fluid%dm_Xh, fieldi,&
                                         ignore_existing=.true.)
    end do

    if (present(filename)) then
       if (present(precision)) then
          call this%output%init(precision, filename, size(fields))
       else
          call this%output%init(sp, filename, size(fields))
       end if
       do i=1, size(fields)
          fieldi = trim(fields(i))
          call this%output%fields%assign(i, neko_field_registry%get_field(fieldi))
       end do

       call this%case%output_controller%add(this%output, &
            this%output_controller%control_value, &
            this%output_controller%control_mode)
    else
      do i=1, size(fields)
         fieldi = trim(fields(i))
         call this%case%f_out%fluid%append(neko_field_registry%get_field(fieldi))
      end do
    end if

  end subroutine field_writer_init_from_attributes

  !> Destructor.
  subroutine field_writer_free(this)
    class(field_writer_t), intent(inout) :: this
    call this%free_base()
  end subroutine field_writer_free

  !> Here to comply with the interface, does nothing.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine field_writer_compute(this, t, tstep)
    class(field_writer_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

  end subroutine field_writer_compute

end module field_writer
