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
!> Implements the `weak_grad_t` type.

module weak_grad
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : opgrad
  use case, only : case_t
  use fld_file_output, only : fld_file_output_t
  use json_utils, only : json_get, json_get_or_default
  use field_writer, only : field_writer_t
  implicit none
  private

  !> A simulation component that computes the weak gradient of a field.
  !! Wraps the `opgrad` operator.
  type, public, extends(simulation_component_t) :: weak_grad_t
     !> The scalar field to compute the weak gradient of.
     type(field_t), pointer :: u
     !> X weak grad component.
     type(field_t), pointer :: grad_x
     !> Y weak grad component.
     type(field_t), pointer :: grad_y
     !> Z weak grad component.
     type(field_t), pointer :: grad_z
     !> Output writer.
     type(field_writer_t) :: writer

   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => weak_grad_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
        weak_grad_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => weak_grad_free
     !> Compute the weak_grad field.
     procedure, pass(this) :: compute_ => weak_grad_compute
  end type weak_grad_t

contains

  !> Constructor from json.
  subroutine weak_grad_init_from_json(this, json, case)
    class(weak_grad_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: fieldname
    character(len=20) :: fields(3)

    ! Add fields keyword to the json so that the field_writer picks it up.
    ! Will also add fields to the registry.
    call json_get(json, "field", fieldname)

    fields(1) = "weak_grad_" // trim(fieldname) // "_x"
    fields(2) = "weak_grad_" // trim(fieldname) // "_y"
    fields(3) = "weak_grad_" // trim(fieldname) // "_z"

    call json%add("fields", fields)

    call this%init_base(json, case)
    call this%writer%init(json, case)

    call weak_grad_init_from_attributes(this, fieldname)
  end subroutine weak_grad_init_from_json

  !> Actual constructor.
  subroutine weak_grad_init_from_attributes(this, fieldname)
    class(weak_grad_t), intent(inout) :: this
    character(len=*) :: fieldname

    this%u => neko_field_registry%get_field_by_name(trim(fieldname))

    this%grad_x => neko_field_registry%get_field_by_name(&
                        "weak_grad_" // fieldname // "_x")
    this%grad_y => neko_field_registry%get_field_by_name(&
                        "weak_grad_" // fieldname // "_y")
    this%grad_z => neko_field_registry%get_field_by_name(&
                        "weak_grad_" // fieldname // "_z")


  end subroutine weak_grad_init_from_attributes

  !> Destructor.
  subroutine weak_grad_free(this)
    class(weak_grad_t), intent(inout) :: this
    call this%free_base()
    call this%writer%free()
    nullify(this%grad_x)
    nullify(this%grad_y)
    nullify(this%grad_z)
    nullify(this%u)
  end subroutine weak_grad_free

  !> Compute the weak_grad field.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine weak_grad_compute(this, t, tstep)
    class(weak_grad_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    call opgrad(this%grad_x%x, this%grad_y%x, this%grad_z%x, this%u%x,&
                this%case%fluid%c_Xh)
  end subroutine weak_grad_compute

end module weak_grad
