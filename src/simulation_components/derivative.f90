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
!> Implements the `derivative_t` type.

module derivative
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : dudxyz
  use case, only : case_t
  use fld_file_output, only : fld_file_output_t
  use json_utils, only : json_get, json_get_or_default
  use field_writer, only : field_writer_t
  use utils, only : neko_error
  implicit none
  private

  !> A simulation component that computes a derivative of a field.
  !! Wraps the `duxyz` operator.
  type, public, extends(simulation_component_t) :: derivative_t
     !> The scalar field to compute the weak gradient of.
     type(field_t), pointer :: u
     !> The derivative field
     type(field_t), pointer :: du
     !> Derivatives of r with respect to the direction of derivation.
     real(kind=rp), pointer :: dr(:,:,:,:)
     !> Derivatives of s with respect to the direction of derivation.
     real(kind=rp), pointer :: ds(:,:,:,:)
     !> Derivatives of t with respect to the direction of derivation.
     real(kind=rp), pointer :: dt(:,:,:,:)
     !> Output writer.
     type(field_writer_t) :: writer

   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => derivative_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
        derivative_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => derivative_free
     !> Compute the derivative field.
     procedure, pass(this) :: compute_ => derivative_compute
  end type derivative_t

contains

  !> Constructor from json.
  subroutine derivative_init_from_json(this, json, case)
    class(derivative_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: fieldname
    character(len=:), allocatable :: direction
    character(len=20) :: fields(1)

    ! Add fields keyword to the json so that the field_writer picks it up.
    ! Will also add fields to the registry.
    call json_get(json, "field", fieldname)
    call json_get(json, "direction", direction)

    fields(1) = "d" // trim(fieldname) // "_d" // direction
    call json%add("fields", fields)

    call this%init_base(json, case)
    call this%writer%init(json, case)

    call derivative_init_from_attributes(this, fieldname, direction)
  end subroutine derivative_init_from_json

  !> Actual constructor.
  subroutine derivative_init_from_attributes(this, fieldname, direction)
    class(derivative_t), intent(inout) :: this
    character(len=*) :: fieldname
    character(len=*) :: direction

    this%u => neko_field_registry%get_field_by_name(trim(fieldname))

    this%du => neko_field_registry%get_field_by_name(&
                        "d" // fieldname // "_d" // direction)

    if (direction .eq. "x") then
       this%dr => this%case%fluid%c_Xh%drdx
       this%ds => this%case%fluid%c_Xh%dsdx
       this%dt => this%case%fluid%c_Xh%dtdx
    else if (direction .eq. "y") then
       this%dr => this%case%fluid%c_Xh%drdy
       this%ds => this%case%fluid%c_Xh%dsdy
       this%dt => this%case%fluid%c_Xh%dtdy
    else if (direction .eq. "z") then
       this%dr => this%case%fluid%c_Xh%drdz
       this%ds => this%case%fluid%c_Xh%dsdz
       this%dt => this%case%fluid%c_Xh%dtdz
    else
        call neko_error("The direction of the derivative must be x, y or z")
    end if
  end subroutine derivative_init_from_attributes

  !> Destructor.
  subroutine derivative_free(this)
    class(derivative_t), intent(inout) :: this
    call this%free_base()
    call this%writer%free()
    nullify(this%du)
    nullify(this%u)
    nullify(this%dr)
    nullify(this%ds)
    nullify(this%dt)
  end subroutine derivative_free

  !> Compute the derivative field.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine derivative_compute(this, t, tstep)
    class(derivative_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    call dudxyz(this%du%x, this%u%x, this%dr, this%dr, this%dr,&
                this%case%fluid%c_Xh)
  end subroutine derivative_compute

end module derivative
