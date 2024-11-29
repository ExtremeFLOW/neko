! Copyright (c) 2023, The Neko Authors
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
!> A simulation component that computes veldiv
!! The values are stored in the field registry under the name 'veldiv'

module veldiv
  use num_types, only : rp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : div
  use case, only : case_t
  use field_writer, only : field_writer_t
  use device
  implicit none
  private

  type, public, extends(simulation_component_t) :: veldiv_t
     !> X velocity component.
     type(field_t), pointer :: u
     !> Y velocity component.
     type(field_t), pointer :: v
     !> Z velocity component.
     type(field_t), pointer :: w

     !> X veldiv component.
     type(field_t), pointer :: veldiv

     !> Work arrays.
     type(field_t) :: temp1
     type(field_t) :: temp2

     !> Output writer.
     type(field_writer_t) :: writer

   contains
     !> Constructor from json.
     procedure, pass(this) :: init => veldiv_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
          veldiv_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => veldiv_free
     !> Compute the veldiv field
     procedure, pass(this) :: compute_ => veldiv_compute
  end type veldiv_t

contains

  !> Constructor from json.
  subroutine veldiv_init_from_json(this, json, case)
    class(veldiv_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target ::case
    character(len=20) :: fields(1)
    type(field_t), pointer :: u, v, w, veldiv

    ! Add fields keyword to the json so that the field_writer picks it up.
    ! Will also add fields to the registry.
    fields(1) = "veldiv"
    call json%add("fields", fields)

    call this%init_base(json, case)
    call this%writer%init(json, case)
    u => neko_field_registry%get_field("u")
    v => neko_field_registry%get_field("v")
    w => neko_field_registry%get_field("w")
    veldiv => neko_field_registry%get_field("veldiv")

    call veldiv_init_from_attributes(this, u, v, w, veldiv)
  end subroutine veldiv_init_from_json

  !> Actual constructor.
  subroutine veldiv_init_from_attributes(this, u, v, w, veldiv)
    class(veldiv_t), intent(inout) :: this
    type(field_t), pointer, intent(inout) :: u, v, w, veldiv

    this%u => u
    this%v => v
    this%w => w
    this%veldiv => veldiv

  end subroutine veldiv_init_from_attributes

  !> Destructor.
  subroutine veldiv_free(this)
    class(veldiv_t), intent(inout) :: this
    call this%free_base()
  end subroutine veldiv_free

  !> Compute the veldiv field.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine veldiv_compute(this, t, tstep)
    class(veldiv_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    call div(this%veldiv, this%u, this%v, this%w, this%case%fluid%c_Xh)

  end subroutine veldiv_compute

end module veldiv
