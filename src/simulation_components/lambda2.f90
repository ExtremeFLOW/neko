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
!> A simulation component that computes lambda2
!! The values are stored in the field registry under the name 'lambda2'

module lambda2
  use num_types, only : rp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : lambda2op
  use case, only : case_t
  use field_writer, only : field_writer_t
  use device
  implicit none
  private

  type, public, extends(simulation_component_t) :: lambda2_t
     !> X velocity component.
     type(field_t), pointer :: u
     !> Y velocity component.
     type(field_t), pointer :: v
     !> Z velocity component.
     type(field_t), pointer :: w

     !> X lambda2 component.
     type(field_t), pointer :: lambda2

     !> Work arrays.
     type(field_t) :: temp1
     type(field_t) :: temp2

     !> Output writer.
     type(field_writer_t) :: writer

   contains
     !> Constructor from json.
     procedure, pass(this) :: init => lambda2_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
          lambda2_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => lambda2_free
     !> Compute the lambda2 field
     procedure, pass(this) :: compute_ => lambda2_compute
  end type lambda2_t

contains

  !> Constructor from json.
  subroutine lambda2_init_from_json(this, json, case)
    class(lambda2_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target ::case
    character(len=20) :: fields(1)
    type(field_t), pointer :: u, v, w, lambda2

    ! Add fields keyword to the json so that the field_writer picks it up.
    ! Will also add fields to the registry.
    fields(1) = "lambda2"
    call json%add("fields", fields)

    call this%init_base(json, case)
    call this%writer%init(json, case)
    u => neko_field_registry%get_field("u")
    v => neko_field_registry%get_field("v")
    w => neko_field_registry%get_field("w")
    lambda2 => neko_field_registry%get_field("lambda2")

    call lambda2_init_from_attributes(this, u, v, w, lambda2)
  end subroutine lambda2_init_from_json

  !> Actual constructor.
  subroutine lambda2_init_from_attributes(this, u, v, w, lambda2)
    class(lambda2_t), intent(inout) :: this
    type(field_t), pointer, intent(inout) :: u, v, w, lambda2

    this%u => u
    this%v => v
    this%w => w
    this%lambda2 => lambda2

  end subroutine lambda2_init_from_attributes

  !> Destructor.
  subroutine lambda2_free(this)
    class(lambda2_t), intent(inout) :: this
    call this%free_base()
  end subroutine lambda2_free

  !> Compute the lambda2 field.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine lambda2_compute(this, t, tstep)
    class(lambda2_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    call lambda2op(this%lambda2, this%u, this%v, this%w, this%case%fluid%c_Xh)

  end subroutine lambda2_compute

end module lambda2
