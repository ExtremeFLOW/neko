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
!> Implements the `les_simcomp_t` type.

module les_simcomp
  use num_types, only : rp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use case, only : case_t
  use les_model, only : les_model_t, les_model_factory
  use json_utils, only : json_get, json_get_or_default
  use field_writer, only : field_writer_t
  implicit none
  private

  !> A simulation component that drives the computation of the SGS
  !! viscosity.
   type, public, extends(simulation_component_t) :: les_simcomp_t
     !> The LES model.
     class(les_model_t), allocatable :: les_model
     !> Output writer.
     type(field_writer_t) :: writer
   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => les_simcomp_init_from_json
     !> Destructor.
     procedure, pass(this) :: free => les_simcomp_free
     !> Compute the les_simcomp field.
     procedure, pass(this) :: compute_ => les_simcomp_compute
  end type les_simcomp_t

contains

  !> Constructor from json.
  subroutine les_simcomp_init_from_json(this, json, case)
    class(les_simcomp_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: name
    character(len=:), allocatable :: nut_field
    character(len=20) :: fields(2)

    call this%free()

    ! Add fields keyword to the json so that the field_writer picks it up.
    ! Will also add fields to the registry if missing.
    call json_get_or_default(json, "nut_field", nut_field, "nut")
    fields(1) = "les_delta"
    fields(2) = nut_field

    call json%add("fields", fields)

    call this%init_base(json, case)
    call this%writer%init(json, case)

    call json_get(json, "model", name)

    call les_model_factory(this%les_model, name, case%fluid%dm_Xh,&
                           case%fluid%c_Xh, json)

  end subroutine les_simcomp_init_from_json

  !> Destructor.
  subroutine les_simcomp_free(this)
    class(les_simcomp_t), intent(inout) :: this
    call this%free_base()
    call this%writer%free()

    if (allocated(this%les_model)) then
      call this%les_model%free()
      deallocate(this%les_model)
    end if
  end subroutine les_simcomp_free

  !> Compute the les_simcomp field.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine les_simcomp_compute(this, t, tstep)
    class(les_simcomp_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    call this%les_model%compute(t, tstep)
  end subroutine les_simcomp_compute

end module les_simcomp
