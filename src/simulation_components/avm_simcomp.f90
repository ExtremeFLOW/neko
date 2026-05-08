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
!
!> Implements the `avm_simcomp_t` type.

module avm_simcomp
  use num_types, only : rp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use case, only : case_t
  use time_state, only : time_state_t
  use avm_model, only : avm_model_t, avm_model_factory
  use json_utils, only : json_get, json_get_or_default
  use field_writer, only : field_writer_t
  use utils, only : neko_error, NEKO_VARNAME_LEN
  implicit none
  private

  !> A simulation component that drives the computation of the artificial
  !! viscosity method.
  type, public, extends(simulation_component_t) :: avm_simcomp_t
     !> The avm model.
     class(avm_model_t), allocatable :: avm_model
     !> Output writer.
     type(field_writer_t) :: writer
   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => avm_simcomp_init_from_json
     !> Destructor.
     procedure, pass(this) :: free => avm_simcomp_free
     !> Compute the avm_simcomp field before the time step.
     procedure, pass(this) :: preprocess_ => avm_simcomp_preprocess
     !> Compute the avm_simcomp field.
     procedure, pass(this) :: compute_ => avm_simcomp_compute
     !> Compute the avm_simcomp field when restart.
     procedure, pass(this) :: restart_ => avm_simcomp_restart
  end type avm_simcomp_t

contains

  !> Constructor from json.
  subroutine avm_simcomp_init_from_json(this, json, case)
    class(avm_simcomp_t), intent(inout), target :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: name
    character(len=:), allocatable :: model_name
    character(len=:), allocatable :: reg_coeff
    character(len=NEKO_VARNAME_LEN) :: fields(1)

    call this%free()

    ! Add fields keyword to the json so that the field_writer picks it up.
    ! Will also add fields to the registry if missing.
    call json_get_or_default(json, "name", name, "avm_model")
    call json_get_or_default(json, "reg_coeff", reg_coeff, "reg_coeff")
    fields(1) = reg_coeff

    call json%add("fields", fields)

    call this%init_base(json, case)
    call this%writer%init(json, case)

    call json_get(json, "model", model_name)
    this%name = name

    call avm_model_factory(this%avm_model, model_name, case, json)

    if (allocated(name)) deallocate(name)
    if (allocated(model_name)) deallocate(model_name)
    if (allocated(reg_coeff)) deallocate(reg_coeff)
  end subroutine avm_simcomp_init_from_json

  !> Destructor.
  subroutine avm_simcomp_free(this)
    class(avm_simcomp_t), intent(inout) :: this
    call this%free_base()
    call this%writer%free()

    if (allocated(this%avm_model)) then
       call this%avm_model%free()
       deallocate(this%avm_model)
    end if
  end subroutine avm_simcomp_free

  !> Compute the avm_simcomp field.
  !! @param time The current time info
  subroutine avm_simcomp_preprocess(this, time)
    class(avm_simcomp_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    call this%avm_model%preprocess(time)
  end subroutine avm_simcomp_preprocess

  !> Compute the avm_simcomp field.
  !! @param time The current time info
  subroutine avm_simcomp_compute(this, time)
    class(avm_simcomp_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    call this%avm_model%compute(time)
  end subroutine avm_simcomp_compute

  !> Compute the avm_simcomp field when restart.
  !! @param time The current time info
  subroutine avm_simcomp_restart(this, time)
    class(avm_simcomp_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    call this%avm_model%compute(time)
  end subroutine avm_simcomp_restart

end module avm_simcomp
