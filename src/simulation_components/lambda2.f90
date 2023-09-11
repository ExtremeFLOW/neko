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
  use field, only : field_t, field_free, field_init, field_ptr_t
  use operators, only : lambda2op
  use case, only : case_t
  use neko_config
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

    call this%init_base(json, case)

    call lambda2_init_from_attributes(this)
  end subroutine lambda2_init_from_json

  !> Actual constructor.
  subroutine lambda2_init_from_attributes(this)
    class(lambda2_t), intent(inout) :: this

    this%u => neko_field_registry%get_field_by_name("u")
    this%v => neko_field_registry%get_field_by_name("v")
    this%w => neko_field_registry%get_field_by_name("w")

    if (.not. neko_field_registry%field_exists("lambda2")) then
       call neko_field_registry%add_field(this%u%dof, "lambda2")
    end if
    this%lambda2 => neko_field_registry%get_field_by_name("lambda2")
    !Store lambda2 in the next free field in the fluid output.
    !If running without scalar this means the temperature field.
    call this%case%f_out%fluid%append(this%lambda2)
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
    !$omp parallel if((NEKO_BCKND_DEVICE .eq. 0) .and. (NEKO_BCKND_SX .eq. 0))
    call lambda2op(this%lambda2, this%u, this%v, this%w, this%case%fluid%c_Xh)
    !$omp end parallel
  end subroutine lambda2_compute

end module lambda2
