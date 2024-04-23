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
!> Implements the `power_iterations_t` type.

module power_iterations
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : curl
  use case, only : case_t
  use fld_file_output, only : fld_file_output_t
  use json_utils, only : json_get, json_get_or_default
  use field_writer, only : field_writer_t
  use math, only: add3, copy
  implicit none
  private

  !> A simulation component that computes the power_iterations field.
  !! Added to the field registry as `omega_x`, `omega_y``, and `omega_z`.
  type, public, extends(simulation_component_t) :: power_iterations_t
     !> X velocity component.
     type(field_t), pointer :: u
     !> Y velocity component.
     type(field_t), pointer :: v
     !> Z velocity component.
     type(field_t), pointer :: w

     !> X power_iterations component.
     type(field_t), pointer :: u_b
     !> Y power_iterations component.
     type(field_t), pointer :: v_b
     !> Z power_iterations component.
     type(field_t), pointer :: w_b

     type(field_t), pointer :: u_full
     type(field_t), pointer :: v_full
     type(field_t), pointer :: w_full

   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => power_iterations_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
       power_iterations_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => power_iterations_free
     !> Compute the power_iterations field.
     procedure, pass(this) :: compute_ => power_iterations_compute
  end type power_iterations_t

contains

  !> Constructor from json.
  subroutine power_iterations_init_from_json(this, json, case)
    class(power_iterations_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: filename
    character(len=:), allocatable :: precision
    character(len=20) :: fields(3)


    call neko_field_registry%add_field(case%fluid%dm_Xh, "ub",&
                                       ignore_existing=.true.)
    call neko_field_registry%add_field(case%fluid%dm_Xh, "vb",&
                                       ignore_existing=.true.)
    call neko_field_registry%add_field(case%fluid%dm_Xh, "wb",&
                                       ignore_existing=.true.)

    call neko_field_registry%add_field(case%fluid%dm_Xh, "u_full",&
                                       ignore_existing=.true.)
    call neko_field_registry%add_field(case%fluid%dm_Xh, "v_full",&
                                       ignore_existing=.true.)
    call neko_field_registry%add_field(case%fluid%dm_Xh, "w_full",&
                                       ignore_existing=.true.)

    call power_iterations_init_from_attributes(this)
  end subroutine power_iterations_init_from_json

  !> Actual constructor.
  subroutine power_iterations_init_from_attributes(this)
    class(power_iterations_t), intent(inout) :: this

    this%u => neko_field_registry%get_field("u")
    this%v => neko_field_registry%get_field("v")
    this%w => neko_field_registry%get_field("w")

    this%u_b => neko_field_registry%get_field("ub")
    this%v_b => neko_field_registry%get_field("vb")
    this%w_b => neko_field_registry%get_field("wb")

    this%u_full => neko_field_registry%get_field("u_full")
    this%v_full => neko_field_registry%get_field("v_full")
    this%w_full => neko_field_registry%get_field("w_full")

  end subroutine power_iterations_init_from_attributes

  !> Destructor.
  subroutine power_iterations_free(this)
    class(power_iterations_t), intent(inout) :: this
    call this%free_base()
  end subroutine power_iterations_free

  !> Compute the power_iterations field.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine power_iterations_compute(this, t, tstep)
    class(power_iterations_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    integer :: n

    n = this%u%dof%size()

    call add3(this%u_full%x, this%u_b%x, this%u%x, n)
    call add3(this%v_full%x, this%v_b%x, this%v%x, n)
    call add3(this%w_full%x, this%w_b%x, this%w%x, n)

  end subroutine power_iterations_compute

end module power_iterations
