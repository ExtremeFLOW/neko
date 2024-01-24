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
!> Implements the `vorticity_t` type.

module vorticity
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : curl
  use case, only : case_t
  use fld_output, only : fld_output_t
  use json_utils, only : json_get, json_get_or_default
  implicit none
  private

  !> A simulation component that computes the vorticity field.
  !! Added to the field registry as `omega_x`, `omega_y``, and `omega_z`.
  type, public, extends(simulation_component_t) :: vorticity_t
     !> X velocity component.
     type(field_t), pointer :: u
     !> Y velocity component.
     type(field_t), pointer :: v
     !> Z velocity component.
     type(field_t), pointer :: w

     !> X vorticity component.
     type(field_t), pointer :: omega_x
     !> Y vorticity component.
     type(field_t), pointer :: omega_y
     !> Z vorticity component.
     type(field_t), pointer :: omega_z

     !> Work array.
     type(field_t) :: temp1
     !> Work array.
     type(field_t) :: temp2

     !> Output writer.
     type(fld_output_t), private :: output

   contains
     !> Constructor from json, wrapping the actual constructor.
     procedure, pass(this) :: init => vorticity_init_from_json
     !> Actual constructor.
     procedure, pass(this) :: init_from_attributes => &
        vorticity_init_from_attributes
     !> Destructor.
     procedure, pass(this) :: free => vorticity_free
     !> Compute the vorticity field.
     procedure, pass(this) :: compute_ => vorticity_compute
  end type vorticity_t

contains

  !> Constructor from json.
  subroutine vorticity_init_from_json(this, json, case)
    class(vorticity_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(case_t), intent(inout), target :: case
    character(len=:), allocatable :: filename
    character(len=:), allocatable :: precision

    call this%init_base(json, case)

    if (json%valid_path("output_filename")) then
       call json_get(json, "output_filename", filename)
       if (json%valid_path("output_precision")) then
           call json_get(json, "output_precision", precision)
           if (precision == "double") then
              call vorticity_init_from_attributes(this, filename, dp)
           else
              call vorticity_init_from_attributes(this, filename, sp)
           end if
       else
           call vorticity_init_from_attributes(this, filename)
       end if
    else
       call vorticity_init_from_attributes(this)
    end if
  end subroutine vorticity_init_from_json

  !> Actual constructor.
  subroutine vorticity_init_from_attributes(this, filename, precision)
    class(vorticity_t), intent(inout) :: this
    character(len=*), intent(in), optional :: filename
    integer, intent(in), optional :: precision
    type(fld_output_t) :: output

    this%u => neko_field_registry%get_field_by_name("u")
    this%v => neko_field_registry%get_field_by_name("v")
    this%w => neko_field_registry%get_field_by_name("w")

    if (.not. neko_field_registry%field_exists("omega_x")) then
       call neko_field_registry%add_field(this%u%dof, "omega_x")
    end if
    if (.not. neko_field_registry%field_exists("omega_y")) then
       call neko_field_registry%add_field(this%u%dof, "omega_y")
    end if
    if (.not. neko_field_registry%field_exists("omega_z")) then
       call neko_field_registry%add_field(this%u%dof, "omega_z")
    end if
    this%omega_x => neko_field_registry%get_field_by_name("omega_x")
    this%omega_y => neko_field_registry%get_field_by_name("omega_y")
    this%omega_z => neko_field_registry%get_field_by_name("omega_z")

    call this%temp1%init(this%u%dof)
    call this%temp2%init(this%u%dof)


    if (present(filename)) then
       if (present(precision)) then
          call this%output%init(precision, filename, 3)
       else
          call this%output%init(sp, filename, 3)
       end if
       this%output%fields%fields(1)%f => this%omega_x
       this%output%fields%fields(2)%f => this%omega_y
       this%output%fields%fields(3)%f => this%omega_z
       call this%case%s%add(this%output, this%output_controller%control_value, &
                            this%output_controller%control_mode)
    else
       call this%case%f_out%fluid%append(this%omega_x)
       call this%case%f_out%fluid%append(this%omega_y)
       call this%case%f_out%fluid%append(this%omega_z)
    end if

  end subroutine vorticity_init_from_attributes

  !> Destructor.
  subroutine vorticity_free(this)
    class(vorticity_t), intent(inout) :: this
    call this%free_base()
    call this%temp1%free()
    call this%temp2%free()
  end subroutine vorticity_free

  !> Compute the vorticity field.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine vorticity_compute(this, t, tstep)
    class(vorticity_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep

    call curl(this%omega_x, this%omega_y, this%omega_z, this%u, this%v, &
                 this%w, this%temp1, this%temp2, this%case%fluid%c_Xh)
  end subroutine vorticity_compute

end module vorticity
