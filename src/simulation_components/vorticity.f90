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
!> A simulation component that computes the vorticity field.
!! The values are stored in the field registry under the names omega_x, omega_y
!! and omega_z.

module vorticity
  use num_types, only : rp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t, field_free, field_init
  use operators, only : curl
  use case, only : case_t
  implicit none
  private
  
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

     !> Work arrays.
     type(field_t) :: temp1
     type(field_t) :: temp2

     contains
      !> Constructor from json.
      procedure, pass(this) :: init => vorticity_init_from_json
      !> Actual constructor.
      procedure, pass(this) :: init_from_attributes => &
        vorticity_init_from_attributes
      !> Destructor.
      procedure, pass(this) :: free => vorticity_free
      !> Compute the vorticity field
     procedure, pass(this) :: compute => vorticity_compute
  end type vorticity_t
  
  contains
  
  !> Constructor from json.
  subroutine vorticity_init_from_json(this, json, case)
       class(vorticity_t), intent(inout) :: this
       type(json_file), intent(inout) :: json
       class(case_t), intent(inout), target ::case 
       
       call this%init_base(json, case)

       call vorticity_init_from_attributes(this)
  end subroutine vorticity_init_from_json

  !> Actual constructor.
  subroutine vorticity_init_from_attributes(this)
       class(vorticity_t), intent(inout) :: this

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

       call field_init(this%temp1, this%u%dof)
       call field_init(this%temp2, this%u%dof)
  end subroutine vorticity_init_from_attributes

  !> Destructor.
  subroutine vorticity_free(this)
       class(vorticity_t), intent(inout) :: this
       call field_free(this%temp1)
       call field_free(this%temp2)
  end subroutine vorticity_free

  !> Compute the vorticity field.
  !! @param t The time value.
  !! @param tstep The current time-step
  subroutine vorticity_compute(this, t, tstep)
       class(vorticity_t), intent(inout) :: this
       real(kind=rp), intent(in) :: t
       integer, intent(in) :: tstep
       !$omp parallel
       call curl(this%omega_x, this%omega_y, this%omega_z, this%u, this%v, &
                 this%w, this%temp1, this%temp2, this%case%fluid%c_Xh)
       !$omp end parallel
  end subroutine vorticity_compute
  
end module vorticity
