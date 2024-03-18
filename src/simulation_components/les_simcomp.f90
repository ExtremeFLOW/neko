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
  use num_types, only : rp, dp, sp
  use json_module, only : json_file
  use simulation_component, only : simulation_component_t
  use field_registry, only : neko_field_registry
  use field, only : field_t
  use operators, only : curl
  use case, only : case_t
  use les_model, only : les_model_t
  use les_model_fctry, only : les_model_factory
  use json_utils, only : json_get
  use fluid_pnpn_stress, only : fluid_pnpn_stress_t
  use math, only : addcol3, add2s2, cfill
  use device_math, only : device_addcol3, device_add2s2, device_cfill
  use neko_config, only : NEKO_BCKND_DEVICE
  use fld_file_output, only : fld_file_output_t
  implicit none
  private

  !> A simulation component that drives the computation of the SGS
  !! viscosity.
   type, public, extends(simulation_component_t) :: les_simcomp_t
     !> The LES model.
     class(les_model_t), allocatable :: les_model
     !> Density field
     type(field_t), pointer :: rho => null()
     !> Variable dynamic viscosity fields.
     type(field_t), pointer :: mu_field => null()
     !> Molecular dynamic visocity.
     real(kind=rp) :: mu = 0.0
     !> Output writer.
     type(fld_file_output_t), private :: output
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
    character(len=:), allocatable :: filename
    character(len=:), allocatable :: precision

    call this%free()

    call this%free()

    call json_get(json, "model", name)

    call this%init_base(json, case)

    call les_model_factory(this%les_model, name, case%fluid%dm_Xh,&
                           case%fluid%c_Xh, json)

    select type(f => case%fluid)
    type is(fluid_pnpn_stress_t)
      this%rho => f%rho_field
      this%mu_field => f%mu_field
      this%mu = f%mu
    end select

    ! Configure output for delta and nut
    if (json%valid_path("output_filename")) then
       call json_get(json, "output_filename", filename)
       if (json%valid_path("output_precision")) then
           call json_get(json, "output_precision", precision)
           if (precision == "double") then
              call this%output%init(dp, filename, 2)
           end if
       else
           call this%output%init(sp, filename, 2)
       end if
       this%output%fields%fields(1)%f => this%les_model%nut
       this%output%fields%fields(2)%f => this%les_model%delta
       call this%case%s%add(this%output, this%output_controller%control_value, &
                            this%output_controller%control_mode)
    else
       call this%case%f_out%fluid%append(this%les_model%nut)
       call this%case%f_out%fluid%append(this%les_model%delta)
    end if

  end subroutine les_simcomp_init_from_json

  !> Destructor.
  subroutine les_simcomp_free(this)
    class(les_simcomp_t), intent(inout) :: this
    call this%free_base()
    nullify(this%rho)
    nullify(this%mu_field)

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
    this%mu_field = this%mu
    if (NEKO_BCKND_DEVICE .eq. 1) then
        call device_addcol3(this%mu_field%x_d, this%rho%x_d, &
                            this%les_model%nut%x_d, this%rho%dof%size())
    else
        call addcol3(this%mu_field%x, this%rho%x, this%les_model%nut%x, &
                     this%rho%dof%size())
    end if

    if (allocated(this%case%scalar)) then
      associate (s => this%case%scalar)
        if (NEKO_BCKND_DEVICE .eq. 1) then
           call device_cfill(s%lambda_field%x_d, s%lambda, s%dm_Xh%size())
           call device_add2s2(s%lambda_field%x_d, this%les_model%nut%x_d, &
                              s%rho*s%cp, s%dm_Xh%size())
        else
           call cfill(s%lambda_field%x, s%lambda, s%dm_Xh%size())
           call add2s2(s%lambda_field%x, this%les_model%nut%x, s%rho*s%cp, &
                       s%dm_Xh%size())
        end if
      end associate
    end if

  end subroutine les_simcomp_compute

end module les_simcomp
