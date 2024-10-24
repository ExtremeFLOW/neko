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
!> Implements the `boussinesq_source_term_t` type.
module boussinesq_source_term
  use num_types, only : rp
  use field_list, only : field_list_t
  use field, only : field_t
  use json_module, only : json_file
  use json_utils, only: json_get, json_get_or_default
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  use boussinesq_source_term_cpu, only : boussinesq_source_term_compute_cpu
  use boussinesq_source_term_device, only : &
       boussinesq_source_term_compute_device
  use field_registry, only : neko_field_registry
  implicit none
  private

  !> Bouyancy source term accroding to the Boussinesq approximation.
  !! @details Called "boussinesq" in the JSON.
  !! Controlled by the following parameters:
  !! - "scalar_field": The name of the scalar that drives the source term,
  !!   defaults to "s".
  !! - "ref_value": The reference value of the scalar.
  !! - "g": The gravity vector.
  !! - "beta": The the thermal expansion coefficeint, defaults to `1/ref_value`.
  type, public, extends(source_term_t) :: boussinesq_source_term_t
     !> The scalar field used to drive the source term, typically temperature.
     type(field_t), pointer :: s => null()
     !> The reference value of the scalar field.
     real(kind=rp) :: ref_value = 0
     !> Gravity vector
     real(kind=rp) :: g(3)
     !> Thermal expantion coefficient
     real(kind=rp) :: beta
   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => boussinesq_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_compenents => &
       boussinesq_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => boussinesq_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => boussinesq_source_term_compute
  end type boussinesq_source_term_t

contains
  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  subroutine boussinesq_source_term_init_from_json(this, json, fields, coef)
    class(boussinesq_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(inout), target :: fields
    type(coef_t), intent(inout), target :: coef
    real(kind=rp), allocatable :: values(:)
    real(kind=rp) :: start_time, end_time, ref_value
    character(len=:), allocatable :: scalar_name
    real(kind=rp), allocatable :: g(:)
    real(kind=rp) :: beta

    if (.not. fields%size() == 3) then
       call neko_error("Boussinesq term expects 3 fields to work on.")
    end if

    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))

    call json_get_or_default(json, "scalar_field", scalar_name, "s")
    call json_get(json, "g", g)

    if (.not. size(g) == 3) then
       call neko_error("The gravity vector should have 3 components")
    end if

    call json_get(json, "reference_value", ref_value)
    call json_get_or_default(json, "beta", beta, 1.0_rp/ref_value)

    call boussinesq_source_term_init_from_components(this, fields, scalar_name,&
       ref_value, g, beta, coef, start_time, end_time)

  end subroutine boussinesq_source_term_init_from_json

  !> The constructor from type components.
  !! @param fields A list of fields for adding the source values.
  !! @param scalar_name The name of the scalar field driving the source.
  !! @param ref_value The reference value of the scalar field.
  !! @param g The gravity vector.
  !! @param beta The thermal expansion coefficient.
  !! @param coef The SEM coeffs.
  !! @param start_time When to start adding the source term.
  !! @param end_time When to stop adding the source term.
  subroutine boussinesq_source_term_init_from_components(this, fields, &
    scalar_name, ref_value, g, beta, coef, start_time, end_time)
    class(boussinesq_source_term_t), intent(inout) :: this
    class(field_list_t), intent(inout), target :: fields
    character(len=*), intent(in) :: scalar_name
    real(kind=rp), intent(in) :: ref_value
    real(kind=rp), intent(in) :: g(3)
    real(kind=rp), intent(in) :: beta
    type(coef_t) :: coef
    real(kind=rp), intent(in) :: start_time
    real(kind=rp), intent(in) :: end_time

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)

    if (.not. neko_field_registry%field_exists(scalar_name)) then
       call neko_field_registry%add_field(this%fields%dof(1), "s")
    end if
    this%s => neko_field_registry%get_field("s")

    this%ref_value = ref_value
    this%g = g
    this%beta = beta
  end subroutine boussinesq_source_term_init_from_components

  !> Destructor.
  subroutine boussinesq_source_term_free(this)
    class(boussinesq_source_term_t), intent(inout) :: this

    call this%free_base()
    nullify(this%s)
  end subroutine boussinesq_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine boussinesq_source_term_compute(this, t, tstep)
    class(boussinesq_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: n_fields, i, n

    n_fields = this%fields%size()
    n = this%fields%item_size(1)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call boussinesq_source_term_compute_device(this%fields, this%s,&
         this%ref_value, this%g, this%beta)
    else
       call boussinesq_source_term_compute_cpu(this%fields, this%s,&
         this%ref_value, this%g, this%beta)
    end if
  end subroutine boussinesq_source_term_compute

end module boussinesq_source_term
