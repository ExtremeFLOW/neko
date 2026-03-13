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
!> Implements the `phase_field_sharpening_source_term_t` type.
!!
!! @details
!! The Allen-Cahn reinitialization forcing applied to the phase-field (scalar)
!! equation in a diffuse-interface two-phase flow model.
!! Called "phase_field_sharpening" in the JSON.
!!
!! Computes:
!! \f$ f = \gamma u_\text{max} \nabla \cdot \left( -\phi(1-\phi)
!!         \frac{\nabla\phi}{|\nabla\phi|} \right) \f$
!!
!! Controlled by the following parameters:
!! - "scalar_field": Name of the phase field, defaults to "phase".
!! - "gamma": Sharpening coefficient.
!! - "u_max": Reference velocity scale.
module phase_field_sharpening_source_term
  use num_types, only : rp
  use field, only : field_t
  use field_list, only : field_list_t
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use operators, only : grad, div
  use gather_scatter, only : GS_OP_ADD
  use scratch_registry, only : neko_scratch_registry
  use field_registry, only : neko_field_registry
  use math, only : col2, add2, cmult
  use utils, only : neko_error
  use time_state, only : time_state_t
  implicit none
  private

  !> Phase-field sharpening (Allen-Cahn reinitialization) source term.
  !! @details Called "phase_field_sharpening" in the JSON.
  !! Controlled by the following parameters:
  !! - "scalar_field": The name of the phase field, defaults to "phase".
  !! - "gamma": Sharpening coefficient.
  !! - "u_max": Reference velocity scale.
  type, public, extends(source_term_t) :: phase_field_sharpening_source_term_t
     !> The phase field.
     type(field_t), pointer :: phi => null()
     !> Sharpening coefficient.
     real(kind=rp) :: gamma
     !> Reference velocity scale.
     real(kind=rp) :: u_max
   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => &
          phase_field_sharpening_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_components => &
          phase_field_sharpening_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => phase_field_sharpening_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => &
          phase_field_sharpening_source_term_compute
  end type phase_field_sharpening_source_term_t

contains

  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param variable_name The name of the variable where the source term acts.
  subroutine phase_field_sharpening_source_term_init_from_json(this, json, &
       fields, coef, variable_name)
    class(phase_field_sharpening_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name
    real(kind=rp) :: start_time, end_time, gamma, u_max
    character(len=:), allocatable :: scalar_name

    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))
    call json_get_or_default(json, "scalar_field", scalar_name, "phase")
    call json_get(json, "gamma", gamma)
    call json_get(json, "u_max", u_max)

    call phase_field_sharpening_source_term_init_from_components(this, fields, &
         scalar_name, gamma, u_max, coef, start_time, end_time)

  end subroutine phase_field_sharpening_source_term_init_from_json

  !> The constructor from type components.
  !! @param fields A list of fields for adding the source values.
  !! @param scalar_name The name of the phase field.
  !! @param gamma Sharpening coefficient.
  !! @param u_max Reference velocity scale.
  !! @param coef The SEM coeffs.
  !! @param start_time When to start adding the source term.
  !! @param end_time When to stop adding the source term.
  subroutine phase_field_sharpening_source_term_init_from_components(this, &
       fields, scalar_name, gamma, u_max, coef, start_time, end_time)
    class(phase_field_sharpening_source_term_t), intent(inout) :: this
    type(field_list_t), intent(in), target :: fields
    character(len=*), intent(in) :: scalar_name
    real(kind=rp), intent(in) :: gamma, u_max
    type(coef_t), intent(in), target :: coef
    real(kind=rp), intent(in) :: start_time, end_time

    if (.not. fields%size() == 1) then
       call neko_error("Phase-field sharpening term expects 1 field to work on.")
    end if

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)

    if (.not. neko_field_registry%field_exists(scalar_name)) then
       call neko_field_registry%add_field(this%fields%dof(1), scalar_name)
    end if
    this%phi => neko_field_registry%get_field(scalar_name)

    this%gamma = gamma
    this%u_max = u_max

  end subroutine phase_field_sharpening_source_term_init_from_components

  !> Destructor.
  subroutine phase_field_sharpening_source_term_free(this)
    class(phase_field_sharpening_source_term_t), intent(inout) :: this

    nullify(this%phi)
    call this%free_base()
  end subroutine phase_field_sharpening_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param time The time state.
  subroutine phase_field_sharpening_source_term_compute(this, time)
    class(phase_field_sharpening_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: rhs_s
    type(field_t), pointer :: work1, work2, work3, work4
    integer :: ind(4), i, n
    real(kind=rp) :: absgrad

    rhs_s => this%fields%get(1)
    n = rhs_s%size()

    call neko_scratch_registry%request_field(work1, ind(1))
    call neko_scratch_registry%request_field(work2, ind(2))
    call neko_scratch_registry%request_field(work3, ind(3))
    call neko_scratch_registry%request_field(work4, ind(4))

    ! Compute gradient of phase field
    call grad(work1%x, work2%x, work3%x, this%phi%x, this%coef)

    ! Apply gather-scatter and multiplicity for continuity across elements
    call this%coef%gs_h%op(work1, GS_OP_ADD)
    call this%coef%gs_h%op(work2, GS_OP_ADD)
    call this%coef%gs_h%op(work3, GS_OP_ADD)
    call col2(work1%x, this%coef%mult, n)
    call col2(work2%x, this%coef%mult, n)
    call col2(work3%x, this%coef%mult, n)

    ! Compute F = -phi*(1-phi)*grad(phi)/|grad(phi)|
    do i = 1, n
       absgrad = sqrt(work1%x(i,1,1,1)**2 + work2%x(i,1,1,1)**2 &
            + work3%x(i,1,1,1)**2)
       if (absgrad < 1.0e-12_rp) then
          work1%x(i,1,1,1) = 0.0_rp
          work2%x(i,1,1,1) = 0.0_rp
          work3%x(i,1,1,1) = 0.0_rp
       else
          work1%x(i,1,1,1) = &
               -this%phi%x(i,1,1,1)*(1.0_rp - this%phi%x(i,1,1,1)) &
               * work1%x(i,1,1,1) / absgrad
          work2%x(i,1,1,1) = &
               -this%phi%x(i,1,1,1)*(1.0_rp - this%phi%x(i,1,1,1)) &
               * work2%x(i,1,1,1) / absgrad
          work3%x(i,1,1,1) = &
               -this%phi%x(i,1,1,1)*(1.0_rp - this%phi%x(i,1,1,1)) &
               * work3%x(i,1,1,1) / absgrad
       end if
    end do

    ! Compute divergence of F, scale by gamma*u_max, then add to RHS
    call div(work4%x, work1%x, work2%x, work3%x, this%coef)
    call cmult(work4%x, this%gamma * this%u_max, n)
    call add2(rhs_s%x, work4%x, n)

    call neko_scratch_registry%relinquish_field(ind)

  end subroutine phase_field_sharpening_source_term_compute

end module phase_field_sharpening_source_term
