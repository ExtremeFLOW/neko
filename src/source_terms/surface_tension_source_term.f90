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
!> Implements the `surface_tension_source_term_t` type.
!!
!! @details
!! Continuum Surface Force (CSF) model for surface tension in a
!! phase-field two-phase flow simulation.
!! Called "surface_tension" in the JSON.
!!
!! Computes:
!! \f$ \mathbf{F} = \sigma \kappa \nabla\phi \f$
!! where \f$\kappa = \nabla \cdot (\nabla\phi / |\nabla\phi|)\f$ is the
!! interface curvature and \f$\phi\f$ is the phase field.
!!
!! Controlled by the following parameters:
!! - "phase_field": The name of the phase field, defaults to "phase".
!! - "sigma": Surface tension coefficient.
module surface_tension_source_term
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
  use math, only : col2, add2, copy, cmult
  use utils, only : neko_error
  use time_state, only : time_state_t
  implicit none
  private

  !> Surface tension source term based on the Continuum Surface Force (CSF)
  !! model for phase-field two-phase flow.
  !! @details Called "surface_tension" in the JSON.
  !! Controlled by the following parameters:
  !! - "phase_field": The name of the phase field, defaults to "phase".
  !! - "sigma": Surface tension coefficient.
  type, public, extends(source_term_t) :: surface_tension_source_term_t
     !> The phase field.
     type(field_t), pointer :: phi => null()
     !> Surface tension coefficient.
     real(kind=rp) :: sigma
   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => surface_tension_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_components => &
          surface_tension_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => surface_tension_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => surface_tension_source_term_compute
  end type surface_tension_source_term_t

contains

  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param variable_name The name of the variable where the source term acts.
  subroutine surface_tension_source_term_init_from_json(this, json, fields, &
       coef, variable_name)
    class(surface_tension_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name
    real(kind=rp) :: start_time, end_time, sigma
    character(len=:), allocatable :: phase_name

    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))
    call json_get_or_default(json, "phase_field", phase_name, "phase")
    call json_get(json, "sigma", sigma)

    call surface_tension_source_term_init_from_components(this, fields, &
         phase_name, sigma, coef, start_time, end_time)

  end subroutine surface_tension_source_term_init_from_json

  !> The constructor from type components.
  !! @param fields A list of fields for adding the source values.
  !! @param phase_name The name of the phase field.
  !! @param sigma Surface tension coefficient.
  !! @param coef The SEM coeffs.
  !! @param start_time When to start adding the source term.
  !! @param end_time When to stop adding the source term.
  subroutine surface_tension_source_term_init_from_components(this, fields, &
       phase_name, sigma, coef, start_time, end_time)
    class(surface_tension_source_term_t), intent(inout) :: this
    type(field_list_t), intent(in), target :: fields
    character(len=*), intent(in) :: phase_name
    real(kind=rp), intent(in) :: sigma
    type(coef_t), intent(in), target :: coef
    real(kind=rp), intent(in) :: start_time, end_time

    if (.not. fields%size() == 3) then
       call neko_error("Surface tension term expects 3 fields to work on.")
    end if

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)

    if (.not. neko_field_registry%field_exists(phase_name)) then
       call neko_field_registry%add_field(this%fields%dof(1), phase_name)
    end if
    this%phi => neko_field_registry%get_field(phase_name)

    this%sigma = sigma

  end subroutine surface_tension_source_term_init_from_components

  !> Destructor.
  subroutine surface_tension_source_term_free(this)
    class(surface_tension_source_term_t), intent(inout) :: this

    nullify(this%phi)
    call this%free_base()
  end subroutine surface_tension_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param time The time state.
  subroutine surface_tension_source_term_compute(this, time)
    class(surface_tension_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: rhs_u, rhs_v, rhs_w
    type(field_t), pointer :: dphi_x, dphi_y, dphi_z, kappa, nx, ny, nz
    integer :: ind(7), i, n
    real(kind=rp) :: absgrad

    rhs_u => this%fields%get(1)
    rhs_v => this%fields%get(2)
    rhs_w => this%fields%get(3)
    n = rhs_u%size()

    call neko_scratch_registry%request_field(dphi_x, ind(1))
    call neko_scratch_registry%request_field(dphi_y, ind(2))
    call neko_scratch_registry%request_field(dphi_z, ind(3))
    call neko_scratch_registry%request_field(kappa,  ind(4))
    call neko_scratch_registry%request_field(nx,     ind(5))
    call neko_scratch_registry%request_field(ny,     ind(6))
    call neko_scratch_registry%request_field(nz,     ind(7))

    ! Compute gradient of phase field: grad(phi)
    call grad(dphi_x%x, dphi_y%x, dphi_z%x, this%phi%x, this%coef)

    ! Apply gather-scatter and multiplicity for continuity across elements
    call this%coef%gs_h%op(dphi_x, GS_OP_ADD)
    call this%coef%gs_h%op(dphi_y, GS_OP_ADD)
    call this%coef%gs_h%op(dphi_z, GS_OP_ADD)
    call col2(dphi_x%x, this%coef%mult, n)
    call col2(dphi_y%x, this%coef%mult, n)
    call col2(dphi_z%x, this%coef%mult, n)

    ! Normalize grad(phi) to get interface normal n = grad(phi)/|grad(phi)|
    ! Store grad(phi) in nx/ny/nz (will overwrite with the normalized version)
    call copy(nx%x, dphi_x%x, n)
    call copy(ny%x, dphi_y%x, n)
    call copy(nz%x, dphi_z%x, n)

    do i = 1, n
       absgrad = sqrt(nx%x(i,1,1,1)**2 + ny%x(i,1,1,1)**2 &
            + nz%x(i,1,1,1)**2)
       if (absgrad < 1.0e-12_rp) then
          nx%x(i,1,1,1) = 0.0_rp
          ny%x(i,1,1,1) = 0.0_rp
          nz%x(i,1,1,1) = 0.0_rp
       else
          nx%x(i,1,1,1) = nx%x(i,1,1,1) / absgrad
          ny%x(i,1,1,1) = ny%x(i,1,1,1) / absgrad
          nz%x(i,1,1,1) = nz%x(i,1,1,1) / absgrad
       end if
    end do

    ! Compute curvature: kappa = div(n)
    call div(kappa%x, nx%x, ny%x, nz%x, this%coef)

    ! Add surface tension force F = sigma * kappa * grad(phi) to momentum RHS
    do i = 1, n
       dphi_x%x(i,1,1,1) = this%sigma * kappa%x(i,1,1,1) * dphi_x%x(i,1,1,1)
       dphi_y%x(i,1,1,1) = this%sigma * kappa%x(i,1,1,1) * dphi_y%x(i,1,1,1)
       dphi_z%x(i,1,1,1) = this%sigma * kappa%x(i,1,1,1) * dphi_z%x(i,1,1,1)
    end do

    call add2(rhs_u%x, dphi_x%x, n)
    call add2(rhs_v%x, dphi_y%x, n)
    call add2(rhs_w%x, dphi_z%x, n)

    call neko_scratch_registry%relinquish_field(ind)

  end subroutine surface_tension_source_term_compute

end module surface_tension_source_term
