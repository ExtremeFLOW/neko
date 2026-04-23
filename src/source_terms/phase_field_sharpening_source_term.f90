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
!! - "dealias": Use dealias for the normal/divergence term.
module phase_field_sharpening_source_term
  use num_types, only : rp
  use field, only : field_t
  use field_list, only : field_list_t
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use device_math, only : device_add2, device_cmult, device_col2, &
       device_copy, device_invcol2
  use dofmap, only : dofmap_t
  use gather_scatter, only : GS_OP_ADD, gs_t
  use interpolation, only : interpolator_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use operators, only : grad, div, dudxyz
  use profiler, only : profiler_start_region, profiler_end_region
  use scratch_registry, only : neko_scratch_registry, scratch_registry_t
  use field_registry, only : neko_field_registry
  use math, only : col2, add2, cmult, copy, invcol2
  use space, only : space_t, GL
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
  !! - "dealias": Use dealias for the normal/divergence term.
  type, public, extends(source_term_t) :: phase_field_sharpening_source_term_t
     !> The phase field.
     type(field_t), pointer :: phi => null()
     !> Sharpening coefficient.
     real(kind=rp) :: gamma
     !> Reference velocity scale.
     real(kind=rp) :: u_max
     !> Use dealias for the normal/divergence evaluation.
     logical :: dealias = .false.
     !> Function space used for the dealias evaluation.
     type(space_t) :: Xh_GL
     !> Dofmap associated with Xh_GL.
     type(dofmap_t) :: dm_Xh_GL
     !> Gather-scatter associated with Xh_GL.
     type(gs_t) :: gs_Xh_GL
     !> Coefficients associated with Xh_GL.
     type(coef_t) :: coef_GL
     !> Interpolator between the original and dealias spaces.
     type(interpolator_t) :: GLL_to_GL
     !> Scratch registry on the dealias space.
     type(scratch_registry_t) :: scratch_GL
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
    logical :: dealias

    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))
    call json_get_or_default(json, "scalar_field", scalar_name, "phase")
    call json_get(json, "gamma", gamma)
    call json_get(json, "u_max", u_max)
    call json_get_or_default(json, "dealias", dealias, .false.)

    call phase_field_sharpening_source_term_init_from_components(this, fields, &
         scalar_name, gamma, u_max, coef, start_time, end_time, &
         dealias)

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
       fields, scalar_name, gamma, u_max, coef, start_time, end_time, &
       dealias)
    class(phase_field_sharpening_source_term_t), intent(inout) :: this
    type(field_list_t), intent(in), target :: fields
    character(len=*), intent(in) :: scalar_name
    real(kind=rp), intent(in) :: gamma, u_max
    type(coef_t), intent(in), target :: coef
    real(kind=rp), intent(in) :: start_time, end_time
    logical, intent(in), optional :: dealias
    integer :: lx, lxd

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
    this%dealias = .false.
    if (present(dealias)) then
       this%dealias = dealias
    end if

    if (this%dealias) then
       lx = coef%Xh%lx
       lxd = 2 * lx

       call this%Xh_GL%init(GL, lxd, lxd, lxd)
       call this%dm_Xh_GL%init(this%coef%msh, this%Xh_GL)
       call this%gs_Xh_GL%init(this%dm_Xh_GL)
       call this%coef_GL%init(this%gs_Xh_GL)
       call this%GLL_to_GL%init(this%Xh_GL, this%coef%Xh)
       call this%scratch_GL%init(size=5, expansion_size=2, &
            dof=this%dm_Xh_GL)
    end if

  end subroutine phase_field_sharpening_source_term_init_from_components

  !> Destructor.
  subroutine phase_field_sharpening_source_term_free(this)
    class(phase_field_sharpening_source_term_t), intent(inout) :: this

    nullify(this%phi)
    if (this%dealias) then
       call this%scratch_GL%free()
       call this%GLL_to_GL%free()
       call this%coef_GL%free()
       call this%gs_Xh_GL%free()
       call this%dm_Xh_GL%free()
       call this%Xh_GL%free()
       this%dealias = .false.
    end if
    call this%free_base()
  end subroutine phase_field_sharpening_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param time The time state.
  subroutine phase_field_sharpening_source_term_compute(this, time)
    class(phase_field_sharpening_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: rhs_s

    rhs_s => this%fields%get(1)

    if (this%dealias) then
       call phase_field_sharpening_source_term_compute_dealias(this, rhs_s)
    else
       call phase_field_sharpening_source_term_compute_gll(this, rhs_s)
    end if

  end subroutine phase_field_sharpening_source_term_compute

  !> Computes the source term on the original GLL space.
  subroutine phase_field_sharpening_source_term_compute_gll(this, rhs_s)
    class(phase_field_sharpening_source_term_t), intent(inout) :: this
    type(field_t), pointer, intent(inout) :: rhs_s
    type(field_t), pointer :: work1, work2, work3, work4
    integer :: ind(4), i, n
    real(kind=rp) :: absgrad

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

  end subroutine phase_field_sharpening_source_term_compute_gll

  !> Computes the source term on a dealias GL space.
  subroutine phase_field_sharpening_source_term_compute_dealias(this, rhs_s)
    class(phase_field_sharpening_source_term_t), intent(inout) :: this
    type(field_t), pointer, intent(inout) :: rhs_s
    type(field_t), pointer :: work1_GL, work2_GL, work3_GL, work4_GL
    type(field_t), pointer :: phi_GL, accum_GLL, term_GLL
    integer :: ind_GL(5), ind_GLL(2), i, n_GL, n_GLL, nel
    real(kind=rp) :: absgrad

    call profiler_start_region("phase_sharpen_dealias", 26)

    call this%scratch_GL%request_field(work1_GL, ind_GL(1))
    call this%scratch_GL%request_field(work2_GL, ind_GL(2))
    call this%scratch_GL%request_field(work3_GL, ind_GL(3))
    call this%scratch_GL%request_field(work4_GL, ind_GL(4))
    call this%scratch_GL%request_field(phi_GL, ind_GL(5))
    call neko_scratch_registry%request_field(accum_GLL, ind_GLL(1))
    call neko_scratch_registry%request_field(term_GLL, ind_GLL(2))

    n_GL = work1_GL%size()
    n_GLL = rhs_s%size()
    nel = this%coef%msh%nelv

    call this%GLL_to_GL%map(phi_GL%x, this%phi%x, nel, this%Xh_GL)
    call grad(work1_GL%x, work2_GL%x, work3_GL%x, phi_GL%x, this%coef_GL)

    do i = 1, n_GL
       absgrad = sqrt(work1_GL%x(i,1,1,1)**2 + work2_GL%x(i,1,1,1)**2 &
            + work3_GL%x(i,1,1,1)**2)
       if (absgrad < 1.0e-12_rp) then
          work1_GL%x(i,1,1,1) = 0.0_rp
          work2_GL%x(i,1,1,1) = 0.0_rp
          work3_GL%x(i,1,1,1) = 0.0_rp
       else
          work1_GL%x(i,1,1,1) = &
               -phi_GL%x(i,1,1,1)*(1.0_rp - phi_GL%x(i,1,1,1)) &
               * work1_GL%x(i,1,1,1) / absgrad
          work2_GL%x(i,1,1,1) = &
               -phi_GL%x(i,1,1,1)*(1.0_rp - phi_GL%x(i,1,1,1)) &
               * work2_GL%x(i,1,1,1) / absgrad
          work3_GL%x(i,1,1,1) = &
               -phi_GL%x(i,1,1,1)*(1.0_rp - phi_GL%x(i,1,1,1)) &
               * work3_GL%x(i,1,1,1) / absgrad
       end if
    end do

    call dudxyz(work4_GL%x, work1_GL%x, this%coef_GL%drdx, &
         this%coef_GL%dsdx, this%coef_GL%dtdx, this%coef_GL)
    call phase_field_sharpening_source_term_map_to_gll(this, term_GLL, &
         work4_GL, n_GL, n_GLL, nel)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(accum_GLL%x_d, term_GLL%x_d, n_GLL)
    else
       call copy(accum_GLL%x, term_GLL%x, n_GLL)
    end if

    call dudxyz(work4_GL%x, work2_GL%x, this%coef_GL%drdy, &
         this%coef_GL%dsdy, this%coef_GL%dtdy, this%coef_GL)
    call phase_field_sharpening_source_term_map_to_gll(this, term_GLL, &
         work4_GL, n_GL, n_GLL, nel)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add2(accum_GLL%x_d, term_GLL%x_d, n_GLL)
    else
       call add2(accum_GLL%x, term_GLL%x, n_GLL)
    end if

    call dudxyz(work4_GL%x, work3_GL%x, this%coef_GL%drdz, &
         this%coef_GL%dsdz, this%coef_GL%dtdz, this%coef_GL)
    call phase_field_sharpening_source_term_map_to_gll(this, term_GLL, &
         work4_GL, n_GL, n_GLL, nel)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_add2(accum_GLL%x_d, term_GLL%x_d, n_GLL)
       call device_cmult(accum_GLL%x_d, this%gamma * this%u_max, n_GLL)
       call device_add2(rhs_s%x_d, accum_GLL%x_d, n_GLL)
    else
       call add2(accum_GLL%x, term_GLL%x, n_GLL)
       call cmult(accum_GLL%x, this%gamma * this%u_max, n_GLL)
       call add2(rhs_s%x, accum_GLL%x, n_GLL)
    end if

    call this%scratch_GL%relinquish_field(ind_GL)
    call neko_scratch_registry%relinquish_field(ind_GLL)
    call profiler_end_region("phase_sharpen_dealias", 26)

  end subroutine phase_field_sharpening_source_term_compute_dealias

  !> Maps a GL-space derivative contribution to the GLL RHS space.
  subroutine phase_field_sharpening_source_term_map_to_gll(this, term_GLL, &
       term_GL, n_GL, n_GLL, nel)
    class(phase_field_sharpening_source_term_t), intent(inout) :: this
    type(field_t), pointer, intent(inout) :: term_GLL
    type(field_t), pointer, intent(inout) :: term_GL
    integer, intent(in) :: n_GL, n_GLL, nel

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(term_GL%x_d, this%coef_GL%B_d, n_GL)
       call this%GLL_to_GL%map(term_GLL%x, term_GL%x, nel, this%coef%Xh)
       call device_invcol2(term_GLL%x_d, this%coef%B_d, n_GLL)
    else
       call col2(term_GL%x, this%coef_GL%B, n_GL)
       call this%GLL_to_GL%map(term_GLL%x, term_GL%x, nel, this%coef%Xh)
       call invcol2(term_GLL%x, this%coef%B, n_GLL)
    end if

  end subroutine phase_field_sharpening_source_term_map_to_gll

end module phase_field_sharpening_source_term
