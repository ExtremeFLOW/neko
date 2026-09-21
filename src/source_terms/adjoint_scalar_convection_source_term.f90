!> @file adjoint_scalar_convection_source_term.f90
!! @copyright
!! Copyright (c) 2025-2026, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!
!> Implements the `adjoint_scalar_convection_source_term` type.
! this is a such a dumb name
module adjoint_scalar_convection_source_term
  use num_types, only: rp
  use field_list, only: field_list_t
  use field, only: field_t
  use json_module, only: json_file
  use time_state, only: time_state_t
  use source_term, only: source_term_t
  use interpolation, only: interpolator_t
  use space, only: space_t, GL
  use coefs, only: coef_t
  use field_math, only: field_subcol3, field_sub2, field_col3
  use operators, only: grad, dudxyz
  use utils, only: neko_error
  use gather_scatter, only: GS_OP_ADD
  use scratch_registry, only: neko_scratch_registry, scratch_registry_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use math, only: col2
  use device_math, only: device_col2
  implicit none
  private
  public :: adjoint_scalar_convection_source_term_allocate

  ! I don't know how to name this term, but when you have a passive
  ! scalar you get an extra term in the adjoint velocity equation, which comes
  ! from the convective term in the passive scalar equation.
  ! Maybe this should be called just `adjoint_scalar_convection`?
  ! but it does come in a source term...

  ! In any case,
  ! it's a source term acting on the adjoint velocity equations, of the form:
  ! \f$\nabla s s_adj\f$
  type, public, extends(source_term_t) :: &
       adjoint_scalar_convection_source_term_t
     !> adjoint passive scalar
     type(field_t), pointer :: s_adj => null()
     !> forward passive scalar
     type(field_t), pointer :: s => null()
     ! --- for over-integration
     !> The original space used in the simulation
     type(space_t), pointer :: Xh_GLL
     !> The additional higher-order space used in dealiasing
     type(space_t), pointer :: Xh_GL
     !> cfs of the higher-order space
     type(coef_t), pointer :: c_Xh_GL
     !> Interpolator between the original and higher-order spaces
     type(interpolator_t), pointer :: GLL_to_GL
     !> if dealiasing should be applied
     logical :: dealias
     !> GL scratch registry
     type(scratch_registry_t), pointer :: scratch_GL

   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => &
          adjoint_scalar_convection_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_components => &
          adjoint_scalar_convection_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => adjoint_scalar_convection_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => &
          adjoint_scalar_convection_source_term_compute
  end type adjoint_scalar_convection_source_term_t

contains

  !> Allocator for the adjoint scalar convection source term.
  subroutine adjoint_scalar_convection_source_term_allocate(obj)
    class(source_term_t), allocatable, intent(inout) :: obj
    allocate(adjoint_scalar_convection_source_term_t::obj)
  end subroutine adjoint_scalar_convection_source_term_allocate

  !> The common constructor using a JSON object.
  !! @param this The object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param variable_name The name of the variable where the source term acts.
  subroutine adjoint_scalar_convection_source_term_init_from_json(this, &
       json, fields, coef, variable_name)
    class(adjoint_scalar_convection_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name

    ! this is a bit weird... because I don't think this should come from the
    ! JSON.
    ! Maybe we should think of all these source terms as only "appendable"
    !
    ! Because we'll never have the whole case here, so we'll never be able
    ! init from components anyway...


  end subroutine adjoint_scalar_convection_source_term_init_from_json

  !> The constructor from type components.
  !! @param this The source term.
  !! @param f_x, f_y, f_z the RHS of the equation (either primal or adjoint).
  !! @param s the primal scalar
  !! @param s_adj the primal scalar
  !! @param coef The SEM coeffs.
  !! @param c_Xh_GL The SEM coeffs on the over integration mesh.
  !! @param GLL_to_GL Interpolator between GLL and GL.
  !! @param dealias if dealiasing should be applied.
  !! @param scratch_GL A scratch registry on the GL space.
  subroutine adjoint_scalar_convection_source_term_init_from_components(this,&
       f_x, f_y, f_z, s, s_adj, coef, c_Xh_GL, GLL_to_GL, dealias, scratch_GL)
    class(adjoint_scalar_convection_source_term_t), intent(inout) :: this
    type(field_t), pointer, intent(in) :: f_x, f_y, f_z
    type(field_t), intent(in), target :: s, s_adj
    type(coef_t), intent(in), target :: coef
    type(coef_t), intent(in), target :: c_Xh_GL
    type(interpolator_t), intent(in), target :: GLL_to_GL
    logical, intent(in) :: dealias
    type(scratch_registry_t), intent(in), target :: scratch_GL

    type(field_list_t) :: fields
    real(kind=rp) :: start_time
    real(kind=rp) :: end_time

    ! I wish you didn't need a start time and end time...
    ! but I'm just going to set a super big number...
    start_time = 0.0_rp
    end_time = 100000000.0_rp

    call this%free()

    ! this is copying the fluid source term init
    ! We package the fields for the source term to operate on in a field list.
    call fields%init(3)
    call fields%assign(1, f_x)
    call fields%assign(2, f_y)
    call fields%assign(3, f_z)

    call this%init_base(fields, coef, start_time, end_time)
    call fields%free()

    ! point everything in the correct places
    this%s_adj => s_adj
    this%s => s

    ! for over integration
    this%dealias = dealias
    this%c_Xh_GL => c_Xh_GL
    this%Xh_GL => this%c_Xh_GL%Xh
    this%Xh_GLL => this%coef%Xh
    this%GLL_to_GL => GLL_to_GL
    this%scratch_GL => scratch_GL

  end subroutine adjoint_scalar_convection_source_term_init_from_components

  !> Destructor.
  subroutine adjoint_scalar_convection_source_term_free(this)
    class(adjoint_scalar_convection_source_term_t), intent(inout) :: this

    call this%free_base()
    nullify(this%s_adj)
    nullify(this%s)
    nullify(this%c_Xh_GL)
    nullify(this%Xh_GL)
    nullify(this%Xh_GLL)
    nullify(this%GLL_to_GL)
    nullify(this%scratch_GL)

  end subroutine adjoint_scalar_convection_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param this The object.
  !! @param time The time state.
  subroutine adjoint_scalar_convection_source_term_compute(this, time)
    class(adjoint_scalar_convection_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: fu, fv, fw
    integer :: temp_indices(4)
    type(field_t), pointer :: dsdx, dsdy, dsdz, work
    type(field_t), pointer :: accumulate, fld_GL, s_GL, s_adj_GL
    integer :: temp_indices_GL(4)
    integer :: n_GL, nel


    call neko_scratch_registry%request_field(dsdx, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(dsdy, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(dsdz, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(work, temp_indices(4), .false.)

    fu => this%fields%get(1)
    fv => this%fields%get(2)
    fw => this%fields%get(3)

    ! we need the term \f$\nabla s s_adj\f$
    if (this%dealias) then
       nel = this%coef%msh%nelv
       n_GL = nel * this%Xh_GL%lxyz
       call this%scratch_GL%request_field(accumulate, temp_indices_GL(1), .false.)
       call this%scratch_GL%request_field(fld_GL, temp_indices_GL(2), .false.)
       call this%scratch_GL%request_field(s_GL, temp_indices_GL(3), .false.)
       call this%scratch_GL%request_field(s_adj_GL, temp_indices_GL(4), .false.)

       call this%GLL_to_GL%map(s_GL%x, this%s%x, nel, this%Xh_GL)
       call this%GLL_to_GL%map(s_adj_GL%x, this%s_adj%x, nel, this%Xh_GL)

       ! u
       call dudxyz(fld_GL%x, s_GL%x, this%c_Xh_GL%drdx, &
            this%c_Xh_GL%dsdx, this%c_Xh_GL%dtdx, this%c_Xh_GL)
       call field_col3(accumulate, s_adj_GL, fld_GL)
       ! Evaluate term on GL and preempt the GLL premultiplication
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_col2(accumulate%x_d, this%c_Xh_GL%B_d, n_GL)
          call this%GLL_to_GL%map(work%x, accumulate%x, nel, this%Xh_GLL)
          ! Sum contributions from adjoining elements at shared dofs before
          ! normalizing.
          call this%coef%gs_h%op(work, GS_OP_ADD)
          call device_col2(work%x_d, this%coef%Binv_d, work%size())
       else
          call col2(accumulate%x, this%c_Xh_GL%B, n_GL)
          call this%GLL_to_GL%map(work%x, accumulate%x, nel, this%Xh_GLL)
          ! See comment in the device branch above.
          call this%coef%gs_h%op(work, GS_OP_ADD)
          call col2(work%x, this%coef%Binv, work%size())
       end if
       call field_sub2(fu, work)

       ! v
       call dudxyz(fld_GL%x, s_GL%x, this%c_Xh_GL%drdy, &
            this%c_Xh_GL%dsdy, this%c_Xh_GL%dtdy, this%c_Xh_GL)
       call field_col3(accumulate, s_adj_GL, fld_GL)
       ! Evaluate term on GL and preempt the GLL premultiplication
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_col2(accumulate%x_d, this%c_Xh_GL%B_d, n_GL)
          call this%GLL_to_GL%map(work%x, accumulate%x, nel, this%Xh_GLL)
          ! Sum contributions from adjoining elements at shared dofs before
          ! normalizing.
          call this%coef%gs_h%op(work, GS_OP_ADD)
          call device_col2(work%x_d, this%coef%Binv_d, work%size())
       else
          call col2(accumulate%x, this%c_Xh_GL%B, n_GL)
          call this%GLL_to_GL%map(work%x, accumulate%x, nel, this%Xh_GLL)
          ! See comment in the device branch above.
          call this%coef%gs_h%op(work, GS_OP_ADD)
          call col2(work%x, this%coef%Binv, work%size())
       end if
       call field_sub2(fv, work)

       ! w
       call dudxyz(fld_GL%x, s_GL%x, this%c_Xh_GL%drdz, &
            this%c_Xh_GL%dsdz, this%c_Xh_GL%dtdz, this%c_Xh_GL)
       call field_col3(accumulate, s_adj_GL, fld_GL)
       ! Evaluate term on GL and preempt the GLL premultiplication
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_col2(accumulate%x_d, this%c_Xh_GL%B_d, n_GL)
          call this%GLL_to_GL%map(work%x, accumulate%x, nel, this%Xh_GLL)
          ! Sum contributions from adjoining elements at shared dofs before
          ! normalizing.
          call this%coef%gs_h%op(work, GS_OP_ADD)
          call device_col2(work%x_d, this%coef%Binv_d, work%size())
       else
          call col2(accumulate%x, this%c_Xh_GL%B, n_GL)
          call this%GLL_to_GL%map(work%x, accumulate%x, nel, this%Xh_GLL)
          ! See comment in the device branch above.
          call this%coef%gs_h%op(work, GS_OP_ADD)
          call col2(work%x, this%coef%Binv, work%size())
       end if
       call field_sub2(fw, work)

       call this%scratch_GL%relinquish_field(temp_indices_GL)

    else
       call grad(dsdx%x, dsdy%x, dsdz%x, this%s%x, this%coef)
       call field_subcol3(fu, this%s_adj, dsdx)
       call field_subcol3(fv, this%s_adj, dsdy)
       call field_subcol3(fw, this%s_adj, dsdz)
    end if

    ! free the scratch
    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine adjoint_scalar_convection_source_term_compute

end module adjoint_scalar_convection_source_term
