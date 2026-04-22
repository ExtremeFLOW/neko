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
!> Implements the `hpfrt_source_term_t` type.
module hpfrt_source_term
  use num_types, only : rp
  use field, only : field_t
  use field_list, only : field_list_t
  use field_math, only : field_add2, field_cmult, field_sub3
  use json_module, only : json_file
  use json_utils, only : json_get_or_default, json_get_or_lookup
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use elementwise_filter, only : elementwise_filter_t
  use registry, only : neko_registry
  use scratch_registry, only : neko_scratch_registry
  use time_state, only : time_state_t
  use utils, only : neko_error
  implicit none
  private

  !> High-pass filter relaxation source term.
  !! @details
  !! The low-pass filter is built from a one-dimensional modal transfer
  !! function. For `filter_modes` = \f$k_f\f$ and \f$n = lx\f$, let
  !! \f$k_0 = n - k_f\f$. The modal transfer function is
  !! \f[
  !! \sigma_k =
  !! \begin{cases}
  !! 1, & k \leq k_0,\\
  !! 1 - \left(\frac{k-k_0}{k_f}\right)^2, & k_0 < k \leq n.
  !! \end{cases}
  !! \f]
  !! The source term adds \f$\chi (I - F)\phi\f$, where
  !! \f$\chi = -|\texttt{filter_weight}|\f$.
  type, public, extends(source_term_t) :: hpfrt_source_term_t
     !> Number of high modes affected by the filter.
     integer :: filter_modes = 0
     !> Damping coefficient, stored as a non-positive value.
     real(kind=rp) :: chi = 0.0_rp
     !> Low-pass elementwise filter.
     type(elementwise_filter_t) :: filter
     !> Fields to filter before adding to the source RHS.
     type(field_list_t) :: source_fields
   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => hpfrt_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_components => &
          hpfrt_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => hpfrt_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => hpfrt_source_term_compute
  end type hpfrt_source_term_t

contains
  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param variable_name The name of the variable where the source term acts.
  subroutine hpfrt_source_term_init_from_json(this, json, fields, coef, &
       variable_name)
    class(hpfrt_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name
    real(kind=rp) :: filter_weight, start_time, end_time
    integer :: filter_modes

    call json_get_or_lookup(json, "filter_weight", filter_weight)
    call json_get_or_lookup(json, "filter_modes", filter_modes)
    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))

    call hpfrt_source_term_init_from_components(this, fields, coef, &
         filter_weight, filter_modes, start_time, end_time, variable_name)

  end subroutine hpfrt_source_term_init_from_json

  !> The constructor from type components.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param filter_weight The positive damping coefficient magnitude.
  !! @param filter_modes Number of high modes affected by the filter.
  !! @param start_time When to start adding the source term.
  !! @param end_time When to stop adding the source term.
  !! @param field_name Name of the scalar field for scalar source terms.
  subroutine hpfrt_source_term_init_from_components(this, fields, coef, &
       filter_weight, filter_modes, start_time, end_time, field_name)
    class(hpfrt_source_term_t), intent(inout) :: this
    class(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    real(kind=rp), intent(in) :: filter_weight
    integer, intent(in) :: filter_modes
    real(kind=rp), intent(in) :: start_time
    real(kind=rp), intent(in) :: end_time
    character(len=*), intent(in) :: field_name
    integer :: nx, i, k0
    real(kind=rp) :: amp
    real(kind=rp), allocatable :: transfer(:)
    type(field_t), pointer :: source_field

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)

    if (fields%size() .ne. 1 .and. fields%size() .ne. 3) then
       call neko_error("HPFRT source term expects either 1 or 3 fields.")
    end if

    nx = coef%dof%xh%lx
    if (filter_modes .lt. 1 .or. filter_modes .gt. nx) then
       call neko_error("HPFRT filter_modes must be between 1 and lx.")
    end if

    this%filter_modes = filter_modes
    this%chi = -abs(filter_weight)

    allocate(transfer(nx))
    transfer = 1.0_rp
    k0 = nx - filter_modes
    do i = k0 + 1, nx
       amp = real((i - k0) * (i - k0), kind=rp) / &
            real(filter_modes * filter_modes, kind=rp)
       transfer(i) = 1.0_rp - amp
    end do

    call this%filter%init_from_components(coef, "nonBoyd", transfer)

    call this%source_fields%init(fields%size())

    if (fields%size() .eq. 3) then
       call this%source_fields%assign(1, neko_registry%get_field("u"))
       call this%source_fields%assign(2, neko_registry%get_field("v"))
       call this%source_fields%assign(3, neko_registry%get_field("w"))
    else
       if (.not. neko_registry%field_exists(field_name)) then
          call neko_error("HPFRT source field does not exist: " // &
               trim(field_name))
       end if
       source_field => neko_registry%get_field(field_name)
       call this%source_fields%assign(1, source_field)
    end if

  end subroutine hpfrt_source_term_init_from_components

  !> Destructor.
  subroutine hpfrt_source_term_free(this)
    class(hpfrt_source_term_t), intent(inout) :: this

    call this%free_base()
    call this%filter%free()
    call this%source_fields%free()

    this%filter_modes = 0
    this%chi = 0.0_rp

  end subroutine hpfrt_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param time The time state.
  subroutine hpfrt_source_term_compute(this, time)
    class(hpfrt_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: rhs, source, work
    integer :: i, work_idx

    call neko_scratch_registry%request_field(work, work_idx, .false.)

    do i = 1, this%fields%size()
       rhs => this%fields%get(i)
       source => this%source_fields%get(i)

       call this%filter%apply(work, source)
       call field_sub3(work, source, work)
       call field_cmult(work, this%chi)
       call field_add2(rhs, work)
    end do

    call neko_scratch_registry%relinquish_field(work_idx)

  end subroutine hpfrt_source_term_compute

end module hpfrt_source_term
