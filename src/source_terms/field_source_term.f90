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
!> Implements the `field_source_term_t` type.
module field_source_term
  use num_types, only : rp
  use field_list, only : field_list_t
  use json_module, only : json_file
  use json_utils, only : json_get_or_default, json_get_or_lookup, json_get
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use utils, only : neko_error, NEKO_VARNAME_LEN
  use time_state, only : time_state_t
  use registry, only : neko_registry
  use field_math, only : field_add2
  implicit none
  private

  !> A source term that grabs the values from fields in the registry.
  !! The fields are specified with the `field_names` keyword, which should be an
  !! array, with a value for each component of the source.
  type, public, extends(source_term_t) :: field_source_term_t
     !> The values for the source term, one for each field.
     type(field_list_t) :: registry_fields
   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => field_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_compenents => &
          field_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => field_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => field_source_term_compute
  end type field_source_term_t

contains
  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  !! @param variable_name The name of the variable where the source term
  !! acts.
  subroutine field_source_term_init_from_json(this, json, fields, coef, &
       variable_name)
    class(field_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name
    real(kind=rp) :: start_time, end_time
    character(len=NEKO_VARNAME_LEN), allocatable :: field_names(:)

    call json_get(json, "field_names", field_names)
    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))


    call field_source_term_init_from_components(this, fields, field_names, &
         coef, start_time, end_time)

  end subroutine field_source_term_init_from_json

  !> The constructor from type components.
  !! @param fields A list of fields for adding the source values.
  !! @param field_names The name of the fields that define the source term
  !! strength.
  !! @param coef The SEM coeffs.
  !! @param start_time When to start adding the source term.
  !! @param end_time When to stop adding the source term.
  subroutine field_source_term_init_from_components(this, fields, field_names, &
       coef, start_time, end_time)
    class(field_source_term_t), intent(inout) :: this
    class(field_list_t), intent(in), target :: fields
    character(len=*), intent(in) :: field_names(:)
    type(coef_t), target :: coef
    real(kind=rp), intent(in) :: start_time
    real(kind=rp), intent(in) :: end_time

    integer :: i

    call this%free()
    call this%init_base(fields, coef, start_time, end_time)

    if (size(field_names) .ne. fields%size()) then
       call neko_error("Number of fields and field names inconsistent.")
    end if

    call this%registry_fields%init(size(field_names))

    do i = 1, size(field_names)
       ! Add zero-valued if doesn't exist.
       ! May occur due to initialization order
       call neko_registry%add_field(this%coef%dof, field_names(i), &
            ignore_existing = .true.)
       call this%registry_fields%assign(i, &
            neko_registry%get_field(field_names(i)))
    end do

  end subroutine field_source_term_init_from_components

  !> Destructor.
  subroutine field_source_term_free(this)
    class(field_source_term_t), intent(inout) :: this

    call this%free_base()
    call this%registry_fields%free()

  end subroutine field_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param time The time state.
  subroutine field_source_term_compute(this, time)
    class(field_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer :: n_fields, i, n

    n_fields = this%fields%size()

    do i = 1, n_fields
       call field_add2(this%fields%get(i), this%registry_fields%get(i))
    end do
  end subroutine field_source_term_compute

end module field_source_term
