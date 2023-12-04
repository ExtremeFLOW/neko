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
!> Implements the `const_source_term_t` type.
module const_source_term
  use num_types, only : rp
  use field_list, only : field_list_t
  use json_module, only : json_file
  use json_utils, only: json_get
  use source_term, only : source_term_t
  use coefs, only : coef_t
  use neko_config, only : NEKO_BCKND_DEVICE
  use utils, only : neko_error
  use const_source_term_cpu, only : const_source_term_compute_cpu
  use const_source_term_device, only : const_source_term_compute_device
  implicit none
  private

  !> A constant source term.
  !! The strength is specified with the `values` keyword, which should be an
  !! array, with a value for each component of the source.
  type, public, extends(source_term_t) :: const_source_term_t
     !> The value of the source term.
     real(kind=rp), allocatable :: values(:)
   contains
     !> The common constructor using a JSON object.
     procedure, pass(this) :: init => const_source_term_init_from_json
     !> The constructor from type components.
     procedure, pass(this) :: init_from_compenents => & 
       const_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => const_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute => const_source_term_compute
  end type const_source_term_t

contains
  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  subroutine const_source_term_init_from_json(this, json, fields, coef) 
    class(const_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(inout), target :: fields
    type(coef_t), intent(inout) :: coef
    real(kind=rp), allocatable :: values(:)

    call json_get(json, "values", values)
    call const_source_term_init_from_components(this, fields, values, coef)

  end subroutine const_source_term_init_from_json

  !> The constructor from type components.
  !! @param fields A list of fields for adding the source values.
  !! @param values The array of values, one for each field.
  !! @param coef The SEM coeffs.
  subroutine const_source_term_init_from_components(this, fields, values, &
                                                    coef) 
    class(const_source_term_t), intent(inout) :: this
    class(field_list_t), intent(inout), target :: fields
    real(kind=rp), intent(in) :: values(:)
    type(coef_t) :: coef

    call this%free()
    call this%init_base(fields, coef)

    if (size(values) .ne. size(fields%fields)) then
       call neko_error("Number of fields and values inconsistent.")
    end if

    this%values = values
  end subroutine const_source_term_init_from_components

  !> Destructor.
  subroutine const_source_term_free(this) 
    class(const_source_term_t), intent(inout) :: this

    call this%free_base()
  end subroutine const_source_term_free

  !> Computes the source term and adds the result to `fields`.
  !! @param t The time value.
  !! @param tstep The current time-step.
  subroutine const_source_term_compute(this, t, tstep) 
    class(const_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: n_fields, i, n

    n_fields = size(this%fields%fields)
    n = this%fields%fields(1)%f%dof%size()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call const_source_term_compute_device(this%fields, this%values)
    else
       call const_source_term_compute_cpu(this%fields, this%values)
    end if
  end subroutine const_source_term_compute
  
end module const_source_term
