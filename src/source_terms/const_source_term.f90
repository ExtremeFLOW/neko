
! Copyright (c) 2020-2021, The Neko Authors
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
  implicit none
  private

  type, public, extends(source_term_t) :: const_source_term_t
     !> The value of the source term.
     real(kind=rp) :: value
   contains
     !> The common constructor using a JSON dictionary.
     procedure, pass(this) :: init => const_source_term_init_from_json
     !> The construct from type components.
     procedure, pass(this) :: init_from_compenents => & 
       const_source_term_init_from_components
     !> Destructor.
     procedure, pass(this) :: free => const_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute => const_source_term_compute
  end type const_source_term_t

contains
  subroutine const_source_term_init_from_json(this, json, fields, coef) 
    class(const_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(field_list_t), intent(inout), target :: fields
    type(coef_t), intent(inout) :: coef
    real(kind=rp) :: value

    call json_get(json, "value", value)
    call const_source_term_init_from_components(this, json, fields, value, coef)

  end subroutine const_source_term_init_from_json

  subroutine const_source_term_init_from_components(this, json, fields, value, &
                                                    coef) 
    class(const_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(field_list_t), intent(inout), target :: fields
    real(kind=rp), intent(in) :: value
    type(coef_t) :: coef

    call this%init_base(fields, coef)
  end subroutine const_source_term_init_from_components

  subroutine const_source_term_free(this) 
    class(const_source_term_t), intent(inout) :: this

    call this%free_base()
  end subroutine const_source_term_free

  subroutine const_source_term_compute(this, t, tstep) 
    class(const_source_term_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: n_fields, i

    n_fields = size(this%fields%fields)

    do i=1, n_fields
       this%fields%fields(i)%f%x  = this%fields%fields(i)%f%x + this%value
    end do
  end subroutine const_source_term_compute
  
end module const_source_term
