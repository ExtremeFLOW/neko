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
!> Implements the `source_term_t` type and a wrapper `source_term_wrapper_t`.
module source_term
  use neko_config
  use num_types, only : rp
  use coefs, only : coef_t
  use field_list, only : field_list_t
  use json_module, only : json_file
  implicit none
  private

  !> Base abstract type for source terms selected at run time.
  !! @note The user-provided source term is treated seperately in the type
  !! `user_source_term.`
  type, abstract, public:: source_term_t
     !> The fields to be updated with the source term values
     type(field_list_t) :: fields
     !> Coefficients for the SEM.
     type(coef_t), pointer :: coef => null()

   contains
     !> Constructor for the source_term_t (base) type.
     procedure, pass(this) :: init_base => source_term_init_base
     !> Destructor for the source_term_t (base) type.
     procedure, pass(this) :: free_base => source_term_free_base
     !> The common constructor using a JSON object.
     procedure(source_term_init), pass(this), deferred :: init
     !> Destructor.
     procedure(source_term_free), pass(this), deferred :: free
     !> Computes the source term and adds the result to `fields`.
     procedure(source_term_compute), pass(this), deferred :: compute
  end type source_term_t

  !> A helper type that is needed to have an array of polymorphic objects
  type, public :: source_term_wrapper_t
     !> Wrapped polymorphic source term.
     class(source_term_t), allocatable :: source_term
   contains
     !> Destructor.
     procedure, pass(this) :: free => source_term_wrapper_free 
  end type source_term_wrapper_t

  abstract interface
    !> The common constructor using a JSON object.
    !! @param json The JSON object for the source.
    !! @param fields A list of fields for adding the source values.
    !! @param coef The SEM coeffs.
     subroutine source_term_init(this, json, fields, coef)  
       import source_term_t, json_file, field_list_t, coef_t
       class(source_term_t), intent(inout) :: this
       type(json_file), intent(inout) :: json
       class(field_list_t), intent(inout), target :: fields
       type(coef_t), intent(inout) :: coef
     end subroutine
  end interface

  abstract interface
     !> Destructor.
     subroutine source_term_free(this)  
       import source_term_t
       class(source_term_t), intent(inout) :: this
     end subroutine
  end interface

  abstract interface
     !> Computes the source term and adds the result to `fields`.
     !! @param t The time value.
     !! @param tstep The current time-step
     subroutine source_term_compute(this, t, tstep)  
       import source_term_t, rp
       class(source_term_t), intent(inout) :: this
       real(kind=rp), intent(in) :: t
       integer, intent(in) :: tstep
     end subroutine
  end interface
contains

  !> Constructor for the `source_term_t` (base) type.
  !> @param fields A list of pointers to fields to be updated by the source 
  !! term.
  !> @param coef SEM coefs.
  subroutine source_term_init_base(this, fields, coef) 
    class(source_term_t), intent(inout) :: this
    type(field_list_t) :: fields
    type(coef_t), intent(inout), target :: coef
    integer :: n_fields, i

    this%coef => coef
    n_fields = size(fields%fields)
    allocate(this%fields%fields(n_fields))

    ! A lot of attribute nesting here due to Fortran needing wrapper types
    ! but this is just pointer assignement for the fields.
    do i=1, n_fields
       this%fields%fields(i)%f => fields%fields(i)%f
    end do
  end subroutine source_term_init_base

  !> Destructor for the `source_term_t` (base) type.
  subroutine source_term_free_base(this) 
    class(source_term_t), intent(inout) :: this
    integer :: n_fields, i

    n_fields = size(this%fields%fields)

    do i=1, n_fields
       nullify(this%fields%fields(i)%f)
    end do
    deallocate(this%fields%fields)
    nullify(this%coef)
  end subroutine source_term_free_base

  !> Destructor for the `source_term_wrapper_t` type.
  subroutine source_term_wrapper_free(this) 
    class(source_term_wrapper_t), intent(inout) :: this
    integer :: n_fields, i

    if (allocated(this%source_term)) then
       call this%source_term%free()
       deallocate(this%source_term)
    end if
  end subroutine source_term_wrapper_free
  
end module source_term
