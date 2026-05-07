! Copyright (c) 2025, The Neko Authors
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
!> In this module we implement a custom source term, `my_source_term` and prep
!! it for being recognized by Neko at run time. The source term doesn't really
!! do anything, we just show how such a thing could be done.

!! After compiling this file with makeneko, you can add the following to your
!! JSON file to use this custom source term:
!!        "source_terms": [
!!          {
!!            "type": "my_source_term",
!!            "greeting": "Hello there!"
!!
!!          }
!!        ],
!! NOTE: the module name must be the same as the file name sans the extension.
module custom_types
  use num_types, only : rp
  use json_module, only : json_file
  use json_utils, only: json_get, json_get_or_default
  use field_list, only : field_list_t
  use time_state, only : time_state_t

  ! These imports are needed for registering our new type with Neko
  use source_term, only : source_term_t, register_source_term, &
       source_term_allocate

  use coefs, only : coef_t
  implicit none
  private

  type, public, extends(source_term_t) :: my_source_term_t
     ! We will read this greeting from the JSON file.
     character(len=:), allocatable :: greeting
   contains
     !> The common constructor using a JSON object. This is where we parse
     !! the JSON object and initialize the source term.
     procedure, pass(this) :: init => my_source_term_init_from_json
     !> Destructor.
     procedure, pass(this) :: free => my_source_term_free
     !> Computes the source term and adds the result to `fields`.
     procedure, pass(this) :: compute_ => my_source_term_compute
  end type my_source_term_t

  ! Have to make this public!
  public :: custom_types_register_types

contains
  !> The common constructor using a JSON object.
  !! @param json The JSON object for the source.
  !! @param fields A list of fields for adding the source values.
  !! @param coef The SEM coeffs.
  subroutine my_source_term_init_from_json(this, json, fields, coef, &
       variable_name)
    class(my_source_term_t), intent(inout) :: this
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(in), target :: fields
    type(coef_t), intent(in), target :: coef
    character(len=*), intent(in) :: variable_name
    real(kind=rp) :: start_time, end_time

    call json_get_or_default(json, "start_time", start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", end_time, huge(0.0_rp))

    ! Start fresh
    call this%free()
    ! Initialize the base class source_term_t
    call this%init_base(fields, coef, start_time, end_time)

    call json_get(json, "greeting", this%greeting)

  end subroutine my_source_term_init_from_json

  !> Destructor.
  subroutine my_source_term_free(this)
    class(my_source_term_t), intent(inout) :: this

    if (allocated(this%greeting)) then
       deallocate(this%greeting)
    end if
    call this%free_base()
  end subroutine my_source_term_free

  !> Will just bring our greeting to the console.
  !! @param time The time state.
  subroutine my_source_term_compute(this, time)
    class(my_source_term_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time

    write(*,*) this%greeting

  end subroutine my_source_term_compute

  ! The name of this routine follows a convention: the name of the module +
  ! register_types. This is where we call the register_source_term routine, or
  ! a similar routine for other types. You can register as much types as you
  ! want in this routine, in case you have several types in the module.
  subroutine custom_types_register_types()
    ! Just a helper variable
    procedure(source_term_allocate), pointer :: allocator

    allocator => my_source_term_allocate

    ! Based on this the name of the source term will be "my_source_term",
    ! This is what you set in the JSON file, in the type keyword
    call register_source_term("my_source_term", allocator)
  end subroutine custom_types_register_types

  ! This thing does nothing except allocate a polymorphic object to our new
  ! custom type.
  subroutine my_source_term_allocate(obj)
    class(source_term_t), allocatable, intent(inout) :: obj
    allocate(my_source_term_t::obj)
  end subroutine my_source_term_allocate
end module custom_types
