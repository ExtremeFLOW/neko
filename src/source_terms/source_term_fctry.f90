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
!
!> Defines a factory subroutine for source terms.
module source_term_fctry
  use source_term, only : source_term_t
  use const_source_term, only : const_source_term_t
  use boussinesq_source_term, only : boussinesq_source_term_t
  use brinkman_source_term, only: brinkman_source_term_t
  use json_module, only : json_file
  use json_utils, only : json_get
  use field_list, only : field_list_t
  use utils, only : concat_string_array, neko_error
  use coefs, only : coef_t
  implicit none
  private

  public :: source_term_factory

  ! List of all possible types created by the factory routine
  character(len=20) :: KNOWN_TYPES(3) = [character(len=20) :: &
     "constant", &
     "boussinesq", &
     "brinkman"]

contains

  !> Source term factory. Both constructs and initializes the object.
  !! @param json JSON object initializing the source term.
  !! @param fields The list of fields updated by the source term.
  !! @param coef The SEM coefficients.
  subroutine source_term_factory(object, json, fields, coef)
    class(source_term_t), allocatable, intent(inout) :: object
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(inout) :: fields
    type(coef_t), intent(inout) :: coef
    character(len=:), allocatable :: type_name
    character(len=:), allocatable :: type_string

    type_string =  concat_string_array(KNOWN_TYPES, NEW_LINE('A') // "-  ", &
                                       .true.)

    call json_get(json, "type", type_name)

    if (trim(type_name) .eq. "constant") then
       allocate(const_source_term_t::object)
    else if (trim(type_name) .eq. "boussinesq") then
       allocate(boussinesq_source_term_t::object)
    else if (trim(type_name) .eq. "brinkman") then
       allocate(brinkman_source_term_t::object)
    else
       call neko_error("Unknown source term type: " &
                       // trim(type_name) // ".  Known types are: " &
                       // type_string)
    end if

    ! Initialize
    call object%init(json, fields, coef)

  end subroutine source_term_factory

end module source_term_fctry
