! Copyright (c) 2023-2024, The Neko Authors
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
submodule (source_term) source_term_fctry
  use const_source_term, only : const_source_term_t
  use boussinesq_source_term, only : boussinesq_source_term_t
  use brinkman_source_term, only: brinkman_source_term_t
  use coriolis_source_term, only : coriolis_source_term_t
  use centrifugal_source_term, only : centrifugal_source_term_t
  use gradient_jump_penalty, only : gradient_jump_penalty_t
  use json_utils, only : json_get
  use utils, only : neko_type_error, neko_type_registration_error
  implicit none

  ! List of all possible types created by the factory routine
  character(len=20) :: SOURCE_KNOWN_TYPES(6) = [character(len=20) :: &
       "constant", &
       "boussinesq", &
       "coriolis", &
       "centrifugal", &
       "gradient_jump_penalty", &
       "brinkman"]

contains

  !> Source term factory. Both constructs and initializes the object.
  !! @param json JSON object initializing the source term.
  !! @param fields The list of fields updated by the source term.
  !! @param coef The SEM coefficients.
  module subroutine source_term_factory(object, json, fields, coef, &
       variable_name)
    class(source_term_t), allocatable, intent(inout) :: object
    type(json_file), intent(inout) :: json
    type(field_list_t), intent(inout) :: fields
    type(coef_t), intent(inout) :: coef
    character(len=*), intent(in) :: variable_name
    character(len=:), allocatable :: type_name
    character(len=:), allocatable :: type_string

    call json_get(json, "type", type_name)

    ! Allocate
    call source_term_allocator(object, type_name)

    ! Initialize
    call object%init(json, fields, coef, variable_name)

  end subroutine source_term_factory

  !> Source term allocator.
  !! @param object The object to be allocated.
  !! @param type_name The name of the type to allocate.
  module subroutine source_term_allocator(object, type_name)
    class(source_term_t), allocatable, intent(inout) :: object
    character(len=:), allocatable, intent(in) :: type_name
    integer :: i

    select case (trim(type_name))
    case ("constant")
       allocate(const_source_term_t::object)
    case ("boussinesq")
       allocate(boussinesq_source_term_t::object)
    case ("coriolis")
       allocate(coriolis_source_term_t::object)
    case ("centrifugal")
       allocate(centrifugal_source_term_t::object)
    case ("brinkman")
       allocate(brinkman_source_term_t::object)
    case ("gradient_jump_penalty")
       allocate(gradient_jump_penalty_t::object)
    case default
       do i = 1, source_term_registry_size
          if (trim(type_name) .eq. trim(source_term_registry(i)%type_name)) then
             call source_term_registry(i)%allocator(object)
             return
          end if
       end do

       call neko_type_error("source term", type_name, SOURCE_KNOWN_TYPES)
    end select
  end subroutine source_term_allocator

  !> Register a custom source term allocator.
  !! Called in custom user modules inside the `module_name_register_types`
  !! routine to add a custom type allocator to the registry.
  !! @param type_name The name of the type to allocate.
  !! @param allocator The allocator for the custom user type.
  module subroutine register_source_term(type_name, allocator)
    character(len=*), intent(in) :: type_name
    procedure(source_term_allocate), pointer, intent(in) :: allocator
    type(allocator_entry), allocatable :: temp(:)
    integer :: i

    do i = 1, size(SOURCE_KNOWN_TYPES)
       if (trim(type_name) .eq. trim(SOURCE_KNOWN_TYPES(i))) then
          call neko_type_registration_error("source term", type_name, .true.)
       end if
    end do

    do i = 1, source_term_registry_size
       if (trim(type_name) .eq. trim(source_term_registry(i)%type_name)) then
          call neko_type_registration_error("source term", type_name, .false.)
       end if
    end do

    ! Expand registry
    if (source_term_registry_size .eq. 0) then
       allocate(source_term_registry(1))
    else
       allocate(temp(source_term_registry_size + 1))
       temp(1:source_term_registry_size) = source_term_registry
       call move_alloc(temp, source_term_registry)
    end if

    source_term_registry_size = source_term_registry_size + 1
    source_term_registry(source_term_registry_size)%type_name = type_name
    source_term_registry(source_term_registry_size)%allocator => allocator
  end subroutine register_source_term

end submodule source_term_fctry
