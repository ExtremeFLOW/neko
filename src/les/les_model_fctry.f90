! Copyright (c) 2021-2025, The Neko Authors
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
submodule (les_model) les_model_fctry
  use vreman, only : vreman_t
  use smagorinsky, only : smagorinsky_t
  use dynamic_smagorinsky, only : dynamic_smagorinsky_t
  use sigma, only : sigma_t
  use fluid_scheme_base, only : fluid_scheme_base_t
  use wale, only : wale_t
  use utils, only : neko_type_registration_error
  implicit none

  ! List of all possible types created by the factory routine
  character(len=20) :: LES_KNOWN_TYPES(5) = [character(len=20) :: &
       "vreman", &
       "smagorinsky", &
       "dymamic_smagorinsky", &
       "sigma", &
       "wale"]

contains
  !> LES model factory.
  !! @param object The object to be allocated.
  !! @param type_name The name of the LES model.
  !! @param fluid The fluid scheme base type pointer.
  !! @param dofmap SEM map of degrees of freedom.
  !! @param coef SEM coefficients.
  !! @param json A dictionary with parameters.
  module subroutine les_model_factory(object, type_name, fluid, json)
    class(les_model_t), allocatable, intent(inout) :: object
    character(len=*), intent(in) :: type_name
    class(fluid_scheme_base_t), intent(inout) :: fluid
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: type_string

    call les_model_allocator(object, type_name)
    call object%init(fluid, json)
  end subroutine les_model_factory

  !> LES model allocator.
  !! @param object The object to be allocated.
  !! @param type_name The name of the LES model.
  module subroutine les_model_allocator(object, type_name)
    class(les_model_t), allocatable, intent(inout) :: object
    character(len=*), intent(in) :: type_name
    integer :: i

    if (allocated(object)) deallocate(object)

    select case (trim(type_name))
    case ('vreman')
       allocate(vreman_t::object)
    case ('smagorinsky')
       allocate(smagorinsky_t::object)
    case ('dynamic_smagorinsky')
       allocate(dynamic_smagorinsky_t::object)
    case ('sigma')
       allocate(sigma_t::object)
    case ('wale')
       allocate(wale_t::object)
    case default
       do i = 1, les_model_registry_size
          if (trim(type_name) == trim(les_model_registry(i)%type_name)) then
             call les_model_registry(i)%allocator(object)
             return
          end if
       end do

       call neko_type_error("LES model", type_name, LES_KNOWN_TYPES)
    end select

  end subroutine les_model_allocator

  !> Register a custom LES model allocator.
  !! Called in custom user modules inside the `module_name_register_types`
  !! routine to add a custom type allocator to the registry.
  !! @param type_name The name of the type to allocate.
  !! @param allocator The allocator for the custom user type.
  module subroutine register_les_model(type_name, allocator)
    character(len=*), intent(in) :: type_name
    procedure(les_model_allocate), pointer, intent(in) :: allocator
    type(allocator_entry), allocatable :: temp(:)
    integer :: i

    do i = 1, size(LES_KNOWN_TYPES)
       if (trim(type_name) .eq. trim(LES_KNOWN_TYPES(i))) then
          call neko_type_registration_error("LES model", type_name, .true.)
       end if
    end do

    do i = 1, les_model_registry_size
       if (trim(type_name) .eq. trim(les_model_registry(i)%type_name)) then
          call neko_type_registration_error("LES model", type_name, .false.)
       end if
    end do

    ! Expand registry
    if (les_model_registry_size == 0) then
       allocate(les_model_registry(1))
    else
       allocate(temp(les_model_registry_size + 1))
       temp(1:les_model_registry_size) = les_model_registry
       call move_alloc(temp, les_model_registry)
    end if

    les_model_registry_size = les_model_registry_size + 1
    les_model_registry(les_model_registry_size)%type_name = type_name
    les_model_registry(les_model_registry_size)%allocator => allocator
  end subroutine register_les_model

end submodule les_model_fctry
