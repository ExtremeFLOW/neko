! Copyright (c) 2021-2024, The Neko Authors
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
submodule (wall_model) wall_model_fctry
  use vreman, only : vreman_t
  use spalding, only : spalding_t
  use rough_log_law, only : rough_log_law_t
  use utils, only : neko_type_error
  use utils, only : neko_type_registration_error
  implicit none

  ! List of all possible types created by the factory routine
  character(len=20) :: WALLM_KNOWN_TYPES(2) = [character(len=20) :: &
       "spalding", &
       "rough_log_law"]

contains

  !> Wall model factory. Both constructs and initializes the object.
  !! @param object The object to be allocated.
  !! @param scheme_name The name of the scheme for which the wall model is used.
  !! @param coef SEM coefficients.
  !! @param msk The boundary mask.
  !! @param facet The boundary facets.
  !! @param h_index The off-wall index of the sampling cell.
  !! @param json A dictionary with parameters.
  module subroutine wall_model_factory(object, scheme_name, coef, msk, facet, &
       json)
    class(wall_model_t), allocatable, intent(inout) :: object
    character(len=*), intent(in) :: scheme_name
    type(coef_t), intent(in) :: coef
    integer, intent(in) :: msk(:)
    integer, intent(in) :: facet(:)
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: type_name
    character(len=:), allocatable :: type_string
    integer :: h_index

    call json_get(json, "model", type_name)
    call json_get(json, "h_index", h_index)

    call wall_model_allocator(object, type_name)

    ! Initialize
    call object%init(scheme_name, coef, msk, facet, h_index, json)

  end subroutine wall_model_factory

  !> Wall model allocator.
  !! @param object The object to be allocated.
  !! @param type_name The name of the type to allocate.
  module subroutine wall_model_allocator(object, type_name)
    class(wall_model_t), allocatable, intent(inout) :: object
    character(len=:), allocatable, intent(in) :: type_name
    integer :: i

    select case (trim(type_name) )
    case ("spalding")
       allocate(spalding_t::object)
    case ("rough_log_law")
       allocate(rough_log_law_t::object)
    case default
       do i = 1, wall_model_registry_size
          if (trim(type_name) .eq. trim(wall_model_registry(i)%type_name)) then
             call wall_model_registry(i)%allocator(object)
             return
          end if
       end do
       call neko_type_error("wall model", trim(type_name), WALLM_KNOWN_TYPES)
    end select
  end subroutine wall_model_allocator

  !> Register a custom wall model allocator.
  !! Called in custom user modules inside the `module_name_register_types`
  !! routine to add a custom type allocator to the registry.
  !! @param type_name The name of the type to allocate.
  !! @param allocator The allocator for the custom user type.
  module subroutine register_wall_model(type_name, allocator)
    character(len=*), intent(in) :: type_name
    procedure(wall_model_allocate), pointer, intent(in) :: allocator
    type(allocator_entry), allocatable :: temp(:)
    integer :: i

    do i = 1, size(WALLM_KNOWN_TYPES)
       if (trim(type_name) .eq. trim(WALLM_KNOWN_TYPES(i))) then
          call neko_type_registration_error("wall model", type_name, .true.)
       end if
    end do

    do i = 1, wall_model_registry_size
       if (trim(type_name) .eq. trim(wall_model_registry(i)%type_name)) then
          call neko_type_registration_error("wall model", type_name, .false.)
       end if
    end do

    ! Expand registry
    if (wall_model_registry_size .eq. 0) then
       allocate(wall_model_registry(1))
    else
       allocate(temp(wall_model_registry_size + 1))
       temp(1:wall_model_registry_size) = wall_model_registry
       call move_alloc(temp, wall_model_registry)
    end if

    wall_model_registry_size = wall_model_registry_size + 1
    wall_model_registry(wall_model_registry_size)%type_name = type_name
    wall_model_registry(wall_model_registry_size)%allocator => allocator
  end subroutine register_wall_model

end submodule wall_model_fctry
