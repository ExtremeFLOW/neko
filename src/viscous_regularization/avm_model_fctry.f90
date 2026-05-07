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
submodule (avm_model) avm_model_fctry
  use case, only : case_t
  use entropy_viscosity, only : entropy_viscosity_t
  use utils, only : neko_type_error, neko_type_registration_error
  implicit none

  ! List of all possible types created by the factory routine
  character(len=20) :: AVM_KNOWN_TYPES(1) = [character(len=20) :: &
       "entropy_viscosity"]

contains
  !> avm model factory.
  !! @param object The object to be allocated.
  !! @param type_name The name of the avm model.
  !! @param case The case_t object.
  !! @param json A dictionary with parameters.
  module subroutine avm_model_factory(object, type_name, case, json)
    class(avm_model_t), allocatable, intent(inout) :: object
    character(len=*), intent(in) :: type_name
    class(case_t), intent(inout), target :: case
    type(json_file), intent(inout) :: json
    character(len=:), allocatable :: type_string

    call avm_model_allocator(object, type_name)
    call object%init(case, json)
  end subroutine avm_model_factory

  !> avm model allocator.
  !! @param object The object to be allocated.
  !! @param type_name The name of the avm model.
  module subroutine avm_model_allocator(object, type_name)
    class(avm_model_t), allocatable, intent(inout) :: object
    character(len=*), intent(in) :: type_name
    integer :: i

    if (allocated(object)) then
       call object%free()
       deallocate(object)
    end if

    select case (trim(type_name))
    case ('entropy_viscosity')
       allocate(entropy_viscosity_t::object)
    case default
       do i = 1, avm_model_registry_size
          if (trim(type_name) == trim(avm_model_registry(i)%type_name)) then
             call avm_model_registry(i)%allocator(object)
             return
          end if
       end do

       call neko_type_error("avm model", type_name, AVM_KNOWN_TYPES)
    end select

  end subroutine avm_model_allocator

  !> Register a custom avm model allocator.
  !! Called in custom user moduavm inside the `module_name_register_types`
  !! routine to add a custom type allocator to the registry.
  !! @param type_name The name of the type to allocate.
  !! @param allocator The allocator for the custom user type.
  module subroutine register_avm_model(type_name, allocator)
    character(len=*), intent(in) :: type_name
    procedure(avm_model_allocate), pointer, intent(in) :: allocator
    type(allocator_entry), allocatable :: temp(:)
    integer :: i

    do i = 1, size(AVM_KNOWN_TYPES)
       if (trim(type_name) .eq. trim(AVM_KNOWN_TYPES(i))) then
          call neko_type_registration_error("avm model", type_name, .true.)
       end if
    end do

    do i = 1, avm_model_registry_size
       if (trim(type_name) .eq. trim(avm_model_registry(i)%type_name)) then
          call neko_type_registration_error("avm model", type_name, .false.)
       end if
    end do

    ! Expand registry
    if (avm_model_registry_size == 0) then
       allocate(avm_model_registry(1))
    else
       allocate(temp(avm_model_registry_size + 1))
       temp(1:avm_model_registry_size) = avm_model_registry
       call move_alloc(temp, avm_model_registry)
    end if

    avm_model_registry_size = avm_model_registry_size + 1
    avm_model_registry(avm_model_registry_size)%type_name = type_name
    avm_model_registry(avm_model_registry_size)%allocator => allocator
  end subroutine register_avm_model

end submodule avm_model_fctry
