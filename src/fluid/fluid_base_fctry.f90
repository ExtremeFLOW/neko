! Copyright (c) 2025-2026, The Neko Authors
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
!> Factory for all fluid schemes
submodule (fluid_scheme_base) fluid_base_fctry
  use fluid_pnpn, only : fluid_pnpn_t
  use fluid_scheme_compressible_ns, only : fluid_scheme_compressible_ns_t
  use utils, only : neko_type_error, neko_type_registration_error
  implicit none

  ! List of all possible types created by the factory routine
  character(len=20) :: FLUID_KNOWN_TYPES(2) = [character(len=20) :: &
       "pnpn", "compressible"]

contains

  !> Fluid scheme factory. Both constructs and initializes the object.
  !! @param object The object to be created and initialized.
  !! @param msh The mesh.
  !! @param lx The number of GLL points in each direction.
  !! @param params The fluid parameters, i.e. `case.fluid`.
  !! @param numerics_params The numerical parameters, i.e. `case.numerics`.
  !! @param user Type with user-defined procedures.
  !! @param chkp Checkpoint for restarts.
  !! @param output_boundary Whether to write the boundary condition markers
  !! to file.
  module subroutine fluid_scheme_factory(object, msh, lx, params, &
       numerics_params, user, chkp, output_boundary)
    class(fluid_scheme_base_t), allocatable, intent(inout) :: object
    type(mesh_t), target, intent(inout) :: msh
    integer, intent(in) :: lx
    type(json_file), target, intent(inout) :: params
    type(json_file), target, intent(inout) :: numerics_params
    type(user_t), target, intent(in) :: user
    type(chkp_t), target, intent(inout) :: chkp
    logical, intent(in) :: output_boundary
    character(len=:), allocatable :: type_name

    call json_get(params, "scheme", type_name)

    call fluid_scheme_allocator(object, type_name)

    call object%init(msh, lx, params, numerics_params, user, chkp, &
         output_boundary)

  end subroutine fluid_scheme_factory

  !> Fluid scheme allocator.
  !! @param object The object to be allocated.
  !! @param type_name The name of the fluid scheme type.
  module subroutine fluid_scheme_allocator(object, type_name)
    class(fluid_scheme_base_t), allocatable, intent(inout) :: object
    character(len=*), intent(in) :: type_name
    integer :: i

    if (allocated(object)) then
       call object%free()
       deallocate(object)
    end if

    select case (trim(type_name))
    case ("pnpn")
       allocate(fluid_pnpn_t::object)
    case ("compressible")
       allocate(fluid_scheme_compressible_ns_t::object)
    case default
       do i = 1, fluid_scheme_registry_size
          if (trim(type_name) == &
               trim(fluid_scheme_registry(i)%type_name)) then
             call fluid_scheme_registry(i)%allocator(object)
             return
          end if
       end do
       call neko_type_error("fluid scheme", trim(type_name), &
            FLUID_KNOWN_TYPES)
    end select

  end subroutine fluid_scheme_allocator

  !> Register a custom fluid scheme allocator.
  !! Called in custom user modules inside the `module_name_register_types`
  !! routine to add a custom type allocator to the registry.
  !! @param type_name The name of the type to allocate.
  !! @param allocator The allocator for the custom user type.
  module subroutine register_fluid_scheme(type_name, allocator)
    character(len=*), intent(in) :: type_name
    procedure(fluid_scheme_allocate), pointer, intent(in) :: allocator
    type(fluid_scheme_allocator_entry), allocatable :: temp(:)
    integer :: i

    do i = 1, size(FLUID_KNOWN_TYPES)
       if (trim(type_name) .eq. trim(FLUID_KNOWN_TYPES(i))) then
          call neko_type_registration_error("fluid scheme", type_name, &
               .true.)
       end if
    end do

    do i = 1, fluid_scheme_registry_size
       if (trim(type_name) .eq. &
            trim(fluid_scheme_registry(i)%type_name)) then
          call neko_type_registration_error("fluid scheme", type_name, &
               .false.)
       end if
    end do

    ! Expand registry
    if (fluid_scheme_registry_size == 0) then
       allocate(fluid_scheme_registry(1))
    else
       allocate(temp(fluid_scheme_registry_size + 1))
       temp(1:fluid_scheme_registry_size) = fluid_scheme_registry
       call move_alloc(temp, fluid_scheme_registry)
    end if

    fluid_scheme_registry_size = fluid_scheme_registry_size + 1
    fluid_scheme_registry(fluid_scheme_registry_size)%type_name = type_name
    fluid_scheme_registry(fluid_scheme_registry_size)%allocator => allocator
  end subroutine register_fluid_scheme

end submodule fluid_base_fctry
