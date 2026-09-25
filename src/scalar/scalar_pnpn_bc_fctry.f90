! Copyright (c) 2024-2026, The Neko Authors
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
!> Defines boundary condition factory and allocator routines for
!! `scalar_pnpn_t`.
submodule(scalar_pnpn) scalar_pnpn_bc_fctry
  use dirichlet, only : dirichlet_t
  use neumann, only : neumann_t
  use user_intf, only : user_t
  use utils, only : neko_type_error, neko_type_registration_error
  use field_dirichlet, only : field_dirichlet_t
  use field_neumann, only : field_neumann_t
  use expression_dirichlet, only : expression_dirichlet_t
  use overset_interface, only : overset_interface_t
  implicit none

  ! List of all possible types created by the boundary condition factories
  character(len=25), parameter :: SCALAR_PNPN_KNOWN_BCS(6) = [ &
       character(len=25) :: &
       "dirichlet", &
       "expression_dirichlet", &
       "user_dirichlet", &
       "user_neumann", &
       "neumann", &
       "overset_interface"]

contains

  !> Scalar Pn/Pn boundary condition allocator.
  !! @param[inout] object The object to be allocated.
  !! @param[in] type_name The name of the boundary condition type.
  module subroutine scalar_pnpn_bc_allocator(object, type_name)
    class(bc_t), pointer, intent(inout) :: object
    character(len=*), intent(in) :: type_name
    integer :: i

    if (associated(object)) then
       call object%free()
       deallocate(object)
    end if

    select case (trim(type_name))
    case ("user_dirichlet")
       allocate(field_dirichlet_t::object)
    case ("dirichlet")
       allocate(dirichlet_t::object)
    case ("expression_dirichlet")
       allocate(expression_dirichlet_t::object)
    case ("user_neumann")
       allocate(field_neumann_t::object)
    case ("neumann")
       allocate(neumann_t::object)
    case ("overset_interface")
       allocate(overset_interface_t::object)
    case default
       do i = 1, scalar_pnpn_bc_registry_size
          if (trim(type_name) .eq. &
               trim(scalar_pnpn_bc_registry(i)%type_name)) then
             call scalar_pnpn_bc_registry(i)%allocator(object)
             return
          end if
       end do

       call neko_type_error("scalar_pnpn boundary conditions", type_name, &
            SCALAR_PNPN_KNOWN_BCS)
    end select
  end subroutine scalar_pnpn_bc_allocator

  !> Register a custom scalar Pn/Pn boundary condition allocator.
  !! Called in custom user modules inside the `module_name_register_types`
  !! routine to add a custom type allocator to the registry.
  !! @param[in] type_name The name of the boundary condition type.
  !! @param[in] allocator The allocator for the custom user type.
  module subroutine register_scalar_pnpn_bc(type_name, allocator)
    character(len=*), intent(in) :: type_name
    procedure(scalar_pnpn_bc_allocate), pointer, intent(in) :: allocator
    type(scalar_pnpn_bc_allocator_entry), allocatable :: temp(:)
    integer :: i

    do i = 1, size(SCALAR_PNPN_KNOWN_BCS)
       if (trim(type_name) .eq. trim(SCALAR_PNPN_KNOWN_BCS(i))) then
          call neko_type_registration_error("scalar_pnpn boundary condition", &
               type_name, .true.)
       end if
    end do

    do i = 1, scalar_pnpn_bc_registry_size
       if (trim(type_name) .eq. &
            trim(scalar_pnpn_bc_registry(i)%type_name)) then
          call neko_type_registration_error("scalar_pnpn boundary condition", &
               type_name, .false.)
       end if
    end do

    if (scalar_pnpn_bc_registry_size .eq. 0) then
       allocate(scalar_pnpn_bc_registry(1))
    else
       allocate(temp(scalar_pnpn_bc_registry_size + 1))
       temp(1:scalar_pnpn_bc_registry_size) = scalar_pnpn_bc_registry
       call move_alloc(temp, scalar_pnpn_bc_registry)
    end if

    scalar_pnpn_bc_registry_size = scalar_pnpn_bc_registry_size + 1
    scalar_pnpn_bc_registry(scalar_pnpn_bc_registry_size)%type_name = type_name
    scalar_pnpn_bc_registry(scalar_pnpn_bc_registry_size)%allocator => allocator
  end subroutine register_scalar_pnpn_bc

  !> Boundary condition factory. Both constructs and initializes the object.
  !! Will mark a mesh zone for the bc and finalize.
  !! @param[inout] object The boundary condition to be allocated.
  !! @param[in] scheme The `scalar_pnpn` scheme.
  !! @param[inout] json JSON object for initializing the bc.
  !! @param[in] coef SEM coefficients.
  !! @param[in] user The user interface.
  module subroutine scalar_pnpn_bc_factory(object, scheme, json, coef, user)
    class(bc_t), pointer, intent(inout) :: object
    type(scalar_pnpn_t), intent(in) :: scheme
    type(json_file), intent(inout) :: json
    type(coef_t), target, intent(in) :: coef
    type(user_t), intent(in) :: user
    character(len=:), allocatable :: type
    integer :: i
    integer, allocatable :: zone_indices(:)
    character(len=:), allocatable :: default_name
    character(len=64) :: buf

    call json_get(json, "type", type)
    call scalar_pnpn_bc_allocator(object, type)

    select case (trim(type))
    case ("user_dirichlet")
       select type (obj => object)
       type is (field_dirichlet_t)
          obj%update => user%dirichlet_conditions
          ! Add the name of the dummy field in the bc, matching the scalar
          ! solved for.
          call json%add("field_name", scheme%s%name)
       end select
    case ("user_neumann")
       select type (obj => object)
       type is (field_neumann_t)
          obj%update => user%neumann_conditions
          ! Add the name of the dummy field in the bc, matching the scalar
          ! solved for.
          call json%add("field_name", scheme%s%name)
       end select
    case ("overset_interface")
       select type (obj => object)
       type is (overset_interface_t)
          call json%add("field_name", scheme%s%name)
       end select
    end select

    call json_get(json, "zone_indices", zone_indices)
    call object%init(coef, json)
    do i = 1, size(zone_indices)
       call object%mark_labeled_zone(zone_indices(i))
    end do

    write(buf,'("scalar_bc_",I0)') zone_indices(1)
    default_name = trim(buf)
    call json_get_or_default(json, "name", object%name, default_name)
    object%zone_indices = zone_indices
    call object%finalize()

    deallocate(type)
    deallocate(zone_indices)
    deallocate(default_name)

  end subroutine scalar_pnpn_bc_factory


end submodule scalar_pnpn_bc_fctry
