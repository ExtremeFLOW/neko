! Copyright (c) 2024, The Neko Authors
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
!> Defines a factory subroutine for `scalar_pnpn_t`.
submodule(scalar_pnpn) scalar_pnpn_bc_fctry
  use dirichlet, only : dirichlet_t
  use neumann, only : neumann_t
  use user_intf, only : user_t
  use utils, only : neko_type_error
  use field_dirichlet, only : field_dirichlet_t
  implicit none

  ! List of all possible types created by the boundary condition factories
  character(len=25) :: SCALAR_PNPN_KNOWN_BCS(3) = [character(len=25) :: &
       "dirichlet", &
       "user", &
       "neumann"]

contains

  !> Boundary condition factory. Both constructs and initializes the object.
  !! Will mark a mesh zone for the bc and finalize.
  !! @param[object] object The boundary condition to be allocated.
  !! @param[in] scheme The `scalar_pnpn` scheme.
  !! @param[inout] json JSON object for initializing the bc.
  !! @param[in] coef SEM coefficients.
  !! @param[in] user The user interface.
  module subroutine bc_factory(object, scheme, json, coef, user)
    class(bc_t), pointer, intent(inout) :: object
    type(scalar_pnpn_t), intent(in) :: scheme
    type(json_file), intent(inout) :: json
    type(coef_t), intent(in) :: coef
    type(user_t), intent(in) :: user
    character(len=:), allocatable :: type
    integer :: i
    integer, allocatable :: zone_indices(:)

    call json_get(json, "type", type)

    select case (trim(type))
    case ("user")
       allocate(field_dirichlet_t::object)
       select type (obj => object)
       type is (field_dirichlet_t)
          obj%update => user%dirichlet_conditions
          ! Add the name of the dummy field in the bc, matching the scalar
          ! solved for.
          call json%add("field_name", scheme%s%name)
       end select
    case ("dirichlet")
       allocate(dirichlet_t::object)
    case ("neumann")
       allocate(neumann_t::object)
    case default
       call neko_type_error("scalar_pnpn boundary conditions", type, &
            SCALAR_PNPN_KNOWN_BCS)
    end select

    call json_get(json, "zone_indices", zone_indices)
    call object%init(coef, json)
    do i = 1, size(zone_indices)
       call object%mark_zone(coef%msh%labeled_zones(zone_indices(i)))
    end do
    call object%finalize()

  end subroutine bc_factory


end submodule scalar_pnpn_bc_fctry
