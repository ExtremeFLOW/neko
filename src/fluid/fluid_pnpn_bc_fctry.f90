
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
!> Defines factory subroutines for `fluid_pnpn_t`.
submodule(fluid_pnpn) fluid_pnpn_bc_fctry
  use user_intf, only : user_t
  use utils, only : neko_type_error
  use field_dirichlet, only : field_dirichlet_t
  use inflow, only : inflow_t
  use blasius, only : blasius_t
  use dirichlet, only : dirichlet_t
  use dong_outflow, only : dong_outflow_t
  use symmetry, only : symmetry_t
  use non_normal, only : non_normal_t
  use field_dirichlet_vector, only : field_dirichlet_vector_t
  implicit none

  ! List of all possible types created by the boundary condition factories
  character(len=25) :: FLUID_PNPN_KNOWN_BCS(14) = [character(len=25) :: &
       "symmetry", &
       "velocity_value", &
       "no_slip", &
       "outflow", &
       "normal_outflow", &
       "outflow+user", &
       "normal_outflow+user", &
       "outflow+dong", &
       "normal_outflow+dong", &
       "shear_stress", &
       "user_velocity", &
       "user_pressure", &
       "blasius_profile", &
       "wall_model"]

contains

  !> Factory routine for pressure boundary conditions.
  !! @param object The boundary condition to be allocated.
  !! @param scheme The `fluid_pnpn_t`  scheme.
  !! @param json The parameter dictionary for the boundary.
  !! @param coef The SEM coeffcients.
  !! @param user The user interface.
  module subroutine pressure_bc_factory(object, scheme, json, coef, user)
    class(bc_t), pointer, intent(inout) :: object
    type(fluid_pnpn_t), intent(in) :: scheme
    type(json_file), intent(inout) :: json
    type(coef_t), target, intent(in) :: coef
    type(user_t), intent(in) :: user
    character(len=:), allocatable :: type
    integer :: i, j, k
    integer, allocatable :: zone_indices(:)

    call json_get(json, "type", type)

    select case (trim(type))
    case ("outflow", "normal_outflow")
       allocate(zero_dirichlet_t::object)

    case ("outflow+dong", "normal_outflow+dong")
       allocate(dong_outflow_t::object)

    case ("user_pressure", "outflow+user", "normal_outflow+user")
       allocate(field_dirichlet_t::object)
       select type (obj => object)
       type is (field_dirichlet_t)
          obj%update => user%dirichlet_conditions
          call json%add("field_name", scheme%p%name)
       end select

    case default
       do i = 1, size(FLUID_PNPN_KNOWN_BCS)
          if (trim(type) .eq. trim(FLUID_PNPN_KNOWN_BCS(i))) return
       end do
       call neko_type_error("fluid_pnpn boundary conditions", type, &
            FLUID_PNPN_KNOWN_BCS)
    end select

    call json_get(json, "zone_indices", zone_indices)
    call object%init(coef, json)

    do i = 1, size(zone_indices)
       call object%mark_zone(coef%msh%labeled_zones(zone_indices(i)))
    end do
    call object%finalize()

    ! All pressure bcs are currently strong, so for all of them we
    ! mark with value 1 in the mesh
    do i = 1, size(zone_indices)
       do j = 1, scheme%msh%nelv
          do k = 1, 2 * scheme%msh%gdim
             if (scheme%msh%facet_type(k,j) .eq. -zone_indices(i)) then
                scheme%msh%facet_type(k, j) = 1
             end if
          end do
       end do
    end do
  end subroutine pressure_bc_factory

  !> Factory routine for velocity boundary conditions.
  !! @param object The boundary condition to be allocated.
  !! @param scheme The `fluid_pnpn_t`  scheme.
  !! @param json The parameter dictionary for the boundary.
  !! @param coef The SEM coeffcients.
  !! @param user The user interface.
  module subroutine velocity_bc_factory(object, scheme, json, coef, user)
    class(bc_t), pointer, intent(inout) :: object
    type(fluid_pnpn_t), intent(inout) :: scheme
    type(json_file), intent(inout) :: json
    type(coef_t), target, intent(in) :: coef
    type(user_t), intent(in) :: user
    character(len=:), allocatable :: type
    integer :: i, j, k
    integer, allocatable :: zone_indices(:)

    call json_get(json, "type", type)

    select case (trim(type))
    case ("symmetry")
       allocate(symmetry_t::object)
    case ("velocity_value")
       allocate(inflow_t::object)
    case ("no_slip")
       allocate(zero_dirichlet_t::object)
    case ("normal_outflow", "normal_outflow+dong", "normal_outflow+user")
       allocate(non_normal_t::object)
    case ("blasius_profile")
       allocate(blasius_t::object)
    case ("shear_stress")
       allocate(shear_stress_t::object)
    case ("wall_model")
       allocate(wall_model_bc_t::object)
       ! Kind of hack, but  OK for now
       call json%add("scheme_name", scheme%name)

    case ("user_velocity")
       allocate(field_dirichlet_vector_t::object)
       select type (obj => object)
       type is (field_dirichlet_vector_t)
          obj%update => user%dirichlet_conditions
       end select

    case default
       do i = 1, size(FLUID_PNPN_KNOWN_BCS)
          if (trim(type) .eq. trim(FLUID_PNPN_KNOWN_BCS(i))) return
       end do
       call neko_type_error("fluid_pnpn boundary conditions", type, &
            FLUID_PNPN_KNOWN_BCS)
    end select

    call json_get(json, "zone_indices", zone_indices)
    call object%init(coef, json)
    do i = 1, size(zone_indices)
       call object%mark_zone(coef%msh%labeled_zones(zone_indices(i)))
    end do
    call object%finalize()

    ! Exclude these two because they are bcs for the residual, not velocity
    if (trim(type) .ne. "normal_outflow" .and. &
         trim(type) .ne. "normal_outflow+dong") then
       do i = 1, size(zone_indices)
          do j = 1, scheme%msh%nelv
             do k = 1, 2 * scheme%msh%gdim
                if (scheme%msh%facet_type(k,j) .eq. -zone_indices(i)) then
                   scheme%msh%facet_type(k, j) = 2
                end if
             end do
          end do
       end do
    end if
  end subroutine velocity_bc_factory

end submodule fluid_pnpn_bc_fctry
