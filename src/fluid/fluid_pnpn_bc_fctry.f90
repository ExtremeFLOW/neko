
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
  use usr_inflow, only : usr_inflow_t, usr_inflow_eval
  use blasius, only : blasius_t
  use dirichlet, only : dirichlet_t
  use dong_outflow, only : dong_outflow_t
  use symmetry, only : symmetry_t
  use non_normal, only : non_normal_t
  use field_dirichlet_vector, only : field_dirichlet_vector_t
  implicit none

  ! List of all possible types created by the boundary condition factories
  character(len=25) :: FLUID_PNPN_KNOWN_BCS(13) = [character(len=25) :: &
     "symmetry", &
     "velocity_value", &
     "no_slip", &
     "outflow", &
     "normal_outflow", &
     "outflow+dong", &
     "normal_outflow+dong", &
     "shear_stress", &
     "user_velocity", &
     "user_pressure", &
     "blasius_profile", &
     "user_velocity_pointwise", &
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
    type(coef_t), intent(in) :: coef
    type(user_t), intent(in) :: user
    character(len=:), allocatable :: type
    integer :: zone_index, i, j, k

    call json_get(json, "type", type)

    if ( (trim(type) .eq. "outflow") .or. &
         (trim(type) .eq. "normal_outflow"))  then
       allocate(zero_dirichlet_t::object)
    else if ((trim(type) .eq. "outflow+dong") .or. &
            (trim(type) .eq. "normal_outflow+dong")) then
       allocate(dong_outflow_t::object)
    else if (trim(type) .eq. "user_pressure") then
       allocate(field_dirichlet_t::object)
       select type(obj => object)
       type is(field_dirichlet_t)
          obj%update => user%user_dirichlet_update
          call json%add("field_name", scheme%p%name)
       end select
    else
      do i=1, size(FLUID_PNPN_KNOWN_BCS)
         if (trim(type) .eq. trim(FLUID_PNPN_KNOWN_BCS(i))) return
      end do
      call neko_type_error("fluid_pnpn boundary conditions", type, &
           FLUID_PNPN_KNOWN_BCS)
    end if

    call json_get(json, "zone_index", zone_index)
    call object%init(coef, json)
    call object%mark_zone(coef%msh%labeled_zones(zone_index))
    call object%finalize()

    ! All pressure bcs are currently strong, so for all of them we
    ! mark with value 1 in the mesh
    do j = 1, scheme%msh%nelv
          do k = 1, 2 * scheme%msh%gdim
             if (scheme%msh%facet_type(k,j) .eq. -zone_index) then
                scheme%msh%facet_type(k, j) = 1
          end if
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
    type(fluid_pnpn_t), intent(in) :: scheme
    type(json_file), intent(inout) :: json
    type(coef_t), intent(in) :: coef
    type(user_t), intent(in) :: user
    character(len=:), allocatable :: type
    integer :: zone_index, i, j, k

    call json_get(json, "type", type)

    if (trim(type) .eq. "symmetry") then
       allocate(symmetry_t::object)
    else if (trim(type) .eq. "velocity_value") then
       allocate(inflow_t::object)
    else if (trim(type) .eq. "no_slip") then
       allocate(zero_dirichlet_t::object)
    else if (trim(type) .eq. "normal_outflow" .or. &
             trim(type) .eq. "normal_outflow+dong") then
       allocate(non_normal_t::object)
    else if (trim(type) .eq. "blasius_profile") then
       allocate(blasius_t::object)
    else if (trim(type) .eq. "shear_stress") then
       allocate(shear_stress_t::object)
    else if (trim(type) .eq. "wall_model") then
       allocate(wall_model_bc_t::object)
       ! Kind of hack, but maybe OK? The thing is, we need the nu for
       ! initing the wall model, and forcing the user duplicate that there
       ! would be a nightmare.
       call json%add("nu", scheme%mu / scheme%rho)
    else if (trim(type) .eq. "user_velocity") then
       allocate(field_dirichlet_vector_t::object)
       select type(obj => object)
       type is(field_dirichlet_vector_t)
          obj%update => user%user_dirichlet_update
       end select
    else if (trim(type) .eq. "user_velocity_pointwise") then
       allocate(usr_inflow_t::object)
       select type(obj => object)
       type is(usr_inflow_t)
          call obj%set_eval(user%fluid_user_if)
       end select
    else
      do i=1, size(FLUID_PNPN_KNOWN_BCS)
         if (trim(type) .eq. trim(FLUID_PNPN_KNOWN_BCS(i))) return
      end do
      call neko_type_error("fluid_pnpn boundary conditions", type, &
           FLUID_PNPN_KNOWN_BCS)
    end if

    call json_get(json, "zone_index", zone_index)
    call object%init(coef, json)
    call object%mark_zone(coef%msh%labeled_zones(zone_index))
    call object%finalize()

    if (type .ne. "normal_outflow" .and. type .ne. "normal_outflow+dong") then
       do j = 1, scheme%msh%nelv
          do k = 1, 2 * scheme%msh%gdim
             if (scheme%msh%facet_type(k,j) .eq. -zone_index) then
                 scheme%msh%facet_type(k, j) = 2
             end if
          end do
       end do
    end if
  end subroutine velocity_bc_factory

end submodule fluid_pnpn_bc_fctry
