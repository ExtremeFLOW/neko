! Copyright (c) 2025, The Neko Authors
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
submodule(fluid_scheme_compressible_euler) euler_bc_fctry
  use dirichlet, only : dirichlet_t
  use inflow, only : inflow_t
  use zero_dirichlet, only : zero_dirichlet_t
  use symmetry, only : symmetry_t
  implicit none

  ! List of all possible types created by the boundary condition factories
  character(len=25) :: EULER_KNOWN_BCS(7) = [character(len=25) :: &
       "velocity_value", &
       "density_value", &
       "pressure_value", &
       "no_slip", &
       "symmetry", &
       "outflow", &
       "normal_outflow"]

contains
  !> Factory routine for pressure boundary conditions.
  !! @param object The boundary condition to be allocated.
  !! @param scheme The `fluid_scheme_compressible_euler_t`  scheme.
  !! @param json The parameter dictionary for the boundary.
  !! @param coef The SEM coeffcients.
  !! @param user The user interface.
  module subroutine density_bc_factory(object, scheme, json, coef, user)
    class(bc_t), pointer, intent(inout) :: object
    type(fluid_scheme_compressible_euler_t), intent(in) :: scheme
    type(json_file), intent(inout) :: json
    type(coef_t), intent(in) :: coef
    type(user_t), intent(in) :: user
    character(len=:), allocatable :: type
    integer :: i, j, k
    integer, allocatable :: zone_indices(:)

    call json_get(json, "type", type)

    select case (trim(type))
    case ("density_value")
       allocate(dirichlet_t::object)
    case default
       do i = 1, size(EULER_KNOWN_BCS)
          if (trim(type) .eq. trim(EULER_KNOWN_BCS(i))) return
       end do
       call neko_type_error("compressible_euler boundary conditions", type, &
            EULER_KNOWN_BCS)
    end select

    call json_get(json, "zone_indices", zone_indices)
    call object%init(coef, json)

    do i = 1, size(zone_indices)
       call object%mark_zone(coef%msh%labeled_zones(zone_indices(i)))
    end do
    call object%finalize()
  end subroutine density_bc_factory

  !> Factory routine for pressure boundary conditions.
  !! @param object The boundary condition to be allocated.
  !! @param scheme The `fluid_scheme_compressible_euler_t`  scheme.
  !! @param json The parameter dictionary for the boundary.
  !! @param coef The SEM coeffcients.
  !! @param user The user interface.
  module subroutine pressure_bc_factory(object, scheme, json, coef, user)
    class(bc_t), pointer, intent(inout) :: object
    type(fluid_scheme_compressible_euler_t), intent(inout) :: scheme
    type(json_file), intent(inout) :: json
    type(coef_t), intent(in) :: coef
    type(user_t), intent(in) :: user
    character(len=:), allocatable :: type
    integer :: i, j, k
    integer, allocatable :: zone_indices(:)

    call json_get(json, "type", type)

    select case (trim(type))
    case ("outflow", "normal_outflow")
       allocate(zero_dirichlet_t::object)
    case ("pressure_value")
       allocate(dirichlet_t::object)
    case default
       do i = 1, size(EULER_KNOWN_BCS)
          if (trim(type) .eq. trim(EULER_KNOWN_BCS(i))) return
       end do
       call neko_type_error("compressible_euler boundary conditions", type, &
            EULER_KNOWN_BCS)
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
  !! @param scheme The `fluid_scheme_compressible_euler_t` scheme.
  !! @param json The parameter dictionary for the boundary.
  !! @param coef The SEM coeffcients.
  !! @param user The user interface.
  module subroutine velocity_bc_factory(object, scheme, json, coef, user)
    class(bc_t), pointer, intent(inout) :: object
    type(fluid_scheme_compressible_euler_t), intent(in) :: scheme
    type(json_file), intent(inout) :: json
    type(coef_t), intent(in) :: coef
    type(user_t), intent(in) :: user
    character(len=:), allocatable :: type
    integer :: i, j, k
    integer, allocatable :: zone_indices(:)

    call json_get(json, "type", type)

    select case (trim(type))
    case ("symmetry")
       allocate(symmetry_t::object)
    case ("no_slip")
       allocate(zero_dirichlet_t::object)
    case ("velocity_value")
       allocate(inflow_t::object)
    case default
       do i = 1, size(EULER_KNOWN_BCS)
          if (trim(type) .eq. trim(EULER_KNOWN_BCS(i))) return
       end do
       call neko_type_error("compressible_euler boundary conditions", type, &
            EULER_KNOWN_BCS)
    end select

    call json_get(json, "zone_indices", zone_indices)
    call object%init(coef, json)
    do i = 1, size(zone_indices)
       call object%mark_zone(coef%msh%labeled_zones(zone_indices(i)))
    end do
    call object%finalize()

  end subroutine velocity_bc_factory

end submodule euler_bc_fctry
