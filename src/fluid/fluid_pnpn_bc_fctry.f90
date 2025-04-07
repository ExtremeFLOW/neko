
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
  use logger, only: neko_log, LOG_SIZE, NEKO_LOG_INFO
  implicit none

  ! List of all possible types created by the boundary condition factories
  character(len=25) :: FLUID_PNPN_KNOWN_BCS(12) = [character(len=25) :: &
       "symmetry", &
       "velocity_value", &
       "no_slip", &
       "outflow", &
       "normal_outflow", &
       "outflow+dong", &
       "normal_outflow+dong", &
       "shear_stress", &
       "user_dirichlet", &
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
    integer :: i, j, k
    integer, allocatable :: zone_indices(:)

    ! For user_dirichlet
    character(len=20), allocatable :: user_dirichlet_components(:)
    logical :: pressure_found = .false.

    call json_get(json, "type", type)

    select case (trim(type))
    case ("outflow", "normal_outflow")
       allocate(zero_dirichlet_t::object)

    case ("outflow+dong", "normal_outflow+dong")
       allocate(dong_outflow_t::object)

    case ("user_dirichlet")
       ! Check if "p" is specified in the dirichlet components
       call json_get(json, "fields", user_dirichlet_components)
       do i = 1, size(user_dirichlet_components)
          if (trim(user_dirichlet_components(i)) .eq. "p") then
             print *, "found p"
             allocate(field_dirichlet_t::object)
             select type (obj => object)
             type is (field_dirichlet_t)
                obj%update => user%user_dirichlet_update
                call json%add("field_name", scheme%p%name)
             end select
             pressure_found = .true.
             exit
          end if

          if (.not. pressure_found) return

       end do

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
    type(fluid_pnpn_t), intent(in) :: scheme
    type(json_file), intent(inout) :: json
    type(coef_t), intent(in) :: coef
    type(user_t), intent(in) :: user
    character(len=:), allocatable :: type
    integer :: i, j, k
    integer, allocatable :: zone_indices(:)

    ! For user_dirichlet
    character(len=20), allocatable :: user_dirichlet_fields(:)
    integer :: n_user_dirichlet_comps, idx
    character(len=20), allocatable :: user_dirichlet_vel_comps(:)
    character(len=LOG_SIZE) :: log_buf

    call json_get(json, "type", type)

    call neko_log%message("BC is " // trim(type))

    select case (trim(type))
    case ("symmetry")
       allocate(symmetry_t::object)
    case ("velocity_value")
       allocate(inflow_t::object)
    case ("no_slip")
       allocate(zero_dirichlet_t::object)
    case ("normal_outflow", "normal_outflow+dong")
       allocate(non_normal_t::object)
    case ("blasius_profile")
       allocate(blasius_t::object)
    case ("shear_stress")
       allocate(shear_stress_t::object)
    case ("wall_model")
       allocate(wall_model_bc_t::object)
       ! Kind of hack, but maybe OK? The thing is, we need the nu for
       ! initing the wall model, and forcing the user duplicate that there
       ! would be a nightmare.
       call json%add("nu", scheme%mu / scheme%rho)

    case ("user_dirichlet")

       call json_get(json, "fields", user_dirichlet_fields)

       ! Check if any velocity fields are specified and if so, how many
       n_user_dirichlet_comps = 0
       do i = 1, size(user_dirichlet_fields)
          if ( trim(user_dirichlet_fields(i)) .eq. "u" .or. &
               trim(user_dirichlet_fields(i)) .eq. "v" .or. &
               trim(user_dirichlet_fields(i)) .eq. "w") &
               n_user_dirichlet_comps = n_user_dirichlet_comps + 1
       end do

       call neko_log%section("(F-D, pnpn)", lvl=NEKO_LOG_INFO)
       write (log_buf, *) "Found", n_user_dirichlet_comps, "components"
       call neko_log%message(log_buf, lvl=NEKO_LOG_INFO)

       ! If there are velocity components...
       if (n_user_dirichlet_comps .gt. 0) then

          ! ... duplicate only those velocity components, and inject the
          ! modified array in the json, which will be picked
          ! up by the field_dirichlet_vector_t constructor.
          allocate(user_dirichlet_vel_comps(n_user_dirichlet_comps))

          idx = 1
          do i = 1, size(user_dirichlet_fields)
             if ( trim(user_dirichlet_fields(i)) .eq. "u" .or. &
                  trim(user_dirichlet_fields(i)) .eq. "v" .or. &
                  trim(user_dirichlet_fields(i)) .eq. "w") then
                call neko_log%message("found component " // trim(user_dirichlet_fields(i)), lvl=NEKO_LOG_INFO)
                user_dirichlet_vel_comps(idx) = user_dirichlet_fields(i)
                idx = idx + 1
             end if
          end do

          ! the constructor of field_dirichlet_vector_t will search for
          ! the "velocity_components" array
          call neko_log%message("Adding new json array", lvl=NEKO_LOG_INFO)
          call json%add("velocity_components", user_dirichlet_vel_comps)

          allocate(field_dirichlet_vector_t::object)
          select type (obj => object)
          type is (field_dirichlet_vector_t)
             obj%update => user%user_dirichlet_update
          end select
       else
          return ! is this wonky? This simulates the fact that if
                 ! we don't specify any velocity components we assume
                 ! we want a "do nothing". In the case default this is
                 ! what is done so...
       end if

    call neko_log%end_section(lvl=NEKO_LOG_INFO)

    case ("user_velocity_pointwise")
       allocate(usr_inflow_t::object)
       select type (obj => object)
       type is (usr_inflow_t)
          call obj%set_eval(user%fluid_user_if)
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
