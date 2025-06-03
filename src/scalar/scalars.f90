! Copyright (c) 2022-2025, The Neko Authors
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
!> Contains the scalars_t type that manages multiple scalar fields.

module scalars
  use num_types, only: rp
  use scalar_pnpn, only: scalar_pnpn_t
  use scalar_scheme, only: scalar_scheme_t
  use mesh, only: mesh_t
  use space, only: space_t
  use gather_scatter, only: gs_t
  use time_scheme_controller, only: time_scheme_controller_t
  use time_step_controller, only: time_step_controller_t
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default, json_extract_object, &
                        json_extract_item
  use field, only: field_t
  use field_series, only: field_series_t
  use field_registry, only: neko_field_registry
  use checkpoint, only: chkp_t
  use krylov, only: ksp_t
  use logger, only: neko_log, LOG_SIZE, NEKO_LOG_VERBOSE
  use user_intf, only: user_t
  use utils, only: neko_error
  use coefs, only : coef_t
  use time_state, only : time_state_t
  implicit none
  private

  !> Type to manage multiple scalar transport equations
  type, public :: scalars_t
     !> The scalar fields
     class(scalar_scheme_t), allocatable :: scalar_fields(:)
     !> Shared KSP solver for all scalar fields
     class(ksp_t), allocatable :: shared_ksp
   contains
     !> Initialize the scalars container
     generic :: init => scalars_init, scalars_init_single
     procedure, private :: scalars_init
     procedure, private :: scalars_init_single
     !> Perform a time step for all scalar fields
     procedure :: step => scalars_step
     !> Update the material properties for all scalar fields
     procedure :: update_material_properties => &
          scalars_update_material_properties
     !> Restart from checkpoint data
     procedure :: restart => scalars_restart
     !> Check if the configuration is valid
     procedure :: validate => scalars_validate
     !> Clean up all resources
     procedure :: free => scalars_free
  end type scalars_t

contains

  !> Initialize the scalars container
  subroutine scalars_init(this, n_scalars, msh, coef, gs, params, &
       numerics_params, user, chkp, ulag, vlag, wlag, time_scheme, rho)
    class(scalars_t), intent(inout) :: this
    integer, intent(in) :: n_scalars
    type(mesh_t), target, intent(in) :: msh
    type(coef_t), target, intent(in) :: coef
    type(gs_t), target, intent(inout) :: gs
    type(json_file), target, intent(inout) :: params
    type(json_file), target, intent(inout) :: numerics_params
    type(user_t), target, intent(in) :: user
    type(field_series_t), target, intent(in) :: ulag, vlag, wlag
    type(time_scheme_controller_t), target, intent(in) :: time_scheme
    TYPE(field_t), TARGET, INTENT(IN) :: rho
    type(chkp_t), target, intent(inout) :: chkp
    type(json_file) :: json_subdict
    integer :: i, j
    character(len=:), allocatable :: field_name
    character(len=:), allocatable :: field_names(:)
    character(len=256) :: error_msg

    ! Allocate the scalar fields
    ! If there are more scalar_scheme_t types, add a factory function here
    allocate(scalar_pnpn_t::this%scalar_fields(n_scalars))

    ! For multiple scalars, collect and validate field names
    if (n_scalars > 1) then
       allocate(character(len=256) :: field_names(n_scalars))

       do i = 1, n_scalars
          ! Extract element i from the "scalars" array
          call json_extract_item(params, "", i, json_subdict)

          ! Try to get name from JSON, generate one if not found or empty
          if (json_subdict%valid_path('name')) then
             call json_get(json_subdict, 'name', field_name)
          else
             field_name = ''
          end if

          ! If name is empty or not provided, generate a default one
          if (len_trim(field_name) == 0) then
             write(field_name, '(A,I0)') 'scalar_', i
          end if

          field_names(i) = trim(field_name)

          ! If there's a duplicate, append a number until unique
          j = 1
          do while (j < i)
             if (trim(field_names(i)) == trim(field_names(j))) then
                write(field_name, '(A,I0)') trim(field_names(i))//'_', j
                field_names(i) = trim(field_name)
                j = 1  ! Start over to check if new name is unique
             else
                j = j + 1
             end if
          end do
       end do
    end if

    do i = 1, n_scalars
       call json_extract_item(params, "", i, json_subdict)
       call this%scalar_fields(i)%init(msh, coef, gs, json_subdict, numerics_params, &
            user, chkp, ulag, vlag, wlag, time_scheme, rho)
    end do
  end subroutine scalars_init

  subroutine scalars_init_single(this, msh, coef, gs, params, numerics_params, &
       user, chkp, ulag, vlag, wlag, time_scheme, rho)
    class(scalars_t), intent(inout) :: this
    type(mesh_t), target, intent(in) :: msh
    type(coef_t), target, intent(in) :: coef
    type(gs_t), target, intent(inout) :: gs
    type(json_file), target, intent(inout) :: params
    type(json_file), target, intent(inout) :: numerics_params
    type(user_t), target, intent(in) :: user
    type(chkp_t), target, intent(inout) :: chkp
    type(field_series_t), target, intent(in) :: ulag, vlag, wlag
    type(time_scheme_controller_t), target, intent(in) :: time_scheme
    TYPE(field_t), TARGET, INTENT(IN) :: rho

    ! Allocate a single scalar field
    allocate(scalar_pnpn_t::this%scalar_fields(1))

    ! Initialize it directly with the params
    call this%scalar_fields(1)%init(msh, coef, gs, params, numerics_params, user, &
         chkp, ulag, vlag, wlag, time_scheme, rho)
  end subroutine scalars_init_single

  !> Perform a time step for all scalar fields
  subroutine scalars_step(this, time, ext_bdf, dt_controller)
    class(scalars_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_scheme_controller_t), intent(inout) :: ext_bdf
    type(time_step_controller_t), intent(inout) :: dt_controller
    integer :: i

    ! Iterate through all scalar fields
    do i = 1, size(this%scalar_fields)
       call this%scalar_fields(i)%step(time, ext_bdf, dt_controller)
    end do
  end subroutine scalars_step

  !> Update the material properties for all scalar fields
  subroutine scalars_update_material_properties(this, t, tstep)
    class(scalars_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    integer :: i

    ! Iterate through all scalar fields
    do i = 1, size(this%scalar_fields)
       this%scalar_fields(i)%cp = 1.0_rp
       this%scalar_fields(i)%lambda = 1e-16_rp
       call this%scalar_fields(i)%update_material_properties(t, tstep)
    end do
  end subroutine scalars_update_material_properties

  !> Restart from checkpoint data
  subroutine scalars_restart(this, chkp)
    class(scalars_t), intent(inout) :: this
    type(chkp_t), intent(inout) :: chkp
    integer :: i
    ! Iterate through all scalar fields
    do i = 1, size(this%scalar_fields)
       call this%scalar_fields(i)%restart(chkp)
    end do
  end subroutine scalars_restart

  !> Check if the configuration is valid
  subroutine scalars_validate(this)
    class(scalars_t), intent(inout) :: this
    integer :: i
    ! Iterate through all scalar fields
    do i = 1, size(this%scalar_fields)
       call this%scalar_fields(i)%slag%set(this%scalar_fields(i)%s)
       call this%scalar_fields(i)%validate()
    end do
  end subroutine scalars_validate

  !> Clean up all resources
  subroutine scalars_free(this)
    class(scalars_t), intent(inout) :: this
    integer :: i

    ! Iterate through all scalar fields
    if (allocated(this%scalar_fields)) then
       do i = 1, size(this%scalar_fields)
          call this%scalar_fields(i)%free()
       end do
       deallocate(this%scalar_fields)
    end if
  end subroutine scalars_free

end module scalars
