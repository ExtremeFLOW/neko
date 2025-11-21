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
  use scalar_aux, only: scalar_step_info
  use mesh, only: mesh_t
  use space, only: space_t
  use gather_scatter, only: gs_t
  use time_scheme_controller, only: time_scheme_controller_t
  use time_step_controller, only: time_step_controller_t
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default, json_extract_item
  use field, only: field_t
  use field_list, only: field_list_t
  use field_series, only: field_series_t
  use field_registry, only: neko_field_registry
  use checkpoint, only: chkp_t
  use krylov, only: ksp_t, ksp_monitor_t
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
     !> Restart from checkpoint data
     procedure :: restart => scalars_restart
     !> Check if the configuration is valid
     procedure :: validate => scalars_validate
     !> Clean up all resources
     procedure :: free => scalars_free
     !> Register scalar lag fields with checkpoint
     procedure, private :: register_lags_with_checkpoint
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

    ! Collect and validate field names for all scalars
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
          if (n_scalars == 1) then
             field_name = 's' ! Single scalar gets default name 's'
          else
             write(field_name, '(A,I0)') 's_', i
          end if
       end if

       field_names(i) = trim(field_name)

       ! If there's a duplicate, append a number until unique
       if (n_scalars > 1) then
          j = 1
          do while (j < i)
             if (trim(field_names(i)) == trim(field_names(j))) then
                write(field_name, '(A,I0)') trim(field_names(i))//'_', j
                field_names(i) = trim(field_name)
                j = 1 ! Start over to check if new name is unique
             else
                j = j + 1
             end if
          end do
       end if
    end do

    do i = 1, n_scalars
       call json_extract_item(params, "", i, json_subdict)

       ! Use the processed field names for all scalars
       call json_subdict%add('name', trim(field_names(i)))

       call this%scalar_fields(i)%init(msh, coef, gs, json_subdict, &
            numerics_params, user, chkp, ulag, vlag, wlag, time_scheme, rho)
    end do

    ! Register all scalar lag fields with checkpoint using scalable approach
    if (n_scalars > 1) then
       call this%register_lags_with_checkpoint(chkp)
    else
       ! For single scalar, use legacy interface
       select type(scalar => this%scalar_fields(1))
       type is (scalar_pnpn_t)
          call chkp%add_scalar(scalar%s, scalar%slag, scalar%abx1, scalar%abx2)
       end select
    end if
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

    ! Set the scalar name to "s"
    if (.not. params%valid_path('name')) then
       call params%add('name', 's')
    end if

    ! Initialize it directly with the params
    call this%scalar_fields(1)%init(msh, coef, gs, params, numerics_params, &
         user, chkp, ulag, vlag, wlag, time_scheme, rho)

    ! Register single scalar with checkpoint
    select type(scalar => this%scalar_fields(1))
    type is (scalar_pnpn_t)
       call chkp%add_scalar(scalar%s, scalar%slag, scalar%abx1, scalar%abx2)
    end select
  end subroutine scalars_init_single

  !> Perform a time step for all scalar fields
  subroutine scalars_step(this, time, ext_bdf, dt_controller)
    class(scalars_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_scheme_controller_t), intent(inout) :: ext_bdf
    type(time_step_controller_t), intent(inout) :: dt_controller
    integer :: i
    type(ksp_monitor_t), dimension(size(this%scalar_fields)) :: ksp_results

    ! Iterate through all scalar fields
    do i = 1, size(this%scalar_fields)
       call this%scalar_fields(i)%step(time, ext_bdf, dt_controller, &
            ksp_results(i))
    end do

    call scalar_step_info(time, ksp_results)
  end subroutine scalars_step

  !> Restart from checkpoint data
  subroutine scalars_restart(this, chkp)
    class(scalars_t), intent(inout) :: this
    type(chkp_t), intent(inout) :: chkp
    integer :: i, n_scalars

    n_scalars = size(this%scalar_fields)
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

    if (allocated(this%shared_ksp)) then
       call this%shared_ksp%free()
       deallocate(this%shared_ksp)
    end if
  end subroutine scalars_free

  !> Register scalar lag fields with checkpoint
  subroutine register_lags_with_checkpoint(this, chkp)
    class(scalars_t), intent(inout) :: this
    type(chkp_t), intent(inout) :: chkp
    integer :: i, n_scalars

    n_scalars = size(this%scalar_fields)

    ! Allocate ABX field arrays
    allocate(chkp%scalar_abx1(n_scalars))
    allocate(chkp%scalar_abx2(n_scalars))

    ! Add all scalar lag fields to the checkpoint list and populate ABX fields
    do i = 1, n_scalars
       call chkp%scalar_lags%append(this%scalar_fields(i)%slag)

       ! Cast to scalar_pnpn_t to access ABX fields
       select type(scalar_field => this%scalar_fields(i))
       type is(scalar_pnpn_t)
          call associate_scalar_abx_fields(chkp, i, scalar_field)
       end select
    end do

  end subroutine register_lags_with_checkpoint

  !> Helper subroutine to associate ABX field pointers with proper TARGET attribute
  subroutine associate_scalar_abx_fields(chkp, index, scalar_field)
    type(chkp_t), intent(inout) :: chkp
    integer, intent(in) :: index
    type(scalar_pnpn_t), target, intent(in) :: scalar_field

    chkp%scalar_abx1(index)%ptr => scalar_field%abx1
    chkp%scalar_abx2(index)%ptr => scalar_field%abx2
  end subroutine associate_scalar_abx_fields

end module scalars
