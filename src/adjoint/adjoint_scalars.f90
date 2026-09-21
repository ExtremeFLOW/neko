!> @file adjoint_scalars.f90
!! @copyright
!! Copyright (c) 2025, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!
!> Contains the adjoint_scalars_t type that manages multiple scalar fields.

module adjoint_scalars
  use num_types, only: rp
  use adjoint_scalar_pnpn, only: adjoint_scalar_pnpn_t
  use adjoint_scalar_scheme, only: adjoint_scalar_scheme_t
  use mesh, only: mesh_t
  use space, only: space_t
  use gather_scatter, only: gs_t
  use time_scheme_controller, only: time_scheme_controller_t
  use time_step_controller, only: time_step_controller_t
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default, json_extract_item
  use field, only: field_t
  use field_series, only: field_series_t
  use scalar_aux, only : scalar_step_info
  use krylov, only : ksp_monitor_t
  use registry, only: neko_registry
  use checkpoint, only: chkp_t
  use krylov, only: ksp_t
  use user_intf, only: user_t
  use utils, only: neko_error
  use coefs, only : coef_t
  use time_state, only : time_state_t
  implicit none
  private

  !> Type to manage multiple adjoint scalar transport equations
  type, public :: adjoint_scalars_t
     !> The scalar fields
     class(adjoint_scalar_scheme_t), allocatable :: adjoint_scalar_fields(:)
     !> Shared KSP solver for all scalar fields
     class(ksp_t), allocatable :: shared_ksp
   contains
     !> Initialize the adjoint_scalars container
     generic :: init => adjoint_scalars_init, adjoint_scalars_init_single
     procedure, private :: adjoint_scalars_init
     procedure, private :: adjoint_scalars_init_single
     !> Perform a time step for all scalar fields
     procedure :: step => adjoint_scalars_step
     !> Restart from checkpoint data
     procedure :: restart => adjoint_scalars_restart
     !> Check if the configuration is valid
     procedure :: validate => adjoint_scalars_validate
     !> Clean up all resources
     procedure :: free => adjoint_scalars_free
  end type adjoint_scalars_t

contains

  !> Constructor.
  !! @details Initialize the adjoint_scalars container.
  !! @param[inout] this The adjoint_scalars container.
  !! @param[in] n_adjoint_scalars The number of adjoint scalars.
  !! @param[in] n_primal_scalars The number of primal scalars.
  !! @param[in] msh The mesh structure used to define the field topology.
  !! @param coef The SEM coefficients.
  !! @param[inout] gs Gather scatter object.
  !! @param[inout] params_adjoint JSON parameters specific to the
  !!     adjoint scalars.
  !! @param[inout] params_primal JSON parameters specific to the primal scalars.
  !! @param[inout] numerics_params JSON structure containing
  !!     numerics parameters.
  !! @param[in] user User-defined interface.
  !! @param[inout] chkp Checkpointing structure.
  !! @param[in] ulag, vlag, wlag Field history of the primal velocity fields.
  !! @param[in] time_scheme Time scheme controller.
  !! @param[in] rho Density field.
  subroutine adjoint_scalars_init(this, n_adjoint_scalars, n_primal_scalars, &
       msh, coef, gs, params_adjoint, params_primal, &
       numerics_params, user, chkp, ulag, vlag, wlag, time_scheme, rho)
    class(adjoint_scalars_t), intent(inout) :: this
    integer, intent(in) :: n_adjoint_scalars, n_primal_scalars
    type(mesh_t), target, intent(in) :: msh
    type(coef_t), target, intent(in) :: coef
    type(gs_t), target, intent(inout) :: gs
    type(json_file), target, intent(inout) :: params_adjoint, params_primal
    type(json_file), target, intent(inout) :: numerics_params
    type(user_t), target, intent(in) :: user
    type(field_series_t), target, intent(in) :: ulag, vlag, wlag
    type(time_scheme_controller_t), target, intent(in) :: time_scheme
    TYPE(field_t), TARGET, INTENT(IN) :: rho
    type(chkp_t), target, intent(inout) :: chkp
    type(json_file) :: json_subdict_adjoint, json_subdict_primal
    integer :: i, j
    character(len=:), allocatable :: field_name
    character(len=256), allocatable :: field_names(:)
    character(len=:), allocatable :: primal_name, trial_primal_name
    logical :: found_primal

    if (n_adjoint_scalars .gt. n_primal_scalars) then
       call neko_error("More adjoint scalars than forward scalars")
       ! Note. This assumes every adjoint scalar must have a corresponding
       ! primal scalar, but not the other way around.
    end if

    ! Allocate the scalar fields
    ! If there are more adjoint_scalar_scheme_t types, add a factory function
    ! here
    allocate( &
         adjoint_scalar_pnpn_t::this%adjoint_scalar_fields(n_adjoint_scalars))

    allocate(field_names(n_adjoint_scalars))

    ! For multiple adjoint_scalars, collect and validate field names
    if (n_adjoint_scalars > 1) then
       do i = 1, n_adjoint_scalars
          ! Extract element i from the "adjoint_scalars" array
          call json_extract_item(params_adjoint, "", i, json_subdict_adjoint)

          ! Try to get name from JSON, generate one if not found or empty
          if (json_subdict_adjoint%valid_path('name')) then
             call json_get(json_subdict_adjoint, 'name', field_name)
          else
             field_name = ''
          end if

          ! If name is empty or not provided, generate a default one
          if (len_trim(field_name) == 0) then
             write(field_name, '(A,I0,A)') 's_', i, '_adj'
          end if

          field_names(i) = trim(field_name)

          ! If there's a duplicate, append a number until unique
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
       end do
    end if

    do i = 1, n_adjoint_scalars
       call json_extract_item(params_adjoint, "", i, json_subdict_adjoint)

       ! Use the processed field names for multiple adjoint_scalars
       if (n_adjoint_scalars > 1) then
          call json_subdict_adjoint%add('name', trim(field_names(i)))
       else
          call json_subdict_adjoint%add('name', 's_adj')
       end if

       ! Before initializing, find the corresponding primal scalar, and find
       ! its subdict
       ! note: the default primal name is s.
       ! So we search for it if not specified
       call json_get_or_default(params_adjoint, 'primal_name', primal_name, 's')

       found_primal = .false.
       do j = 1, n_primal_scalars
          ! Extract element i from the "adjoint_scalars" array
          call json_extract_item(params_primal, "", j, json_subdict_primal)

          ! Try to get name from JSON
          if (json_subdict_primal%valid_path('name')) then
             call json_get(json_subdict_primal, 'name', trial_primal_name)
             if (trim(trial_primal_name) == trim(primal_name)) then
                found_primal = .true.
                exit
             end if
          end if
       end do

       if (.not. found_primal) then
          call neko_error('Could not find corresponding primal scalar name')
       else
          ! now we modify the json in case a name wasn't specified
          call json_subdict_adjoint%add('primal_name', primal_name)
       end if

       call this%adjoint_scalar_fields(i)%init(msh, coef, gs, &
            json_subdict_adjoint, json_subdict_primal, numerics_params, &
            user, chkp, ulag, vlag, wlag, time_scheme, rho)
    end do
  end subroutine adjoint_scalars_init

  !> Constructor.
  !! @details Initialize a single adjoint_scalar for backwards compatibility.
  !! @param[inout] this The adjoint_scalars container.
  !! @param[in] msh The mesh structure used to define the field topology.
  !! @param coef The SEM coefficients.
  !! @param[inout] gs Gather scatter object.
  !! @param[inout] params_adjoint JSON parameters specific to the adjoint
  !!     scalars.
  !! @param[inout] params_primal JSON parameters specific to the primal scalars.
  !! @param[inout] numerics_params JSON structure containing numerics
  !!     parameters.
  !! @param[in] user User-defined interface.
  !! @param[inout] chkp Checkpointing structure.
  !! @param[in] ulag, vlag, wlag Field history of the primal velocity fields.
  !! @param[in] time_scheme Time scheme controller.
  !! @param[in] rho Density field.
  subroutine adjoint_scalars_init_single(this, msh, coef, gs, params_adjoint, &
       params_primal, numerics_params, user, chkp, ulag, vlag, wlag, &
       time_scheme, rho)
    class(adjoint_scalars_t), intent(inout) :: this
    type(mesh_t), target, intent(in) :: msh
    type(coef_t), target, intent(in) :: coef
    type(gs_t), target, intent(inout) :: gs
    type(json_file), target, intent(inout) :: params_adjoint
    type(json_file), target, intent(inout) :: params_primal
    type(json_file), target, intent(inout) :: numerics_params
    type(user_t), target, intent(in) :: user
    type(chkp_t), target, intent(inout) :: chkp
    type(field_series_t), target, intent(in) :: ulag, vlag, wlag
    type(time_scheme_controller_t), target, intent(in) :: time_scheme
    TYPE(field_t), TARGET, INTENT(IN) :: rho

    ! Allocate a single scalar field
    allocate(adjoint_scalar_pnpn_t::this%adjoint_scalar_fields(1))

    ! TODO, there may be a corner case with multiple primal and a single adjoint
    ! we'll need to catch that before entering this subroutine.

    ! Initialize it directly with the params
    call this%adjoint_scalar_fields(1)%init(msh, coef, gs, params_adjoint, &
         params_primal, numerics_params, user, chkp, ulag, vlag, wlag, &
         time_scheme, rho)
  end subroutine adjoint_scalars_init_single

  !> Perform a time step for all scalar fields
  subroutine adjoint_scalars_step(this, time, ext_bdf, dt_controller)
    class(adjoint_scalars_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    type(time_scheme_controller_t), intent(inout) :: ext_bdf
    type(time_step_controller_t), intent(inout) :: dt_controller
    type(ksp_monitor_t), dimension(size(this%adjoint_scalar_fields)) :: ksp_results
    integer :: i
    logical :: all_frozen

    all_frozen = .true.

    ! Iterate through all scalar fields
    do i = 1, size(this%adjoint_scalar_fields)
       all_frozen = all_frozen .and. this%adjoint_scalar_fields(i)%freeze
       call this%adjoint_scalar_fields(i)%step(time, ext_bdf, dt_controller, &
            ksp_results(i))
    end do

    if (.not. all_frozen) then
       call ksp_results(i)%print_header()
    end if

    do i = 1, size(this%adjoint_scalar_fields)
       if (this%adjoint_scalar_fields(i)%freeze) cycle
       call scalar_step_info(time, ksp_results(i))
    end do
  end subroutine adjoint_scalars_step

  !> Restart from checkpoint data
  subroutine adjoint_scalars_restart(this, chkp)
    class(adjoint_scalars_t), intent(inout) :: this
    type(chkp_t), intent(inout) :: chkp
    integer :: i
    ! Iterate through all scalar fields
    do i = 1, size(this%adjoint_scalar_fields)
       call this%adjoint_scalar_fields(i)%restart(chkp)
    end do
  end subroutine adjoint_scalars_restart

  !> Check if the configuration is valid
  subroutine adjoint_scalars_validate(this)
    class(adjoint_scalars_t), intent(inout) :: this
    integer :: i
    ! Iterate through all scalar fields
    do i = 1, size(this%adjoint_scalar_fields)
       call this%adjoint_scalar_fields(i)%s_adj_lag%set(&
            this%adjoint_scalar_fields(i)%s_adj)
       call this%adjoint_scalar_fields(i)%validate()
    end do
  end subroutine adjoint_scalars_validate

  !> Clean up all resources
  subroutine adjoint_scalars_free(this)
    class(adjoint_scalars_t), intent(inout) :: this
    integer :: i

    ! Iterate through all scalar fields
    if (allocated(this%adjoint_scalar_fields)) then
       do i = 1, size(this%adjoint_scalar_fields)
          call this%adjoint_scalar_fields(i)%free()
       end do
       deallocate(this%adjoint_scalar_fields)
    end if
  end subroutine adjoint_scalars_free

end module adjoint_scalars
