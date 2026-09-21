!> @file adjoint_case.f90
!! @copyright
!! Copyright (c) 2024-2025, The Neko-TOP Authors
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

! Implements the `adjoint_case_t` type.
module adjoint_case
  use num_types, only: rp, dp, sp
  use case, only: case_t
  use adjoint_fluid_scheme, only: adjoint_fluid_scheme_t
  use adjoint_fluid_fctry, only: adjoint_fluid_scheme_factory
  use adjoint_fluid_pnpn, only: adjoint_fluid_pnpn_t
  use adjoint_output, only: adjoint_output_t
  use scalar_ic, only: set_scalar_ic
  use checkpoint, only: chkp_t
  use chkp_output, only: chkp_output_t
  use flow_ic, only: set_flow_ic
  use output_controller, only: output_controller_t
  use time_based_controller, only: time_based_controller_t
  use file, only: file_t
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default, json_extract_item, &
       json_get_or_lookup
  use adjoint_scalar_scheme, only: adjoint_scalar_scheme_t
  use adjoint_scalar_pnpn, only : adjoint_scalar_pnpn_t
  use logger, only : neko_log
  use time_state, only : time_state_t
  use utils, only: neko_error
  use adjoint_scalar_convection_source_term, only: &
       adjoint_scalar_convection_source_term_t
  use json_utils_ext, only: json_key_fallback
  use adjoint_scalars, only: adjoint_scalars_t
  implicit none
  private
  public :: adjoint_case_t

  !> Adjoint case type.
  !! Todo: This should Ideally be a subclass of case_t, however, that is not yet
  !! supported by Neko.
  type :: adjoint_case_t
     !> Adjoint fluid
     class(adjoint_fluid_scheme_t), allocatable :: fluid_adj
     !> Adjoint scalar
     type(adjoint_scalars_t), allocatable :: adjoint_scalars
     !> Source term coupling the adjoint scalar to the adjoint fluid
     type(adjoint_scalar_convection_source_term_t), allocatable :: &
          adjoint_convection_term
     type(case_t), pointer :: case
     type(time_state_t) :: time
     type(chkp_t) :: chkp
     type(chkp_output_t) :: chkp_out

     ! Fields
     type(adjoint_output_t) :: f_out
     type(output_controller_t) :: output_controller
     type(time_based_controller_t) :: norm_output_ctrl
     type(file_t) :: norm_output_file
     logical :: norm_output_enabled = .false.

     logical :: have_scalar = .false.

   contains
     procedure, pass(this) :: init => adjoint_init_from_json
     procedure, pass(this) :: free => adjoint_free
  end type adjoint_case_t

contains

  ! Constructor from json.
  subroutine adjoint_init_from_json(this, neko_case)
    class(adjoint_case_t), intent(inout) :: this
    type(case_t), target, intent(inout) :: neko_case

    this%case => neko_case
    call adjoint_case_init_common(this, neko_case)

  end subroutine adjoint_init_from_json

  !> Initialize a neko_case from its (loaded) params object
  subroutine adjoint_case_init_common(this, neko_case)
    class(adjoint_case_t), target, intent(inout) :: this
    type(case_t), intent(inout) :: neko_case
    integer :: lx = 0
    real(kind=rp) :: real_val = 0.0_rp
    character(len=:), allocatable :: string_val, file_format, name
    character(len=:), allocatable :: norm_control, norm_file
    integer :: precision, integer_val, layout
    integer :: n_scalars_primal, n_scalars_adjoint, i
    logical :: scalar = .false.
    logical :: temperature_found = .false.
    logical :: logical_val

    ! extra things for json
    type(json_file) :: ic_json, numerics_params
    type(json_file) :: scalar_params_primal, scalar_params_adjoint, json_subdict
    character(len=:), allocatable :: json_key
    logical :: dealias_adjoint_scalar_convection

    !
    ! Setup adjoint fluid
    !
    call json_get(neko_case%params, 'case.fluid.scheme', string_val)
    call adjoint_fluid_scheme_factory(this%fluid_adj, trim(string_val))

    call json_get(neko_case%params, 'case.numerics.polynomial_order', lx)
    lx = lx + 1 ! add 1 to get number of gll points

    call this%chkp%init()
    this%chkp%tlag => this%time%tlag
    this%chkp%dtlag => this%time%dtlag

    select type (f => this%fluid_adj)
    type is (adjoint_fluid_pnpn_t)
       call f%init(neko_case%msh, lx, neko_case%params, &
            neko_case%user, this%chkp)
    end select
    !
    ! Setup adjoint scalar
    !
    ! @todo no adjoint_scalars factory for now, probably not needed

    ! check how many adjoint scalars
    scalar = .false.
    n_scalars_adjoint = 0
    if (neko_case%params%valid_path('case.adjoint_scalar')) then
       call json_get_or_default(neko_case%params, &
            'case.adjoint_scalar.enabled', scalar, .true.)
       n_scalars_adjoint = 1
       n_scalars_primal = 1
    else if (neko_case%params%valid_path('case.adjoint_scalars')) then
       call neko_case%params%info('case.adjoint_scalars', &
            n_children = n_scalars_adjoint)
       call neko_case%params%info('case.scalars', n_children = n_scalars_primal)
       if (n_scalars_adjoint > 0) then
          scalar = .true.
       end if
    end if

    this%have_scalar = scalar






    if (this%have_scalar) then
       allocate(this%adjoint_scalars)
       call json_get(neko_case%params, 'case.numerics', &
            numerics_params)
       if (neko_case%params%valid_path('case.adjoint_scalar')) then
          ! For backward compatibility
          call json_get(neko_case%params, 'case.adjoint_scalar', &
               scalar_params_adjoint)
          call json_get(neko_case%params, 'case.scalar', &
               scalar_params_primal)
          call this%adjoint_scalars%init(neko_case%msh, neko_case%fluid%c_Xh, &
               neko_case%fluid%gs_Xh, scalar_params_adjoint, &
               scalar_params_primal, numerics_params, neko_case%user, &
               neko_case%chkp, neko_case%fluid%ulag, neko_case%fluid%vlag, &
               neko_case%fluid%wlag, neko_case%fluid%ext_bdf, &
               neko_case%fluid%rho)
          ! allocate the coupling term
          allocate(this%adjoint_convection_term)
          ! initialize the coupling term
          call json_get_or_default(neko_case%params, &
               'case.adjoint_scalar.dealias_coupling_term', &
               dealias_adjoint_scalar_convection, .true.)
          call this%adjoint_convection_term%init_from_components( &
               this%fluid_adj%f_adj_x, this%fluid_adj%f_adj_y, &
               this%fluid_adj%f_adj_z, this%case%scalars%scalar_fields(1)%scalar%s, &
               this%adjoint_scalars%adjoint_scalar_fields(1)%s_adj, &
               this%fluid_adj%c_Xh, this%fluid_adj%c_Xh_GL, &
               this%fluid_adj%GLL_to_GL, &
               dealias_adjoint_scalar_convection, this%fluid_adj%scratch_GL)

          select type (f => this%fluid_adj)
          type is (adjoint_fluid_pnpn_t)
             ! append the coupling term to the adjoint velocity equation
             call f%source_term%add(this%adjoint_convection_term)
          end select
       else
          ! Multiple scalars

          call json_get(this%case%params, &
               'case.adjoint_scalars', scalar_params_adjoint)
          call json_get(this%case%params, &
               'case.scalars', scalar_params_primal)
          call this%adjoint_scalars%init(n_scalars_adjoint, n_scalars_primal, &
               neko_case%msh, neko_case%fluid%c_Xh, neko_case%fluid%gs_Xh, &
               scalar_params_adjoint, scalar_params_primal, numerics_params, &
               neko_case%user, neko_case%chkp, neko_case%fluid%ulag, &
               neko_case%fluid%vlag, neko_case%fluid%wlag, &
               neko_case%fluid%ext_bdf, neko_case%fluid%rho)
          call neko_error('The adjoint scaling coupling term have not yet' // &
               'been implemented for multiple scalars')
       end if
    end if

    !
    ! Time control
    !
    call json_get(this%case%params, 'case.time', json_subdict)
    call this%time%init(json_subdict)

    !
    ! Setup user defined conditions
    !
    ! if (neko_case%params%valid_path('case.fluid.inflow_condition')) then
    !    call json_get(neko_case%params, 'case.fluid.inflow_condition.type', &
    !         string_val)
    !    if (trim(string_val) .eq. 'user') then
    !       call neko_case%fluid%set_usr_inflow(neko_case%user%fluid_user_if)
    !    end if
    ! end if

    ! Setup user boundary conditions for the scalar.
    ! if (adjoint_scalars) then
    !    call neko_case%adjoint_scalars%set_user_bc(&
    !         neko_case%user%scalar_user_bc)
    ! end if

    !
    ! Setup initial conditions
    !

    call neko_log%section("Adjoint initial condition")
    json_key = json_key_fallback(neko_case%params, &
         'case.adjoint_fluid.initial_condition', 'case.fluid.initial_condition')

    call json_get(neko_case%params, json_key, ic_json)
    call json_get(ic_json, 'type', string_val)

    if (trim(string_val) .ne. 'user') then
       call set_flow_ic( &
            this%fluid_adj%u_adj, this%fluid_adj%v_adj, this%fluid_adj%w_adj, &
            this%fluid_adj%p_adj, this%fluid_adj%c_Xh, this%fluid_adj%gs_Xh, &
            string_val, ic_json)
    else
       call set_flow_ic( &
            this%fluid_adj%u_adj, this%fluid_adj%v_adj, this%fluid_adj%w_adj, &
            this%fluid_adj%p_adj, this%fluid_adj%c_Xh, this%fluid_adj%gs_Xh, &
            neko_case%user%initial_conditions, neko_case%fluid%name)
    end if

    call neko_log%end_section()

    if (this%have_scalar) then

       if (neko_case%params%valid_path('case.adjoint_scalar')) then
          ! we shouldn't fallback to the primal here.
          call json_get(neko_case%params, &
               'case.adjoint_scalar.initial_condition.type', string_val)
          call json_get(neko_case%params, &
               'case.adjoint_scalar.initial_condition', ic_json)

          !call neko_log%section("Adjoint scalar initial condition ")

          if (trim(string_val) .ne. 'user') then
             if (trim(neko_case%scalars%scalar_fields(1)%scalar%name) .eq. &
                  'temperature') then
                call set_scalar_ic(&
                     this%adjoint_scalars%adjoint_scalar_fields(1)%s_adj, &
                     this%adjoint_scalars%adjoint_scalar_fields(1)%c_Xh, &
                     this%adjoint_scalars%adjoint_scalar_fields(1)%gs_Xh, &
                     string_val, ic_json, 0)
             else
                call set_scalar_ic(&
                     this%adjoint_scalars%adjoint_scalar_fields(1)%s_adj, &
                     this%adjoint_scalars%adjoint_scalar_fields(1)%c_Xh, &
                     this%adjoint_scalars%adjoint_scalar_fields(1)%gs_Xh, &
                     string_val, ic_json, 1)
             end if
          else
             call neko_error("user ICs not implemented for adjoint scalar")
             ! call set_scalar_ic(this%adjoint_scalars%s_adj, &
             !      this%adjoint_scalars%c_Xh, this%adjoint_scalars%gs_Xh, &
             !      this%usr%scalar_user_ic, neko_case%params)
          end if

          ! call neko_log%end_section()
       else

          ! Handle multiple scalars
          do i = 1, n_scalars_adjoint
             call json_extract_item(neko_case%params, 'case.adjoint_scalars', &
                  i, scalar_params_adjoint)
             call json_get(scalar_params_adjoint, &
                  'initial_condition.type', string_val)
             call json_get(scalar_params_adjoint, &
                  'initial_condition', json_subdict)

             if (trim(string_val) .ne. 'user') then
                if (trim(neko_case%scalars%scalar_fields(i)%scalar%name) .eq. &
                     'temperature') then
                   call set_scalar_ic(&
                        this%adjoint_scalars%adjoint_scalar_fields(i)%s_adj, &
                        this%adjoint_scalars%adjoint_scalar_fields(i)%c_Xh, &
                        this%adjoint_scalars%adjoint_scalar_fields(i)%gs_Xh, &
                        string_val, json_subdict, 0)
                   temperature_found = .true.
                else
                   if (temperature_found) then
                      ! if temperature is found, scalars start from index 1
                      call set_scalar_ic(&
                           this%adjoint_scalars%adjoint_scalar_fields(i)%s_adj, &
                           this%adjoint_scalars%adjoint_scalar_fields(i)%c_Xh, &
                           this%adjoint_scalars%adjoint_scalar_fields(i)%gs_Xh, &
                           string_val, json_subdict, i - 1)
                   else
                      ! if temperature is not found, scalars start from index 0
                      call set_scalar_ic(&
                           this%adjoint_scalars%adjoint_scalar_fields(i)%s_adj, &
                           this%adjoint_scalars%adjoint_scalar_fields(i)%c_Xh, &
                           this%adjoint_scalars%adjoint_scalar_fields(i)%gs_Xh, &
                           string_val, json_subdict, i)
                   end if
                end if
             else
                call neko_error("user ICs not implemented for adjoint scalar")
             end if
          end do
       end if
    end if

    ! Add initial conditions to BDF fluid_adj (if present)
    select type (f => this%fluid_adj)
    type is (adjoint_fluid_pnpn_t)
       call f%ulag%set(f%u_adj)
       call f%vlag%set(f%v_adj)
       call f%wlag%set(f%w_adj)
    end select

    !
    ! Validate that the neko_case is properly setup for time-stepping
    !
    call this%fluid_adj%validate

    if (this%have_scalar) then
       call this%adjoint_scalars%validate()
    end if

    !
    ! Setup output precision of the field files
    !
    call json_get_or_default(neko_case%params, 'case.output_precision', &
         string_val, 'single')

    if (trim(string_val) .eq. 'double') then
       precision = dp
    else
       precision = sp
    end if

    !
    ! Setup output_controller
    !
    call json_get_or_default(neko_case%params, &
         'case.adjoint_fluid.output_filename', name, "adjoint")
    call json_get_or_default(neko_case%params, &
         'case.adjoint_fluid.output_format', file_format, 'fld')
    call json_get_or_default(neko_case%params, &
         'case.adjoint_fluid.output_mesh_in_all_files', &
         logical_val, .false.)
    call this%output_controller%init(this%time%end_time)
    if (scalar) then
       call this%f_out%init(precision, this%fluid_adj, &
            this%adjoint_scalars, name = name, &
            path = trim(neko_case%output_directory), &
            fmt = trim(file_format), layout = layout, &
            always_write_mesh = logical_val)
    else
       call this%f_out%init(precision, this%fluid_adj, name = name, &
            path = trim(neko_case%output_directory), &
            fmt = trim(file_format), layout = layout, &
            always_write_mesh = logical_val)
    end if

    call json_get_or_default(neko_case%params, &
         'case.adjoint_fluid.output_subdivide', logical_val, .false.)
    call this%f_out%file_%set_subdivide(logical_val)

    call json_get_or_default(neko_case%params, &
         'case.adjoint_fluid.output_control', string_val, 'never')

    if (trim(string_val) .eq. 'org') then
       ! yes, it should be real_val below for type compatibility
       call json_get_or_lookup(neko_case%params, 'case.nsamples', integer_val)
       real_val = real(integer_val, kind=rp)
       call this%output_controller%add(this%f_out, real_val, 'nsamples')
    else if (trim(string_val) .eq. 'never') then
       call this%output_controller%add(this%f_out, 0.0_rp, 'never')
    else if (trim(string_val) .eq. 'tsteps' .or. &
         trim(string_val) .eq. 'nsamples') then
       call json_get_or_lookup(neko_case%params, &
            'case.adjoint_fluid.output_value', integer_val)
       real_val = real(integer_val, kind=rp)
       call this%output_controller%add(this%f_out, real_val, string_val)
    else if (trim(string_val) .eq. 'simulationtime') then
       call json_get_or_lookup(neko_case%params, &
            'case.adjoint_fluid.output_value', real_val)
       call this%output_controller%add(this%f_out, real_val, string_val)
    else
       call neko_log%error('Unknown output control type for the fluid: ' // &
            trim(string_val))
    end if

    !
    ! Setup adjoint norm output
    !
    call json_get_or_default(neko_case%params, &
         'case.adjoint_fluid.norm_output_control', norm_control, 'never')
    if (trim(norm_control) .ne. 'never') then
       call json_get_or_default(neko_case%params, &
            'case.adjoint_fluid.norm_output_value', real_val, 1.0_rp)
       call json_get_or_default(neko_case%params, &
            'case.adjoint_fluid.norm_output_file', norm_file, &
            'adjoint_norm.csv')
       call this%norm_output_file%init(trim(neko_case%output_directory) // &
            trim(norm_file))
       call this%norm_output_file%set_header('Time, Norm')
       call this%norm_output_file%set_overwrite(.true.)
       call this%norm_output_ctrl%init(this%time%start_time, &
            this%time%end_time, trim(norm_control), real_val)
       this%norm_output_enabled = .true.
    end if

    ! !
    ! ! Save checkpoints (if nothing specified, default to saving at end of sim)
    ! !
    ! call json_get_or_default(neko_case%params, 'case.output_checkpoints',&
    !      logical_val, .true.)
    ! if (logical_val) then
    !    call json_get_or_default(neko_case%params, 'case.checkpoint_format', &
    !         string_val, "chkp")
    !   neko_case%f_chkp = chkp_output_t(this%fluid_adj%chkp, &
    ! path = output_directory, &
    !         ! fmt = trim(string_val))
    !    call json_get_or_default(neko_case%params, 'case.checkpoint_control', &
    !         string_val, "simulationtime")
    !    call json_get_or_default(neko_case%params, 'case.checkpoint_value', &
    ! real_val,&
    !         1e10_rp)
    !   call this%output_controller%add(this%f_chkp, real_val, string_val)
    ! end if

    !
    ! Initialize time and step
    !
    this%time%t = 0d0
    this%time%tstep = 0

  end subroutine adjoint_case_init_common

  ! Destructor.
  subroutine adjoint_free(this)
    class(adjoint_case_t), intent(inout) :: this

    if (allocated(this%fluid_adj)) then
       call this%fluid_adj%free()
       deallocate(this%fluid_adj)
    end if

    if (allocated(this%adjoint_scalars)) then
       call this%adjoint_scalars%free()
       deallocate(this%adjoint_scalars)
    end if

    if (allocated(this%adjoint_convection_term)) then
       call this%adjoint_convection_term%free()
       deallocate(this%adjoint_convection_term)
    end if

    nullify(this%case)

    ! call this%time%free()
    call this%chkp%free()
    call this%chkp_out%free()

    ! Fields
    call this%f_out%free()
    call this%output_controller%free()
    call this%norm_output_ctrl%free()
    call this%norm_output_file%free()

    this%have_scalar = .false.

  end subroutine adjoint_free

end module adjoint_case
