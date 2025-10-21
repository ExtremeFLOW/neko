! Copyright (c) 2020-2025, The Neko Authors
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
!> Defines a simulation case
module case
  use num_types, only : rp, sp, dp
  use fluid_pnpn, only : fluid_pnpn_t
  use fluid_scheme_incompressible, only : fluid_scheme_incompressible_t
  use fluid_scheme_base, only: fluid_scheme_base_t, fluid_scheme_base_factory
  use fluid_output, only : fluid_output_t
  use chkp_output, only : chkp_output_t
  use mesh_field, only : mesh_fld_t
  use parmetis, only : parmetis_partmeshkway
  use redist, only : redist_mesh
  use output_controller, only : output_controller_t
  use flow_ic, only : set_flow_ic
  use scalar_ic, only : set_scalar_ic
  use file, only : file_t
  use utils, only : neko_error
  use mesh, only : mesh_t
  use math, only : NEKO_EPS
  use checkpoint, only: chkp_t
  use time_scheme_controller, only : time_scheme_controller_t
  use logger, only : neko_log, NEKO_LOG_QUIET
  use jobctrl, only : jobctrl_set_time_limit
  use user_intf, only : user_t
  use scalar_pnpn, only : scalar_pnpn_t
  use scalar_scheme, only : scalar_scheme_t
  use time_state, only : time_state_t
  use json_module, only : json_file
  use json_utils, only : json_get, json_get_or_default, json_extract_item, json_no_defaults
  use scratch_registry, only : scratch_registry_t, neko_scratch_registry
  use point_zone_registry, only: neko_point_zone_registry
  use scalars, only : scalars_t
  use comm, only : NEKO_COMM, pe_rank, pe_size
  use mpi_f08, only : MPI_Bcast, MPI_CHARACTER, MPI_INTEGER

  implicit none
  private

  type, public :: case_t
     type(mesh_t) :: msh
     type(json_file) :: params
     character(len=:), allocatable :: output_directory
     type(output_controller_t) :: output_controller
     type(fluid_output_t) :: f_out
     type(time_state_t) :: time
     type(chkp_output_t) :: chkp_out
     type(chkp_t) :: chkp
     type(user_t) :: user
     class(fluid_scheme_base_t), allocatable :: fluid
     type(scalars_t), allocatable :: scalars
  end type case_t

  interface case_init
     module procedure case_init_from_file, case_init_from_json
  end interface case_init

  public :: case_init, case_free

contains

  !> Initialize a case from an input file @a case_file
  subroutine case_init_from_file(this, case_file)
    type(case_t), target, intent(inout) :: this
    character(len=*), intent(in) :: case_file
    integer :: ierr, integer_val
    character(len=:), allocatable :: json_buffer
    logical :: exist

    ! Check if the file exists
    inquire(file = trim(case_file), exist = exist)
    if (.not. exist) then
       call neko_error('The case file '//trim(case_file)//' does not exist.')
    end if

    call neko_log%section('Case')
    call neko_log%message('Reading case file ' // trim(case_file), &
         NEKO_LOG_QUIET)

    if (pe_rank .eq. 0) then
       call this%params%load_file(filename = trim(case_file))
       call this%params%print_to_string(json_buffer)
       integer_val = len(json_buffer)
    end if

    call MPI_Bcast(integer_val, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)
    if (pe_rank .ne. 0) allocate(character(len = integer_val) :: json_buffer)
    call MPI_Bcast(json_buffer, integer_val, MPI_CHARACTER, 0, NEKO_COMM, ierr)
    call this%params%load_from_string(json_buffer)

    deallocate(json_buffer)

    call case_init_common(this)

  end subroutine case_init_from_file

  !> Initialize a case from a JSON object describing a case
  subroutine case_init_from_json(this, case_json)
    type(case_t), target, intent(inout) :: this
    type(json_file), intent(in) :: case_json

    call neko_log%section('Case')
    call neko_log%message('Creating case from JSON object', NEKO_LOG_QUIET)

    this%params = case_json

    call case_init_common(this)

  end subroutine case_init_from_json

  !> Initialize a case from its (loaded) params object
  subroutine case_init_common(this)
    type(case_t), target, intent(inout) :: this
    integer :: lx = 0
    logical :: scalar = .false.
    type(file_t) :: msh_file, bdry_file, part_file
    type(mesh_fld_t) :: msh_part, parts
    logical :: found, logical_val
    logical :: temperature_found = .false.
    integer :: integer_val
    real(kind=rp) :: real_val
    character(len = :), allocatable :: string_val, name, file_format
    integer :: output_dir_len
    integer :: precision, layout
    type(json_file) :: scalar_params, numerics_params
    type(json_file) :: json_subdict
    integer :: n_scalars, i

    !
    ! Setup user defined functions
    !
    call this%user%init()

    ! Run user startup routine
    call this%user%startup(this%params)

    ! Check if default value fill-in is allowed
    if (this%params%valid_path('case.no_defaults')) then
       call json_get(this%params, 'case.no_defaults', json_no_defaults)
    end if

    !
    ! Load mesh
    !
    call json_get_or_default(this%params, 'case.mesh_file', string_val, &
         'no mesh')
    if (trim(string_val) .eq. 'no mesh') then
       call neko_error('The mesh_file keyword could not be found in the .' // &
            'case file. Often caused by incorrectly formatted json.')
    end if
    call msh_file%init(string_val)

    call msh_file%read(this%msh)

    !
    ! Load Balancing
    !
    call json_get_or_default(this%params, 'case.load_balancing', logical_val,&
         .false.)

    if (pe_size .gt. 1 .and. logical_val) then
       call neko_log%section('Load Balancing')
       call parmetis_partmeshkway(this%msh, parts)
       call redist_mesh(this%msh, parts)

       ! store the balanced mesh (for e.g. restarts)
       string_val = trim(string_val(1:scan(trim(string_val), &
            '.', back = .true.) - 1)) // '_lb.nmsh'
       call msh_file%init(string_val)
       call msh_file%write(this%msh)

       call neko_log%end_section()
    end if

    !
    ! Time control
    !
    call json_get(this%params, 'case.time', json_subdict)
    call this%time%init(json_subdict)

    !
    ! Initialize point_zones registry
    !
    call neko_point_zone_registry%init(this%params, this%msh)

    ! Run user mesh motion routine
    call this%user%mesh_setup(this%msh, this%time)

    call json_get(this%params, 'case.numerics', numerics_params)

    !
    ! Setup fluid scheme
    !
    call json_get(this%params, 'case.fluid.scheme', string_val)
    call fluid_scheme_base_factory(this%fluid, trim(string_val))

    call json_get(this%params, 'case.numerics.polynomial_order', lx)
    lx = lx + 1 ! add 1 to get number of gll points
    ! Set time lags in chkp
    this%chkp%tlag => this%time%tlag
    this%chkp%dtlag => this%time%dtlag
    call this%fluid%init(this%msh, lx, this%params, this%user, this%chkp)


    !
    ! Setup scratch registry
    !
    call neko_scratch_registry%init(this%fluid%dm_Xh, 10, 10)

    !
    ! Setup scalar scheme
    !
    ! @todo no scalar factory for now, probably not needed
    scalar = .false.
    n_scalars = 0
    if (this%params%valid_path('case.scalar')) then
       call json_get_or_default(this%params, 'case.scalar.enabled', scalar, &
            .true.)
       n_scalars = 1
    else if (this%params%valid_path('case.scalars')) then
       call this%params%info('case.scalars', n_children = n_scalars)
       if (n_scalars > 0) then
          scalar = .true.
       end if
    end if

    if (scalar) then
       allocate(this%scalars)
       if (this%params%valid_path('case.scalar')) then
          ! For backward compatibility
          call json_get(this%params, 'case.scalar', scalar_params)
          call this%scalars%init(this%msh, this%fluid%c_Xh, this%fluid%gs_Xh, &
               scalar_params, numerics_params, this%user, this%chkp, this%fluid%ulag, &
               this%fluid%vlag, this%fluid%wlag, this%fluid%ext_bdf, &
               this%fluid%rho)
       else
          ! Multiple scalars
          call json_get(this%params, 'case.scalars', json_subdict)
          call this%scalars%init(n_scalars, this%msh, this%fluid%c_Xh, this%fluid%gs_Xh, &
               json_subdict, numerics_params, this%user, this%chkp, this%fluid%ulag, &
               this%fluid%vlag, this%fluid%wlag, this%fluid%ext_bdf, &
               this%fluid%rho)
       end if
    end if

    !
    ! Setup initial conditions
    !
    call json_get(this%params, 'case.fluid.initial_condition.type', &
         string_val)
    call json_get(this%params, 'case.fluid.initial_condition', &
         json_subdict)

    call neko_log%section("Fluid initial condition ")

    if (this%params%valid_path('case.restart_file')) then
       call neko_log%message("Restart file specified, " // &
            "initial conditions ignored")
    else if (trim(string_val) .ne. 'user') then
       call set_flow_ic(this%fluid%u, this%fluid%v, this%fluid%w, &
            this%fluid%p, this%fluid%c_Xh, this%fluid%gs_Xh, string_val, &
            json_subdict)
    else
       call json_get(this%params, 'case.fluid.scheme', string_val)
       if (trim(string_val) .eq. 'compressible') then
          call set_flow_ic(this%fluid%rho, &
               this%fluid%u, this%fluid%v, this%fluid%w, this%fluid%p, &
               this%fluid%c_Xh, this%fluid%gs_Xh, &
               this%user%initial_conditions, this%fluid%name)
       else
          call set_flow_ic(this%fluid%u, this%fluid%v, this%fluid%w, &
               this%fluid%p, this%fluid%c_Xh, this%fluid%gs_Xh, &
               this%user%initial_conditions, this%fluid%name)
       end if
    end if

    call neko_log%end_section()

    if (scalar) then
       call neko_log%section("Scalar initial condition ")

       if (this%params%valid_path('case.restart_file')) then
          call neko_log%message("Restart file specified, " // &
               "initial conditions ignored")
       else if (this%params%valid_path('case.scalar')) then
          ! For backward compatibility with single scalar
          call json_get(this%params, 'case.scalar.initial_condition.type', &
               string_val)
          call json_get(this%params, &
               'case.scalar.initial_condition', json_subdict)

          if (trim(string_val) .ne. 'user') then
             if (trim(this%scalars%scalar_fields(1)%name) .eq. 'temperature') then
                call set_scalar_ic(this%scalars%scalar_fields(1)%s, &
                     this%scalars%scalar_fields(1)%c_Xh, &
                     this%scalars%scalar_fields(1)%gs_Xh, &
                     string_val, json_subdict, 0)
             else
                call set_scalar_ic(this%scalars%scalar_fields(1)%s, &
                     this%scalars%scalar_fields(1)%c_Xh, &
                     this%scalars%scalar_fields(1)%gs_Xh, &
                     string_val, json_subdict, 1)
             end if
          else
             call set_scalar_ic(this%scalars%scalar_fields(1)%name, &
                  this%scalars%scalar_fields(1)%s, &
                  this%scalars%scalar_fields(1)%c_Xh, &
                  this%scalars%scalar_fields(1)%gs_Xh, &
                  this%user%initial_conditions)
          end if

       else
          ! Handle multiple scalars
          do i = 1, n_scalars
             call json_extract_item(this%params, 'case.scalars', i, &
                  scalar_params)
             call json_get(scalar_params, 'initial_condition.type', string_val)
             call json_get(scalar_params, 'initial_condition', &
                  json_subdict)

             if (trim(string_val) .ne. 'user') then
                if (trim(this%scalars%scalar_fields(i)%name) .eq. 'temperature') then
                   call set_scalar_ic(this%scalars%scalar_fields(i)%s, &
                        this%scalars%scalar_fields(i)%c_Xh, &
                        this%scalars%scalar_fields(i)%gs_Xh, &
                        string_val, json_subdict, 0)
                   temperature_found = .true.
                else
                   if (temperature_found) then
                      ! if temperature is found, other scalars start from index 1
                      call set_scalar_ic(this%scalars%scalar_fields(i)%s, &
                           this%scalars%scalar_fields(i)%c_Xh, &
                           this%scalars%scalar_fields(i)%gs_Xh, &
                           string_val, json_subdict, i - 1)
                   else
                      ! if temperature is not found, other scalars start from index 0
                      call set_scalar_ic(this%scalars%scalar_fields(i)%s, &
                           this%scalars%scalar_fields(i)%c_Xh, &
                           this%scalars%scalar_fields(i)%gs_Xh, &
                           string_val, json_subdict, i)
                   end if
                end if
             else
                call set_scalar_ic(this%scalars%scalar_fields(i)%name,&
                     this%scalars%scalar_fields(i)%s, &
                     this%scalars%scalar_fields(i)%c_Xh, &
                     this%scalars%scalar_fields(i)%gs_Xh, &
                     this%user%initial_conditions)
             end if
          end do
       end if

       call neko_log%end_section()
    end if

    ! Add initial conditions to BDF scheme (if present)
    select type (f => this%fluid)
    type is (fluid_pnpn_t)
       call f%ulag%set(f%u)
       call f%vlag%set(f%v)
       call f%wlag%set(f%w)
    end select

    !
    ! Validate that the case is properly setup for time-stepping
    !
    call this%fluid%validate

    if (scalar) then
       call this%scalars%validate()
    end if

    !
    ! Get and process output directory
    !
    call json_get_or_default(this%params, 'case.output_directory',&
         this%output_directory, '')

    output_dir_len = len(trim(this%output_directory))
    if (output_dir_len .gt. 0) then
       if (this%output_directory(output_dir_len:output_dir_len) .ne. "/") then
          this%output_directory = trim(this%output_directory)//"/"
          if (pe_rank .eq. 0) then
             call execute_command_line('mkdir -p '//this%output_directory)
          end if
       end if
    end if

    !
    ! Save mesh partitions (if requested)
    !
    call json_get_or_default(this%params, 'case.output_partitions',&
         logical_val, .false.)
    if (logical_val) then
       call msh_part%init(this%msh, 'MPI_Rank')
       msh_part%data = pe_rank
       call part_file%init(trim(this%output_directory)//'partitions.vtk')
       call part_file%write(msh_part)
       call msh_part%free()
    end if

    !
    ! Setup output precision of the field files
    !
    call json_get_or_default(this%params, 'case.output_precision', string_val,&
         'single')

    if (trim(string_val) .eq. 'double') then
       precision = dp
    else
       precision = sp
    end if

    !
    ! Setup output layout of the field bp file
    !
    call json_get_or_default(this%params, 'case.output_layout', layout, 1)

    !
    ! Setup output_controller
    !
    call json_get_or_default(this%params, 'case.fluid.output_filename', &
         name, "field")
    call json_get_or_default(this%params, 'case.fluid.output_format', &
         file_format, 'fld')
    call this%output_controller%init(this%time%end_time)
    if (scalar) then
       call this%f_out%init(precision, this%fluid, this%scalars, name = name, &
            path = trim(this%output_directory), &
            fmt = trim(file_format), layout = layout)
    else
       call this%f_out%init(precision, this%fluid, name = name, &
            path = trim(this%output_directory), &
            fmt = trim(file_format), layout = layout)
    end if

    call json_get_or_default(this%params, 'case.fluid.output_control',&
         string_val, 'org')

    if (trim(string_val) .eq. 'org') then
       ! yes, it should be real_val below for type compatibility
       call json_get(this%params, 'case.nsamples', real_val)
       call this%output_controller%add(this%f_out, real_val, 'nsamples')
    else if (trim(string_val) .eq. 'never') then
       ! Fix a dummy 0.0 output_value
       call json_get_or_default(this%params, 'case.fluid.output_value', &
            real_val, 0.0_rp)
       call this%output_controller%add(this%f_out, 0.0_rp, string_val)
    else
       call json_get(this%params, 'case.fluid.output_value', real_val)
       call this%output_controller%add(this%f_out, real_val, string_val)
    end if

    !
    ! Save checkpoints (if nothing specified, default to saving at end of sim)
    !
    call json_get_or_default(this%params, 'case.output_checkpoints',&
         logical_val, .true.)
    if (logical_val) then
       call json_get_or_default(this%params, 'case.checkpoint_filename', &
            name, "fluid")
       call json_get_or_default(this%params, 'case.checkpoint_format', &
            string_val, "chkp")
       call this%chkp_out%init(this%chkp, name = name,&
            path = this%output_directory, fmt = trim(string_val))
       call json_get_or_default(this%params, 'case.checkpoint_control', &
            string_val, "simulationtime")
       call json_get_or_default(this%params, 'case.checkpoint_value', &
            real_val, 1e10_rp)
       call this%output_controller%add(this%chkp_out, real_val, string_val, &
            NEKO_EPS)
    end if

    !
    ! Setup joblimit
    !
    if (this%params%valid_path('case.job_timelimit')) then
       call json_get(this%params, 'case.job_timelimit', string_val)
       call jobctrl_set_time_limit(string_val)
    end if

    call neko_log%end_section()

  end subroutine case_init_common

  !> Deallocate a case
  subroutine case_free(this)
    type(case_t), intent(inout) :: this

    if (allocated(this%fluid)) then
       call this%fluid%free()
       deallocate(this%fluid)
    end if

    if (allocated(this%scalars)) then
       call this%scalars%free()
       deallocate(this%scalars)
    end if

    call this%msh%free()

    call this%f_out%free()

    call this%output_controller%free()

  end subroutine case_free

end module case
