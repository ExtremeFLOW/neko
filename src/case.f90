! Copyright (c) 2020-2023, The Neko Authors
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
  use fluid_scheme, only : fluid_scheme_t, fluid_scheme_factory
  use fluid_output, only : fluid_output_t
  use chkp_output, only : chkp_output_t
  use mesh_field, only : mesh_fld_t, mesh_field_init, mesh_field_free
  use parmetis, only : parmetis_partmeshkway
  use redist, only : redist_mesh
  use sampler, only : sampler_t
  use flow_ic, only : set_flow_ic
  use scalar_ic, only : set_scalar_ic
  use stats, only : stats_t
  use file, only : file_t
  use utils, only : neko_error
  use mesh, only : mesh_t
  use comm
  use time_scheme_controller, only : time_scheme_controller_t
  use logger, only : neko_log, NEKO_LOG_QUIET, LOG_SIZE
  use jobctrl, only : jobctrl_set_time_limit
  use user_intf, only : user_t
  use scalar_pnpn, only : scalar_pnpn_t
  use json_module, only : json_file, json_core, json_value
  use json_utils, only : json_get, json_get_or_default
  use scratch_registry, only : scratch_registry_t, neko_scratch_registry
  use point_zone_registry, only: neko_point_zone_registry
  implicit none
  private

  type, public :: case_t
     type(mesh_t) :: msh
     type(json_file) :: params
     type(time_scheme_controller_t) :: ext_bdf
     real(kind=rp), dimension(10) :: tlag
     real(kind=rp), dimension(10) :: dtlag
     real(kind=rp) :: dt
     real(kind=rp) :: end_time
     character(len=:), allocatable :: output_directory
     type(sampler_t) :: s
     type(fluid_output_t) :: f_out
     type(chkp_output_t) :: f_chkp
     type(user_t) :: usr
     class(fluid_scheme_t), allocatable :: fluid
     type(scalar_pnpn_t), allocatable :: scalar
  end type case_t

  interface case_init
     module procedure case_init_from_file, case_init_from_json
  end interface case_init

  public :: case_init, case_free

contains

  !> Initialize a case from an input file @a case_file
  subroutine case_init_from_file(C, case_file)
    type(case_t), target, intent(inout) :: C
    character(len=*), intent(in) :: case_file
    integer :: ierr, integer_val
    character(len=:), allocatable :: json_buffer

    call neko_log%section('Case')
    call neko_log%message('Reading case file ' // trim(case_file), &
                          NEKO_LOG_QUIET)

    if (pe_rank .eq. 0) then
       call C%params%load_file(filename = trim(case_file))
       call C%params%print_to_string(json_buffer)
       integer_val = len(json_buffer)
    end if

    call MPI_Bcast(integer_val, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)
    if (pe_rank .ne. 0) allocate(character(len = integer_val) :: json_buffer)
    call MPI_Bcast(json_buffer, integer_val, MPI_CHARACTER, 0, NEKO_COMM, ierr)
    call C%params%load_from_string(json_buffer)

    deallocate(json_buffer)

    call case_init_common(C)

  end subroutine case_init_from_file

  !> Initialize a case from a JSON object describing a case
  subroutine case_init_from_json(C, case_json)
    type(case_t), target, intent(inout) :: C
    type(json_file), intent(in) :: case_json

    call neko_log%section('Case')
    call neko_log%message('Creating case from JSON object', NEKO_LOG_QUIET)

    C%params = case_json

    call case_init_common(C)

  end subroutine case_init_from_json

  !> Initialize a case from its (loaded) params object
  subroutine case_init_common(C)
    type(case_t), target, intent(inout) :: C
    integer :: lx = 0
    logical :: scalar = .false.
    type(file_t) :: msh_file, bdry_file, part_file
    type(mesh_fld_t) :: msh_part, parts
    logical :: found, logical_val
    integer :: integer_val
    real(kind=rp) :: real_val
    character(len = :), allocatable :: string_val
    real(kind=rp) :: stats_start_time, stats_output_val
    integer :: stats_sampling_interval
    integer :: output_dir_len
    integer :: precision

    !
    ! Load mesh
    !
    call json_get_or_default(C%params, 'case.mesh_file', string_val, 'no mesh')
    if (trim(string_val) .eq. 'no mesh') then
       call neko_error('The mesh_file keyword could not be found in the .' // &
                       'case file. Often caused by incorrectly formatted json.')
    end if
    msh_file = file_t(string_val)

    call msh_file%read(C%msh)

    !
    ! Load Balancing
    !
    call json_get_or_default(C%params, 'case.load_balancing', logical_val,&
                             .false.)

    if (pe_size .gt. 1 .and. logical_val) then
       call neko_log%section('Load Balancing')
       call parmetis_partmeshkway(C%msh, parts)
       call redist_mesh(C%msh, parts)
       call neko_log%end_section()
    end if

    !
    ! Time step
    !
    call C%params%get('case.variable_timestep', logical_val, found)
    if (.not. logical_val) then
       call json_get(C%params, 'case.timestep', C%dt)
    else
       ! randomly set an initial dt to get cfl when dt is variable
       C%dt = 1.0_rp
    end if

    !
    ! End time
    !
    call json_get(C%params, 'case.end_time', C%end_time)

    !
    ! Initialize point_zones registry
    !
    call neko_point_zone_registry%init(C%params, C%msh)

    !
    ! Setup user defined functions
    !
    call C%usr%init()
    call C%usr%user_mesh_setup(C%msh)

    !
    ! Set order of timestepper
    !
    call json_get(C%params, 'case.numerics.time_order', integer_val)
    call C%ext_bdf%init(integer_val)

    !
    ! Setup fluid scheme
    !
    call json_get(C%params, 'case.fluid.scheme', string_val)
    call fluid_scheme_factory(C%fluid, trim(string_val))

    call json_get(C%params, 'case.numerics.polynomial_order', lx)
    lx = lx + 1 ! add 1 to get number of gll points
    C%fluid%chkp%tlag => C%tlag
    C%fluid%chkp%dtlag => C%dtlag
    call C%fluid%init(C%msh, lx, C%params, C%usr, C%ext_bdf)
    select type (f => C%fluid)
    type is (fluid_pnpn_t)
       f%chkp%abx1 => f%abx1
       f%chkp%abx2 => f%abx2
       f%chkp%aby1 => f%aby1
       f%chkp%aby2 => f%aby2
       f%chkp%abz1 => f%abz1
       f%chkp%abz2 => f%abz2
    end select


    !
    ! Setup scratch registry
    !
    neko_scratch_registry = scratch_registry_t(C%fluid%dm_Xh, 10, 10)

    !
    ! Setup scalar scheme
    !
    ! @todo no scalar factory for now, probably not needed
    if (C%params%valid_path('case.scalar')) then
       call json_get_or_default(C%params, 'case.scalar.enabled', scalar,&
                                .true.)
    end if

    if (scalar) then
       allocate(C%scalar)
       C%scalar%chkp%tlag => C%tlag
       C%scalar%chkp%dtlag => C%dtlag
       call C%scalar%init(C%msh, C%fluid%c_Xh, C%fluid%gs_Xh, C%params, C%usr,&
            C%fluid%ulag, C%fluid%vlag, C%fluid%wlag, C%ext_bdf, C%fluid%rho)
       call C%fluid%chkp%add_scalar(C%scalar%s)
       C%fluid%chkp%abs1 => C%scalar%abx1
       C%fluid%chkp%abs2 => C%scalar%abx2
       C%fluid%chkp%slag => C%scalar%slag
    end if

    !
    ! Setup user defined conditions
    !
    if (C%params%valid_path('case.fluid.inflow_condition')) then
       call json_get(C%params, 'case.fluid.inflow_condition.type',&
                     string_val)
       if (trim(string_val) .eq. 'user') then
          call C%fluid%set_usr_inflow(C%usr%fluid_user_if)
       end if
    end if

    ! Setup user boundary conditions for the scalar.
    if (scalar) then
       call C%scalar%set_user_bc(C%usr%scalar_user_bc)
    end if

    !
    ! Setup initial conditions
    !
    call json_get(C%params, 'case.fluid.initial_condition.type',&
                  string_val)

    call neko_log%section("Fluid initial condition ")
    call neko_log%message("Type: " // trim(string_val))

    if (trim(string_val) .ne. 'user') then
       call set_flow_ic(C%fluid%u, C%fluid%v, C%fluid%w, C%fluid%p, &
            C%fluid%c_Xh, C%fluid%gs_Xh, string_val, C%params)
    else
       call set_flow_ic(C%fluid%u, C%fluid%v, C%fluid%w, C%fluid%p, &
            C%fluid%c_Xh, C%fluid%gs_Xh, C%usr%fluid_user_ic, C%params)
    end if

    call neko_log%end_section()

    if (scalar) then

       call json_get(C%params, 'case.scalar.initial_condition.type', string_val)

       call neko_log%section("Scalar initial condition ")
       call neko_log%message("Type: " // trim(string_val))

       if (trim(string_val) .ne. 'user') then
          call set_scalar_ic(C%scalar%s, &
            C%scalar%c_Xh, C%scalar%gs_Xh, string_val, C%params)
       else
          call set_scalar_ic(C%scalar%s, &
            C%scalar%c_Xh, C%scalar%gs_Xh, C%usr%scalar_user_ic, C%params)
       end if

       call neko_log%end_section()

    end if

    ! Add initial conditions to BDF scheme (if present)
    select type (f => C%fluid)
    type is (fluid_pnpn_t)
       call f%ulag%set(f%u)
       call f%vlag%set(f%v)
       call f%wlag%set(f%w)
    end select

    !
    ! Validate that the case is properly setup for time-stepping
    !
    call C%fluid%validate

    if (scalar) then
       call C%scalar%slag%set(C%scalar%s)
       call C%scalar%validate
    end if

    !
    ! Get and process output directory
    !
    call json_get_or_default(C%params, 'case.output_directory',&
                             C%output_directory, '')

    output_dir_len = len(trim(C%output_directory))
    if (output_dir_len .gt. 0) then
       if (C%output_directory(output_dir_len:output_dir_len) .ne. "/") then
          C%output_directory = trim(C%output_directory)//"/"
          if (pe_rank .eq. 0) then
             call execute_command_line('mkdir -p '//C%output_directory)
          end if
       end if
    end if

    !
    ! Save boundary markings for fluid (if requested)
    !
    call json_get_or_default(C%params, 'case.output_boundary',&
                             logical_val, .false.)
    if (logical_val) then
       bdry_file = file_t(trim(C%output_directory)//'bdry.fld')
       call bdry_file%write(C%fluid%bdry)
    end if

    !
    ! Save mesh partitions (if requested)
    !
    call json_get_or_default(C%params, 'case.output_partitions',&
                             logical_val, .false.)
    if (logical_val) then
       call mesh_field_init(msh_part, C%msh, 'MPI_Rank')
       msh_part%data = pe_rank
       part_file = file_t(trim(C%output_directory)//'partitions.vtk')
       call part_file%write(msh_part)
       call mesh_field_free(msh_part)
    end if

    !
    ! Setup output precision of the field files
    !
    call json_get_or_default(C%params, 'case.output_precision', string_val,&
         'single')

    if (trim(string_val) .eq. 'double') then
       precision = dp
    else
       precision = sp
    end if

    !
    ! Setup sampler
    !
    call C%s%init(C%end_time)
    if (scalar) then
       C%f_out = fluid_output_t(precision, C%fluid, C%scalar, &
            path = trim(C%output_directory))
    else
       C%f_out = fluid_output_t(precision, C%fluid, &
            path = trim(C%output_directory))
    end if

    call json_get_or_default(C%params, 'case.fluid.output_control',&
                             string_val, 'org')

    if (trim(string_val) .eq. 'org') then
       ! yes, it should be real_val below for type compatibility
       call json_get(C%params, 'case.nsamples', real_val)
       call C%s%add(C%f_out, real_val, 'nsamples')
    else if (trim(string_val) .eq. 'never') then
       ! Fix a dummy 0.0 output_value
       call json_get_or_default(C%params, 'case.fluid.output_value', real_val, &
                                0.0_rp)
       call C%s%add(C%f_out, 0.0_rp, string_val)
    else
       call json_get(C%params, 'case.fluid.output_value', real_val)
       call C%s%add(C%f_out, real_val, string_val)
    end if

    !
    ! Save checkpoints (if nothing specified, default to saving at end of sim)
    !
    call json_get_or_default(C%params, 'case.output_checkpoints',&
                             logical_val, .true.)
    if (logical_val) then
       call json_get_or_default(C%params, 'case.checkpoint_format', &
            string_val, "chkp")
       C%f_chkp = chkp_output_t(C%fluid%chkp, path = C%output_directory, &
            fmt = trim(string_val))
       call json_get_or_default(C%params, 'case.checkpoint_control', &
            string_val, "simulationtime")
       call json_get_or_default(C%params, 'case.checkpoint_value', real_val,&
            1e10_rp)
       call C%s%add(C%f_chkp, real_val, string_val)
    end if

    !
    ! Setup joblimit
    !
    if (C%params%valid_path('case.job_timelimit')) then
       call json_get(C%params, 'case.job_timelimit', string_val)
       call jobctrl_set_time_limit(string_val)
    end if

    call neko_log%end_section()

  end subroutine case_init_common

  !> Deallocate a case
  subroutine case_free(C)
    type(case_t), intent(inout) :: C

    if (allocated(C%fluid)) then
       call C%fluid%free()
       deallocate(C%fluid)
    end if

    if (allocated(C%scalar)) then
       call C%scalar%free()
       deallocate(C%scalar)
    end if

    call C%msh%free()

    call C%s%free()

  end subroutine case_free

end module case
