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
  use output_controller, only : output_controller_t
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
     type(output_controller_t) :: output_controller
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
  subroutine case_init_from_file(this, case_file)
    type(case_t), target, intent(inout) :: this
    character(len=*), intent(in) :: case_file
    integer :: ierr, integer_val
    character(len=:), allocatable :: json_buffer

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
    call json_get_or_default(this%params, 'case.mesh_file', string_val, &
         'no mesh')
    if (trim(string_val) .eq. 'no mesh') then
       call neko_error('The mesh_file keyword could not be found in the .' // &
                       'case file. Often caused by incorrectly formatted json.')
    end if
    msh_file = file_t(string_val)

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
       call neko_log%end_section()
    end if

    !
    ! Time step
    !
    call this%params%get('case.variable_timestep', logical_val, found)
    if (.not. logical_val) then
       call json_get(this%params, 'case.timestep', this%dt)
    else
       ! randomly set an initial dt to get cfl when dt is variable
       this%dt = 1.0_rp
    end if

    !
    ! End time
    !
    call json_get(this%params, 'case.end_time', this%end_time)

    !
    ! Initialize point_zones registry
    !
    call neko_point_zone_registry%init(this%params, this%msh)

    !
    ! Setup user defined functions
    !
    call this%usr%init()
    call this%usr%user_mesh_setup(this%msh)

    !
    ! Set order of timestepper
    !
    call json_get(this%params, 'case.numerics.time_order', integer_val)
    call this%ext_bdf%init(integer_val)

    !
    ! Setup fluid scheme
    !
    call json_get(this%params, 'case.fluid.scheme', string_val)
    call fluid_scheme_factory(this%fluid, trim(string_val))

    call json_get(this%params, 'case.numerics.polynomial_order', lx)
    lx = lx + 1 ! add 1 to get number of gll points
    this%fluid%chkp%tlag => this%tlag
    this%fluid%chkp%dtlag => this%dtlag
    call this%fluid%init(this%msh, lx, this%params, this%usr, this%ext_bdf)
    select type (f => this%fluid)
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
    neko_scratch_registry = scratch_registry_t(this%fluid%dm_Xh, 10, 10)

    !
    ! Setup scalar scheme
    !
    ! @todo no scalar factory for now, probably not needed
    if (this%params%valid_path('case.scalar')) then
       call json_get_or_default(this%params, 'case.scalar.enabled', scalar,&
                                .true.)
    end if

    if (scalar) then
       allocate(this%scalar)
       this%scalar%chkp%tlag => this%tlag
       this%scalar%chkp%dtlag => this%dtlag
       call this%scalar%init(this%msh, this%fluid%c_Xh, this%fluid%gs_Xh, &
            this%params, this%usr, this%fluid%ulag, this%fluid%vlag, &
            this%fluid%wlag, this%ext_bdf, this%fluid%rho)

       call this%fluid%chkp%add_scalar(this%scalar%s)

       this%fluid%chkp%abs1 => this%scalar%abx1
       this%fluid%chkp%abs2 => this%scalar%abx2
       this%fluid%chkp%slag => this%scalar%slag
    end if

    !
    ! Setup user defined conditions
    !
    if (this%params%valid_path('case.fluid.inflow_condition')) then
       call json_get(this%params, 'case.fluid.inflow_condition.type',&
                     string_val)
       if (trim(string_val) .eq. 'user') then
          call this%fluid%set_usr_inflow(this%usr%fluid_user_if)
       end if
    end if

    ! Setup user boundary conditions for the scalar.
    if (scalar) then
       call this%scalar%set_user_bc(this%usr%scalar_user_bc)
    end if

    !
    ! Setup initial conditions
    !
    call json_get(this%params, 'case.fluid.initial_condition.type',&
                  string_val)


    call neko_log%section("Fluid initial condition ")

    if (trim(string_val) .ne. 'user') then
       call set_flow_ic(this%fluid%u, this%fluid%v, this%fluid%w, &
            this%fluid%p, this%fluid%c_Xh, this%fluid%gs_Xh, string_val, &
            this%params)

    else
       call set_flow_ic(this%fluid%u, this%fluid%v, this%fluid%w, this%fluid%p,&
            this%fluid%c_Xh, this%fluid%gs_Xh, this%usr%fluid_user_ic, &
            this%params)
    end if

    call neko_log%end_section()

    if (scalar) then

       call json_get(this%params, 'case.scalar.initial_condition.type', &
            string_val)

       call neko_log%section("Scalar initial condition ")

       if (trim(string_val) .ne. 'user') then
          call set_scalar_ic(this%scalar%s, &
            this%scalar%c_Xh, this%scalar%gs_Xh, string_val, this%params)
       else
          call set_scalar_ic(this%scalar%s, &
            this%scalar%c_Xh, this%scalar%gs_Xh, this%usr%scalar_user_ic, &
            this%params)
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
       call this%scalar%slag%set(this%scalar%s)
       call this%scalar%validate
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
    ! Save boundary markings for fluid (if requested)
    !
    call json_get_or_default(this%params, 'case.output_boundary',&
                             logical_val, .false.)
    if (logical_val) then
       bdry_file = file_t(trim(this%output_directory)//'bdry.fld')
       call bdry_file%write(this%fluid%bdry)
    end if

    !
    ! Save mesh partitions (if requested)
    !
    call json_get_or_default(this%params, 'case.output_partitions',&
                             logical_val, .false.)
    if (logical_val) then
       call mesh_field_init(msh_part, this%msh, 'MPI_Rank')
       msh_part%data = pe_rank
       part_file = file_t(trim(this%output_directory)//'partitions.vtk')
       call part_file%write(msh_part)
       call mesh_field_free(msh_part)
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
    ! Setup output_controller
    !
    call this%output_controller%init(this%end_time)
    if (scalar) then
       this%f_out = fluid_output_t(precision, this%fluid, this%scalar, &
            path = trim(this%output_directory))
    else
       this%f_out = fluid_output_t(precision, this%fluid, &
            path = trim(this%output_directory))
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
       call json_get_or_default(this%params, 'case.checkpoint_format', &
            string_val, "chkp")
       this%f_chkp = chkp_output_t(this%fluid%chkp, &
            path = this%output_directory, fmt = trim(string_val))
       call json_get_or_default(this%params, 'case.checkpoint_control', &
            string_val, "simulationtime")
       call json_get_or_default(this%params, 'case.checkpoint_value', real_val,&
            1e10_rp)
       call this%output_controller%add(this%f_chkp, real_val, string_val)
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

    if (allocated(this%scalar)) then
       call this%scalar%free()
       deallocate(this%scalar)
    end if

    call this%msh%free()

    call this%output_controller%free()

  end subroutine case_free

end module case
