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
  use num_types
  use fluid_fctry
  use fluid_output
  use chkp_output
  use mean_sqr_flow_output
  use mean_flow_output
  use fluid_stats_output
  use mpi_types
  use mpi_f08
  use mesh_field
  use parmetis
  use redist
  use sampler
  use flow_ic    
  use stats
  use file
  use utils
  use mesh
  use comm
  use time_scheme_controller, only : time_scheme_controller_t
  use logger
  use jobctrl
  use user_intf  
  use scalar_pnpn ! todo directly load the pnpn? can we have other
  use json_module, only : json_file, json_core, json_value
  use json_utils, only : json_get, json_get_or_default
  use scratch_registry, only : scratch_registry_t, neko_scratch_registry
  use point_zone_registry, only: neko_point_zone_registry
  use material_properties, only : material_properties_t
  implicit none
  
  type :: case_t
     type(mesh_t) :: msh
     type(json_file) :: params
     type(time_scheme_controller_t) :: ext_bdf
     real(kind=rp), dimension(10) :: tlag
     real(kind=rp), dimension(10) :: dtlag
     real(kind=rp) :: dt
     real(kind=rp) :: end_time
     type(sampler_t) :: s
     type(fluid_output_t) :: f_out
     type(fluid_stats_output_t) :: f_stats_output
     type(chkp_output_t) :: f_chkp
     type(mean_flow_output_t) :: f_mf
     type(mean_sqr_flow_output_t) :: f_msqrf
     type(stats_t) :: q   
     type(user_t) :: usr
     class(fluid_scheme_t), allocatable :: fluid
     type(scalar_pnpn_t), allocatable :: scalar 
     type(material_properties_t):: material_properties
  end type case_t

  interface case_init
     module procedure case_init_from_file
  end interface case_init

  private :: case_init_from_file, case_init_from_json, case_init_common

contains

  !> Initialize a case from an input file @a case_file
  subroutine case_init_from_file(C, case_file)
    type(case_t), target, intent(inout) :: C
    character(len=*), intent(in) :: case_file
    integer :: ierr, integer_val
    character(len=:), allocatable :: json_buffer
   
    call neko_log%section('Case')
    call neko_log%message('Reading case file ' // trim(case_file))
    
    if (pe_rank .eq. 0) then
      call C%params%load_file(filename=trim(case_file))
      call C%params%print_to_string(json_buffer)
      integer_val = len(json_buffer)
    end if

    call MPI_Bcast(integer_val, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)
    if (pe_rank .ne. 0) allocate(character(len=integer_val)::json_buffer)
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
    call neko_log%message('Creating case from JSON object')

    C%params = case_json

    call case_init_common(C)
    
  end subroutine case_init_from_json

  !> Initialize a case from its (loaded) params object
  subroutine case_init_common(C)
    type(case_t), target, intent(inout) :: C
    character(len=:), allocatable :: output_directory
    integer :: lx = 0
    logical :: scalar = .false.
    type(file_t) :: msh_file, bdry_file, part_file
    type(mesh_fld_t) :: msh_part, parts
    logical :: found, logical_val
    integer :: integer_val
    real(kind=rp) :: real_val
    character(len=:), allocatable :: string_val
    real(kind=rp) :: stats_start_time, stats_output_val
    integer ::  stats_sampling_interval 
    integer :: output_dir_len
    integer :: n_simcomps, i
    type(json_core) :: core
    type(json_value), pointer :: json_val1, json_val2 
    type(json_file) :: json_subdict

    !
    ! Load mesh
    !
    call json_get(C%params, 'case.mesh_file', string_val)
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
    call json_get(C%params, 'case.timestep', C%dt)

    !
    ! End time
    !
    call json_get(C%params, 'case.end_time', C%end_time)

    !
    ! Setup user defined functions
    !
    call C%usr%init()
    call C%usr%user_mesh_setup(C%msh)


    !
    ! Material properties
    !
    call C%material_properties%init(C%params, C%usr)
    
    !
    ! Setup fluid scheme
    !
    call json_get(C%params, 'case.fluid.scheme', string_val)
    call fluid_scheme_factory(C%fluid, trim(string_val))

    call json_get(C%params, 'case.numerics.polynomial_order', lx)
    lx = lx + 1 ! add 1 to get poly order
    call C%fluid%init(C%msh, lx, C%params, C%usr, C%material_properties)

    
    !
    ! Setup scratch registry
    !
    neko_scratch_registry = scratch_registry_t(C%fluid%dm_Xh, 10, 10)

    !
    ! Initialize point_zones registry
    !
    call neko_point_zone_registry%init(C%params, C%fluid%u%dof)

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
       call C%scalar%init(C%msh, C%fluid%c_Xh, C%fluid%gs_Xh, C%params, C%usr,&
                          C%material_properties)
       call C%fluid%chkp%add_scalar(C%scalar%s)
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
    
    ! Setup source term for the scalar
    ! @todo should be expanded for user sources etc. Now copies the fluid one
    if (scalar) then
       logical_val = C%params%valid_path('case.scalar.source_term')
       call json_get_or_default(C%params, 'case.scalar.source_term.type',&
                                string_val, 'noforce')
       if (trim(string_val) .eq. 'user') then
          call C%scalar%set_source(trim(string_val), &
               usr_f=C%usr%scalar_user_f)
       else if (trim(string_val) .eq. 'user_vector') then
          call C%scalar%set_source(trim(string_val), &
               usr_f_vec=C%usr%scalar_user_f_vector)
       else
          call C%scalar%set_source(trim(string_val))
       end if

       call C%scalar%set_user_bc(C%usr%scalar_user_bc)
    end if

    !
    ! Setup initial conditions
    ! 
    call json_get(C%params, 'case.fluid.initial_condition.type',&
                  string_val)
    if (trim(string_val) .ne. 'user') then
       call set_flow_ic(C%fluid%u, C%fluid%v, C%fluid%w, C%fluid%p, &
            C%fluid%c_Xh, C%fluid%gs_Xh, string_val, C%params)
    else
       call set_flow_ic(C%fluid%u, C%fluid%v, C%fluid%w, C%fluid%p, &
            C%fluid%c_Xh, C%fluid%gs_Xh, C%usr%fluid_user_ic, C%params)
    end if

    ! Add initial conditions to BDF scheme (if present)
    select type(f => C%fluid)
    type is(fluid_pnpn_t)
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
    ! Set order of timestepper
    !
    call json_get(C%params, 'case.numerics.time_order', integer_val)
    call C%ext_bdf%init(integer_val)

    !
    ! Get and process output directory
    !
    call json_get_or_default(C%params, 'case.output_directory',&
                             output_directory, '')

    output_dir_len = len(trim(output_directory))
    if (output_dir_len .gt. 0) then
       if (output_directory(output_dir_len:output_dir_len) .ne. "/") then
          output_directory = trim(output_directory)//"/"
          if (pe_rank .eq. 0) then
             call execute_command_line('mkdir -p '//output_directory)
          end if
       end if
    end if
    
    !
    ! Save boundary markings for fluid (if requested)
    ! 
    call json_get_or_default(C%params, 'case.output_boundary',&
                             logical_val, .false.)
    if (logical_val) then
       bdry_file = file_t(trim(output_directory)//'bdry.fld')
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
       part_file = file_t(trim(output_directory)//'partitions.vtk')
       call part_file%write(msh_part)
       call mesh_field_free(msh_part)
    end if

    !
    ! Setup sampler
    !
    call C%s%init(C%end_time)
    if (scalar) then
       C%f_out = fluid_output_t(C%fluid, C%scalar, path=trim(output_directory))
    else
       C%f_out = fluid_output_t(C%fluid, path=trim(output_directory))
    end if

    call json_get_or_default(C%params, 'case.fluid.output_control',&
                             string_val, 'org')

    if (trim(string_val) .eq. 'org') then
       ! yes, it should be real_val below for type compatibility
       call json_get(C%params, 'case.nsamples', real_val)
       call C%s%add(C%f_out, real_val, 'nsamples')
    else if (trim(string_val) .eq. 'never') then
       ! Fix a dummy 0.0 output_value
       call C%s%add(C%f_out, 0.0_rp, string_val)
    else 
       call json_get(C%params, 'case.fluid.output_value', real_val)
       call C%s%add(C%f_out, real_val, string_val)
    end if   

    !
    ! Save checkpoints (if nothing specified, default to saving at end of sim)
    !
    call json_get_or_default(C%params, 'case.output_checkpoints',&
                             logical_val, .false.)
    if (logical_val) then
       C%f_chkp = chkp_output_t(C%fluid%chkp, path=output_directory)
       call json_get(C%params, 'case.checkpoint_control', string_val)
       call json_get(C%params, 'case.checkpoint_value', real_val)
       call C%s%add(C%f_chkp, real_val, string_val)
    end if

    !
    ! Setup statistics
    !
    
    ! Always init, so that we can call eval in simulation.f90 with no if.
    ! Note, don't use json_get_or_default here, because that will break the
    ! valid_path if statement below (the path will become valid always).
    call C%params%get('case.statistics.start_time', stats_start_time,&
                           found)
    if (.not. found) stats_start_time = 0.0_rp

    call C%params%get('case.statistics.sampling_interval', &
                           stats_sampling_interval, found)
    if (.not. found) stats_sampling_interval = 10

    call C%q%init(stats_start_time, stats_sampling_interval)

    found = C%params%valid_path('case.statistics')
    if (found) then
       call json_get_or_default(C%params, 'case.statistics.enabled',&
                                logical_val, .true.)
       if (logical_val) then
          call C%q%add(C%fluid%mean%u)
          call C%q%add(C%fluid%mean%v)
          call C%q%add(C%fluid%mean%w)
          call C%q%add(C%fluid%mean%p)

          C%f_mf = mean_flow_output_t(C%fluid%mean, stats_start_time, &
                                      path=output_directory)

          call json_get(C%params, 'case.statistics.output_control', &
                        string_val)
          call json_get(C%params, 'case.statistics.output_value', &
                        stats_output_val)
       
          call C%s%add(C%f_mf, stats_output_val, string_val)
          call C%q%add(C%fluid%stats)

          C%f_stats_output = fluid_stats_output_t(C%fluid%stats, &
            stats_start_time, path=output_directory)
          call C%s%add(C%f_stats_output, stats_output_val, string_val)
       end if
    end if

!    if (C%params%stats_mean_sqr_flow) then
!       call C%q%add(C%fluid%mean_sqr%uu)
!       call C%q%add(C%fluid%mean_sqr%vv)
!       call C%q%add(C%fluid%mean_sqr%ww)
!       call C%q%add(C%fluid%mean_sqr%pp)

!       if (C%params%output_mean_sqr_flow) then
!          C%f_msqrf = mean_sqr_flow_output_t(C%fluid%mean_sqr, &
!                                             C%params%stats_begin, &
!                                             path=output_directory)
!          call C%s%add(C%f_msqrf, C%params%stats_write_par, &
!               C%params%stats_write_control)
!       end if
!    end if

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

    call C%q%free()
    
  end subroutine case_free
  
end module case
