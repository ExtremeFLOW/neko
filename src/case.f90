! Copyright (c) 2020-2021, The Neko Authors
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
  use scalar_output
  use chkp_output
  use mean_sqr_flow_output
  use mean_flow_output
  use fluid_stats_output
  use field_list_output
  use parameters
  use mpi_types
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
  use time_scheme_controller
  use logger
  use jobctrl
  use user_intf  
  use scalar_pnpn ! todo directly load the pnpn? can we have other
  use json_module, only : json_file_t => json_file, json_value_t => json_value,&
                          json_core_t => json_core
  implicit none

  type :: case_t
     type(mesh_t) :: msh
     type(param_t) :: params
     type(json_file_t) :: json_params
     type(time_scheme_controller_t) :: ext_bdf
     real(kind=rp), dimension(10) :: tlag
     real(kind=rp), dimension(10) :: dtlag
     real(kind=rp) :: dt
     real(kind=rp) :: end_time
     type(sampler_t) :: s
     type(fluid_output_t) :: f_out
     type(fluid_stats_output_t) :: f_stats_output
     type(scalar_output_t) :: s_out
     type(chkp_output_t) :: f_chkp
     type(mean_flow_output_t) :: f_mf
     type(mean_sqr_flow_output_t) :: f_msqrf
     type(stats_t) :: q   
     type(user_t) :: usr
     class(fluid_scheme_t), allocatable :: fluid
     type(scalar_pnpn_t), allocatable :: scalar 
  end type case_t

contains

  !> Initialize a case from an input file @a case_file
  subroutine case_init(C, case_file)
    implicit none
    type(case_t), target, intent(inout) :: C
    character(len=*), intent(in) :: case_file

    ! Namelist for case description
    character(len=:), allocatable :: mesh_file_name
    character(len=:), allocatable :: output_directory
    character(len=80) :: fluid_scheme  = ''
    character(len=80) :: mesh_file  = ''
    character(len=80) :: source_term = ''
    character(len=80) :: initial_condition = ''
    integer :: lx = 0
    logical :: scalar = .false.
    character(len=80) :: scalar_source_term = ''
    type(param_io_t) :: params
    namelist /NEKO_CASE/ mesh_file, fluid_scheme, lx,  &
         source_term, initial_condition, scalar, scalar_source_term
    
    integer :: ierr
    type(file_t) :: msh_file, bdry_file, part_file
    type(mesh_fld_t) :: msh_part
    integer, parameter :: nbytes = NEKO_FNAME_LEN + (4 * 80) + 4 + 4
    character buffer(nbytes)
    integer :: pack_index
    type(mesh_fld_t) :: parts

    type(json_file_t) :: json_file
    type(json_file_t) :: json_file2
    type(json_value_t), pointer :: json_value
    type(json_core_t) :: json_core
    logical found, logical_val
    integer integer_val
    real(kind=rp) real_val
    character(len=:), allocatable :: string_val
    character(len=:), allocatable :: json_buffer
    real(kind=rp) :: stats_start_time, stats_output_val
    integer ::  stats_sampling_interval 
   
    call neko_log%section('Case')
    call neko_log%message('Reading case file ' // trim(case_file))
    
    call json_file%Initialize()
    call json_file%load(filename='hemi.json')
    call C%json_params%load(filename='hemi.json')

!    call json_core%serialize(json_value, json_buffer)
    !call json_core%create_object(json_value, '')
    !call json_core%print(json_value,'test.json')
    !json_file2 = json_file
!    call json_file2%deserialize(json_buffer)
!    call json_file2%print()

!    write(*,*) "JSON TEST", json_buffer
!    write(*,*) "JSON TEST", int_val
!    write(*,*) "JSON TEST", bool_val
    
    
    !
    ! Read case description
    !
    if (pe_rank .eq. 0) then
       open(10, file=trim(case_file))
       read(10, nml=NEKO_CASE)
       read(10, *) params
       close(10)
       
       pack_index = 0
       call MPI_Pack(mesh_file, NEKO_FNAME_LEN, MPI_CHARACTER, &
            buffer, nbytes, pack_index, NEKO_COMM, ierr)
       call MPI_Pack(fluid_scheme, 80, MPI_CHARACTER, &
            buffer, nbytes, pack_index, NEKO_COMM, ierr)
       call MPI_Pack(source_term, 80, MPI_CHARACTER, &
            buffer, nbytes, pack_index, NEKO_COMM, ierr)
       call MPI_Pack(initial_condition, 80, MPI_CHARACTER, &
            buffer, nbytes, pack_index, NEKO_COMM, ierr)
       call MPI_Pack(scalar_source_term, 80, MPI_CHARACTER, &
            buffer, nbytes, pack_index, NEKO_COMM, ierr)
       call MPI_Pack(lx, 1, MPI_INTEGER, &
            buffer, nbytes, pack_index, NEKO_COMM, ierr)
       call MPI_Pack(scalar, 1, MPI_LOGICAL, &
            buffer, nbytes, pack_index, NEKO_COMM, ierr)
       call MPI_Bcast(buffer, nbytes, MPI_PACKED, 0, NEKO_COMM, ierr)
       call MPI_Bcast(params%p, 1, MPI_NEKO_PARAMS, 0, NEKO_COMM, ierr)
    else
       call MPI_Bcast(buffer, nbytes, MPI_PACKED, 0, NEKO_COMM, ierr)
       pack_index = 0

       call MPI_Unpack(buffer, nbytes, pack_index, &
            mesh_file, NEKO_FNAME_LEN, MPI_CHARACTER, NEKO_COMM, ierr)
       call MPI_Unpack(buffer, nbytes, pack_index, &
            fluid_scheme, 80, MPI_CHARACTER, NEKO_COMM, ierr)
       call MPI_Unpack(buffer, nbytes, pack_index, &
            source_term, 80, MPI_CHARACTER, NEKO_COMM, ierr)
       call MPI_Unpack(buffer, nbytes, pack_index, &
            initial_condition, 80, MPI_CHARACTER, NEKO_COMM, ierr)
       call MPI_Unpack(buffer, nbytes, pack_index, &
            scalar_source_term, 80, MPI_CHARACTER, NEKO_COMM, ierr)
       call MPI_Unpack(buffer, nbytes, pack_index, &
            lx, 1, MPI_INTEGER, NEKO_COMM, ierr)
       call MPI_Unpack(buffer, nbytes, pack_index, &
            scalar, 1, MPI_LOGICAL, NEKO_COMM, ierr)
       call MPI_Bcast(params%p, 1, MPI_NEKO_PARAMS, 0, NEKO_COMM, ierr)
    end if
    

    call C%json_params%get('case.mesh_file', string_val, found)
    if (.not. found) then
       call neko_error("Parameter mesh_file missing in the case file")
    end if
    msh_file = file_t(string_val)
    
    call msh_file%read(C%msh)
    C%params = params%p
    
    write(*,*) "BC LABELS", C%params%bc_labels(1)

    !
    ! Load Balancing
    !
    call C%json_params%get('case.load_balancing', logical_val, found)
    if (.not. found) logical_val = .false.

    if (pe_size .gt. 1 .and. logical_val) then
       call neko_log%section('Load Balancing')
       call parmetis_partmeshkway(C%msh, parts)
       call redist_mesh(C%msh, parts)
       call neko_log%end_section()       
    end if

    !
    ! Time step
    !
    call C%json_params%get('case.timestep', real_val, found)
    if (.not. found) then
       call neko_error("Parameter case.timestep missing in the case file")
    end if
    C%dt = real_val

    !
    ! End time
    !
    call C%json_params%get('case.end_time', real_val, found)
    if (.not. found) then
       call neko_error("Parameter case.end_time missing in the case file")
    end if
    C%end_time = real_val

    !
    ! Setup user defined functions
    !
    call C%usr%init()
    call C%usr%user_mesh_setup(C%msh)
    
    !
    ! Setup fluid scheme
    !
    call C%json_params%get('case.fluid.scheme', string_val, found)
    if (.not. found) then
       call neko_error("Parameter fluid.scheme missing in the case file")
    end if
    call fluid_scheme_factory(C%fluid, trim(string_val))

    ! Switch to JSON
    call C%json_params%get('case.numerics.polynomial_order', lx, found)
    if (.not. found) then
       call neko_error( &
         "Parameter numerics.polynomial_order missing in the case file")
    else 
      lx = lx + 1 ! add 1 to poly order
    end if
    call C%fluid%init(C%msh, lx, C%json_params)

    !
    ! Setup scalar scheme
    !
    ! @todo no scalar factroy for now, probably not needed
    call json_file%get('case.scalar', json_value, found)
    if (found) then
       call C%json_params%get('case.scalar.enabled', scalar, found)
       if (found .and. scalar) then
          allocate(C%scalar)
          ! Switch to json
          call C%scalar%init(C%msh, C%fluid%c_Xh, C%fluid%gs_Xh, C%params)
          call C%fluid%chkp%add_scalar(C%scalar%s)
       end if
    end if
    !
    ! Setup user defined conditions    
    !
    call C%json_params%get('case.fluid.inflow_condition.type', string_val, found)
    if (found .and. trim(string_val) .eq. 'user') then
       call C%fluid%set_usr_inflow(C%usr%fluid_user_if)
    end if
    
    !
    ! Setup source term
    ! 
    call C%json_params%get('case.fluid.source_term.type', string_val, found)
    if (found .and. trim(string_val) .eq. 'user') then
       call C%fluid%set_source(trim(string_val), usr_f=C%usr%fluid_user_f)
    else if (found .and. trim(string_val) .eq. 'user_vector') then
       call C%fluid%set_source(trim(string_val), &
            usr_f_vec=C%usr%fluid_user_f_vector)
    else if (found) then
       call C%fluid%set_source(trim(string_val))
    else
       if (.not. found) then
          call neko_error("Parameter end_time missing in the case file")
       end if
    end if

    ! Setup source term for the scalar
    ! @todo should be expanded for user sources etc. Now copies the fluid one
    call C%json_params%get('case.scalar.source_term.type', string_val, found)
    if (scalar) then
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
    call C%json_params%get('case.fluid.intial_condition.type', string_val, found)
    if (found) then
       if (trim(string_val) .ne. 'user') then
          ! Switch to json
          call set_flow_ic(C%fluid%u, C%fluid%v, C%fluid%w, C%fluid%p, &
               C%fluid%c_Xh, C%fluid%gs_Xh, initial_condition, C%params)
       else
          call set_flow_ic(C%fluid%u, C%fluid%v, C%fluid%w, C%fluid%p, &
               C%fluid%c_Xh, C%fluid%gs_Xh, C%usr%fluid_user_ic, C%params)
       end if
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
    call C%json_params%get('case.numerics.time_order', integer_val, found)
    call C%ext_bdf%init(integer_val)


    !
    ! Get output directory
    !
    call C%json_params%get('case.output_directory', output_directory, found)
    if (.not. found) then
       output_directory = ""
    else if (len(trim(output_directory)) .ne. 0) then
       ! append forward-slash
       output_directory = trim(output_directory)//"/"
    end if
    
    !
    ! Save boundary markings for fluid (if requested)
    ! 
    call C%json_params%get('case.output_boundary', logical_val, found)
    if (found .and. logical_val) then
       bdry_file = file_t(trim(output_directory)//'bdry.fld')
       call bdry_file%write(C%fluid%bdry)
    end if

    !
    ! Save mesh partitions (if requested)
    !
    call C%json_params%get('case.output_partitions', logical_val, found)
    if (found .and. logical_val) then
       call mesh_field_init(msh_part, C%msh, 'MPI_Rank')
       msh_part%data = pe_rank
       part_file = file_t(trim(output_directory)//'partitions.vtk')
       call part_file%write(msh_part)
       call mesh_field_free(msh_part)
    end if

    !
    ! Setup sampler
    !
    call C%json_params%get('case.end_time', real_val, found)
    if (.not. found) then
       call neko_error("Parameter end_time missing in the case file")
    end if
    call C%s%init(real_val)

    C%f_out = fluid_output_t(C%fluid, path=output_directory)

    call C%json_params%get('case.fluid.output_control', string_val, found)
    if (.not. found) string_val = 'org'

    if (trim(string_val) .eq. 'org') then
       call C%json_params%get('case.nsamples', real_val, found)
       if (.not. found) then
          call neko_error("Parameter nsamples missing in the case file")
       end if
       call C%s%add(C%f_out, real_val, 'nsamples')
    else 
       call C%json_params%get('case.fluid.output_value', real_val, found)
       if (.not. found) then
          call neko_error("Parameter fuid.output_value missing in the case file")
       end if
       call C%s%add(C%f_out, real_val, string_val)
    end if
    
    if (scalar) then
       C%s_out = scalar_output_t(C%scalar, path=output_directory)

    end if

    !
    ! Save checkpoints (if nothing specified, default to saving at end of sim)
    !
    call C%json_params%get('case.output_checkpoints', logical_val, found)
    if (.not. found .or. logical_val) then
       C%f_chkp = chkp_output_t(C%fluid%chkp, path=output_directory)

       call C%json_params%get('case.checkpoint_control', string_val, found)
       if (.not. found) then
          call neko_error(&
            "Parameter checkpoint_control missing in the case file")
       end if

       call C%json_params%get('case.checkpoint_value', real_val, found)
       if (.not. found) then
          call neko_error("Parameter checkpoint_value missing in the case file")
       end if

       call C%s%add(C%f_chkp, real_val, string_val)
    end if

    !
    ! Setup statistics
    !
    call C%json_params%get('case.statistics', json_value, found)
    if (found) then
       call C%json_params%get('case.statistics.enabled', logical_val, found)
       if (.not. found .or. (found .and. logical_val)) then

          call C%json_params%get('case.statistics.start_time', stats_start_time,&
                                  found)
          if (.not. found) then
             call neko_error( &
               "Parameter statistics.start_time missing in the case file")
          end if

          call C%json_params%get('case.statistics.sampling_interval', &
                                 stats_sampling_interval, found)
          if (.not. found) then
             call neko_error( &
               "Parameter statistics.sampling_interval missing in the case file")
          end if
          call C%q%init(stats_start_time, stats_sampling_interval)
          

          call C%q%add(C%fluid%mean%u)
          call C%q%add(C%fluid%mean%v)
          call C%q%add(C%fluid%mean%w)
          call C%q%add(C%fluid%mean%p)

          C%f_mf = mean_flow_output_t(C%fluid%mean, stats_start_time, &
                                      path=output_directory)

          call C%json_params%get('case.statistics.output_control', string_val,&
                                 found)
          if (.not. found) then
             call neko_error(&
               "Parameter statistics.output_control missing in the case file")
          end if

          call C%json_params%get('case.statistics.output_value', & 
                                 stats_output_val, found)
          if (.not. found) then
             call neko_error(&
               "Parameter statistics.output_value missing in the case file")
          end if
       
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
    call C%json_params%get('case.job_timelimit', string_val, found)
    if (found) call jobctrl_set_time_limit(string_val)

    call neko_log%end_section()
    
  end subroutine case_init
  
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

    call mesh_free(C%msh)

    call C%s%free()

    call C%q%free()
    
  end subroutine case_free
  
end module case
