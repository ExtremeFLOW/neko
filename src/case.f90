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
     type(time_scheme_controller_t) :: ext_bdf
     real(kind=rp), dimension(10) :: tlag
     real(kind=rp), dimension(10) :: dtlag
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
    character(len=NEKO_FNAME_LEN) :: mesh_file = ''
    character(len=80) :: fluid_scheme  = ''
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
    logical found
    integer int_val
    logical bool_val
    character(len=:), allocatable :: string_val
    character(len=:), allocatable :: json_buffer
   
    call neko_log%section('Case')
    call neko_log%message('Reading case file ' // trim(case_file))
    
    call json_file%Initialize()
    call json_file%load(filename='hemi.json')
    call json_file%get('number', int_val, found) 
    call json_file%get('boolean', bool_val, found) 
    call json_file%get('object.a', string_val, found) 
    
    
    call json_file%get('object', json_value, found)
    call json_core%serialize(json_value, json_buffer)

    !call json_core%create_object(json_value, '')
    !call json_core%print(json_value,'test.json')
    
    !json_file2 = json_file
    call json_file2%deserialize(json_buffer)
    call json_file2%print()

    write(*,*) "JSON TEST", json_buffer
    write(*,*) "JSON TEST", int_val
    write(*,*) "JSON TEST", bool_val
    write(*,*) "JSON TEST", string_val
    
    
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

    msh_file = file_t(trim(mesh_file))
    call msh_file%read(C%msh)
    C%params = params%p

    !
    ! Load Balancing
    !
    if (pe_size .gt. 1 .and. C%params%loadb) then
       call neko_log%section('Load Balancing')
       call parmetis_partmeshkway(C%msh, parts)
       call redist_mesh(C%msh, parts)
       call neko_log%end_section()       
    end if


    !
    ! Setup user defined functions
    !
    call C%usr%init()
    call C%usr%user_mesh_setup(C%msh)
    
    !
    ! Setup fluid scheme
    !
    call fluid_scheme_factory(C%fluid, trim(fluid_scheme))
    call C%fluid%init(C%msh, lx, C%params)

    !
    ! Setup scalar scheme
    !
    ! @todo no scalar factroy for now, probably not needed
    if (scalar) then
       allocate(C%scalar)
       call C%scalar%init(C%msh, C%fluid%c_Xh, C%fluid%gs_Xh, C%params)
       call C%fluid%chkp%add_scalar(C%scalar%s)
    end if
    !
    ! Setup user defined conditions    
    !
    if (trim(C%params%fluid_inflow) .eq. 'user') then
       call C%fluid%set_usr_inflow(C%usr%fluid_user_if)
    end if
    
    !
    ! Setup source term
    ! 
    if (trim(source_term) .eq. 'user') then
       call C%fluid%set_source(trim(source_term), usr_f=C%usr%fluid_user_f)
    else if (trim(source_term) .eq. 'user_vector') then
       call C%fluid%set_source(trim(source_term), &
            usr_f_vec=C%usr%fluid_user_f_vector)
    else
       call C%fluid%set_source(trim(source_term))
    end if

    ! Setup source term for the scalar
    ! @todo should be expanded for user sources etc. Now copies the fluid one
    if (scalar) then
       if (trim(scalar_source_term) .eq. 'user') then
          call C%scalar%set_source(trim(scalar_source_term), &
               usr_f=C%usr%scalar_user_f)
       else if (trim(scalar_source_term) .eq. 'user_vector') then
          call C%scalar%set_source(trim(scalar_source_term), &
               usr_f_vec=C%usr%scalar_user_f_vector)
       else
          call C%scalar%set_source(trim(scalar_source_term))
       end if
    end if

    !
    ! Setup initial conditions
    ! 
    if (len_trim(initial_condition) .gt. 0) then
       if (trim(initial_condition) .ne. 'user') then
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
    call C%ext_bdf%init(C%params%time_order)

    !
    ! Save boundary markings for fluid (if requested)
    ! 
    if (C%params%output_bdry) then
       bdry_file = file_t(trim(C%params%output_dir)//'bdry.fld')
       call bdry_file%write(C%fluid%bdry)
    end if

    !
    ! Save mesh partitions (if requested)
    !
    if (C%params%output_part) then
       call mesh_field_init(msh_part, C%msh, 'MPI_Rank')
       msh_part%data = pe_rank
       part_file = file_t(trim(C%params%output_dir)//'partitions.vtk')
       call part_file%write(msh_part)
       call mesh_field_free(msh_part)
    end if

    !
    ! Setup sampler
    !
    call C%s%init(C%params%T_end)
    C%f_out = fluid_output_t(C%fluid, path=C%params%output_dir)
    if (trim(C%params%fluid_write_control) .eq. 'org') then
       call C%s%add(C%f_out, real(C%params%nsamples,rp), 'nsamples')
    else 
       call C%s%add(C%f_out, C%params%fluid_write_par, &
            C%params%fluid_write_control)
    end if
    
    if (scalar) then
       C%s_out = scalar_output_t(C%scalar, path=C%params%output_dir)
       call C%s%add(C%s_out, C%params%fluid_write_par, &
            C%params%fluid_write_control)
    end if

    !
    ! Save checkpoints (if nothing specified, default to saving at end of sim)
    !
    if (C%params%output_chkp) then
       C%f_chkp = chkp_output_t(C%fluid%chkp, path=C%params%output_dir)
       call C%s%add(C%f_chkp, C%params%chkp_write_par, &
            C%params%chkp_write_control)
    end if

    !
    ! Setup statistics
    !
    call C%q%init(C%params%stats_begin, C%params%stats_sample_nstep)

    if (C%params%stats_mean_flow .or. C%params%stats_fluid) then
       call C%q%add(C%fluid%mean%u)
       call C%q%add(C%fluid%mean%v)
       call C%q%add(C%fluid%mean%w)
       call C%q%add(C%fluid%mean%p)

       C%f_mf = mean_flow_output_t(C%fluid%mean, C%params%stats_begin, &
                                   path=C%params%output_dir)
       call C%s%add(C%f_mf, C%params%stats_write_par, &
                    C%params%stats_write_control)
    end if
    
    if (C%params%stats_fluid) then
       call C%q%add(C%fluid%stats)
       C%f_stats_output = fluid_stats_output_t(C%fluid%stats, &
            C%params%stats_begin, path=C%params%output_dir)
       call C%s%add(C%f_stats_output, C%params%stats_write_par, &
            C%params%stats_write_control)
    end if

    if (C%params%stats_mean_sqr_flow) then
       call C%q%add(C%fluid%mean_sqr%uu)
       call C%q%add(C%fluid%mean_sqr%vv)
       call C%q%add(C%fluid%mean_sqr%ww)
       call C%q%add(C%fluid%mean_sqr%pp)

       if (C%params%output_mean_sqr_flow) then
          C%f_msqrf = mean_sqr_flow_output_t(C%fluid%mean_sqr, &
                                             C%params%stats_begin, &
                                             path=C%params%output_dir)
          call C%s%add(C%f_msqrf, C%params%stats_write_par, &
               C%params%stats_write_control)
       end if
    end if

    !
    ! Setup joblimit
    !
    call jobctrl_set_time_limit(C%params%jlimit)

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
