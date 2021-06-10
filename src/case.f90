!> Defines a simulation case
module case
  use num_types
  use fluid_schemes
  use fluid_output
  use parameters
  use mpi_types
  use mesh_field
  use parmetis
  use redist
  use sampler
  use file
  use utils
  use mesh
  use comm
  use abbdf
  use log    
  use user_intf  
  implicit none

  type :: case_t
     type(mesh_t) :: msh
     type(param_t) :: params
     type(abbdf_t) :: ab_bdf
     real(kind=rp), dimension(10) :: tlag
     real(kind=rp), dimension(10) :: dtlag
     type(sampler_t) :: s
     type(fluid_output_t) :: f_out
     type(user_t) :: usr
     class(fluid_scheme_t), allocatable :: fluid
  end type case_t

contains

  !> Initialize a case from an input file @a case_file
  subroutine case_init(C, case_file)
    type(case_t), intent(inout) :: C
    character(len=*), intent(in) :: case_file

    ! Namelist for case description
    character(len=NEKO_FNAME_LEN) :: mesh_file = ''
    character(len=80) :: fluid_scheme  = ''
    character(len=80) :: source_term = ''
    character(len=80) :: initial_condition = ''
    integer :: lx = 0
    type(param_io_t) :: params
    namelist /NEKO_CASE/ mesh_file, fluid_scheme, lx,  &
         source_term, initial_condition
    
    integer :: ierr
    type(file_t) :: msh_file, bdry_file, part_file
    type(mesh_fld_t) :: msh_part
    integer, parameter :: nbytes = NEKO_FNAME_LEN + 240 + 8
    character buffer(nbytes)
    integer :: pack_index, temp, i
    type(mesh_fld_t) :: parts
    real(kind=rp) :: eps, uvw(3)
    
    call neko_log%section('Case')
    call neko_log%message('Reading case file ' // trim(case_file))
    
    !
    ! Read case description
    !
    
    if (pe_rank .eq. 0) then
       open(10, file=trim(case_file))
       read(10, nml=NEKO_CASE)
       read(10, *) params
       close(10)
       
       pack_index = 1
       call MPI_Pack(mesh_file, NEKO_FNAME_LEN, MPI_CHARACTER, &
            buffer, nbytes, pack_index, NEKO_COMM, ierr)
       call MPI_Pack(fluid_scheme, 80, MPI_CHARACTER, &
            buffer, nbytes, pack_index, NEKO_COMM, ierr)
       call MPI_Pack(source_term, 80, MPI_CHARACTER, &
            buffer, nbytes, pack_index, NEKO_COMM, ierr)
       call MPI_Pack(initial_condition, 80, MPI_CHARACTER, &
            buffer, nbytes, pack_index, NEKO_COMM, ierr)
       call MPI_Pack(lx, 1, MPI_INTEGER, &
            buffer, nbytes, pack_index, NEKO_COMM, ierr)
       call MPI_Bcast(buffer, nbytes, MPI_PACKED, 0, NEKO_COMM, ierr)
       call MPI_Bcast(params%p, 1, MPI_NEKO_PARAMS, 0, NEKO_COMM, ierr)
    else
       call MPI_Bcast(buffer, nbytes, MPI_PACKED, 0, NEKO_COMM, ierr)
       pack_index = 1

       call MPI_Unpack(buffer, nbytes, pack_index, &
            mesh_file, NEKO_FNAME_LEN, MPI_CHARACTER, NEKO_COMM, ierr)
       call MPI_Unpack(buffer, nbytes, pack_index, &
            fluid_scheme, 80, MPI_CHARACTER, NEKO_COMM, ierr)
       call MPI_Unpack(buffer, nbytes, pack_index, &
            source_term, 80, MPI_CHARACTER, NEKO_COMM, ierr)
       call MPI_Unpack(buffer, nbytes, pack_index, &
            initial_condition, 80, MPI_CHARACTER, NEKO_COMM, ierr)
       call MPI_Unpack(buffer, nbytes, pack_index, &
            lx, 1, MPI_INTEGER, NEKO_COMM, ierr)
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
    call C%usr%usr_msh_setup(C%msh)
    
    !
    ! Setup fluid scheme
    !

    if (trim(fluid_scheme) .eq. 'plan1') then
       allocate(fluid_plan1_t::C%fluid)
    else if (trim(fluid_scheme) .eq. 'plan4') then
       allocate(fluid_plan4_t::C%fluid)
    else
       call neko_error('Invalid fluid scheme')
    end if
  
    call C%fluid%init(C%msh, lx, C%params)

    !
    ! Setup user defined conditions    
    !
    if (trim(C%params%fluid_inflow) .eq. 'user') then
       call C%fluid%set_usr_inflow(C%usr%fluid_usr_if)
    end if
    
    !
    ! Setup source term
    ! 
    
    !> @todo We shouldn't really mess with other type's datatypes
    if (trim(source_term) .eq. 'noforce') then
       call source_set_type(C%fluid%f_Xh, source_eval_noforce)
    else if (trim(source_term) .eq. 'user') then
       call source_set_pw_type(C%fluid%f_Xh, C%usr%fluid_usr_f)
    else if (trim(source_term) .eq. '') then
       if (pe_rank .eq. 0) then
          call neko_warning('No source term defined, using default (noforce)')
       end if
       call source_set_type(C%fluid%f_Xh, source_eval_noforce)
    else
       call neko_error('Invalid source term')
    end if

    !
    ! Setup initial conditions
    ! 
    
    !> @todo We shouldn't really mess with other type's datatypes
    if (len_trim(initial_condition) .gt. 0) then
       if (trim(initial_condition) .eq. 'uniform') then
          C%fluid%u = C%params%uinf(1)
          C%fluid%v = C%params%uinf(2)
          C%fluid%w = C%params%uinf(3)
       else if (trim(initial_condition) .eq. 'user') then
          call C%usr%fluid_usr_ic(C%fluid%u, C%fluid%v, &
               C%fluid%w, C%fluid%p, C%params)
       else
          call neko_error('Invalid initial condition')
       end if
    end if

    ! Ensure continuity across elements for initial conditions
    call gs_op_vector(C%fluid%gs_Xh, C%fluid%u%x, C%fluid%dm_Xh%n_dofs, GS_OP_ADD) 
    call col2(C%fluid%u%x, C%fluid%c_Xh%mult,C%fluid%dm_Xh%n_dofs) 
    call gs_op_vector(C%fluid%gs_Xh, C%fluid%v%x, C%fluid%dm_Xh%n_dofs, GS_OP_ADD) 
    call col2(C%fluid%v%x, C%fluid%c_Xh%mult,C%fluid%dm_Xh%n_dofs) 
    call gs_op_vector(C%fluid%gs_Xh, C%fluid%w%x, C%fluid%dm_Xh%n_dofs, GS_OP_ADD) 
    call col2(C%fluid%w%x, C%fluid%c_Xh%mult,C%fluid%dm_Xh%n_dofs) 

    ! Add initial conditions to BDF scheme (if present)
    select type(f => C%fluid)
    type is(fluid_plan4_t)
       call copy(f%ulag, C%fluid%u%x, C%fluid%dm_Xh%n_dofs)
       call copy(f%vlag, C%fluid%v%x, C%fluid%dm_Xh%n_dofs)
       call copy(f%wlag, C%fluid%w%x, C%fluid%dm_Xh%n_dofs)
    end select

    !
    ! Validate that the case is properly setup for time-stepping
    !
    call C%fluid%validate

    !
    ! Set order of timestepper
    !
    call C%ab_bdf%set_time_order(C%params%time_order)

    !
    ! Save boundary markings for fluid (if requested)
    ! 
    if (C%params%output_bdry) then
       bdry_file = file_t('bdry.fld')
       call bdry_file%write(C%fluid%bdry)
    end if

    !
    ! Save mesh partitions (if requested)
    !
    if (C%params%output_part) then
       call mesh_field_init(msh_part, C%msh, 'MPI_Rank')
       msh_part%data = pe_rank
       part_file = file_t('partitions.vtk')
       call part_file%write(msh_part)
       call mesh_field_free(msh_part)
    end if

    !
    ! Setup sampler
    !
    call C%s%init(C%params%nsamples, C%params%T_end)
    C%f_out = fluid_output_t(C%fluid)
    call C%s%add(C%f_out)

    call neko_log%end_section()
    
  end subroutine case_init
  !> Deallocate a case 
  subroutine case_free(C)
    type(case_t), intent(inout) :: C

    if (allocated(C%fluid)) then
       call C%fluid%free()
       deallocate(C%fluid)
    end if

    call mesh_free(C%msh)

    call C%s%free()
    
  end subroutine case_free
  
end module case
