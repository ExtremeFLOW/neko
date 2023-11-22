program interpolate_to_unstructured_mesh
  use neko
  use json_utils
  use fld_io_controller, only : fld_io_controller_t

  implicit none

  !> Declare objects
  type(json_file) :: params
  type(fld_io_controller_t) :: fld_ioc
  type(probes_t) :: pb

  !> declare variables
  integer :: i,j,k

  !> Initialize neko 
  call neko_init 

  !> Initialize the parameters
  call init_params(params)

  !> Initialize the file interpolator object
  call fld_ioc%init(params)

  !> Initialize the standalone simulation components
  call init_simcomps_standalone(params,pb, fld_ioc)

  ! Interpolate the first field
  call pb%compute_(fld_ioc%t, fld_ioc%file_counter)

  do i=1, fld_ioc%number_of_files -1

     call fld_ioc%step()
       
     ! Interpolate
     call pb%compute_(fld_ioc%t, fld_ioc%file_counter)

  end do

  !> Finalize the probes
  call pb%free()

  !> Finalize neko
  call neko_finalize

end program interpolate_to_unstructured_mesh

subroutine init_params(params)
  use neko
  use json_utils
  
  type(json_file), intent(inout) :: params
  character(20) :: case_file = "inputs.json"
  integer :: ierr, integer_val
  character(len=:), allocatable :: json_buffer
  
  !> Read input file
    if (pe_rank .eq. 0) then
      write(*,*)  trim(case_file)
      call params%load_file(filename=trim(case_file))
      call params%print_to_string(json_buffer)
      integer_val = len(json_buffer)
    end if

    call MPI_Bcast(integer_val, 1, MPI_INTEGER, 0, NEKO_COMM, ierr)
    if (pe_rank .ne. 0) allocate(character(len=integer_val)::json_buffer)
    call MPI_Bcast(json_buffer, integer_val, MPI_CHARACTER, 0, NEKO_COMM, ierr)
    call params%load_from_string(json_buffer)

    deallocate(json_buffer)

end subroutine init_params


subroutine init_simcomps_standalone(params,pb, fld_ioc)
  use neko
  use fld_io_controller, only : fld_io_controller_t
  type(json_file), intent(inout) :: params
  type(probes_t), intent(inout) :: pb
  type(fld_io_controller_t), intent(inout) :: fld_ioc

  integer :: n_, i
  type(json_core) :: core
  type(json_value), pointer :: simcomp_object, comp_pointer 
  type(json_file) :: comp_subdict
  character(len=:), allocatable :: buffer
  logical :: found

  if (params%valid_path('case.simulation_components')) then
    
     if (pe_rank .eq. 0) then
        write(*,*)  "Initializing simulation components"
     end if

     call params%info('case.simulation_components', n_children=n_simcomps)

     call params%get_core(core)
     call params%get('case.simulation_components', simcomp_object, found)
     do i=1, n_simcomps
        ! Create a new json containing just the subdict for this simcomp
        call core%get_child(simcomp_object, i, comp_pointer, found)
        call core%print_to_string(comp_pointer, buffer)
        call comp_subdict%load_from_string(buffer)
        call pb%init_standalone(comp_subdict, fld_ioc%dof)
     end do
  end if

end subroutine init_simcomps_standalone
