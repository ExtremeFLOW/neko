program interpolate_to_unstructured_mesh
  use neko
  use json_utils, only : json_get

  implicit none

  !> Variables to read the input parameters
  type(json_file) :: params
  character(20) :: case_file = "inputs.json"
  integer :: ierr, integer_val
  character(len=:), allocatable :: json_buffer
 

  !> Test
  real(kind=rp) :: timestep

  !> Initialize neko 
  call neko_init 
 
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

   

    call json_get(params, 'case.test.test', timestep)
    write(*,*) timestep

    !call json_get_or_default(params, 'case.scalar.enabled', timestep,&
    !                            10)


  !> Finalize neko
  call neko_finalize

end program interpolate_to_unstructured_mesh
