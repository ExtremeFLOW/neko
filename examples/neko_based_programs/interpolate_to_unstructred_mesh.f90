program interpolate_to_unstructured_mesh
  use neko
  use json_utils
  use fld_io_controller, only : fld_io_controller_t

  implicit none

  !> Declare objects
  type(json_file) :: params
  type(fld_io_controller_t) :: fld_ioc

  !> declare variables
  integer :: i,j,k

  !> Initialize neko 
  call neko_init 

  !> Initialize the parameters
  call init_params(params)

  !> Initialize the file interpolator object
  call fld_ioc%init(params)

  do i=1, fld_ioc%number_of_files -1

     call fld_ioc%step()

  end do


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
