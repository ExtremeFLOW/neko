module data_processor
  use neko
  use json_utils
  use fld_io_controller, only : fld_io_controller_t
  use post_processing_library

  implicit none

  !> Declare type
  type, public :: data_processor_t     
     !> Types
     type(json_file) :: params
     type(fld_io_controller_t), allocatable :: fld_ioc(:)
     type(probes_t), allocatable :: pb(:)
  end type data_processor_t
    
contains
  
  !> Initialize reading the case object
  subroutine init_params(params)
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
     
     if (pe_rank.eq.0) then
        write(*,*) "contents of json file"     
        write(*,*) json_buffer
     end if

     deallocate(json_buffer)

  end subroutine init_params

  !> Initialize the reader arrays
  subroutine init_readers(this)
     class(data_processor_t), intent(inout) :: this
     !> Internal variables
     integer :: n_readers, i
     type(json_core) :: core
     type(json_value), pointer :: reader_object, reader_pointer 
     type(json_file) :: reader_subdict
     character(len=:), allocatable :: buffer
     logical :: found

     if (this%params%valid_path('case.field_reader')) then
    
        if (pe_rank .eq. 0) then
           write(*,*)  "Initializing data readers"
        end if

        call this%params%info('case.field_reader', n_children=n_readers)

        call this%params%get_core(core)
        call this%params%get('case.field_reader', reader_object, found)

        !> Allocating reader object array
        allocate(this%fld_ioc(n_readers))

        !> Initialize all the reader objects
        do i=1, n_readers

           ! Create a new json containing just the subdict for this simcomp
           call core%get_child(reader_object, i, reader_pointer, found)
           call core%print_to_string(reader_pointer, buffer)
           call reader_subdict%load_from_string(buffer)
           
           !> Initialize the object
           call this%fld_ioc(i)%init(reader_subdict)
           !> Initialize the data processors not related to reading io
           call this%fld_ioc(i)%init_data_processors(reader_subdict)

        end do
     end if
     
  end subroutine init_readers

  !> Initialize probes_standalone
  subroutine init_probes(this)
     class(data_processor_t), intent(inout) :: this
     !> Internal variables
     integer :: n_simcomps, i
     type(json_core) :: core
     type(json_value), pointer :: simcomp_object, comp_pointer 
     type(json_file) :: comp_subdict
     character(len=:), allocatable :: buffer
     logical :: found

     if (this%params%valid_path('case.simulation_components')) then
    
        if (pe_rank .eq. 0) then
           write(*,*)  "Initializing simulation components"
        end if

        call this%params%info('case.simulation_components', n_children=n_simcomps)

        call this%params%get_core(core)
        call this%params%get('case.simulation_components', simcomp_object, found)
     
        !> Allocating reader object array
        allocate(this%pb(n_simcomps))
        
        do i=1, n_simcomps

           ! Create a new json containing just the subdict for this simcomp
           call core%get_child(simcomp_object, i, comp_pointer, found)
           call core%print_to_string(comp_pointer, buffer)
           call comp_subdict%load_from_string(buffer)
           
           !> Initialize the object
           call this%pb(i)%init_standalone(comp_subdict, this%fld_ioc(1)%dof)

        end do
     end if

  end subroutine init_probes

end module data_processor
