module comm
  use mpi
  implicit none
  
  !> MPI communicator
  integer :: NEKO_COMM

  !> MPI rank
  integer :: pe_rank

  !> MPI size of communicator
  integer :: pe_size

  !> I/O node
  logical :: nio

contains
  subroutine comm_init
    integer :: ierr, neko_group
    logical :: initialized

    pe_rank = -1
    pe_size = 0
    nio = .false.

    call MPI_Initialized(initialized, ierr)
    
    if (.not.initialized) then
       call MPI_Init(ierr)       
    end if

    !> @todo Why not use MPI_Comm_dup?
    call MPI_Comm_group(MPI_COMM_WORLD, neko_group, ierr)
    call MPI_Comm_create(MPI_COMM_WORLD, neko_group, NEKO_COMM, ierr)
    call MPI_Group_free(neko_group, ierr)

    call MPI_Comm_rank(NEKO_COMM, pe_rank, ierr)
    call MPI_Comm_size(NEKO_COMM, pe_size, ierr)


  end subroutine comm_init

  subroutine comm_free
    integer :: ierr

    call MPI_Comm_free(NEKO_COMM, ierr)
    call MPI_Finalize(ierr)

  end subroutine comm_free

end module comm
