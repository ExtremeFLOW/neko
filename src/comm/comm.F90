module comm
  use mpi_f08
  implicit none
  
  !> MPI communicator
  type(MPI_Comm) :: NEKO_COMM

  !> MPI type for working precision of REAL types
#ifdef HAVE_MPI_PARAM_DTYPE
  type(MPI_Datatype), parameter :: MPI_REAL_PRECISION = MPI_DOUBLE_PRECISION
#else
  type(MPI_Datatype) :: MPI_REAL_PRECISION
#endif
  
  !> MPI rank
  integer :: pe_rank

  !> MPI size of communicator
  integer :: pe_size

  !> I/O node
  logical :: nio

contains
  subroutine comm_init
    integer :: ierr
    logical :: initialized
    integer :: global_rank

    pe_rank = -1
    pe_size = 0
    nio = .false.

    call MPI_Initialized(initialized, ierr)
    
    if (.not.initialized) then
       call MPI_Init(ierr)       
    end if

#ifndef HAVE_MPI_PARAM_DTYPE
    MPI_REAL_PRECISION = MPI_DOUBLE_PRECISION
#endif

    !call MPI_Comm_dup(MPI_COMM_WORLD, NEKO_COMM, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, global_rank, ierr)
    call MPI_Comm_split(MPI_COMM_WORLD,0, global_rank,NEKO_COMM,ierr)
    call MPI_Comm_rank(NEKO_COMM, pe_rank, ierr)
    call MPI_Comm_size(NEKO_COMM, pe_size, ierr)

  end subroutine comm_init

  subroutine comm_free
    integer :: ierr

    call MPI_Barrier(NEKO_COMM, ierr)
    call MPI_Comm_free(NEKO_COMM, ierr)
    call MPI_Finalize(ierr)

  end subroutine comm_free

end module comm
