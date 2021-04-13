module comm
  use num_types
  use mpi_f08
  implicit none
  
  !> MPI communicator
  type(MPI_Comm) :: NEKO_COMM

  !> Real precision number mpi type, standard double precision
  type(MPI_Datatype) :: MPI_REAL_PRECISION = MPI_DOUBLE_PRECISION

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
    type(MPI_Group) :: neko_group
        

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
    if (rp .eq. sp) then
       MPI_REAL_PRECISION = MPI_REAL
    else if (rp .eq. dp) then
       MPI_REAL_PRECISION = MPI_DOUBLE_PRECISION
    else if (rp .eq. qp) then
       MPI_REAL_PRECISION = MPI_REAL16
    else 
       call neko_error('Chosen real precision (rp) not supported (yet)!')
    end if

  end subroutine comm_init

  subroutine comm_free
    integer :: ierr

    call MPI_Barrier(NEKO_COMM, ierr)
    call MPI_Comm_free(NEKO_COMM, ierr)
    call MPI_Finalize(ierr)

  end subroutine comm_free

end module comm
