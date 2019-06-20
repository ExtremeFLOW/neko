!> Defines practical data distributions
module datadist
  implicit none
  private

  type dist_t
     integer :: comm
     integer :: pe_rank
     integer :: pe_size
     integer :: L
     integer :: R
     integer :: M
     integer :: Ip
  end type dist_t

  !> Load-balanced linear distribution \f$ M = PL + R \f$
  type, extends(dist_t) :: linear_dist_t
   contains 
     procedure :: num_local => linear_dist_Ip
     procedure :: start_idx => linear_dist_start
     procedure :: end_idx => linear_dist_end
  end type linear_dist_t

  interface linear_dist_t
     module procedure linear_dist_init
  end interface linear_dist_t

  public :: linear_dist_t

contains
  
  function linear_dist_init(n, comm) result(this)
    integer, intent(in) :: n    !< Total size
    integer :: comm             !< comm. to define the dist. over
    type(linear_dist_t), target :: this
    integer :: ierr

    this%M = n
    this%comm = comm
    call MPI_Comm_size(this%comm, this%pe_size, ierr)
    call MPI_Comm_rank(this%comm, this%pe_rank, ierr)
    
    this%L = floor(dble(this%M) / dble(this%pe_size))
    this%R = modulo(this%M, this%pe_size)
    this%Ip = floor((dble(this%M) + dble(this%pe_size) - &
         dble(this%pe_rank) - 1d0) / dble(this%pe_size))    
  end function linear_dist_init

  pure function linear_dist_Ip(this) result(n)
    class(linear_dist_t), intent(in) :: this
    integer, intent(out) :: n
    n = this%Ip
  end function linear_dist_Ip

  pure function linear_dist_start(this) result(start)
    class(linear_dist_t), intent(in) :: this
    integer, intent(out) :: start
    start = this%pe_rank * this%L + min(this%pe_rank, this%R)
  end function linear_dist_start

  function linear_dist_end(this) result(end)
    class(linear_dist_t), intent(inout) :: this
    integer, intent(out) :: end
    end = linear_dist_start(this) + (this%Ip - 1)
  end function linear_dist_end
end module datadist
