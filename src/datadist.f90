!> Defines practical data distributions
module datadist
  implicit none
  private

  type dist_t
     integer :: comm            !< Communicator on which the dist. is defined
     integer :: pe_rank         !< Pe's rank in the given distribution
     integer :: pe_size         !< Size of communicator in the given dist.
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
  
  function linear_dist_init(n, rank, size, comm) result(this)
    integer, intent(in) :: n    !< Total size
    integer :: rank             !< PE's rank to define the dist. over
    integer :: size             !< Size of comm where the dist. is def. on
    integer :: comm             !< comm. to define the dist. over
    type(linear_dist_t), target :: this
    integer :: ierr

    this%M = n
    this%comm = comm
    this%pe_rank = rank
    this%pe_size = size
    
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
