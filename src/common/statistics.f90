!> Defines a container for all statistics
module stats
  use num_types
  use stats_quant
  implicit none

  !> Pointer to an arbitrary quantitiy
  type, private :: quantp_t
     class(stats_quant_t), pointer :: quantp
  end type quantp_t

  !> Statistics backend
  type :: stats_t
     type(quantp_t), allocatable :: quant_list(:)
     integer :: n
     integer :: size
     real(kind=rp) :: T_begin
   contains
     procedure, pass(this) :: init => stats_init
     procedure, pass(this) :: free => stats_free
     procedure, pass(this) :: add => stats_add
     procedure, pass(this) :: eval => stats_eval
  end type stats_t

contains

  !> Initialize statistics, computed after @a T_begin
  subroutine stats_init(this, T_begin, size)
    class(stats_t), intent(inout) :: this
    real(kind=rp), intent(in) :: T_begin
    integer, intent(inout), optional ::size
    integer :: n, i
    
    call this%free()

    if (present(size)) then
       n = size
    else
       n = 1
    end if

    allocate(this%quant_list(n))

    do i = 1, n
       this%quant_list(i)%quantp => null()
    end do

    this%n = 0
    this%size = n
    this%T_begin = T_begin
    
  end subroutine stats_init

  !> Deallocate
  subroutine stats_free(this)
    class(stats_t), intent(inout) :: this

    if (allocated(this%quant_list)) then
       deallocate(this%quant_list)
    end if

    this%n = 0
    this%size = 0    
  end subroutine stats_free

  !> Add a statistic quantitiy @a quant to the backend
  subroutine stats_add(this, quant)
    class(stats_t), intent(inout) :: this
    class(stats_quant_t), intent(inout), target :: quant
    type(quantp_t), allocatable :: tmp(:)

    if (this%n .ge. this%size) then
       allocate(tmp(this%size * 2))
       tmp(1:this%size) = this%quant_list
       call move_alloc(tmp, this%quant_list)
       this%size = this%size * 2
    end if

    this%n = this%n + 1
    this%quant_list(this%n)%quantp => quant
  end subroutine stats_add

  !> Evaluated all statistical quantities
  subroutine stats_eval(this, t, k)
    class(stats_t), intent(inout) :: this
    real(kind=dp), intent(in) :: t
    real(kind=dp), intent(in) :: k
    integer :: i

    if (t .ge. this%T_begin) then
       do i = 1, this%n
          call this%quant_list(i)%quantp%update(k)
       end do
    end if
    
  end subroutine stats_eval

end module stats
