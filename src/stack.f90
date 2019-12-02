!> Implements a dynamic stack ADT
!! @details a stack storing values @a data of a arbitrary type
module stack
  use num_types
  use utils
  use math, only : NEKO_M_LN2
  implicit none
  private

  !> Base type for a stack
  type, private :: stack_t
     class(*), allocatable :: data(:)
     integer :: top_
     integer :: size_
   contains
     procedure, pass(this) :: stack_init
     procedure, public, pass(this) :: free => stack_free
     procedure, public, pass(this) :: clear => stack_clear
     procedure, public, pass(this) :: size => stack_size
  end type stack_t

  !> Integer based stack
  type, public, extends(stack_t) :: stack_i4_t
   contains
     procedure, pass(this) :: init => stack_i4_init
  end type stack_i4_t

  !> Double precision based stack
  type, public, extends(stack_t) :: stack_r8_t
   contains
     procedure, pass(this) :: init => stack_r8_init
  end type stack_r8_t
  
contains

  !> Initialize a stack of type @a data
  subroutine stack_init(this, size, data)
    class(stack_t), intent(inout) :: this 
    integer, value :: size       !< Initial size of the stack
    class(*), intent(in) :: data !< Type of data


    if (size .lt. 4) then
       size =4
    end if

    this%size_ = ishft(4, ceiling(log(dble(size)) / NEKO_M_LN2))
    this%top_ = 0

    allocate(this%data(this%size_), source=data)

  end subroutine stack_init
  
  !> Destroy a stack
  subroutine stack_free(this)
    class(stack_t), intent(inout) :: this
    
    if (allocated(this%data)) then
       deallocate(this%data)
       this%size_ = 0 
       this%top_ = 0
    end if    

  end subroutine stack_free

  !> Clear all entries of a stack
  subroutine stack_clear(this)
    class(stack_t), intent(inout) :: this
    this%top_ = 0
  end subroutine stack_clear

  !> Return number of entries in the stack
  pure function stack_size(this) result(size)
    class(stack_t), intent(in) :: this
    integer :: size
    size = this%size_
  end function stack_size

  !> Initialize an integer based stack
  subroutine stack_i4_init(this, size)
    class(stack_i4_t), intent(inout) :: this
    integer, value, optional :: size
    integer :: data

    if (present(size)) then
       call stack_init(this, size, data)
    else
       call stack_init(this, 32, data)
    end if

  end subroutine stack_i4_init

  !> Initialize a double precision based stack
  subroutine stack_r8_init(this, size)
    class(stack_r8_t), intent(inout) :: this
    integer, value, optional :: size
    real(kind=dp) :: data

    if (present(size)) then
       call stack_init(this, size, data)
    else
       call stack_init(this, 32, data)
    end if

  end subroutine stack_r8_init

  
end module stack
