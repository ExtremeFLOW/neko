!> Stores a series fields
module field_series
  use num_types
  use field
  implicit none
  private
  
  type, public :: field_series_t
     type(field_t), pointer :: f => null()
     type(field_t), allocatable :: lf(:) 
     integer, private :: len = 0
   contains
     procedure, pass(this) :: init => field_series_init
     procedure, pass(this) :: free => field_series_free
     procedure, pass(this) :: update => field_series_update
     procedure, pass(this) :: set => field_series_set
     procedure, pass(this) :: size => field_series_size
  end type field_series_t

contains

  !> Initialize a field series of length @a len for a field @a f
  subroutine field_series_init(this, f, len)
    class(field_series_t), intent(inout) :: this
    type(field_t), intent(inout), target :: f
    integer :: len
    character(len=80) :: name
    integer :: i

    call this%free()

    this%f => f
    this%len = len

    allocate(this%lf(len))

    do i = 1, this%len
       name = trim(f%name)//'_lag'//char(i)
       call field_init(this%lf(i), f%dof, name)
    end do

  end subroutine field_series_init

  !> Deallocates a field series
  subroutine field_series_free(this)
    class(field_series_t), intent(inout) :: this
    integer :: i

    if (associated(this%f)) then
       nullify(this%f)
    end if

    do i = 1, this%len
       call field_free(this%lf(i))
    end do
    
  end subroutine field_series_free

  !> Return the size of the field series
  function field_series_size(this) result(len)
    class(field_series_t), intent(in) :: this
    integer :: len
    len = this%len
  end function field_series_size

  !> Update a field series (evict oldest entry)
  subroutine field_series_update(this)
    class(field_series_t), intent(inout) :: this
    integer :: i

    do i = this%len, 2, -1
       this%lf(i) = this%lf(i-1)
    end do

    this%lf(1) = this%f
    
  end subroutine field_series_update

  !> Set all fields in a series to @a g
  subroutine field_series_set(this, g)
    class(field_series_t), intent(inout) :: this
    type(field_t), intent(in) :: g
    integer :: i

    do i = 1, this%len
       this%lf(i) = g
    end do
    
  end subroutine field_series_set
  
end module field_series
