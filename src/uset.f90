!> Implements an unordered set ADT
!! @details A unordered set storing a fixed data-type @a data 
module uset
  use num_types
  use htable
  implicit none
  private

  !> Integer based unordered set
  type, public :: uset_i4_t
     type(htable_i4_t) :: t
   contains
     procedure, pass(this) :: init => uset_i4_init
     procedure, pass(this) :: free => uset_i4_free
     procedure, pass(this) :: size => uset_i4_size
     procedure, pass(this) :: clear => uset_i4_clear
     procedure, pass(this) :: element => uset_i4_element
     procedure, pass(this) :: add => uset_i4_add
  end type uset_i4_t

  !> Double precision unordered set
  type, public :: uset_r8_t
     type(htable_r8_t) :: t
   contains
     procedure, pass(this) :: init => uset_r8_init
     procedure, pass(this) :: free => uset_r8_free
     procedure, pass(this) :: size => uset_r8_size
     procedure, pass(this) :: clear => uset_r8_clear
     procedure, pass(this) :: element => uset_r8_element
     procedure, pass(this) :: add => uset_r8_add
  end type uset_r8_t

contains

  !> Initialize an empty integer based unordered set
  subroutine uset_i4_init(this, n)
    class(uset_i4_t), intent(inout) :: this
    integer, value, optional :: n
    integer :: key

    if (present(n)) then
       call this%t%init(n)
    else
       call this%t%init(64)
    end if    
  end subroutine uset_i4_init
  
  !> Destroy an integer based unordered set
  subroutine uset_i4_free(this)
    class(uset_i4_t), intent(inout) :: this

    call this%t%free()
    
  end subroutine uset_i4_free

  !> Return the cardinality of an integer based unordered set
  function uset_i4_size(this) result(entries)
    class(uset_i4_t), intent(in) :: this
    integer :: entries

    entries = this%t%num_entries()
    
  end function uset_i4_size

  !> Clear an integer based unordered set
  subroutine uset_i4_clear(this)
    class(uset_i4_t), intent(inout) :: this

    call this%t%clear()
  end subroutine uset_i4_clear

  !> Check if an integer @a key is an element of the set
  function uset_i4_element(this, key) result(res)
    class(uset_i4_t), intent(inout) :: this
    integer, intent(inout) :: key
    integer :: data
    logical :: res

    res = (this%t%get(key, data) .eq. 0)
    
  end function uset_i4_element
  
  !> Add an integer @a key to the set
  subroutine uset_i4_add(this, key)
    class(uset_i4_t), intent(inout) :: this
    integer :: key

    call this%t%set(key, 1)
  end subroutine uset_i4_add

  !> Initialize an empty double precision based unordered set
  subroutine uset_r8_init(this, n)
    class(uset_r8_t), intent(inout) :: this
    integer, value, optional :: n
    integer :: data
    if (present(n)) then
       call this%t%init(n, data)
    else
       call this%t%init(64, data)
    end if
  end subroutine uset_r8_init
  
  !> Destroy a double precision based unordered set
  subroutine uset_r8_free(this)
    class(uset_r8_t), intent(inout) :: this
    
    call this%t%free()
    
  end subroutine uset_r8_free

  !> Return the cardinality of a double precision based unordered set
  function uset_r8_size(this)  result(entries)
    class(uset_r8_t), intent(in) :: this
    integer :: entries

    entries = this%t%num_entries()

  end function uset_r8_size

  !> Clear a double precision based unordered set
  subroutine uset_r8_clear(this)
    class(uset_r8_t), intent(inout) :: this

    call this%t%clear()
  end subroutine uset_r8_clear
  
  !> Check if a double precision @a key is an element of the set
  function uset_r8_element(this, key) result(res)
    class(uset_r8_t), intent(inout) :: this
    real(kind=dp), intent(inout) :: key
    integer :: data
    logical :: res

    res = (this%t%get(key, data) .eq. 0)
    
  end function uset_r8_element

    !> Add a double precision @a key to the set
  subroutine uset_r8_add(this, key)
    class(uset_r8_t), intent(inout) :: this
    real(kind=dp) :: key

    call this%t%set(key, 1)
  end subroutine uset_r8_add


end module uset
