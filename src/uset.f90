!> Implements an unordered set ADT
!! @details A unordered set storing a fixed data-type @a data 
module uset
  use num_types
  use htable
  implicit none
  private

  !> Base type for an unordered set
  type, private, abstract :: uset_t
   contains
     procedure(uset_init), pass(this), deferred :: init
     procedure(uset_free), pass(this), deferred :: free
     procedure(uset_size), pass(this), deferred :: size
     procedure(uset_clear), pass(this), deferred :: clear
     procedure(uset_element), pass(this), deferred :: element
     procedure(uset_add), pass(this), deferred :: add
  end type uset_t

  !> Integer based unordered set
  type, extends(uset_t), public :: uset_i4_t
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
  type, extends(uset_t), public :: uset_r8_t
     type(htable_r8_t) :: t
   contains
     procedure, pass(this) :: init => uset_r8_init
     procedure, pass(this) :: free => uset_r8_free
     procedure, pass(this) :: size => uset_r8_size
     procedure, pass(this) :: clear => uset_r8_clear
     procedure, pass(this) :: element => uset_r8_element
     procedure, pass(this) :: add => uset_r8_add
  end type uset_r8_t

  !> Interface for initializing an unordered set
  abstract interface
     subroutine uset_init(this, n) 
       import uset_t
       class(uset_t), intent(inout) :: this
       integer, value, optional :: n
     end subroutine uset_init
  end interface

  !> Interface for destroying an unordered set
  abstract interface
     subroutine uset_free(this) 
       import uset_t
       class(uset_t), intent(inout) :: this
     end subroutine uset_free
  end interface

  !> Interface for getting the cardinality of an unordered set
  abstract interface
     pure function uset_size(this) result(entries)
       import uset_t
       class(uset_t), intent(in) :: this
       integer :: entries
     end function uset_size
  end interface

  !> Interface for clearing an unordered set
  abstract interface
     subroutine uset_clear(this) 
       import uset_t
       class(uset_t), intent(inout) :: this
     end subroutine uset_clear
  end interface

  !> Interface for checking if @a key is an element of an unordered set
  abstract interface
     function uset_element(this, key) result(res)
       import uset_t
       class(uset_t), intent(inout) :: this
       class(*), intent(inout) :: key
       logical :: res
     end function uset_element
  end interface

  !> Inteface for adding @a key to an unorderd set
  abstract interface     
     subroutine uset_add(this, key)
       import uset_t
       class(uset_t), intent(inout) :: this
       class(*), intent(inout) :: key
     end subroutine uset_add
 end interface

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
  pure function uset_i4_size(this) result(entries)
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
    class(*), intent(inout) :: key
    integer :: data
    logical :: res

    select type(key)
    type is (integer)
       res = (this%t%get(key, data) .eq. 0)
    end select    
  end function uset_i4_element
  
  !> Add an integer @a key to the set
  subroutine uset_i4_add(this, key)
    class(uset_i4_t), intent(inout) :: this
    class(*), intent(inout) :: key
    integer :: data
    data = 1

    select type(key)
    type is (integer)
       call this%t%set(key, data)
    end select
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
  pure function uset_r8_size(this)  result(entries)
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
    class(*), intent(inout) :: key
    logical :: res
    integer :: data


    select type(key)
    type is (double precision)
       res = (this%t%get(key, data) .eq. 0)
    end select
    
  end function uset_r8_element

    !> Add a double precision @a key to the set
  subroutine uset_r8_add(this, key)
    class(uset_r8_t), intent(inout) :: this
    class(*), intent(inout) :: key
    integer :: data
    data = 1

    select type(key)
    type is (double precision)
       call this%t%set(key, data)
    end select
  end subroutine uset_r8_add


end module uset
