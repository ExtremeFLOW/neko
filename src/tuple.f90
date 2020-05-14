!> Implements a tuple
module tuple
  use math  
  use num_types
  implicit none
  private

  !> Base type for an n-tuple 
  type, private, abstract :: tuple_t
   contains
     procedure(tuple_assign), pass(this), deferred :: assign !< Assignment intf.
     procedure(tuple_equal), pass(this) , deferred :: equal  !< Equal operator
  end type tuple_t

  !> Integer based tuple 
  type, extends(tuple_t), public :: tuple_i4_t
     integer :: i, j
   contains
     procedure, pass(this) :: assign => tuple_i4_assign
     procedure, pass(this) :: equal => tuple_i4_equal
     generic :: operator(.eq.) => equal
     generic :: assignment(=) => assign
  end type tuple_i4_t

  !> Double precision based tuple 
  type, extends(tuple_t), public :: tuple_r8_t
     real(kind=dp) :: i, j
   contains
     procedure, pass(this) :: assign => tuple_r8_assign
     procedure, pass(this) :: equal => tuple_r8_equal
     generic :: operator(.eq.) => equal
     generic :: assignment(=) => assign
  end type tuple_r8_t

  abstract interface
     subroutine tuple_assign(this, other)
       import :: tuple_t
       class(tuple_t), intent(inout) :: this
       class(tuple_t), intent(in) :: other
     end subroutine tuple_assign
  end interface
  
  abstract interface
     pure function tuple_equal(this, other) result(res)
       import :: tuple_t
       class(tuple_t), intent(in) :: this
       class(tuple_t), intent(in) :: other
       logical :: res
     end function tuple_equal
  end interface

contains

  subroutine tuple_i4_assign(this, other)
    class(tuple_i4_t), intent(inout) :: this
    type(tuple_i4_t), intent(in) :: other

    this%i = other%i
    this%j = other%j
    
  end subroutine tuple_i4_assign

  !> Check if two integer based tuples are equal
  pure function tuple_i4_equal(this, other) result(res)
    class(tuple_i4_t), intent(in) :: this
    type(tuple_i4_t), intent(in) :: other
    logical :: res

    if ((this%i .eq. other%i) .and. &
         (this%j .eq. other%j)) then
       res = .true.
    else
       res = .false.
    end if
    
  end function tuple_i4_equal

  subroutine tuple_r8_assign(this, other)
    class(tuple_r8_t), intent(inout) :: this
    type(tuple_r8_t), intent(in) :: other

    this%i = other%i
    this%j = other%j
    
  end subroutine tuple_r8_assign

  !> Check if two double precision tuples are equal
  pure function tuple_r8_equal(this, other) result(res)
    class(tuple_r8_t), intent(in) :: this
    type(tuple_r8_t), intent(in) :: other
    logical :: res

    if (abscmp(this%i, other%i) .and. &
         abscmp(this%j, other%j)) then
       res = .true.
    else
       res = .false.
    end if
    
  end function tuple_r8_equal

  
  
end module tuple
