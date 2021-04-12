!> Implements a n-tuple
module tuple
  use math  
  use num_types
  implicit none
  private

  !> Base type for an n-tuple 
  type, public, abstract :: tuple_t
   contains
     procedure(tuple_assign_tuple), pass(this), deferred :: assign_tuple
     procedure(tuple_assign_vector), pass(this), deferred :: assign_vector
     procedure(tuple_equal), pass(this) , deferred :: equal
     generic :: operator(.eq.) => equal
     generic :: assignment(=) => assign_tuple, assign_vector
  end type tuple_t

  !> Integer based 2-tuple 
  type, extends(tuple_t), public :: tuple_i4_t
     integer :: x(2) = (/0, 0/)
   contains
     procedure, pass(this) :: assign_tuple => tuple_i4_assign_tuple
     procedure, pass(this) :: assign_vector => tuple_i4_assign_vector
     procedure, pass(this) :: equal => tuple_i4_equal
  end type tuple_i4_t

    !> Integer based 4-tuple 
  type, extends(tuple_t), public :: tuple4_i4_t
     integer :: x(4) = (/0, 0, 0, 0/)
   contains
     procedure, pass(this) :: assign_tuple => tuple4_i4_assign_tuple
     procedure, pass(this) :: assign_vector => tuple4_i4_assign_vector
     procedure, pass(this) :: equal => tuple4_i4_equal
  end type tuple4_i4_t

  !> Double precision based 2-tuple 
  type, extends(tuple_t), public :: tuple_r8_t
     real(kind=dp) :: x(2) = (/0d0, 0d0/)
   contains
     procedure, pass(this) :: assign_tuple => tuple_r8_assign_tuple
     procedure, pass(this) :: assign_vector => tuple_r8_assign_vector
     procedure, pass(this) :: equal => tuple_r8_equal
  end type tuple_r8_t

  !> Abstract intf. for assigning a tuple to a tuple
  abstract interface
     subroutine tuple_assign_tuple(this, other)
       import :: tuple_t
       class(tuple_t), intent(inout) :: this
       class(tuple_t), intent(in) :: other
     end subroutine tuple_assign_tuple
  end interface

  !> Abstract intf. for assigning a vector to a n-tuple
  abstract interface
     subroutine tuple_assign_vector(this, x)
       import :: tuple_t
       class(tuple_t), intent(inout) :: this
       class(*), dimension(:), intent(in) :: x
     end subroutine tuple_assign_vector
  end interface
  
  !> Abstract intf. for tuple comparison
  abstract interface
     pure function tuple_equal(this, other) result(res)
       import :: tuple_t
       class(tuple_t), intent(in) :: this
       class(tuple_t), intent(in) :: other
       logical :: res
     end function tuple_equal
  end interface

contains
  
  !> Assign an integer 2-tuple to a tuple
  subroutine tuple_i4_assign_tuple(this, other)
    class(tuple_i4_t), intent(inout) :: this
    class(tuple_t), intent(in) :: other

    select type(other)
    type is (tuple_i4_t)
       this%x = other%x
    end select
  end subroutine tuple_i4_assign_tuple

  !> Assign an integer vector to a tuple
  subroutine tuple_i4_assign_vector(this, x)
    class(tuple_i4_t), intent(inout) :: this
    class(*), dimension(:), intent(in) :: x

    select type(x)
    type is (integer)
       this%x = x
    end select    
  end subroutine tuple_i4_assign_vector

  !> Check if two integer based tuples are equal
  pure function tuple_i4_equal(this, other) result(res)
    class(tuple_i4_t), intent(in) :: this
    class(tuple_t), intent(in) :: other
    logical :: res

    res = .false.
    select type(other)
    type is(tuple_i4_t) 
       if ((this%x(1) .eq. other%x(1)) .and. &
            (this%x(2) .eq. other%x(2))) then
          res = .true.
       end if
    end select
  end function tuple_i4_equal

  !> Assign an integer 4-tuple to a tuple
  subroutine tuple4_i4_assign_tuple(this, other)
    class(tuple4_i4_t), intent(inout) :: this
    class(tuple_t), intent(in) :: other

    select type(other)
    type is(tuple4_i4_t)
       this%x = other%x
    end select
  end subroutine tuple4_i4_assign_tuple

  !> Assign an integer vector to a tuple
  subroutine tuple4_i4_assign_vector(this, x)
    class(tuple4_i4_t), intent(inout) :: this
    class(*), dimension(:), intent(in) :: x

    select type(x)
    type is (integer)
       this%x = x
    end select    
  end subroutine tuple4_i4_assign_vector

  !> Check if two integer based tuples are equal
  pure function tuple4_i4_equal(this, other) result(res)
    class(tuple4_i4_t), intent(in) :: this
    class(tuple_t), intent(in) :: other
    logical :: res    

    res = .false.
    select type(other)
    type is(tuple4_i4_t)
       if ((this%x(1) .eq. other%x(1)) .and. &
            (this%x(2) .eq. other%x(2)) .and. &
            (this%x(3) .eq. other%x(3)) .and. &
            (this%x(4) .eq. other%x(4))) then
          res = .true.
       end if
    end select
  end function tuple4_i4_equal

  !> Assign a double precision 2-tuple to a tuple
  subroutine tuple_r8_assign_tuple(this, other)
    class(tuple_r8_t), intent(inout) :: this
    class(tuple_t), intent(in) :: other

    select type(other)
    type is(tuple_r8_t)       
       this%x = other%x
    end select
  end subroutine tuple_r8_assign_tuple

  !> Assign a double precision vector to a tuple
  subroutine tuple_r8_assign_vector(this, x)
    class(tuple_r8_t), intent(inout) :: this
    class(*), dimension(:), intent(in) :: x

    select type(x)
    type is (double precision)
       this%x = x
    end select    
  end subroutine tuple_r8_assign_vector

  !> Check if two double precision tuples are equal
  pure function tuple_r8_equal(this, other) result(res)
    class(tuple_r8_t), intent(in) :: this
    class(tuple_t), intent(in) :: other
    logical :: res

    res = .false.
    select type(other)
    type is(tuple_r8_t)
       if (dabscmp(this%x(1), other%x(1)) .and. &
            dabscmp(this%x(2), other%x(2))) then
          res = .true.          
       end if
    end select
  end function tuple_r8_equal
   
end module tuple
