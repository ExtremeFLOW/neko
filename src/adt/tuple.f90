! Copyright (c) 2020-2022, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
!
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

  !> Integer based 3-tuple 
  type, extends(tuple_t), public :: tuple3_i4_t
     integer :: x(3) = (/0, 0, 0/)
   contains
     procedure, pass(this) :: assign_tuple => tuple3_i4_assign_tuple
     procedure, pass(this) :: assign_vector => tuple3_i4_assign_vector
     procedure, pass(this) :: equal => tuple3_i4_equal
  end type tuple3_i4_t

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

  !> Mixed integer (\f$ x \f$) double precision (\f$ y \f$) 2-tuple \f$(x, y)\f$
  type, extends(tuple_t), public :: tuple_i4r8_t
     integer :: x
     real(kind=dp) :: y
   contains
     procedure, pass(this) :: assign_tuple => tuple_i4r8_assign_tuple
     procedure, pass(this) :: assign_vector => tuple_i4r8_assign_vector
     procedure, pass(this) :: equal => tuple_i4r8_equal
  end type tuple_i4r8_t

  !> Mixed integer (\f$ x, y \f$) double precision (\f$ z \f$) 3-tuple 
  type, extends(tuple_t), public :: tuple_2i4r8_t
     integer :: x, y     
     real(kind=dp) :: z
   contains
     procedure, pass(this) :: assign_tuple => tuple_2i4r8_assign_tuple
     procedure, pass(this) :: assign_vector => tuple_2i4r8_assign_vector
     procedure, pass(this) :: equal => tuple_2i4r8_equal
  end type tuple_2i4r8_t

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
    implicit none
    class(tuple_i4_t), intent(inout) :: this
    class(tuple_t), intent(in) :: other

    select type(other)
    type is (tuple_i4_t)
       this%x = other%x
    end select
  end subroutine tuple_i4_assign_tuple

  !> Assign an integer vector to a tuple
  subroutine tuple_i4_assign_vector(this, x)
    implicit none
    class(tuple_i4_t), intent(inout) :: this
    class(*), dimension(:), intent(in) :: x

    select type(x)
    type is (integer)
       this%x = x
    end select    
  end subroutine tuple_i4_assign_vector

  !> Check if two integer based tuples are equal
  pure function tuple_i4_equal(this, other) result(res)
    implicit none
    class(tuple_i4_t), intent(in) :: this
    class(tuple_t), intent(in) :: other
    logical :: res

    res = .false.
    select type(other)
    type is(tuple_i4_t)
       res = all(this%x .eq. other%x)
    end select
  end function tuple_i4_equal

  !> Assign an integer 3-tuple to a tuple
  subroutine tuple3_i4_assign_tuple(this, other)
    implicit none
    class(tuple3_i4_t), intent(inout) :: this
    class(tuple_t), intent(in) :: other

    select type(other)
    type is(tuple3_i4_t)
       this%x = other%x
    end select
  end subroutine tuple3_i4_assign_tuple

  !> Assign an integer vector to a tuple
  subroutine tuple3_i4_assign_vector(this, x)
    implicit none
    class(tuple3_i4_t), intent(inout) :: this
    class(*), dimension(:), intent(in) :: x

    select type(x)
    type is (integer)
       this%x = x
    end select    
  end subroutine tuple3_i4_assign_vector

  !> Check if two integer based tuples are equal
  pure function tuple3_i4_equal(this, other) result(res)
    implicit none
    class(tuple3_i4_t), intent(in) :: this
    class(tuple_t), intent(in) :: other
    logical :: res    

    res = .false.
    select type(other)
    type is(tuple3_i4_t)
       res = all(this%x .eq. other%x)
    end select
  end function tuple3_i4_equal
  
  !> Assign an integer 4-tuple to a tuple
  subroutine tuple4_i4_assign_tuple(this, other)
    implicit none
    class(tuple4_i4_t), intent(inout) :: this
    class(tuple_t), intent(in) :: other

    select type(other)
    type is(tuple4_i4_t)
       this%x = other%x
    end select
  end subroutine tuple4_i4_assign_tuple

  !> Assign an integer vector to a tuple
  subroutine tuple4_i4_assign_vector(this, x)
    implicit none
    class(tuple4_i4_t), intent(inout) :: this
    class(*), dimension(:), intent(in) :: x

    select type(x)
    type is (integer)
       this%x = x
    end select    
  end subroutine tuple4_i4_assign_vector

  !> Check if two integer based tuples are equal
  pure function tuple4_i4_equal(this, other) result(res)
    implicit none
    class(tuple4_i4_t), intent(in) :: this
    class(tuple_t), intent(in) :: other
    logical :: res    

    res = .false.
    select type(other)
    type is(tuple4_i4_t)
       res = all(this%x .eq. other%x)
    end select
  end function tuple4_i4_equal

  !> Assign a double precision 2-tuple to a tuple
  subroutine tuple_r8_assign_tuple(this, other)
    implicit none
    class(tuple_r8_t), intent(inout) :: this
    class(tuple_t), intent(in) :: other

    select type(other)
    type is(tuple_r8_t)       
       this%x = other%x
    end select
  end subroutine tuple_r8_assign_tuple

  !> Assign a double precision vector to a tuple
  subroutine tuple_r8_assign_vector(this, x)
    implicit none
    class(tuple_r8_t), intent(inout) :: this
    class(*), dimension(:), intent(in) :: x

    select type(x)
    type is (double precision)
       this%x = x
    end select    
  end subroutine tuple_r8_assign_vector

  !> Check if two double precision tuples are equal
  pure function tuple_r8_equal(this, other) result(res)
    implicit none
    class(tuple_r8_t), intent(in) :: this
    class(tuple_t), intent(in) :: other
    logical :: res

    res = .false.
    select type(other)
    type is(tuple_r8_t)
       if (abscmp(this%x(1), other%x(1)) .and. &
            abscmp(this%x(2), other%x(2))) then
          res = .true.          
       end if
    end select
  end function tuple_r8_equal

  !> Assign a mixed integer-double precision 2-tuple to a tuple
  subroutine tuple_i4r8_assign_tuple(this, other)
    implicit none
    class(tuple_i4r8_t), intent(inout) :: this
    class(tuple_t), intent(in) :: other

    select type(other)
    type is(tuple_i4r8_t)       
       this%x = other%x
       this%y = other%y
    end select
  end subroutine tuple_i4r8_assign_tuple

  !> Assign a mixed intreger-double precision vector to a tuple
  subroutine tuple_i4r8_assign_vector(this, x)
    implicit none
    class(tuple_i4r8_t), intent(inout) :: this
    class(*), dimension(:), intent(in) :: x

    select type(x)
    type is (integer)
       this%x = x(1)
       this%y = dble(x(2))
    type is (double precision)
       this%x = int(x(1))
       this%y = x(2)
    end select

  end subroutine tuple_i4r8_assign_vector

  !> Check if two mixed integer-double precision tuples are equal
  pure function tuple_i4r8_equal(this, other) result(res)
    implicit none
    class(tuple_i4r8_t), intent(in) :: this
    class(tuple_t), intent(in) :: other
    logical :: res

    res = .false.
    select type(other)
    type is(tuple_i4r8_t)
       if ((this%x .eq. other%x) .and. &
            abscmp(this%y, other%y)) then
          res = .true.          
       end if
    end select
  end function tuple_i4r8_equal

  !> Assign a mixed integer-double precision 3-tuple to a tuple
  subroutine tuple_2i4r8_assign_tuple(this, other)
    implicit none
    class(tuple_2i4r8_t), intent(inout) :: this
    class(tuple_t), intent(in) :: other

    select type(other)
    type is(tuple_2i4r8_t)       
       this%x = other%x
       this%y = other%y
       this%z = other%z
    end select
  end subroutine tuple_2i4r8_assign_tuple

  !> Assign a mixed intreger-double precision vector to a tuple
  subroutine tuple_2i4r8_assign_vector(this, x)
    implicit none
    class(tuple_2i4r8_t), intent(inout) :: this
    class(*), dimension(:), intent(in) :: x

    select type(x)
    type is (integer)
       this%x = x(1)
       this%y = x(2)
       this%z = dble(x(3))
    type is (double precision)
       this%x = int(x(1))
       this%y = int(x(2))
       this%z = x(3)
    end select

  end subroutine tuple_2i4r8_assign_vector

  !> Check if two mixed integer-double precision tuples are equal
  pure function tuple_2i4r8_equal(this, other) result(res)
    implicit none
    class(tuple_2i4r8_t), intent(in) :: this
    class(tuple_t), intent(in) :: other
    logical :: res

    res = .false.
    select type(other)
    type is(tuple_2i4r8_t)
       if ((this%x .eq. other%x) .and. &
            (this%y .eq. other%y) .and. &
            abscmp(this%z, other%z)) then
          res = .true.          
       end if
    end select
  end function tuple_2i4r8_equal
   
end module tuple
