! Copyright (c) 2019-2023, The Neko Authors
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
!> Implements an unordered set ADT
!! @details A unordered set storing a fixed data-type @a data
module uset
  use utils, only : neko_error
  use num_types, only : i8, dp
  use htable, only : htable_i4_t, htable_iter_i4_t, &
                     htable_i8_t, htable_iter_i8_t, &
                     htable_r8_t, htable_iter_r8_t
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
     procedure(uset_remove), pass(this), deferred :: remove
  end type uset_t

  !> Integer based unordered set
  type, extends(uset_t), public :: uset_i4_t
     type(htable_i4_t) :: t
     type(htable_iter_i4_t) :: it
   contains
     procedure, pass(this) :: init => uset_i4_init
     procedure, pass(this) :: free => uset_i4_free
     procedure, pass(this) :: size => uset_i4_size
     procedure, pass(this) :: clear => uset_i4_clear
     procedure, pass(this) :: element => uset_i4_element
     procedure, pass(this) :: add => uset_i4_add
     procedure, pass(this) :: remove => uset_i4_remove
     procedure, pass(this) :: iter_init => uset_i4_iter_init
     procedure, pass(this) :: iter_next => uset_i4_iter_next
     procedure, pass(this) :: iter_value => uset_i4_iter_value
  end type uset_i4_t

  !> Integer*8 based unordered set
  type, extends(uset_t), public :: uset_i8_t
     type(htable_i8_t) :: t
     type(htable_iter_i8_t) :: it
   contains
     procedure, pass(this) :: init => uset_i8_init
     procedure, pass(this) :: free => uset_i8_free
     procedure, pass(this) :: size => uset_i8_size
     procedure, pass(this) :: clear => uset_i8_clear
     procedure, pass(this) :: element => uset_i8_element
     procedure, pass(this) :: add => uset_i8_add
     procedure, pass(this) :: remove => uset_i8_remove
     procedure, pass(this) :: iter_init => uset_i8_iter_init
     procedure, pass(this) :: iter_next => uset_i8_iter_next
     procedure, pass(this) :: iter_value => uset_i8_iter_value
  end type uset_i8_t

  !> Double precision unordered set
  type, extends(uset_t), public :: uset_r8_t
     type(htable_r8_t) :: t
     type(htable_iter_r8_t) :: it
   contains
     procedure, pass(this) :: init => uset_r8_init
     procedure, pass(this) :: free => uset_r8_free
     procedure, pass(this) :: size => uset_r8_size
     procedure, pass(this) :: clear => uset_r8_clear
     procedure, pass(this) :: element => uset_r8_element
     procedure, pass(this) :: add => uset_r8_add
     procedure, pass(this) :: remove => uset_r8_remove
     procedure, pass(this) :: iter_init => uset_r8_iter_init
     procedure, pass(this) :: iter_next => uset_r8_iter_next
     procedure, pass(this) :: iter_value => uset_r8_iter_value
  end type uset_r8_t

  !> Interface for initializing an unordered set
  abstract interface
     subroutine uset_init(this, n)
       import uset_t
       class(uset_t), intent(inout) :: this
       integer, optional :: n
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

  !> Inteface for removing @a key in an unorderd set
  abstract interface
     subroutine uset_remove(this, key)
       import uset_t
       class(uset_t), intent(inout) :: this
       class(*), intent(inout) :: key
     end subroutine uset_remove
  end interface

contains

  !> Initialize an empty integer based unordered set
  subroutine uset_i4_init(this, n)
    class(uset_i4_t), intent(inout) :: this
    integer, optional :: n

    if (present(n)) then
       call this%t%init(n)
    else
       call this%t%init(64)
    end if
  end subroutine uset_i4_init

  !> Destroy an integer based unordered set
  subroutine uset_i4_free(this)
    class(uset_i4_t), intent(inout) :: this

    nullify(this%it%t)
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
    class default
       res = .false.
    end select
  end function uset_i4_element

  !> Add an integer @a key to the set
  subroutine uset_i4_add(this, key)
    class(uset_i4_t), intent(inout) :: this
    class(*), intent(inout) :: key
    integer :: data

    select type(key)
    type is (integer)
       data = key
       call this%t%set(key, data)
    class default
       call neko_error("Invalid key")
    end select
  end subroutine uset_i4_add

  !> Remove an integer @a key from the set
  subroutine uset_i4_remove(this, key)
    class(uset_i4_t), intent(inout) :: this
    class(*), intent(inout) :: key

    select type(key)
    type is (integer)
       call this%t%remove(key)
    class default
       call neko_error("Invalid key")
    end select
  end subroutine uset_i4_remove

  !> Initialise an integer based set iterator
  subroutine uset_i4_iter_init(this)
    class(uset_i4_t), target, intent(inout) :: this
    call this%it%init(this%t)
  end subroutine uset_i4_iter_init

  !> Advance an integer based set iterator
  function uset_i4_iter_next(this) result(valid)
    class(uset_i4_t), intent(inout) :: this
    logical :: valid
    valid = this%it%next()
  end function uset_i4_iter_next

  !> Return the current value of an integer based set iterator
  function uset_i4_iter_value(this) result(value)
    class(uset_i4_t), target, intent(inout) :: this
    integer, pointer :: value
    value => this%it%value()
  end function uset_i4_iter_value

  !> Initialize an empty integer*8 based unordered set
  subroutine uset_i8_init(this, n)
    class(uset_i8_t), intent(inout) :: this
    integer, optional :: n

    if (present(n)) then
       call this%t%init(n)
    else
       call this%t%init(64)
    end if
  end subroutine uset_i8_init

  !> Destroy an integer*8 based unordered set
  subroutine uset_i8_free(this)
    class(uset_i8_t), intent(inout) :: this

    nullify(this%it%t)
    call this%t%free()

  end subroutine uset_i8_free

  !> Return the cardinality of an integer*8 based unordered set
  pure function uset_i8_size(this) result(entries)
    class(uset_i8_t), intent(in) :: this
    integer :: entries

    entries = this%t%num_entries()

  end function uset_i8_size

  !> Clear an integer*8 based unordered set
  subroutine uset_i8_clear(this)
    class(uset_i8_t), intent(inout) :: this

    call this%t%clear()
  end subroutine uset_i8_clear

  !> Check if an integer*8 @a key is an element of the set
  function uset_i8_element(this, key) result(res)
    class(uset_i8_t), intent(inout) :: this
    class(*), intent(inout) :: key
    integer(kind=i8) :: data
    logical :: res

    select type(key)
    type is (integer(i8))
       res = (this%t%get(key, data) .eq. 0)
    class default
       res = .false.
    end select
  end function uset_i8_element

  !> Add an integer*8 @a key to the set
  subroutine uset_i8_add(this, key)
    class(uset_i8_t), intent(inout) :: this
    class(*), intent(inout) :: key
    integer(kind=i8) :: data

    select type(key)
    type is (integer(i8))
       data = key
       call this%t%set(key, data)
    class default
       call neko_error("Invalid key")
    end select
  end subroutine uset_i8_add

  !> Remove an integer*8 @a key from the set
  subroutine uset_i8_remove(this, key)
    class(uset_i8_t), intent(inout) :: this
    class(*), intent(inout) :: key

    select type(key)
    type is (integer(i8))
       call this%t%remove(key)
    class default
       call neko_error("Invalid key")
    end select
  end subroutine uset_i8_remove

  !> Initialise an integer based set iterator*8
  subroutine uset_i8_iter_init(this)
    class(uset_i8_t), target, intent(inout) :: this
    call this%it%init(this%t)
  end subroutine uset_i8_iter_init

  !> Advance an integer*8 based set iterator
  function uset_i8_iter_next(this) result(valid)
    class(uset_i8_t), intent(inout) :: this
    logical :: valid
    valid = this%it%next()
  end function uset_i8_iter_next

  !> Return the current value of an integer*8 based set iterator
  function uset_i8_iter_value(this) result(value)
    class(uset_i8_t), target, intent(inout) :: this
    integer(kind=i8), pointer :: value

    ! We should not need this extra select block, and it works great
    ! without it for GNU, Intel and NEC, but breaks horribly on Cray
    ! (>11.0.x) when using high opt. levels.
    select type(hp => this)
    type is (uset_i8_t)
       value => hp%it%value()
    class default
       call neko_error('Invalid uset htable iter (i8)')
    end select
  end function uset_i8_iter_value

  !> Initialize an empty double precision based unordered set
  subroutine uset_r8_init(this, n)
    class(uset_r8_t), intent(inout) :: this
    integer, optional :: n

    if (present(n)) then
       call this%t%init(n)
    else
       call this%t%init(64)
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
    real(kind=dp) :: data

    select type(key)
    type is (double precision)
       res = (this%t%get(key, data) .eq. 0)
    class default
       res = .false.
    end select

  end function uset_r8_element

  !> Add a double precision @a key to the set
  subroutine uset_r8_add(this, key)
    class(uset_r8_t), intent(inout) :: this
    class(*), intent(inout) :: key
    real(kind=dp) :: data

    select type(key)
    type is (double precision)
       data = key
       call this%t%set(key, data)
    class default
       call neko_error("Invalid key")
    end select
  end subroutine uset_r8_add

  !> Remove a double precision @a key from the set
  subroutine uset_r8_remove(this, key)
    class(uset_r8_t), intent(inout) :: this
    class(*), intent(inout) :: key

    select type(key)
    type is (double precision)
       call this%t%remove(key)
    class default
       call neko_error("Invalid key")
    end select
  end subroutine uset_r8_remove

  !> Initialise a double precision based set iterator
  subroutine uset_r8_iter_init(this)
    class(uset_r8_t), target, intent(inout) :: this
    call this%it%init(this%t)
  end subroutine uset_r8_iter_init

  !> Advance a double precision based set iterator
  function uset_r8_iter_next(this) result(valid)
    class(uset_r8_t), intent(inout) :: this
    logical :: valid
    valid = this%it%next()
  end function uset_r8_iter_next

  !> Return the current value of a double precision based set iterator
  function uset_r8_iter_value(this) result(value)
    class(uset_r8_t), target, intent(inout) :: this
    real(kind=dp), pointer :: value
    value => this%it%value()
  end function uset_r8_iter_value


end module uset
