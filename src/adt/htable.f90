! Copyright (c) 2019-2025, The Neko Authors
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
!> Implements a hash table ADT
!! @details A hash table storing tuples @a (key, data), with fixed
!! data-type @a key but with arbitrary typed @a data
module htable
  use num_types
  use utils
  use point, only : point_t
  use tuple, only : tuple_i4_t, tuple4_i4_t, tuple_t
  use math, only : NEKO_M_LN2
  use, intrinsic :: iso_c_binding, only : c_ptr, c_associated
  implicit none
  private

  !> Base type for a hash table
  type, public, abstract :: htable_t
     integer, private :: size
     integer, private :: entries
     class(*), private, allocatable :: data(:)
     logical, private, allocatable :: valid(:)
     logical, private, allocatable :: skip(:)
   contains
     procedure(htable_hash), pass(this), deferred :: hash
     procedure(htable_destructor), pass(this), deferred :: free
     procedure, public, pass(this) :: clear => htable_clear
     procedure, public, pass(this) :: num_entries => htable_num_entries
     procedure, public, pass(this) :: get_size => htable_size
     procedure, private, pass(this) :: free_base => htable_free
     procedure, private, nopass :: set_data => htable_set_data
     procedure, private, nopass :: get_data => htable_get_data
  end type htable_t

  abstract interface
     subroutine htable_destructor(this)
       import :: htable_t
       class(htable_t), intent(inout) :: this
     end subroutine htable_destructor
  end interface

  abstract interface
     pure function htable_hash(this, k, c) result(hash)
       import htable_t
       class(htable_t), intent(in) :: this
       class(*), intent(in) :: k
       integer, value :: c
       integer :: hash
     end function htable_hash
  end interface

  !
  ! Type wrappers
  !
  type, public :: h_cptr_t
     type(c_ptr) :: ptr
  end type h_cptr_t


  !
  ! Implementations
  !

  !> Integer based hash table
  type, public, extends(htable_t) :: htable_i4_t
     integer, private, allocatable :: key(:)
   contains
     procedure, pass(this) :: init => htable_i4_init
     procedure, pass(this) :: free => htable_i4_free
     procedure, pass(this) :: set => htable_i4_set
     procedure, pass(this) :: get => htable_i4_get
     procedure, pass(this) :: hash => htable_i4_hash
     procedure, pass(this) :: remove => htable_i4_remove
  end type htable_i4_t

  !> Integer*8 based hash table
  type, public, extends(htable_t) :: htable_i8_t
     integer(kind=i8), private, allocatable :: key(:)
   contains
     procedure, pass(this) :: init => htable_i8_init
     procedure, pass(this) :: free => htable_i8_free
     procedure, pass(this) :: set => htable_i8_set
     procedure, pass(this) :: get => htable_i8_get
     procedure, pass(this) :: hash => htable_i8_hash
     procedure, pass(this) :: remove => htable_i8_remove
  end type htable_i8_t

  !> Double precision based hash table
  type, public, extends(htable_t) :: htable_r8_t
     real(kind=dp), private, allocatable :: key(:)
   contains
     procedure, pass(this) :: init => htable_r8_init
     procedure, pass(this) :: free => htable_r8_free
     procedure, pass(this) :: set => htable_r8_set
     procedure, pass(this) :: get => htable_r8_get
     procedure, pass(this) :: hash => htable_r8_hash
     procedure, pass(this) :: remove => htable_r8_remove
  end type htable_r8_t

  !> Point based hash table
  type, public, extends(htable_t) :: htable_pt_t
     type(point_t), private, allocatable :: key(:)
   contains
     procedure, pass(this) :: init => htable_pt_init
     procedure, pass(this) :: free => htable_pt_free
     procedure, pass(this) :: set => htable_pt_set
     procedure, pass(this) :: get => htable_pt_get
     procedure, pass(this) :: hash => htable_pt_hash
     procedure, pass(this) :: remove => htable_pt_remove
  end type htable_pt_t

  !> Integer 2-tuple based hash table
  type, public, extends(htable_t) :: htable_i4t2_t
     type(tuple_i4_t), private, allocatable :: key(:)
   contains
     procedure, pass(this) :: init => htable_i4t2_init
     procedure, pass(this) :: free => htable_i4t2_free
     procedure, pass(this) :: set => htable_i4t2_set
     procedure, pass(this) :: get => htable_i4t2_get
     procedure, pass(this) :: hash => htable_i4t2_hash
     procedure, pass(this) :: remove => htable_i4t2_remove
  end type htable_i4t2_t

  !> Integer 4-tuple based hash table
  type, public, extends(htable_t) :: htable_i4t4_t
     type(tuple4_i4_t), private, allocatable :: key(:)
   contains
     procedure, pass(this) :: init => htable_i4t4_init
     procedure, pass(this) :: free => htable_i4t4_free
     procedure, pass(this) :: set => htable_i4t4_set
     procedure, pass(this) :: get => htable_i4t4_get
     procedure, pass(this) :: hash => htable_i4t4_hash
     procedure, pass(this) :: remove => htable_i4t4_remove
  end type htable_i4t4_t

  !> C pointer based hash table
  type, public, extends(htable_t) :: htable_cptr_t
     type(h_cptr_t), private, allocatable :: key(:)
   contains
     procedure, pass(this) :: init => htable_cptr_init
     procedure, pass(this) :: free => htable_cptr_free
     procedure, pass(this) :: set => htable_cptr_set
     procedure, pass(this) :: get => htable_cptr_get
     procedure, pass(this) :: hash => htable_cptr_hash
     procedure, pass(this) :: remove => htable_cptr_remove
  end type htable_cptr_t

  !
  ! Iterators
  !

  !> Base type for a hash table iterator
  type, public, abstract :: htable_iter_t
     integer, private :: n
     class(htable_t), pointer :: t
   contains
     procedure, public, pass(this) :: next => htable_iter_next
     procedure, public, pass(this) :: reset => htable_iter_reset
     procedure, public, pass(this) :: data => htable_iter_data
  end type htable_iter_t

  !> Iterator for an integer based hash table
  type, public, extends(htable_iter_t) :: htable_iter_i4_t
   contains
     procedure, pass(this) :: init => htable_iter_i4_init
     procedure, pass(this) :: value => htable_iter_i4_value
     procedure, pass(this) :: key => htable_iter_i4_key
     final :: htable_iter_i4_free
  end type htable_iter_i4_t

  !> Iterator for an integer*8 based hash table
  type, public, extends(htable_iter_t) :: htable_iter_i8_t
   contains
     procedure, pass(this) :: init => htable_iter_i8_init
     procedure, pass(this) :: value => htable_iter_i8_value
     procedure, pass(this) :: key => htable_iter_i8_key
     final :: htable_iter_i8_free
  end type htable_iter_i8_t

  !> Iterator for a double precision based hash table
  type, public, extends(htable_iter_t) :: htable_iter_r8_t
   contains
     procedure, pass(this) :: init => htable_iter_r8_init
     procedure, pass(this) :: value => htable_iter_r8_value
     procedure, pass(this) :: key => htable_iter_r8_key
     final :: htable_iter_r8_free
  end type htable_iter_r8_t

  !> Iterator for a point based hash table
  type, public, extends(htable_iter_t) :: htable_iter_pt_t
   contains
     procedure, pass(this) :: init => htable_iter_pt_init
     procedure, pass(this) :: value => htable_iter_pt_value
     procedure, pass(this) :: key => htable_iter_pt_key
     final :: htable_iter_pt_free
  end type htable_iter_pt_t

  !> Iterator for an integer based 2-tuple hash table
  type, public, extends(htable_iter_t) :: htable_iter_i4t2_t
   contains
     procedure, pass(this) :: init => htable_iter_i4t2_init
     procedure, pass(this) :: value => htable_iter_i4t2_value
     procedure, pass(this) :: key => htable_iter_i4t2_key
     final :: htable_iter_i4t2_free
  end type htable_iter_i4t2_t

  !> Iterator for an integer based 4-tuple hash table
  type, public, extends(htable_iter_t) :: htable_iter_i4t4_t
   contains
     procedure, pass(this) :: init => htable_iter_i4t4_init
     procedure, pass(this) :: value => htable_iter_i4t4_value
     procedure, pass(this) :: key => htable_iter_i4t4_key
     final :: htable_iter_i4t4_free
  end type htable_iter_i4t4_t

  !> Iterator for a C pointer based hash table
  type, public, extends(htable_iter_t) :: htable_iter_cptr_t
   contains
     procedure, pass(this) :: init => htable_iter_cptr_init
     procedure, pass(this) :: value => htable_iter_cptr_value
     procedure, pass(this) :: key => htable_iter_cptr_key
     final :: htable_iter_cptr_free
  end type htable_iter_cptr_t


  !
  ! Integer based implementation
  !
  interface

     !> Initialize an integer based hash table
     module subroutine htable_i4_init(this, size, data)
       class(htable_i4_t), intent(inout) :: this
       integer, value :: size                    !< Initial size of the table
       class(*), intent(inout), optional :: data !< Data to associate with @a key
     end subroutine htable_i4_init

     !> Destroy an integer  hash table
     module subroutine htable_i4_free(this)
       class(htable_i4_t), intent(inout) :: this
     end subroutine htable_i4_free

     !> Insert an integer into the hash table
     recursive module subroutine htable_i4_set(this, key, data)
       class(htable_i4_t), intent(inout) :: this
       integer, intent(inout) :: key   !< Table key
       class(*), intent(inout) :: data !< Data associated with @a key
     end subroutine htable_i4_set

     !> Retrive an integer with key @a key from the hash table
     module function htable_i4_get(this, key, data) result(rcode)
       class(htable_i4_t), intent(inout) :: this
       integer, intent(inout) :: key   !< Key to retrieve
       class(*), intent(inout) :: data !< Retrieved data
       integer :: rcode
     end function htable_i4_get

     !> Hash function for an integer based hash table
     pure module function htable_i4_hash(this, k, c) result(hash)
       class(htable_i4_t), intent(in) :: this
       class(*), intent(in) :: k
       integer, value :: c
       integer :: hash
     end function htable_i4_hash

     !> Remove an integer with key @a key from the hash table
     module subroutine htable_i4_remove(this, key)
       class(htable_i4_t), intent(inout) :: this
       integer, intent(inout) :: key !< Table key
     end subroutine htable_i4_remove

     !> Initialize an integer based hash table iterator
     module subroutine htable_iter_i4_init(this, t)
       class(htable_iter_i4_t), intent(inout) :: this
       type(htable_i4_t), target, intent(inout) :: t
     end subroutine htable_iter_i4_init

     !> Destroy an integer based hash table iterator
     module subroutine htable_iter_i4_free(this)
       type(htable_iter_i4_t), intent(inout) :: this
     end subroutine htable_iter_i4_free

     !> Return the current value of the integer based hash table iterator
     module function htable_iter_i4_value(this) result(value)
       class(htable_iter_i4_t), target, intent(inout) :: this
       integer, pointer :: value
     end function htable_iter_i4_value


     !> Return the current key of the integer based hash table iterator
     module function htable_iter_i4_key(this) result(key)
       class(htable_iter_i4_t), target, intent(inout) :: this
       integer, pointer :: key
     end function htable_iter_i4_key

  end interface

  !
  ! Integer*8 based implementation
  !
  interface

     !> Initialize an integer*8 based hash table
     !! @param size Initial size of the table
     !! @param data to associate with @a key
     module subroutine htable_i8_init(this, size, data)
       class(htable_i8_t), intent(inout) :: this
       integer, value :: size
       class(*), intent(inout), optional :: data
     end subroutine htable_i8_init

     !> Destroy an integer*8 hash table
     module subroutine htable_i8_free(this)
       class(htable_i8_t), intent(inout) :: this
     end subroutine htable_i8_free

     !> Insert an integer*8 into the hash table
     !! @param key in the hash table
     !! @param data to be associated with @a key
     recursive module subroutine htable_i8_set(this, key, data)
       class(htable_i8_t), intent(inout) :: this
       integer(kind=i8), intent(inout) :: key
       class(*), intent(inout) :: data
     end subroutine htable_i8_set

     !> Retrive an integer*8 with key @a key from the hash table
     !! @param key used to retrieve data
     !! @param data returned entry (with @a key) from the hash table
     module function htable_i8_get(this, key, data) result(rcode)
       class(htable_i8_t), intent(inout) :: this
       integer(kind=i8), intent(inout) :: key
       class(*), intent(inout) :: data
       integer :: rcode
     end function htable_i8_get

     !> Hash function for an integer*8 based hash table
     pure module function htable_i8_hash(this, k, c) result(hash)
       class(htable_i8_t), intent(in) :: this
       class(*), intent(in) :: k
       integer, value :: c
       integer :: hash
     end function htable_i8_hash

     !> Remove an integer*8 with key @a key from the hash table
     !! @param key to be removed
     module subroutine htable_i8_remove(this, key)
       class(htable_i8_t), intent(inout) :: this
       integer(kind=i8), intent(inout) :: key
     end subroutine htable_i8_remove

     !> Initialize an integer*8 based hash table iterator
     module subroutine htable_iter_i8_init(this, t)
       class(htable_iter_i8_t), intent(inout) :: this
       type(htable_i8_t), target, intent(inout) :: t
     end subroutine htable_iter_i8_init

     !> Destroy an integer*8 based hash table iterator
     module subroutine htable_iter_i8_free(this)
       type(htable_iter_i8_t), intent(inout) :: this
     end subroutine htable_iter_i8_free

     !> Return the current value of the integer*8 based hash table iterator
     module function htable_iter_i8_value(this) result(value)
       class(htable_iter_i8_t), target, intent(inout) :: this
       integer(kind=i8), pointer :: value
     end function htable_iter_i8_value

     !> Return the current key of the integer*8 based hash table iterator
     module function htable_iter_i8_key(this) result(key)
       class(htable_iter_i8_t), target, intent(inout) :: this
       integer(kind=i8), pointer :: key
     end function htable_iter_i8_key

  end interface

  !
  ! Double precision based implementation
  !
  interface

     !> Initialize a double precision based hash table
     !! @param size Initial size of the table
     !! @param data to associate with @a key
     module subroutine htable_r8_init(this, size, data)
       class(htable_r8_t), intent(inout) :: this
       integer, value :: size
       class(*), intent(inout), optional :: data
     end subroutine htable_r8_init

     !> Destroy a double precision based hash table
     module subroutine htable_r8_free(this)
       class(htable_r8_t), intent(inout) :: this
     end subroutine htable_r8_free

     !> Insert a double precision @a key (with @a data) into the hash table
     !! @param key in the hash table
     !! @param data to be associated with @a key
     recursive module subroutine htable_r8_set(this, key, data)
       class(htable_r8_t), intent(inout) :: this
       real(kind=dp), intent(inout) :: key
       class(*), intent(inout) :: data
     end subroutine htable_r8_set

     !> Retrive a double precision float with key @a key from the hash table
     !! @param key used to retrieve data
     !! @param data returned entry (with @a key) from the hash table
     module function htable_r8_get(this, key, data) result(rcode)
       class(htable_r8_t), intent(inout) :: this
       real(kind=dp), intent(inout) :: key
       class(*), intent(inout) :: data
       integer :: rcode
     end function htable_r8_get

     !> Hash function for a double precision based hash table
     pure module function htable_r8_hash(this, k, c) result(hash)
       class(htable_r8_t), intent(in) :: this
       class(*), intent(in) :: k
       integer, value :: c
       integer :: hash
     end function htable_r8_hash

     !> Remove a double precision key @a key from the hash table
     !! @param key to be removed
     module subroutine htable_r8_remove(this, key)
       class(htable_r8_t), intent(inout) :: this
       real(kind=dp), intent(inout) :: key
     end subroutine htable_r8_remove

     !> Initialize a double precision based hash table iterator
     module subroutine htable_iter_r8_init(this, t)
       class(htable_iter_r8_t), intent(inout) :: this
       type(htable_r8_t), target, intent(inout) :: t
     end subroutine htable_iter_r8_init

     !> Destroy a double precision based hash table iterator
     module subroutine htable_iter_r8_free(this)
       type(htable_iter_r8_t), intent(inout) :: this
     end subroutine htable_iter_r8_free

     !> Return the current value of the double precision based hash table iterator
     module function htable_iter_r8_value(this) result(value)
       class(htable_iter_r8_t), target, intent(inout) :: this
       real(kind=dp), pointer :: value
     end function htable_iter_r8_value

     !> Return the current key of the double precision based hash table iterator
     module function htable_iter_r8_key(this) result(key)
       class(htable_iter_r8_t), target, intent(inout) :: this
       real(kind=dp), pointer :: key
     end function htable_iter_r8_key

  end interface

  !
  ! Point based implementation
  !
  interface

     !> Initialize a point based hash table
     !! @param size Initial size of the table
     !! @param data to associate with @a key
     module subroutine htable_pt_init(this, size, data)
       class(htable_pt_t), intent(inout) :: this
       integer, value :: size
       class(*), intent(inout), optional :: data
     end subroutine htable_pt_init

     !> Destroy a point based hash table
     module subroutine htable_pt_free(this)
       class(htable_pt_t), intent(inout) :: this
     end subroutine htable_pt_free

     !> Insert a point @a key (with @a data) into the hash table
     !! @param key in the hash table
     !! @param data to be associated with @a key
     recursive module subroutine htable_pt_set(this, key, data)
       class(htable_pt_t), intent(inout) :: this
       type(point_t), intent(inout) :: key
       class(*), intent(inout) :: data
     end subroutine htable_pt_set

     !> Retrive a point with key @a key from the hash table
     !! @param key used to retrieve data
     !! @param data returned entry (with @a key) from the hash table
     module function htable_pt_get(this, key, data) result(rcode)
       class(htable_pt_t), intent(inout) :: this
       type(point_t), intent(inout) :: key
       class(*), intent(inout) :: data
       integer :: rcode
     end function htable_pt_get

     !> Hash function for a point based hash table
     pure module function htable_pt_hash(this, k, c) result(hash)
       class(htable_pt_t), intent(in) :: this
       class(*), intent(in) :: k
       integer, value :: c
       integer :: hash
     end function htable_pt_hash

     !> Remove a point with key @a key from the hash table
     !! @param key to be removed
     module subroutine htable_pt_remove(this, key)
       class(htable_pt_t), intent(inout) :: this
       type(point_t), intent(inout) :: key
     end subroutine htable_pt_remove

     !> Initialize a point based hash table iterator
     module subroutine htable_iter_pt_init(this, t)
       class(htable_iter_pt_t), intent(inout) :: this
       type(htable_pt_t), target, intent(inout) :: t
     end subroutine htable_iter_pt_init

     !> Destroy a point based hash table iterator
     module subroutine htable_iter_pt_free(this)
       type(htable_iter_pt_t), intent(inout) :: this
     end subroutine htable_iter_pt_free

     !> Return the current value of the point based hash table iterator
     module function htable_iter_pt_value(this) result(value)
       class(htable_iter_pt_t), target, intent(inout) :: this
       type(point_t), pointer :: value
     end function htable_iter_pt_value

     !> Return the current key of the point based hash table iterator
     module function htable_iter_pt_key(this) result(key)
       class(htable_iter_pt_t), target, intent(inout) :: this
       type(point_t), pointer :: key
     end function htable_iter_pt_key

  end interface

  !
  ! Integer 2-tuple based implementation
  !
  interface

     !> Initialize an integer 2-tuple hash table
     !! @param size Initial size of the table
     !! @param data to associate with @a key
     module subroutine htable_i4t2_init(this, size, data)
       class(htable_i4t2_t), intent(inout) :: this
       integer, value :: size
       class(*), intent(inout), optional :: data
     end subroutine htable_i4t2_init

     !> Destroy an integer 2-tuple based hash table
     module subroutine htable_i4t2_free(this)
       class(htable_i4t2_t), intent(inout) :: this
     end subroutine htable_i4t2_free

     !> Insert an integer 2-tuple into the hash table
     !! @param key in the hash table
     !! @param data to be associated with @a key
     recursive module subroutine htable_i4t2_set(this, key, data)
       class(htable_i4t2_t), intent(inout) :: this
       type(tuple_i4_t), intent(inout) :: key !< Table key
       class(*), intent(inout) :: data !< Data associated with @a key
     end subroutine htable_i4t2_set

     !> Retrive an integer 2-tuple with key @a key from the hash table
     !! @param key used to retrieve data
     !! @param data returned entry (with @a key) from the hash table
     module function htable_i4t2_get(this, key, data) result(rcode)
       class(htable_i4t2_t), intent(inout) :: this
       type(tuple_i4_t), intent(inout) :: key
       class(*), intent(inout) :: data
       integer :: rcode
     end function htable_i4t2_get

     !> Hash function for an integer 2-tuple hash table
     pure module function htable_i4t2_hash(this, k, c) result(hash)
       class(htable_i4t2_t), intent(in) :: this
       class(*), intent(in) :: k
       integer, value :: c
       integer :: hash
     end function htable_i4t2_hash

     !> Remove an integer 2-tuple with key @a key from the hash table
     !! @param key to be removed
     module subroutine htable_i4t2_remove(this, key)
       class(htable_i4t2_t), intent(inout) :: this
       type(tuple_i4_t), intent(inout) :: key !< Table key
     end subroutine htable_i4t2_remove

     !> Initialize an integer 2-tuple based hash table iterator
     module subroutine htable_iter_i4t2_init(this, t)
       class(htable_iter_i4t2_t), intent(inout) :: this
       type(htable_i4t2_t), target, intent(inout) :: t
     end subroutine htable_iter_i4t2_init

     !> Destroy an integer 2-tuple based  hash table iterator
     module subroutine htable_iter_i4t2_free(this)
       type(htable_iter_i4t2_t), intent(inout) :: this
     end subroutine htable_iter_i4t2_free

     !> Return the current value of integer based 2-tuple hash table iterator
     module function htable_iter_i4t2_value(this) result(value)
       class(htable_iter_i4t2_t), target, intent(inout) :: this
       type(tuple_i4_t), pointer :: value
     end function htable_iter_i4t2_value

     !> Return the current key of integer based 2-tuple hash table iterator
     module function htable_iter_i4t2_key(this) result(key)
       class(htable_iter_i4t2_t), target, intent(inout) :: this
       type(tuple_i4_t), pointer :: key
     end function htable_iter_i4t2_key

  end interface

  !
  ! Integer 4-tuple based implementation
  !
  interface

     !> Initialize an integer 4-tuple  hash table
     !! @param size Initial size of the table
     !! @param data to associate with @a key
     module subroutine htable_i4t4_init(this, size, data)
       class(htable_i4t4_t), intent(inout) :: this
       integer, value :: size
       class(*), intent(inout), optional :: data
     end subroutine htable_i4t4_init

     !> Destroy an integer 4-tuple based hash table
     module subroutine htable_i4t4_free(this)
       class(htable_i4t4_t), intent(inout) :: this
     end subroutine htable_i4t4_free

     !> Insert an integer 4-tuple into the hash table
     !! @param key in the hash table
     !! @param data to be associated with @a key
     recursive module subroutine htable_i4t4_set(this, key, data)
       class(htable_i4t4_t), intent(inout) :: this
       type(tuple4_i4_t), intent(inout) :: key
       class(*), intent(inout) :: data
     end subroutine htable_i4t4_set

     !> Retrive an integer 4-tuple with key @a key from the hash table
     !! @param key used to retrieve data
     !! @param data returned entry (with @a key) from the hash table
     module function htable_i4t4_get(this, key, data) result(rcode)
       class(htable_i4t4_t), intent(inout) :: this
       type(tuple4_i4_t), intent(inout) :: key
       class(*), intent(inout) :: data
       integer :: rcode
     end function htable_i4t4_get

     !> Hash function for an integer 4-tuple hash table
     module pure function htable_i4t4_hash(this, k, c) result(hash)
       class(htable_i4t4_t), intent(in) :: this
       class(*), intent(in) :: k
       integer, value :: c
       integer :: hash
     end function htable_i4t4_hash

     !> Remove an integer 4-tuple with key @a key from the hash table
     !! @param key to be removed
     module subroutine htable_i4t4_remove(this, key)
       class(htable_i4t4_t), intent(inout) :: this
       type(tuple4_i4_t), intent(inout) :: key
     end subroutine htable_i4t4_remove

     !> Initialize an integer 4-tuple based hash table iterator
     module subroutine htable_iter_i4t4_init(this, t)
       class(htable_iter_i4t4_t), intent(inout) :: this
       type(htable_i4t4_t), target, intent(inout) :: t
     end subroutine htable_iter_i4t4_init

     !> Destroy an integer 4-tuple based  hash table iterator
     module subroutine htable_iter_i4t4_free(this)
       type(htable_iter_i4t4_t), intent(inout) :: this
     end subroutine htable_iter_i4t4_free

     !> Return the current value of integer based 4-tuple hash table iterator
     module function htable_iter_i4t4_value(this) result(value)
       class(htable_iter_i4t4_t), target, intent(inout) :: this
       type(tuple4_i4_t), pointer :: value
     end function htable_iter_i4t4_value

     !> Return the current key of integer based 4-tuple hash table iterator
     module function htable_iter_i4t4_key(this) result(key)
       class(htable_iter_i4t4_t), target, intent(inout) :: this
       type(tuple4_i4_t), pointer :: key
     end function htable_iter_i4t4_key

  end interface

  !
  ! C pointer based implementation
  !
  interface

     !> Initialize a C pointer based  hash table
     !! @param size Initial size of the table
     !! @param data to associate with @a key
     module subroutine htable_cptr_init(this, size, data)
       class(htable_cptr_t), intent(inout) :: this
       integer, value :: size
       class(*), intent(inout), optional :: data
     end subroutine htable_cptr_init

     !> Destroy a C pointer based hash table
     module subroutine htable_cptr_free(this)
       class(htable_cptr_t), intent(inout) :: this
     end subroutine htable_cptr_free

     !> Insert a C pointer into the hash table
     !! @param key in the hash table
     !! @param data to be associated with @a key
     recursive module subroutine htable_cptr_set(this, key, data)
       class(htable_cptr_t), target, intent(inout) :: this
       type(h_cptr_t), intent(inout) :: key
       class(*), intent(inout) :: data
     end subroutine htable_cptr_set

     !> Retrive a C pointer with key @a key from the hash table
     !! @param key used to retrieve data
     !! @param data returned entry (with @a key) from the hash table
     module function htable_cptr_get(this, key, data) result(rcode)
       class(htable_cptr_t), target, intent(inout) :: this
       type(h_cptr_t), intent(inout) :: key
       class(*), intent(inout) :: data
       integer :: rcode
     end function htable_cptr_get

     !> Hash function for an integer 4-tuple hash table
     module pure function htable_cptr_hash(this, k, c) result(hash)
       class(htable_cptr_t), intent(in) :: this
       class(*), intent(in) :: k
       integer, value :: c
       integer :: hash
     end function htable_cptr_hash

     !> Remove a C pointer with key @a key from the hash table
     !! @param key to be removed
     module subroutine htable_cptr_remove(this, key)
       class(htable_cptr_t), target, intent(inout) :: this
       type(h_cptr_t), intent(inout) :: key
     end subroutine htable_cptr_remove

     !> Initialize a C pointer based hash table iterator
     module subroutine htable_iter_cptr_init(this, t)
       class(htable_iter_cptr_t), intent(inout) :: this
       type(htable_cptr_t), target, intent(inout) :: t
     end subroutine htable_iter_cptr_init

     !> Destroy a C pointer based  hash table iterator
     module subroutine htable_iter_cptr_free(this)
       type(htable_iter_cptr_t), intent(inout) :: this
     end subroutine htable_iter_cptr_free

     !> Return the current value of C pointer based hash table iterator
     module function htable_iter_cptr_value(this) result(value)
       class(htable_iter_cptr_t), target, intent(inout) :: this
       type(h_cptr_t), pointer :: value
     end function htable_iter_cptr_value

     !> Return the current key of a C pointer based hash table iterator
     module function htable_iter_cptr_key(this) result(key)
       class(htable_iter_cptr_t), target, intent(inout) :: this
       type(h_cptr_t), pointer :: key
     end function htable_iter_cptr_key

  end interface

contains

  !> Destroy a hash table
  subroutine htable_free(this)
    class(htable_t), intent(inout) :: this
    integer i

    if (allocated(this%data)) then
       deallocate(this%data)
    end if

    if (allocated(this%valid)) then
       deallocate(this%valid)
    end if

    if (allocated(this%skip)) then
       deallocate(this%skip)
    end if

    this%size = 0
    this%entries = 0

  end subroutine htable_free

  !> Clear all entries in a hash table
  subroutine htable_clear(this)
    class(htable_t), intent(inout) :: this

    if (allocated(this%valid)) then
       this%valid = .false.
       this%entries = 0
    else
       call neko_error("Hash table not allocated")
    end if

  end subroutine htable_clear

  !> Return number of entries in the table
  pure function htable_num_entries(this) result(entries)
    class(htable_t), intent(in) :: this
    integer :: entries
    entries = this%entries
  end function htable_num_entries

  !> Return total size of htable
  pure function htable_size(this) result(size)
    class(htable_t), intent(in) :: this
    integer :: size
    size = this%size
  end function htable_size

  !> Set data at @a idx to @a value
  subroutine htable_set_data(ht_data, data)
    class(*), intent(inout) :: ht_data !< Data entry in table
    class(*), intent(in) :: data       !< Data to set
    
    select type (data)
    type is (integer)
       select type(ht_data)
       type is (integer)
          ht_data = data
       end select
    type is (integer(i8))
       select type(ht_data)
       type is (integer(i8))
          ht_data = data
       end select
    type is (double precision)
       select type(ht_data)
       type is (double precision)
          ht_data = data
       end select
    type is (point_t)
       select type(ht_data)
       type is (point_t)
          ht_data = data
       end select
    class is (tuple_t)
       select type(ht_data)
       type is (tuple_i4_t)
          ht_data = data
       type is (tuple4_i4_t)
          ht_data = data
       end select
    type is (h_cptr_t)
       select type(ht_data)
       type is (h_cptr_t)
          ht_data = data
       end select
    class default
       call neko_error('Invalid htable data (set)')
    end select
  end subroutine htable_set_data

  !> Return data at @a idx in @a value
  subroutine htable_get_data(ht_data, data)
    class(*), intent(in) :: ht_data !< Data entry in table
    class(*), intent(inout) :: data !< Retrieved data

    select type (ht_data)
    type is (integer)
       select type(data)
       type is (integer)
          data = ht_data
       end select
    type is (integer(i8))
       select type(data)
       type is (integer(i8))
          data = ht_data
       end select
    type is (double precision)
       select type(data)
       type is (double precision)
          data = ht_data
       end select
    type is (point_t)
       select type(data)
       type is (point_t)
          data = ht_data
       end select
    class is (tuple_t)
       select type(data)
       type is (tuple_i4_t)
          data = ht_data
       type is (tuple4_i4_t)
          data = ht_data
       end select
    type is (h_cptr_t)
       select type (data)
       type is (h_cptr_t)
          data = ht_data
       end select
    class default
       call neko_error('Invalid htable data (get)')
    end select
  end subroutine htable_get_data

  !> Advance the iterator to the next valid table entry
  function htable_iter_next(this) result(valid)
    class(htable_iter_t), intent(inout) :: this
    logical :: valid

    this%n = this%n + 1
    do while ((.not. this%t%valid(this%n)) .and. (this%n .lt. this%t%size))
       this%n = this%n + 1
    end do

    valid = (this%n .lt. this%t%size)
    if (.not. valid) this%n = -1

  end function htable_iter_next

  !> Reset an iterator
  subroutine htable_iter_reset(this)
    class(htable_iter_t), intent(inout) :: this
    this%n = -1
  end subroutine htable_iter_reset

  !> Return the data at the current iterator position
  !! @attention this will make a deep copy of the data...
  !! @todo Remove once we figure out how to do this with value()
  !! for arbitrary data types
  subroutine htable_iter_data(this, data)
    class(htable_iter_t), target, intent(inout) :: this
    class(*), intent(inout) :: data !< Data to retrieve
    class(*), pointer :: hdp

    hdp => this%t%data(this%n)
    select type(hdp)
    type is (integer)
       select type (data)
       type is (integer)
          data = hdp
       end select
    type is (integer(i8))
       select type (data)
       type is (integer(i8))
          data = hdp
       end select
    type is (double precision)
       select type(data)
       type is (double precision)
          data = hdp
       end select
    type is (point_t)
       select type (data)
       type is (point_t)
          data = hdp
       end select
    class is (tuple_t)
       select type (data)
       type is (tuple_i4_t)
          data = hdp
       type is (tuple4_i4_t)
          data = hdp
       end select
    type is (h_cptr_t)
       select type (data)
       type is (h_cptr_t)
          data = hdp
       end select
    class default
       call neko_error('Invalid htable data (iter)')
    end select

  end subroutine htable_iter_data

end module htable
