! Copyright (c) 2019-2021, The Neko Authors
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
  use point
  use tuple
  use math, only : NEKO_M_LN2
  use, intrinsic :: iso_c_binding
  implicit none
  private

  !> Hash table entry, tuple (key, data)
  type :: h_tuple_t
     logical :: valid = .false.
     logical :: skip = .false.
     class(*), allocatable :: key      
     class(*), allocatable :: data
  end type h_tuple_t

  !> Base type for a hash table
  type, public, abstract :: htable_t
     integer, private :: size
     integer, private :: entries
     type(h_tuple_t), private, allocatable :: t(:)
   contains
     procedure(htable_hash), pass(this), deferred :: hash
     procedure, public, pass(this) :: clear => htable_clear
     procedure, public, pass(this) :: free => htable_free
     procedure, public, pass(this) :: num_entries => htable_num_entries
  end type htable_t

  abstract interface
     pure function htable_hash(this, k) result(hash)
       import htable_t
       class(htable_t), intent(in) :: this
       class(*), intent(in) :: k
       integer :: hash
     end function htable_hash
  end interface

  !
  ! Implementations
  !

  !> Integer based hash table
  type, public, extends(htable_t) :: htable_i4_t
   contains
     procedure, pass(this) :: init => htable_i4_init
     procedure, pass(this) :: set => htable_i4_set
     procedure, pass(this) :: get => htable_i4_get
     procedure, pass(this) :: hash => htable_i4_hash
     procedure, pass(this) :: remove => htable_i4_remove
  end type htable_i4_t

  !> Integer*8 based hash table
  type, public, extends(htable_t) :: htable_i8_t
   contains
     procedure, pass(this) :: init => htable_i8_init
     procedure, pass(this) :: set => htable_i8_set
     procedure, pass(this) :: get => htable_i8_get
     procedure, pass(this) :: hash => htable_i8_hash
     procedure, pass(this) :: remove => htable_i8_remove
  end type htable_i8_t

  !> Double precision based hash table
  type, public, extends(htable_t) :: htable_r8_t
   contains
     procedure, pass(this) :: init => htable_r8_init
     procedure, pass(this) :: set => htable_r8_set
     procedure, pass(this) :: get => htable_r8_get
     procedure, pass(this) :: hash => htable_r8_hash
     procedure, pass(this) :: remove => htable_r8_remove
  end type htable_r8_t

  !> Point based hash table
  type, public, extends(htable_t) :: htable_pt_t
   contains
     procedure, pass(this) :: init => htable_pt_init
     procedure, pass(this) :: set => htable_pt_set
     procedure, pass(this) :: get => htable_pt_get
     procedure, pass(this) :: hash => htable_pt_hash
     procedure, pass(this) :: remove => htable_pt_remove
  end type htable_pt_t

  !> Integer 2-tuple based hash table
  type, public, extends(htable_t) :: htable_i4t2_t
   contains
     procedure, pass(this) :: init => htable_i4t2_init
     procedure, pass(this) :: set => htable_i4t2_set
     procedure, pass(this) :: get => htable_i4t2_get
     procedure, pass(this) :: hash => htable_i4t2_hash
     procedure, pass(this) :: remove => htable_i4t2_remove
  end type htable_i4t2_t

  !> Integer 4-tuple based hash table
  type, public, extends(htable_t) :: htable_i4t4_t
   contains
     procedure, pass(this) :: init => htable_i4t4_init
     procedure, pass(this) :: set => htable_i4t4_set
     procedure, pass(this) :: get => htable_i4t4_get
     procedure, pass(this) :: hash => htable_i4t4_hash
     procedure, pass(this) :: remove => htable_i4t4_remove          
  end type htable_i4t4_t

  !> C pointer based hash table
  type, public, extends(htable_t) :: htable_cptr_t
   contains
     procedure, pass(this) :: init => htable_cptr_init
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
  end type htable_iter_i4_t

  !> Iterator for an integer*8 based hash table
  type, public, extends(htable_iter_t) :: htable_iter_i8_t
   contains
     procedure, pass(this) :: init => htable_iter_i8_init
     procedure, pass(this) :: value => htable_iter_i8_value
     procedure, pass(this) :: key => htable_iter_i8_key
  end type htable_iter_i8_t

  !> Iterator for a double precision based hash table
  type, public, extends(htable_iter_t) :: htable_iter_r8_t
   contains
     procedure, pass(this) :: init => htable_iter_r8_init
     procedure, pass(this) :: value => htable_iter_r8_value
     procedure, pass(this) :: key => htable_iter_r8_key
  end type htable_iter_r8_t

  !> Iterator for a point based hash table
  type, public, extends(htable_iter_t) :: htable_iter_pt_t
   contains
     procedure, pass(this) :: init => htable_iter_pt_init
     procedure, pass(this) :: value => htable_iter_pt_value
     procedure, pass(this) :: key => htable_iter_pt_key
  end type htable_iter_pt_t

  !> Iterator for an integer based 2-tuple hash table
  type, public, extends(htable_iter_t) :: htable_iter_i4t2_t
   contains
     procedure, pass(this) :: init => htable_iter_i4t2_init
     procedure, pass(this) :: value => htable_iter_i4t2_value
     procedure, pass(this) :: key => htable_iter_i4t2_key
  end type htable_iter_i4t2_t

  !> Iterator for an integer based 4-tuple hash table
  type, public, extends(htable_iter_t) :: htable_iter_i4t4_t
   contains
     procedure, pass(this) :: init => htable_iter_i4t4_init
     procedure, pass(this) :: value => htable_iter_i4t4_value
     procedure, pass(this) :: key => htable_iter_i4t4_key
  end type htable_iter_i4t4_t

  !> Iterator for a C pointer based hash table
  type, public, extends(htable_iter_t) :: htable_iter_cptr_t
   contains
     procedure, pass(this) :: init => htable_iter_cptr_init
     procedure, pass(this) :: value => htable_iter_cptr_value
     procedure, pass(this) :: key => htable_iter_cptr_key
  end type htable_iter_cptr_t

  !
  ! Type wrappers
  !
  type, public :: h_cptr_t
     type(c_ptr) :: ptr
  end type h_cptr_t

contains

  !> Initialize a hash table of type @a data
  subroutine htable_init(this, size, key, data)
    class(htable_t), intent(inout) :: this
    integer, value :: size              !< Initial size of the table
    class(*), target, intent(in) :: key            !< Type of key
    class(*), target, intent(in), optional :: data !< Type of data
    class(*), pointer :: dp
    integer :: i
    

    call htable_free(this)
    
    if (size .lt. 4) then
       size = 4
    end if

    size = ishft(1, ceiling(log(dble(size)) / NEKO_M_LN2))

    allocate(this%t(0:size))
    this%t(:)%valid = .false.
    this%size = size
    this%entries = 0

    dp => key
    if (present(data)) then
       dp => data
    end if

    do i = 0, size
       allocate(this%t(i)%key, source=key)
       allocate(this%t(i)%data, source=dp)
    end do
  end subroutine htable_init

  !> Destroy a hash table
  subroutine htable_free(this)
    class(htable_t), intent(inout) :: this
    integer i

    if (allocated(this%t)) then
       do i = 0, this%size !< @todo check range
          deallocate(this%t(i)%key)
          deallocate(this%t(i)%data)
       end do
       deallocate(this%t)
    end if

    this%size = 0
    this%entries = 0
    
  end subroutine htable_free

  !> Clear all entries in a hash table
  subroutine htable_clear(this)
    class(htable_t), intent(inout) :: this

    if (allocated(this%t)) then
       this%t(:)%valid = .false.
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

  !> Insert tuple @a (key, value) into the hash table
  recursive subroutine htable_set(this, key, data) 
    class(htable_t), intent(inout) :: this
    class(*), intent(inout) :: key   !< Table key
    class(*), intent(inout) ::  data !< Data associated with @a key
    class(htable_t), allocatable :: tmp
    integer index, i

    index = this%hash(key)
    if (index .lt. 0) then
       call neko_error("Invalid hash generated")
    end if
    
    i = this%size - 1
    
    do while (i .ge. 0)
       if ((.not. this%t(index)%valid) .or. &
            htable_eq_key(this, index, key)) then
          call htable_set_key(this, index, key)
          call htable_set_data(this, index, data)
          if (.not. this%t(index)%valid) then
             this%entries = this%entries + 1
          end if
          this%t(index)%valid = .true.
          this%t(index)%skip = .false.
          return
       end if
       index = modulo((index + 1), this%size)
       i = i - 1
    end do

    select type(key)
    type is (integer)
       allocate(htable_i4_t::tmp)
    type is (integer(i8))
       allocate(htable_i8_t::tmp)
    type is (double precision)
       allocate(htable_r8_t::tmp)
    type is (point_t)
       allocate(htable_pt_t::tmp)
    type is (tuple_i4_t)
       allocate(htable_i4t2_t::tmp)
    type is (tuple4_i4_t)
       allocate(htable_i4t4_t::tmp)
    type is (h_cptr_t)
       allocate(htable_cptr_t::tmp)
    class default
       call neko_error('Invalid htable key')
    end select

    call htable_init(tmp, ishft(this%size, 1), key, data)

    do i = 0, this%size - 1
       if (this%t(i)%valid) then
          call htable_set(tmp, this%t(i)%key, this%t(i)%data)
       end if
    end do
    this%size = tmp%size
    call move_alloc(tmp%t, this%t)

    call htable_set(this, key, data)

  end subroutine htable_set

  !> Retrieve data associated with @a key into the hash table
  function htable_get(this, key, data) result(rcode)
    class(htable_t), intent(inout) :: this
    class(*), intent(inout) :: key  !< Key to retrieve
    class(*), intent(inout) :: data !< Retrieved data
    integer :: rcode
    integer :: index, i

    index = this%hash(key)
    if (index .lt. 0) then
       call neko_error("Invalid hash generated")
    end if

    i = this%size - 1
    
    do while (i .ge. 0)
       if (.not. this%t(index)%valid .and. &
            .not. this%t(index)%skip) then
          rcode = 1
          return          
       else if ((this%t(index)%valid) .and. &
            htable_eq_key(this, index, key)) then
          call htable_get_data(this, index, data)
          rcode = 0
          return
       end if
       index = modulo((index + 1), this%size)
       i = i - 1 
    end do
    rcode = 1
  end function htable_get

  !> Remove a @a key from the hash table
  subroutine htable_remove(this, key)
    class(htable_t), intent(inout) :: this
    class(*), intent(inout) :: key  !< Key to remove
    integer :: index, i

    index = this%hash(key)
    if (index .lt. 0) then
       call neko_error("Invalid hash generated")
    end if

    i = this%size - 1
    
    do while (i .ge. 0)
       if ((this%t(index)%valid) .and. &
            htable_eq_key(this, index, key)) then
          this%t(index)%valid = .false.
          this%t(index)%skip = .true.
          this%entries = this%entries - 1
          return
       end if
       index = modulo((index + 1), this%size)
       i = i - 1 
    end do
  end subroutine htable_remove

  !> Set data at @a idx to @a value
  subroutine htable_set_data(this, idx, data)
    class(htable_t), target, intent(inout) :: this
    integer, intent(in) :: idx   !< Table index
    class(*), intent(in) :: data !< Data to set at @a idx
    class(*), pointer :: hdp

    hdp => this%t(idx)%data
    select type (data)
    type is (integer)
       select type(hdp)
       type is (integer)
          hdp = data
       end select
    type is (integer(i8))
       select type(hdp)
       type is (integer(i8))
          hdp = data
       end select
    type is (double precision)
       select type(hdp)
       type is (double precision)
          hdp = data
       end select
    type is (point_t)
       select type(hdp)
       type is (point_t)
          hdp = data
       end select
    type is (tuple_i4_t)
       select type(hdp)
       type is (tuple_i4_t)
          hdp = data
       end select
    type is (tuple4_i4_t)
       select type(hdp)
       type is (tuple4_i4_t)
          hdp = data
       end select
    type is (h_cptr_t)
       select type(hdp)
       type is (h_cptr_t)
          hdp = data
       end select
    class default
       call neko_error('Invalid htable data (set)')
    end select
  end subroutine htable_set_data

  !> Return data at @a idx in @a value
  subroutine htable_get_data(this, idx, data)
    class(htable_t), intent(in) :: this
    integer, intent(in) :: idx      !< Table index
    class(*), intent(inout) :: data !< Data to retrieve

    select type (hdp=>this%t(idx)%data)
    type is (integer)
       select type(data)
       type is (integer)
          data = hdp
       end select
    type is (integer(i8))
       select type(data)
       type is (integer(i8))
          data = hdp
       end select
    type is (double precision)
       select type(data)
       type is (double precision)
          data = hdp
       end select
    type is (point_t)
       select type(data)
       type is (point_t)
          data = hdp
       end select
    type is (tuple_i4_t)
       select type(data)
       type is (tuple_i4_t)
          data = hdp
       end select
    type is (tuple4_i4_t)
       select type(data)
       type is (tuple4_i4_t)
          data = hdp
       end select
    type is (h_cptr_t)
       select type (data)
       type is (h_cptr_t)
          data = hdp
       end select
    class default
       call neko_error('Invalid htable data (get)')
    end select
  end subroutine htable_get_data

  !> Compare key at @a idx to @a key
  pure function htable_eq_key(this, idx, key) result(res)
    class(htable_t), intent(in) :: this
    integer, intent(in) :: idx  !< Table index
    class(*), intent(in) :: key !< Key to compare against the key at @a idx
    logical :: res

    res = .true.
    select type (kp=>this%t(idx)%key)
    type is (integer)
       select type(key)
       type is (integer)
          res = (kp .eq. key)
       end select
    type is (integer(i8))
       select type(key)
       type is (integer(i8))
          res = (kp .eq. key)
       end select
    type is (double precision)
       select type(key)
       type is (double precision)
          res = (kp .eq. key)
       end select
    type is (point_t)
       select type (key)
       type is (point_t)
          res = (kp .eq. key)
       end select
    type is (tuple_i4_t)
       select type (key)
       type is (tuple_i4_t)
          res = (kp .eq. key)
       end select
    type is (tuple4_i4_t)
       select type (key)
       type is (tuple4_i4_t)
          res = (kp .eq. key)
       end select
    type is (h_cptr_t)
       select type (key)
       type is (h_cptr_t)
          res = c_associated(kp%ptr, key%ptr)
       end select
    end select
  end function htable_eq_key

  !> Set key at @a idx to @a key
  subroutine htable_set_key(this, idx, key) 
    class(htable_t), target, intent(inout) :: this
    integer, intent(in) :: idx  !< Table index
    class(*), intent(in) :: key !< Key to set at @a idx
    class(*), pointer :: kp
    
    kp => this%t(idx)%key
    select type(key)
    type is (integer)
       select type(kp)
       type is (integer)
          kp = key
       end select
    type is (integer(i8))
       select type(kp)
       type is (integer(i8))
          kp = key
       end select
    type is (double precision)
       select type(kp)
       type is (double precision)
          kp = key
       end select
    type is (point_t)
       select type (kp)
       type is (point_t)
          kp = key
       end select
    type is (tuple_i4_t)
       select type(kp)
       type is (tuple_i4_t)
          kp = key
       end select
    type is (tuple4_i4_t)
       select type(kp)
       type is (tuple4_i4_t)
          kp = key
       end select
    type is (h_cptr_t)
       select type(kp)
       type is (h_cptr_t)
          kp = key
       end select
    class default
       call neko_error('Invalid htable key (set)')
    end select
  end subroutine htable_set_key

  !> Advance the iterator to the next valid table entry
  function htable_iter_next(this) result(valid)
    class(htable_iter_t), intent(inout) :: this
    logical :: valid

    this%n = this%n + 1
    do while ((.not. this%t%t(this%n)%valid) .and. (this%n .lt. this%t%size))
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

    hdp => this%t%t(this%n)%data
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
    type is (tuple_i4_t)
       select type (data)
       type is (tuple_i4_t)
          data = hdp
       end select
    type is (tuple4_i4_t)
       select type (data)
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
  
  !
  ! Integer based implementation
  !
  !> Initialize an integer based hash table
  subroutine htable_i4_init(this, size, data)
    class(htable_i4_t), intent(inout) :: this
    integer, value :: size                    !< Initial size of the table
    class(*), intent(inout), optional :: data !< Data to associate with @a key
    integer :: key

    if (present(data)) then
       call htable_init(this, size, key, data)
    else
       call htable_init(this, size, key)
    end if
    
  end subroutine htable_i4_init

  !> Insert an integer into the hash table
  subroutine htable_i4_set(this, key, data) 
    class(htable_i4_t), intent(inout) :: this
    integer, intent(inout) :: key   !< Table key
    class(*), intent(inout) :: data !< Data associated with @a key

    call htable_set(this, key, data)

  end subroutine htable_i4_set

  !> Retrive an integer with key @a key from the hash table
  function htable_i4_get(this, key, data) result(rcode)
    class(htable_i4_t), intent(inout) :: this
    integer, intent(inout) :: key   !< Key to retrieve
    class(*), intent(inout) :: data !< Retrieved data
    integer :: rcode

    rcode = htable_get(this, key, data)

  end function htable_i4_get

  !> Hash function for an integer based hash table
  pure function htable_i4_hash(this, k) result(hash)
    class(htable_i4_t), intent(in) :: this
    class(*), intent(in) :: k
    integer :: hash
    integer, parameter :: M1 = int(Z'7ed55d1')
    integer, parameter :: M2 = int(Z'c761c23')
    integer, parameter :: M3 = int(Z'165667b')
    integer, parameter :: M4 = int(Z'd3a2646')
    integer, parameter :: M5 = int(Z'fd7046c')
    integer, parameter :: M6 = int(Z'b55a4f0')

    select type(k)
    type is (integer)
       hash = (k + M1) + ishft(k, 12)
       hash = ieor(ieor(hash, M2), ishft(hash, -19))
       hash = (hash + M3) + ishft(hash, 5)
       hash = ieor((hash + M4), ishft(hash, 9))
       hash = (hash + M5) + ishft(hash, 3)
       hash = ieor(ieor(hash, M6), ishft(hash, -16))
       hash = modulo(hash, this%size)
    class default
       hash = -1
    end select
  end function htable_i4_hash

  !> Remove an integer with key @a key from the hash table
  subroutine htable_i4_remove(this, key) 
    class(htable_i4_t), intent(inout) :: this
    integer, intent(inout) :: key   !< Table key

    call htable_remove(this, key)

  end subroutine htable_i4_remove

  !> Initialize an integer based hash table iterator
  subroutine htable_iter_i4_init(this, t)
    class(htable_iter_i4_t), intent(inout) :: this
    type(htable_i4_t), target, intent(inout) :: t

    this%t => t
    this%n = -1    

  end subroutine htable_iter_i4_init

  !> Return the current value of the integer based hash table iterator
  function htable_iter_i4_value(this) result(value)
    class(htable_iter_i4_t), intent(inout) :: this
    integer, pointer :: value

    select type (hdp => this%t%t(this%n)%data)
    type is (integer)
       value => hdp
    class default
       call neko_error('Key and data of different kind (i4)')
    end select
    
  end function htable_iter_i4_value

  !> Return the current key of the integer based hash table iterator
  function htable_iter_i4_key(this) result(key)
    class(htable_iter_i4_t), intent(inout) :: this
    integer, pointer :: key

    select type (kp => this%t%t(this%n)%key)
    type is (integer)
       key => kp
    class default
       call neko_error('Invalid key (i4)')
    end select
    
  end function htable_iter_i4_key

  !
  ! Integer*8 based implementation
  !
  !> Initialize an integer*8 based hash table
  subroutine htable_i8_init(this, size, data)
    class(htable_i8_t), intent(inout) :: this
    integer, value :: size                    !< Initial size of the table
    class(*), intent(inout), optional :: data !< Data to associate with @a key
    integer(kind=i8) :: key

    if (present(data)) then
       call htable_init(this, size, key, data)
    else
       call htable_init(this, size, key)
    end if
    
  end subroutine htable_i8_init

  !> Insert an integer*8 into the hash table
  subroutine htable_i8_set(this, key, data) 
    class(htable_i8_t), intent(inout) :: this
    integer(kind=i8), intent(inout) :: key   !< Table key
    class(*), intent(inout) :: data !< Data associated with @a key

    call htable_set(this, key, data)

  end subroutine htable_i8_set

  !> Retrive an integer*8 with key @a key from the hash table
  function htable_i8_get(this, key, data) result(rcode)
    class(htable_i8_t), intent(inout) :: this
    integer(kind=i8), intent(inout) :: key   !< Key to retrieve
    class(*), intent(inout) :: data !< Retrieved data
    integer :: rcode

    rcode = htable_get(this, key, data)

  end function htable_i8_get

  !> Hash function for an integer*8 based hash table
  pure function htable_i8_hash(this, k) result(hash)
    class(htable_i8_t), intent(in) :: this
    class(*), intent(in) :: k
    integer :: hash
    integer(kind=i8) :: tmp
    integer(kind=i8), parameter :: M1 = int(Z'7ed55d15', i8)
    integer(kind=i8), parameter :: M2 = int(Z'c761c23c', i8)
    integer(kind=i8), parameter :: M3 = int(Z'165667b1', i8)
    integer(kind=i8), parameter :: M4 = int(Z'd3a2646c', i8)
    integer(kind=i8), parameter :: M5 = int(Z'fd7046c5', i8)
    integer(kind=i8), parameter :: M6 = int(Z'b55a4f09', i8)

    select type(k)
    type is (integer(i8))
       tmp = (k + M1) + ishft(k, 12)
       tmp = ieor(ieor(tmp, M2), ishft(tmp, -19))
       tmp = (tmp + M3) + ishft(tmp, 5)
       tmp = ieor((tmp + M4), ishft(tmp, 9))
       tmp = (tmp + M5) + ishft(tmp, 3)
       tmp = ieor(ieor(tmp, M6), ishft(tmp, -16))
       hash = int(modulo(tmp, int(this%size, i8)), i4)
       !> @note I think this hash might be better
       hash = int(modulo(k * 2654435761_i8, int(this%size, i8)), i4)

    class default
       hash = -1
    end select
  end function htable_i8_hash

  !> Remove an integer*8 with key @a key from the hash table
  subroutine htable_i8_remove(this, key) 
    class(htable_i8_t), intent(inout) :: this
    integer(kind=i8), intent(inout) :: key   !< Table key

    call htable_remove(this, key)

  end subroutine htable_i8_remove

  !> Initialize an integer*8 based hash table iterator
  subroutine htable_iter_i8_init(this, t)
    class(htable_iter_i8_t), intent(inout) :: this
    type(htable_i8_t), target, intent(inout) :: t

    this%t => t
    this%n = -1    

  end subroutine htable_iter_i8_init

  !> Return the current value of the integer*8 based hash table iterator
  function htable_iter_i8_value(this) result(value)
    class(htable_iter_i8_t), target, intent(inout) :: this
    integer(kind=i8), pointer :: value


    select type (hdp => this%t%t(this%n)%data)
    type is (integer(i8))
       value => hdp
    class default
       call neko_error('Key and data of different kind (i8)')
    end select
    
  end function htable_iter_i8_value

  !> Return the current key of the integer*8 based hash table iterator
  function htable_iter_i8_key(this) result(key)
    class(htable_iter_i8_t), intent(inout) :: this
    integer(kind=i8), pointer :: key

    ! We should not need this extra select block, and it works great
    ! without it for GNU, Intel and NEC, but breaks horribly on Cray
    ! (>11.0.x) when using high opt. levels.
    select type(hti => this)
    type is(htable_iter_i8_t)
       select type (kp => hti%t%t(this%n)%key)
       type is (integer(i8))
          key => kp
       class default
          call neko_error('Invalid key (i8)')
       end select
    class default
       call neko_error('Corrupt htable iter. (i8)')
    end select
    
  end function htable_iter_i8_key
  

  !
  ! Double precision based implementation
  !
  !> Initialize a double precision based hash table
  subroutine htable_r8_init(this, size, data)
    class(htable_r8_t), intent(inout) :: this
    integer, value :: size                    !< Initial size of the table
    class(*), intent(inout), optional :: data !< Data to associate with @a key
    real(kind=dp) :: key

    if (present(data)) then
       call htable_init(this, size, key, data)
    else
       call htable_init(this, size, key)
    end if
    
  end subroutine htable_r8_init

  !> Insert a double precision @a key (with @a data) into the hash table
  subroutine htable_r8_set(this, key, data) 
    class(htable_r8_t), intent(inout) :: this
    real(kind=dp), intent(inout) :: key !< Table key
    class(*), intent(inout) :: data     !< Data associated with @a key

    call htable_set(this, key, data)

  end subroutine htable_r8_set

  !> Retrive a double precision float with key @a key from the hash table
  function htable_r8_get(this, key, data) result(rcode)
    class(htable_r8_t), intent(inout) :: this
    real(kind=dp), intent(inout) :: key !< Key to retrieve
    class(*), intent(inout) :: data     !< Retrieved data
    integer :: rcode

    rcode = htable_get(this, key, data)

  end function htable_r8_get

  !> Hash function for a double precision based hash table
  pure function htable_r8_hash(this, k) result(hash)
    class(htable_r8_t), intent(in) :: this
    class(*), intent(in) :: k
    integer :: hash
    select type(k)
    type is (double precision)
       hash = modulo(floor((2d0 * abs(fraction(k)) - 1d0) * 2**16), this%size)
    class default
       hash = -1
    end select
  end function htable_r8_hash

  !> Remove a double precision key @a key from the hash table
  subroutine htable_r8_remove(this, key) 
    class(htable_r8_t), intent(inout) :: this
    real(kind=dp), intent(inout) :: key   !< Table key

    call htable_remove(this, key)

  end subroutine htable_r8_remove
  

  !> Initialize a double precision based hash table iterator
  subroutine htable_iter_r8_init(this, t)
    class(htable_iter_r8_t), intent(inout) :: this
    type(htable_r8_t), target, intent(inout) :: t

    this%t => t
    this%n = -1    

  end subroutine htable_iter_r8_init

  !> Return the current value of the double precision based hash table iterator
  function htable_iter_r8_value(this) result(value)
    class(htable_iter_r8_t), intent(inout) :: this
    real(kind=dp), pointer :: value

    select type (hdp => this%t%t(this%n)%data)
    type is (double precision)
       value => hdp
    class default
       call neko_error('Key and data of different kind (r8)')
    end select
    
  end function htable_iter_r8_value

  !> Return the current key of the double precision based hash table iterator
  function htable_iter_r8_key(this) result(key)
    class(htable_iter_r8_t), intent(inout) :: this
    real(kind=dp), pointer :: key

    select type (kp => this%t%t(this%n)%key)
    type is (double precision)
       key => kp
    class default
       call neko_error('Invalid key (r8)')
    end select
    
  end function htable_iter_r8_key
  
  !
  ! Point based implementation
  !
  !> Initialize a point based hash table
  subroutine htable_pt_init(this, size, data)
    class(htable_pt_t), intent(inout) :: this
    integer, value :: size                    !< Initial size of the table
    class(*), intent(inout), optional :: data !< Data to associate with @a key
    type(point_t) :: key

    if (present(data)) then
       call htable_init(this, size, key, data)
    else
       call htable_init(this, size, key)
    end if
    
  end subroutine htable_pt_init

  !> Insert a point @a key (with @a data) into the hash table
  subroutine htable_pt_set(this, key, data) 
    class(htable_pt_t), intent(inout) :: this
    type(point_t), intent(inout) :: key !< Table key
    class(*), intent(inout) :: data     !< Data associated with @a key

    call htable_set(this, key, data)

  end subroutine htable_pt_set

  !> Retrive a point with key @a key from the hash table
  function htable_pt_get(this, key, data) result(rcode)
    class(htable_pt_t), intent(inout) :: this
    type(point_t), intent(inout) :: key !< Key to retrieve
    class(*), intent(inout) :: data     !< Retrieved data
    integer :: rcode

    rcode = htable_get(this, key, data)

  end function htable_pt_get

  !> Hash function for a point based hash table
  pure function htable_pt_hash(this, k) result(hash)
    class(htable_pt_t), intent(in) :: this
    class(*), intent(in) :: k
    integer :: hash
    select type(k)
    type is (point_t)
       hash = modulo(floor((2d0 * &
            abs(fraction(k%x(1)*1d-1 + k%x(2)*1d-2 + k%x(3)*1d-3)) - 1d0) * 2**30), &
            this%size) 
    class default
       hash = -1
    end select

  end function htable_pt_hash

  !> Remove a point with key @a key from the hash table
  subroutine htable_pt_remove(this, key) 
    class(htable_pt_t), intent(inout) :: this
    type(point_t), intent(inout) :: key   !< Table key

    call htable_remove(this, key)

  end subroutine htable_pt_remove


  !> Initialize a point based hash table iterator
  subroutine htable_iter_pt_init(this, t)
    class(htable_iter_pt_t), intent(inout) :: this
    type(htable_pt_t), target, intent(inout) :: t

    this%t => t
    this%n = -1    

  end subroutine htable_iter_pt_init

  !> Return the current value of the point based hash table iterator
  function htable_iter_pt_value(this) result(value)
    class(htable_iter_pt_t), intent(inout) :: this
    type(point_t), pointer :: value

    select type (hdp => this%t%t(this%n)%data)
    type is (point_t)
       value => hdp
    class default
       call neko_error('Key and data of different kind (pt)')
    end select
    
  end function htable_iter_pt_value

  !> Return the current key of the point based hash table iterator
  function htable_iter_pt_key(this) result(key)
    class(htable_iter_pt_t), intent(inout) :: this
    type(point_t), pointer :: key

    select type (kp => this%t%t(this%n)%key)
    type is (point_t)
       key => kp
    class default
       call neko_error('Invalid key (pt)')
    end select
    
  end function htable_iter_pt_key

  !
  ! Integer 2-tuple based implementation
  !
  !> Initialize an integer 2-tuple hash table
  subroutine htable_i4t2_init(this, size, data)
    class(htable_i4t2_t), intent(inout) :: this
    integer, value :: size                    !< Initial size of the table
    class(*), intent(inout), optional :: data !< Data to associate with @a key
    type(tuple_i4_t) :: key

    if (present(data)) then
       call htable_init(this, size, key, data)
    else
       call htable_init(this, size, key)
    end if
    
  end subroutine htable_i4t2_init

  !> Insert an integer 2-tuple into the hash table
  subroutine htable_i4t2_set(this, key, data) 
    class(htable_i4t2_t), intent(inout) :: this
    type(tuple_i4_t), intent(inout) :: key   !< Table key
    class(*), intent(inout) :: data !< Data associated with @a key

    call htable_set(this, key, data)

  end subroutine htable_i4t2_set

  !> Retrive an integer 2-tuple with key @a key from the hash table
  function htable_i4t2_get(this, key, data) result(rcode)
    class(htable_i4t2_t), intent(inout) :: this
    type(tuple_i4_t), intent(inout) :: key   !< Key to retrieve
    class(*), intent(inout) :: data !< Retrieved data
    integer :: rcode

    rcode = htable_get(this, key, data)

  end function htable_i4t2_get

  !> Hash function for an integer 2-tuple hash table
  pure function htable_i4t2_hash(this, k) result(hash)
    class(htable_i4t2_t), intent(in) :: this
    class(*), intent(in) :: k
    integer :: i, tmp, mult, hash
    integer, parameter :: M1 = int(Z'7ed55d1')
    integer, parameter :: M2 = int(Z'c761c23')
    integer, parameter :: M3 = int(Z'165667b')
    integer, parameter :: M4 = int(Z'd3a2646')
    integer, parameter :: M5 = int(Z'fd7046c')
    integer, parameter :: M6 = int(Z'b55a4f0')

    select type(k)
    type is (tuple_i4_t)
       mult = 1000003
       hash = int(Z'345678')
       do i = 1, 2
          tmp = k%x(i)
          tmp = (tmp + M1) + ishft(tmp, 12)
          tmp = ieor(ieor(tmp, M2), ishft(tmp, -19))
          tmp = (tmp + M3) + ishft(tmp, 5)
          tmp = ieor((tmp + M4), ishft(tmp, 9))
          tmp = (tmp + M5) + ishft(tmp, 3)
          tmp = ieor(ieor(tmp, M6), ishft(tmp, -16))
          hash = ieor(hash, tmp) * mult
          mult = mult + 82520 + 4
       end do
       hash = hash + 97531
       hash = modulo(hash, this%size)
    class default
       hash = -1
    end select
  end function htable_i4t2_hash

  !> Remove an integer 2-tuple with key @a key from the hash table
  subroutine htable_i4t2_remove(this, key) 
    class(htable_i4t2_t), intent(inout) :: this
    type(tuple_i4_t), intent(inout) :: key   !< Table key

    call htable_remove(this, key)

  end subroutine htable_i4t2_remove
  
  !> Initialize an integer 2-tuple based hash table iterator
  subroutine htable_iter_i4t2_init(this, t)
    class(htable_iter_i4t2_t), intent(inout) :: this
    type(htable_i4t2_t), target, intent(inout) :: t

    this%t => t
    this%n = -1    

  end subroutine htable_iter_i4t2_init

  !> Return the current value of integer based 2-tuple hash table iterator
  function htable_iter_i4t2_value(this) result(value)
    class(htable_iter_i4t2_t), intent(inout) :: this
    type(tuple_i4_t), pointer :: value

    select type (hdp => this%t%t(this%n)%data)
    type is (tuple_i4_t)
       value => hdp
    class default
       call neko_error('Key and data of different kind (i4t2)')
    end select
    
  end function htable_iter_i4t2_value

  !> Return the current key of integer based 2-tuple hash table iterator
  function htable_iter_i4t2_key(this) result(key)
    class(htable_iter_i4t2_t), intent(inout) :: this
    type(tuple_i4_t), pointer :: key

    select type (kp => this%t%t(this%n)%key)
    type is (tuple_i4_t)
       key => kp
    class default
       call neko_error('Invalid key (i4t2)')
    end select
    
  end function htable_iter_i4t2_key

  !
  ! Integer 4-tuple based implementation
  !
  !> Initialize an integer 4-tuple  hash table
  subroutine htable_i4t4_init(this, size, data)
    class(htable_i4t4_t), intent(inout) :: this
    integer, value :: size                    !< Initial size of the table
    class(*), intent(inout), optional :: data !< Data to associate with @a key
    type(tuple4_i4_t) :: key

    if (present(data)) then
       call htable_init(this, size, key, data)
    else
       call htable_init(this, size, key)
    end if
    
  end subroutine htable_i4t4_init

  !> Insert an integer 4-tuple into the hash table
  subroutine htable_i4t4_set(this, key, data) 
    class(htable_i4t4_t), intent(inout) :: this
    type(tuple4_i4_t), intent(inout) :: key   !< Table key
    class(*), intent(inout) :: data !< Data associated with @a key

    call htable_set(this, key, data)

  end subroutine htable_i4t4_set

  !> Retrive an integer 4-tuple with key @a key from the hash table
  function htable_i4t4_get(this, key, data) result(rcode)
    class(htable_i4t4_t), intent(inout) :: this
    type(tuple4_i4_t), intent(inout) :: key   !< Key to retrieve
    class(*), intent(inout) :: data !< Retrieved data
    integer :: rcode

    rcode = htable_get(this, key, data)

  end function htable_i4t4_get

  !> Hash function for an integer 4-tuple hash table
  pure function htable_i4t4_hash(this, k) result(hash)
    class(htable_i4t4_t), intent(in) :: this
    class(*), intent(in) :: k
    integer :: i, tmp, mult, hash
    integer, parameter :: M1 = int(Z'7ed55d1')
    integer, parameter :: M2 = int(Z'c761c23')
    integer, parameter :: M3 = int(Z'165667b')
    integer, parameter :: M4 = int(Z'd3a2646')
    integer, parameter :: M5 = int(Z'fd7046c')
    integer, parameter :: M6 = int(Z'b55a4f0')

    select type(k)
    type is (tuple4_i4_t)
       mult = 1000003
       hash = int(Z'345678')
       do i = 1, 4
          tmp = k%x(i)
          tmp = (tmp + M1) + ishft(tmp, 12)
          tmp = ieor(ieor(tmp, M2), ishft(tmp, -19))
          tmp = (tmp + M3) + ishft(tmp, 5)
          tmp = ieor((tmp + M4), ishft(tmp, 9))
          tmp = (tmp + M5) + ishft(tmp, 3)
          tmp = ieor(ieor(tmp, M6), ishft(tmp, -16))
          hash = ieor(hash, tmp) * mult
          mult = mult + 82520 + 8
       end do
       hash = hash + 97531
       hash = modulo(hash, this%size)
    class default
       hash = -1
    end select
  end function htable_i4t4_hash

  !> Remove an integer 4-tuple with key @a key from the hash table
  subroutine htable_i4t4_remove(this, key) 
    class(htable_i4t4_t), intent(inout) :: this
    type(tuple4_i4_t), intent(inout) :: key   !< Table key

    call htable_remove(this, key)

  end subroutine htable_i4t4_remove

  !> Initialize an integer 4-tuple based hash table iterator
  subroutine htable_iter_i4t4_init(this, t)
    class(htable_iter_i4t4_t), intent(inout) :: this
    type(htable_i4t4_t), target, intent(inout) :: t

    this%t => t
    this%n = -1    

  end subroutine htable_iter_i4t4_init

  !> Return the current value of integer based 4-tuple hash table iterator
  function htable_iter_i4t4_value(this) result(value)
    class(htable_iter_i4t4_t), intent(inout) :: this
    type(tuple4_i4_t), pointer :: value

    select type (hdp => this%t%t(this%n)%data)
    type is (tuple4_i4_t)
       value => hdp
    class default
       call neko_error('Key and data of different kind (i4t4)')
    end select
    
  end function htable_iter_i4t4_value

  !> Return the current key of integer based 4-tuple hash table iterator
  function htable_iter_i4t4_key(this) result(key)
    class(htable_iter_i4t4_t), target, intent(inout) :: this
    type(tuple4_i4_t), pointer :: key

    ! We should not need this extra select block, and it works great
    ! without it for GNU, Intel and NEC, but breaks horribly on Cray
    ! (>11.0.x) when using high opt. levels.
    select type(hti => this)
    type is(htable_iter_i4t4_t)
       select type (kp => hti%t%t(this%n)%key)
       type is (tuple4_i4_t)
          key => kp
       class default
          call neko_error('Invalid key (i4t4)')
       end select
    class default
       call neko_error('Corrupt htable iter. (i4t4)')
    end select
    
  end function htable_iter_i4t4_key

  !
  ! C pointer based implementation
  !
  !> Initialize a C pointer based  hash table
  subroutine htable_cptr_init(this, size, data)
    class(htable_cptr_t), intent(inout) :: this
    integer, value :: size                    !< Initial size of the table
    class(*), intent(inout), optional :: data !< Data to associate with @a key
    type(h_cptr_t) :: key
    
    if (present(data)) then
       call htable_init(this, size, key, data)
    else
       call htable_init(this, size, key)
    end if
    
  end subroutine htable_cptr_init

  !> Insert a C pointer into the hash table
  subroutine htable_cptr_set(this, key, data) 
    class(htable_cptr_t), target, intent(inout) :: this
    type(h_cptr_t), intent(inout) :: key   !< Table key
    class(*), intent(inout) :: data !< Data associated with @a key

    call htable_set(this, key, data)

  end subroutine htable_cptr_set

  !> Retrive a C pointer with key @a key from the hash table
  function htable_cptr_get(this, key, data) result(rcode)
    class(htable_cptr_t), target, intent(inout) :: this
    type(h_cptr_t), intent(inout) :: key   !< Key to retrieve
    class(*), intent(inout) :: data !< Retrieved data
    integer :: rcode

    rcode = htable_get(this, key, data)

  end function htable_cptr_get

  !> Hash function for an integer 4-tuple hash table
  pure function htable_cptr_hash(this, k) result(hash)
    class(htable_cptr_t), intent(in) :: this
    class(*), intent(in) :: k
    integer :: hash
    integer(kind=i8) :: k_int

    select type(k)
    type is (h_cptr_t)
       k_int = transfer(k%ptr, k_int)
       hash = int(modulo(k_int * 2654435761_i8, int(this%size, i8)), i4)
    class default
       hash = -1
    end select
  end function htable_cptr_hash

  !> Remove a C pointer with key @a key from the hash table
  subroutine htable_cptr_remove(this, key) 
    class(htable_cptr_t), target, intent(inout) :: this
    type(h_cptr_t), intent(inout) :: key   !< Table key

    call htable_remove(this, key)

  end subroutine htable_cptr_remove

  !> Initialize a C pointer based hash table iterator
  subroutine htable_iter_cptr_init(this, t)
    class(htable_iter_cptr_t), intent(inout) :: this
    type(htable_cptr_t), target, intent(inout) :: t

    this%t => t
    this%n = -1    

  end subroutine htable_iter_cptr_init

  !> Return the current value of C pointer based hash table iterator
  function htable_iter_cptr_value(this) result(value)
    class(htable_iter_cptr_t), target, intent(inout) :: this
    class(*), pointer :: hdp
    type(h_cptr_t), pointer :: value

    hdp => this%t%t(this%n)%data
    select type (hdp)
    type is (h_cptr_t)
       value => hdp
    class default
       call neko_error('Key and data of different kind (cptr)')
    end select
    
  end function htable_iter_cptr_value

  !> Return the current key of a C pointer based hash table iterator
  function htable_iter_cptr_key(this) result(key)
    class(htable_iter_cptr_t), target, intent(inout) :: this
    class(*), pointer :: kp
    type(h_cptr_t), pointer :: key

    kp => this%t%t(this%n)%key
    select type (kp)
    type is (h_cptr_t)
       key => kp
    class default
       call neko_error('Invalid key (cptr)')
    end select
    
  end function htable_iter_cptr_key

end module htable
