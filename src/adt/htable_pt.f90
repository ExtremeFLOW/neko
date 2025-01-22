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
submodule (htable) htable_pt
  implicit none

contains
  
  !> Initialize a point based hash table
  module subroutine htable_pt_init(this, size, data)
    class(htable_pt_t), intent(inout) :: this
    integer, value :: size                    !< Initial size of the table
    class(*), intent(inout), optional :: data !< Data to associate with @a key

    call this%free()

    if (size .lt. 4) then
       size = 4
    end if

    size = ishft(1, ceiling(log(dble(size)) / NEKO_M_LN2))

    allocate(this%key(0:size))
    if (present(data)) then
       allocate(this%data(0:size), source=data)       
    else
       allocate(point_t::this%data(0:size))       
    end if

    allocate(this%valid(0:size))
    allocate(this%skip(0:size))
    this%valid = .false.
    this%skip = .false.
    this%size = size
    this%entries = 0

  end subroutine htable_pt_init

  module subroutine htable_pt_free(this)
    class(htable_pt_t), intent(inout) :: this

    if (allocated(this%key)) then
       deallocate(this%key)
    end if

    call this%free_base()
    
  end subroutine htable_pt_free

  !> Insert a point @a key (with @a data) into the hash table
  recursive module subroutine htable_pt_set(this, key, data)
    class(htable_pt_t), intent(inout) :: this
    type(point_t), intent(inout) :: key !< Table key
    class(*), intent(inout) :: data     !< Data associated with @a key
    type(htable_pt_t), allocatable :: tmp
    integer index, i, c

    c = 0
    i = log(1.0/this%size)/log(0.6)    
    index = 0

    do while (i .ge. 0)
       index = this%hash(key, c**2)
       if (index .lt. 0) then
          call neko_error("Invalid hash generated")
       end if
       !> Check if entry at this index is empty or if key matches
       if ((.not. this%valid(index)) .or. &
            this%key(index) .eq. key) then
          this%key(index) = key
          call this%set_data(this%data(index), data)
          if (.not. this%valid(index)) then
             this%entries = this%entries + 1
          end if
          this%valid(index) = .true.
          this%skip(index) = .false.
          return
       end if
       i = i - 1
       c = c + 1
    end do

    allocate(tmp)
    call tmp%init(ishft(this%size, 1), data)

    do i = 0, this%size - 1
       if (this%valid(i)) then
          call tmp%set(this%key(i), this%data(i))
       end if
    end do
    this%size = tmp%size
    this%entries = tmp%entries

    call move_alloc(tmp%key, this%key)
    call move_alloc(tmp%data, this%data)
    call move_alloc(tmp%valid, this%valid)
    call move_alloc(tmp%skip, this%skip)
!    call move_alloc(tmp%t, this%t)

    call this%set(key, data)

  end subroutine htable_pt_set

  !> Retrive a point with key @a key from the hash table
  module function htable_pt_get(this, key, data) result(rcode)
    class(htable_pt_t), intent(inout) :: this
    type(point_t), intent(inout) :: key !< Key to retrieve
    class(*), intent(inout) :: data     !< Retrieved data
    integer :: rcode
    integer :: index, i, c

    c = 0
    i = log(1.0/this%size)/log(0.6)

    do while (i .ge. 0)
       index = this%hash(key, c**2)
       if (index .lt. 0) then
          call neko_error("Invalid hash generated")
       end if

       if (.not. this%valid(index) .and. &
            .not. this%skip(index)) then
          rcode = 1
          return
       else if ((this%valid(index)) .and. &
            this%key(index) .eq. key) then
          call this%get_data(this%data(index), data)
          rcode = 0
          return
       end if
       i = i - 1
       c = c + 1
    end do
    rcode = 1

  end function htable_pt_get

  !> Hash function for a point based hash table
  pure module function htable_pt_hash(this, k, c) result(hash)
    class(htable_pt_t), intent(in) :: this
    class(*), intent(in) :: k
    integer, value :: c
    integer :: hash, i
    integer(kind=i8) :: hash2, tmp, mult
    integer(kind=i8), parameter :: M1 = int(Z'7ed55d15', i8)
    integer(kind=i8), parameter :: M2 = int(Z'c761c23c', i8)
    integer(kind=i8), parameter :: M3 = int(Z'165667b1', i8)
    integer(kind=i8), parameter :: M4 = int(Z'd3a2646c', i8)
    integer(kind=i8), parameter :: M5 = int(Z'fd7046c5', i8)
    integer(kind=i8), parameter :: M6 = int(Z'b55a4f09', i8)

    select type(k)
    type is (point_t)
       mult = 1000003
       hash2 = int(Z'345678')
       do i = 1, 3
          tmp = transfer(k%x(i), tmp)
          tmp = (tmp + M1) + ishft(tmp, 12)
          tmp = ieor(ieor(tmp, M2), ishft(tmp, -19))
          tmp = (tmp + M3) + ishft(tmp, 5)
          tmp = ieor((tmp + M4), ishft(tmp, 9))
          tmp = (tmp + M5) + ishft(tmp, 3)
          tmp = ieor(ieor(tmp, M6), ishft(tmp, -16))
          hash2 = ieor(hash2, tmp) * mult
          mult = mult + 82520 + 8
       end do
       hash2 = hash2 + 97531
       hash2 = modulo(hash2 + int(c, i8), int(this%size,i8))
       hash = int(hash2, i4)
    class default
       hash = -1
    end select

  end function htable_pt_hash

  !> Remove a point with key @a key from the hash table
  module subroutine htable_pt_remove(this, key)
    class(htable_pt_t), intent(inout) :: this
    type(point_t), intent(inout) :: key   !< Table key
        integer :: index, i, c

    c = 0
    i = log(1.0/this%size)/log(0.6)

    do while (i .ge. 0)
       index = this%hash(key, c**2)
       if (index .lt. 0) then
          call neko_error("Invalid hash generated")
       end if

       if ((this%valid(index)) .and. &
            this%key(index) .eq. key) then
          this%valid(index) = .false.
          this%skip(index) = .true.
          this%entries = this%entries - 1
          return
       end if
       i = i - 1
       c = c + 1
    end do

  end subroutine htable_pt_remove


  !> Initialize a point based hash table iterator
  module subroutine htable_iter_pt_init(this, t)
    class(htable_iter_pt_t), intent(inout) :: this
    type(htable_pt_t), target, intent(inout) :: t

    this%t => t
    this%n = -1

  end subroutine htable_iter_pt_init

  !> Destroy a point based hash table iterator
  module subroutine htable_iter_pt_free(this)
    type(htable_iter_pt_t), intent(inout) :: this
    nullify(this%t)
  end subroutine htable_iter_pt_free

  !> Return the current value of the point based hash table iterator
  module function htable_iter_pt_value(this) result(value)
    class(htable_iter_pt_t), target, intent(inout) :: this
    type(point_t), pointer :: value

    select type (hdp => this%t%data)
    type is (point_t)
       value => hdp(this%n)
    class default
       call neko_error('Key and data of different kind (pt)')
    end select

  end function htable_iter_pt_value

  !> Return the current key of the point based hash table iterator
  module function htable_iter_pt_key(this) result(key)
    class(htable_iter_pt_t), target, intent(inout) :: this
    type(point_t), pointer :: key

    select type (kp => this%t)
    type is (htable_pt_t)
       key => kp%key(this%n)
    class default
       call neko_error('Invalid key (pt)')
    end select

  end function htable_iter_pt_key
  
end submodule htable_pt
