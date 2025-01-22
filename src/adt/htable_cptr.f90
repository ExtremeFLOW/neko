submodule (htable) htable_cptr
  implicit none

contains
  
  !> Initialize a C pointer based  hash table
  module subroutine htable_cptr_init(this, size, data)
    class(htable_cptr_t), intent(inout) :: this
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
       allocate(h_cptr_t::this%data(0:size))
    end if

    allocate(this%valid(0:size))
    allocate(this%skip(0:size))
    this%valid = .false.
    this%skip = .false.
    this%size = size
    this%entries = 0

  end subroutine htable_cptr_init

  module subroutine htable_cptr_free(this)
    class(htable_cptr_t), intent(inout) :: this

    if (allocated(this%key)) then
       deallocate(this%key)
    end if

    call this%free_base()
    
  end subroutine htable_cptr_free

  !> Insert a C pointer into the hash table
  recursive module subroutine htable_cptr_set(this, key, data)
    class(htable_cptr_t), target, intent(inout) :: this
    type(h_cptr_t), intent(inout) :: key   !< Table key
    class(*), intent(inout) :: data !< Data associated with @a key
    type(htable_cptr_t), allocatable :: tmp
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
            c_associated(this%key(index)%ptr, key%ptr)) then
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

  end subroutine htable_cptr_set

  !> Retrive a C pointer with key @a key from the hash table
  module function htable_cptr_get(this, key, data) result(rcode)
    class(htable_cptr_t), target, intent(inout) :: this
    type(h_cptr_t), intent(inout) :: key   !< Key to retrieve
    class(*), intent(inout) :: data !< Retrieved data
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
            c_associated(this%key(index)%ptr, key%ptr)) then
          call this%get_data(this%data(index), data)
          rcode = 0
          return
       end if
       i = i - 1
       c = c + 1
    end do
    rcode = 1

  end function htable_cptr_get

  !> Hash function for an integer 4-tuple hash table
  pure module function htable_cptr_hash(this, k, c) result(hash)
    class(htable_cptr_t), intent(in) :: this
    class(*), intent(in) :: k
    integer, value :: c
    integer :: hash
    integer(kind=i8) :: k_int

    select type(k)
    type is (h_cptr_t)
       k_int = transfer(k%ptr, k_int)
       hash = int(modulo(k_int * 2654435761_i8 + int(c, i8),&
            int(this%size, i8)), i4)
    class default
       hash = -1
    end select
  end function htable_cptr_hash

  !> Remove a C pointer with key @a key from the hash table
  module subroutine htable_cptr_remove(this, key)
    class(htable_cptr_t), target, intent(inout) :: this
    type(h_cptr_t), intent(inout) :: key   !< Table key
    integer :: index, i, c

    c = 0
    i = log(1.0/this%size)/log(0.6)

    do while (i .ge. 0)
       index = this%hash(key, c**2)
       if (index .lt. 0) then
          call neko_error("Invalid hash generated")
       end if

       if ((this%valid(index)) .and. &
            c_associated(this%key(index)%ptr, key%ptr)) then
          this%valid(index) = .false.
          this%skip(index) = .true.
          this%entries = this%entries - 1
          return
       end if
       i = i - 1
       c = c + 1
    end do

  end subroutine htable_cptr_remove

  !> Initialize a C pointer based hash table iterator
  module subroutine htable_iter_cptr_init(this, t)
    class(htable_iter_cptr_t), intent(inout) :: this
    type(htable_cptr_t), target, intent(inout) :: t

    this%t => t
    this%n = -1

  end subroutine htable_iter_cptr_init

  !> Destroy a C pointer based  hash table iterator
  module subroutine htable_iter_cptr_free(this)
    type(htable_iter_cptr_t), intent(inout) :: this
    nullify(this%t)
  end subroutine htable_iter_cptr_free

  !> Return the current value of C pointer based hash table iterator
  module function htable_iter_cptr_value(this) result(value)
    class(htable_iter_cptr_t), target, intent(inout) :: this
    type(h_cptr_t), pointer :: value

    select type (hdp => this%t%data)
    type is (h_cptr_t)
       value => hdp(this%n)
    class default
       call neko_error('Key and data of different kind (cptr)')
    end select

  end function htable_iter_cptr_value

  !> Return the current key of a C pointer based hash table iterator
  module function htable_iter_cptr_key(this) result(key)
    class(htable_iter_cptr_t), target, intent(inout) :: this
    type(h_cptr_t), pointer :: key

    select type (kp => this%t)
    type is (htable_cptr_t)
       key => kp%key(this%n)
    class default
       call neko_error('Invalid key (cptr)')
    end select

  end function htable_iter_cptr_key
  
end submodule htable_cptr
