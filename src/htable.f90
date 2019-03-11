!> Implements a hash table ADT
!! @details A hash table storing tuples @a (key, data), with fixed
!! data-type @a key but with arbitrary typed @a data
module htable
  use num_types
  use utils
  use point
  use math, only : NEKO_M_LN2
  implicit none
  private

  !> Hash table entry, tuple (key, data)
  type :: h_tuple_t
     logical :: valid = .false.
     class(*), allocatable :: key      
     class(*), allocatable :: data
  end type h_tuple_t

  !> Base type for a hash table
  type, private, abstract :: htable_t
     integer, private :: size
     integer, private :: entries
     type(h_tuple_t), private, allocatable :: t(:)
   contains
     procedure, pass(this) :: htable_init
     procedure, pass(this) :: htable_set
     procedure, pass(this) :: eq_key => htable_eq_key
     procedure, pass(this) :: set_key => htable_set_key
     procedure, pass(this) :: set_data => htable_set_data
     procedure, pass(this) :: get_data => htable_get_data
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

  !> Integer based hash table
  type, public, extends(htable_t) :: htable_i4_t
   contains
     procedure, pass(this) :: init => htable_i4_init
     procedure, pass(this) :: set => htable_i4_set
     procedure, pass(this) :: get => htable_i4_get
     procedure, pass(this) :: hash => htable_i4_hash
  end type htable_i4_t

  !> Double precision based hash table
  type, public, extends(htable_t) :: htable_r8_t
   contains
     procedure, pass(this) :: init => htable_r8_init
     procedure, pass(this) :: set => htable_r8_set
     procedure, pass(this) :: get => htable_r8_get
     procedure, pass(this) :: hash => htable_r8_hash
  end type htable_r8_t

  !> Point based hash table
  type, public, extends(htable_t) :: htable_pt_t
   contains
     procedure, pass(this) :: init => htable_pt_init
     procedure, pass(this) :: set => htable_pt_set
     procedure, pass(this) :: get => htable_pt_get
     procedure, pass(this) :: hash => htable_pt_hash
  end type htable_pt_t
 
contains

  !> Initialize a hash table of type @a data
  subroutine htable_init(this, size, key, data)
    class(htable_t), intent(inout) :: this
    integer, value :: size              !< Initial size of the table
    class(*), target, intent(in) :: key            !< Type of key
    class(*), target, intent(in), optional :: data !< Type of data
    class(*), pointer :: dp
    integer :: i
    

    call this%free()
    
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
    class(htable_t), target, intent(inout) :: this
    class(*), target, intent(inout) :: key !< Table key
    class(*), intent(inout) ::  data       !< Data associated with @a key
    class(htable_t), allocatable :: tmp
    integer index, i

    index = this%hash(key)
    i = this%size - 1
    
    do while (i .ge. 0)
       if ((.not. this%t(index)%valid) .or. this%eq_key(index, key)) then
          call this%set_key(index, key)
          call this%set_data(index, data)
          if (.not. this%t(index)%valid) then
             this%entries = this%entries + 1
          end if
          this%t(index)%valid = .true.
          return
       end if
       index = modulo((index + 1), this%size)
       i = i - 1
    end do

    select type(key)
    type is (integer)
       allocate(htable_i4_t::tmp)
    type is (double precision)
       allocate(htable_r8_t::tmp)
    type is (point_t)
       allocate(htable_pt_t::tmp)
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
    class(htable_t), target, intent(inout) :: this
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
       if ((this%t(index)%valid) .and. this%eq_key(index, key)) then
          call this%get_data(index, data)          
          rcode = 0
          return
       end if
       index = modulo((index + 1), this%size)
       i = i - 1 
    end do
    rcode = 1
  end function htable_get

  !> Set data at @a idx to @a value
  subroutine htable_set_data(this, idx, data)
    class(htable_t), target, intent(inout) :: this
    integer, intent(in) :: idx   !< Table index
    class(*), intent(in) :: data !< Data to set at @a idx
    class(*), pointer :: dp

    dp => this%t(idx)%data
    select type (data)
    type is (integer)
       select type(dp)
       type is (integer)
          dp = data
       end select
    type is (double precision)
       select type(dp)
       type is (double precision)
          dp = data
       end select
    type is (point_t)
       select type(dp)
       type is (point_t)
          dp = data
       end select
    end select
  end subroutine htable_set_data

  !> Return data at @a idx in @a value
  subroutine htable_get_data(this, idx, data)
    class(htable_t), target, intent(in) :: this
    integer, intent(in) :: idx      !< Table index
    class(*), intent(inout) :: data !< Data to retrieve
    class(*), pointer :: dp

    dp => this%t(idx)%data
    select type (dp)
    type is (integer)
       select type(data)
       type is (integer)
          data = dp
       end select
    type is (double precision)
       select type(data)
       type is (double precision)
          data = dp
       end select
    type is (point_t)
       select type(data)
       type is (point_t)
          data = dp
       end select
    end select    
  end subroutine htable_get_data

  !> Compare key at @a idx to @a key
  function htable_eq_key(this, idx, key) result(res)
    class(htable_t), target, intent(in) :: this
    integer, intent(in) :: idx  !< Table index
    class(*), intent(in) :: key !< Key to compare against the key at @a idx
    class(*), pointer :: kp
    logical :: res

    kp => this%t(idx)%key
    res = .true.
    select type(kp)
    type is (integer)
       select type(key)
       type is (integer)
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
    end select
  end subroutine htable_set_key

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
    class(htable_i4_t), target, intent(inout) :: this
    integer, intent(inout) :: key   !< Table key
    class(*), intent(inout) :: data !< Data associated with @a key

    call htable_set(this, key, data)

  end subroutine htable_i4_set

  !> Retrive an integer with key @a key from the hash table
  function htable_i4_get(this, key, data) result(rcode)
    class(htable_i4_t), target, intent(inout) :: this
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
    select type(k)
    type is (integer)
       hash = modulo(k * (k + 3), this%size)
    class default
       hash = -1
    end select
  end function htable_i4_hash

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

  !> Insert an integer @a key (with @a data) into the hash table
  subroutine htable_r8_set(this, key, data) 
    class(htable_r8_t), target, intent(inout) :: this
    real(kind=dp), intent(inout) :: key !< Table key
    class(*), intent(inout) :: data     !< Data associated with @a key

    call htable_set(this, key, data)

  end subroutine htable_r8_set

  !> Retrive an integer with key @a key from the hash table
  function htable_r8_get(this, key, data) result(rcode)
    class(htable_r8_t), target, intent(inout) :: this
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

  !
  ! Point based implementation
  !
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

  subroutine htable_pt_set(this, key, data) 
    class(htable_pt_t), target, intent(inout) :: this
    type(point_t), intent(inout) :: key !< Table key
    class(*), intent(inout) :: data     !< Data associated with @a key

    call htable_set(this, key, data)

  end subroutine htable_pt_set

  function htable_pt_get(this, key, data) result(rcode)
    class(htable_pt_t), target, intent(inout) :: this
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
            abs(fraction(k%x(1) * k%x(2) * k%x(3))) - 1d0) * 2**16), &
            this%size) 
    class default
       hash = -1
    end select

  end function htable_pt_hash
  
end module htable
