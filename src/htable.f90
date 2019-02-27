!> Implements a hash table ADT
module htable
  use num_types
  use utils
  implicit none
  private

  !> \f$ ln(2) \f$
  real(kind=dp), parameter :: NEKO_M_LN2 = 0.693147180559945d0

  !> Hash table entry, tuple (key, data)
  type :: h_tuple_t
     integer :: key = -1 
     logical :: valid = .false.
     class(*), allocatable :: data
  end type h_tuple_t

  !> Base type for a hash table
  type, private, abstract :: htable_t
     integer, private :: size
     integer, private :: entries
     type(h_tuple_t), private, allocatable :: t(:)
   contains
     procedure, pass(this) :: htable_init
     procedure, pass(this) :: htable_clear
     procedure, pass(this) :: htable_free
     procedure, pass(this) :: htable_set
     procedure(htable_hash), pass(this), deferred :: hash
     generic, public :: clear => htable_clear
     generic, public :: free => htable_free
  end type htable_t

  !> Integer based hash table
  type, public, extends(htable_t) :: htable_i4_t
   contains
     procedure, pass(this) :: htable_i4_init
     procedure, pass(this) :: htable_i4_set
     procedure, pass(this) :: htable_i4_get
     procedure, pass(this) :: hash => htable_i4_hash
     generic, public :: init => htable_i4_init
     generic, public :: set => htable_i4_set
     generic, public :: get => htable_i4_get
  end type htable_i4_t

  !> Double precision based hash table
  type, public, extends(htable_t) :: htable_r8_t
   contains
     procedure, pass(this) :: htable_r8_init
     procedure, pass(this) :: htable_r8_set
     procedure, pass(this) :: htable_r8_get
     procedure, pass(this) :: hash => htable_r8_hash
     generic, public :: init => htable_r8_init
     generic, public :: set => htable_r8_set
     generic, public :: get => htable_r8_get
  end type htable_r8_t
 
  abstract interface
     pure function htable_hash(this, k) result(hash)
       import htable_t
       class(htable_t), intent(in) :: this
       class(*), value, intent(in) :: k
       integer :: hash
     end function htable_hash
  end interface

contains

  !> Initialize a hash table of type @a data
  subroutine htable_init(this, size, data)
    class(htable_t), intent(inout) :: this
    integer, value :: size
    class(*), intent(in) :: data
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

    do i = 0, size
       allocate(this%t(i)%data, source=data)
    end do
  end subroutine htable_init

  !> Destroy a hash table
  subroutine htable_free(this)
    class(htable_t), intent(inout) :: this
    integer i

    if (allocated(this%t)) then
       do i = 0, this%size !< @todo check range
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

  !> Insert tuple @a (key, value) into the hash table
  recursive subroutine htable_set(this, key, value) 
    class(htable_t), target, intent(inout) :: this
    integer, intent(in) :: key
    class(*), intent(in) ::  value
    class(htable_t), allocatable :: tmp
    integer index, i
    class(*), pointer :: dp

    index = this%hash(key)
    i = this%size - 1
    do while (i .ge. 0)
       if ((.not. this%t(index)%valid) .or. &
            (this%t(index)%key .eq. key)) then
          this%t(index)%key = key
          dp => this%t(index)%data
          select type(value)
          type is (integer)
             select type(dp)
             type is (integer)
                dp = value
             end select
          type is (double precision)
             select type(dp)
             type is (double precision)
                dp = value
             end select
          end select
          if (.not. this%t(index)%valid) then
             this%entries = this%entries + 1
          end if
          this%t(index)%valid = .true.
          return
       end if
       index = modulo((index + 1), this%size)
       i = i - 1
    end do  

    select type(value)
    type is (integer)
       allocate(htable_i4_t::tmp)
    type is (double precision)
       allocate(htable_r8_t::tmp)
    end select

    call htable_init(this, ishft(this%size, 1), value)

    do i = 0, this%size - 1
       if (this%t(i)%valid) then
          call htable_set(this, this%t(i)%key, this%t(i)%data)
       end if
    end do
    this%size = tmp%size
    call move_alloc(tmp%t, this%t)

    call htable_set(this, key, value)

  end subroutine htable_set

  !> Retrieve data associated with @a key into the hash table
  function htable_get(this, key, value) result(rcode)
    class(htable_t), target, intent(inout) :: this
    integer, intent(in) :: key
    class(*), intent(inout) :: value
    integer :: rcode
    class(*), pointer :: dp
    integer :: index, i

    index = this%hash(key)
    i = this%size
    
    do while (i .ge. 0)
       if ((this%t(index)%valid) .and. &
            this%t(index)%key .eq. key) then
          dp => this%t(index)%data
          select type (dp)
          type is (integer)
             select type(value)
             type is (integer)
                value = dp
             end select
          type is (double precision)
             select type(value)
             type is (double precision)
                value = dp
             end select
          end select
          rcode = 0
          return
       end if
       index = modulo((index + 1), this%size)
       i = i - 1
    end do
    rcode = 1
  end function htable_get

  !
  ! Integer based implementation
  !
  subroutine htable_i4_init(this, size)
    class(htable_i4_t), intent(inout) :: this
    integer, value :: size
    integer :: data

    call htable_init(this, size, data)
    
  end subroutine htable_i4_init

  subroutine htable_i4_set(this, key, value) 
    class(htable_i4_t), target, intent(inout) :: this
    integer, intent(in) :: key
    integer, intent(in) :: value
    class(*), allocatable :: data
    allocate(data, source=value)

    select type(data)
    type is (integer)
       data = value
       call htable_set(this, key, data)
    end select

  end subroutine htable_i4_set

  function htable_i4_get(this, key, value) result(rcode)
    class(htable_i4_t), target, intent(inout) :: this
    integer, intent(in) :: key
    integer, intent(inout) :: value
    integer :: rcode

    rcode = htable_get(this, key, value)

  end function htable_i4_get

  pure function htable_i4_hash(this, k) result(hash)
    class(htable_i4_t), intent(in) :: this
    integer, value, intent(in) :: k
    integer :: hash
    hash = modulo(k * (k + 3), this%size)
  end function htable_i4_hash

  !
  ! Double precision based implementation
  !
  subroutine htable_r8_init(this, size)
    class(htable_r8_t), intent(inout) :: this
    integer, value :: size
    real(kind=dp) :: data

    call htable_init(this, size, data)
    
  end subroutine htable_r8_init

  subroutine htable_r8_set(this, key, value) 
    class(htable_r8_t), target, intent(inout) :: this
    integer, intent(in) :: key
    real(kind=dp), intent(in) :: value
    class(*), allocatable :: data
    allocate(data, source=value)

    select type(data)
    type is (double precision)
       data = value
       call htable_set(this, key, data)
    end select

  end subroutine htable_r8_set

  function htable_r8_get(this, key, value) result(rcode)
    class(htable_r8_t), target, intent(inout) :: this
    integer, intent(in) :: key
    real(kind=dp), intent(inout) :: value
    integer :: rcode

    rcode = htable_get(this, key, value)

  end function htable_r8_get

  pure function htable_r8_hash(this, k) result(hash)
    class(htable_r8_t), intent(in) :: this
    real(kind=dp), value, intent(in) :: k
    integer :: hash
    hash = modulo(floor((2d0 * abs(fraction(k)) - 1d0) * 2**16), this%size)
  end function htable_r8_hash
  

end module htable
