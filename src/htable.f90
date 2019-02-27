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
  type, public, abstract :: htable_t
     integer, private :: size
     integer, private :: entries
     type(h_tuple_t), allocatable :: t(:)
   contains
     procedure, pass(this) :: htable_clear
     procedure, pass(this) :: htable_free
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
     procedure, pass(this) :: hash =>  htable_i4_hash
     generic, public :: init => htable_i4_init
     generic, public :: set => htable_i4_set
     generic, public :: get => htable_i4_get
  end type htable_i4_t
 
  abstract interface
     pure function htable_hash(this, k) result(hash)
       import htable_t
       class(htable_t), intent(in) :: this
       integer, value, intent(in) :: k
       integer :: hash
     end function htable_hash
  end interface

contains

  !> Destroy hash table
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

  subroutine htable_i4_init(this, size)
    class(htable_i4_t), intent(inout) :: this
    integer, value :: size
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
       allocate(this%t(i)%data, source=i)
    end do
  end subroutine htable_i4_init

  recursive subroutine htable_i4_set(this, key, value) 
    class(htable_i4_t), target, intent(inout) :: this
    integer, intent(in) :: key
    integer, intent(in) :: value
    type(htable_i4_t) :: tmp
    class(*), pointer :: dp
    integer index, i

    index = this%hash(key)
    i = this%size - 1
    do while (i .ge. 0)
       if ((.not. this%t(index)%valid) .or. &
            (this%t(index)%key .eq. key)) then
          this%t(index)%key = key
          dp => this%t(index)%data
          select type(dp)
          type is (integer)
             dp = value
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

    call tmp%init(ishft(this%size, 1))
    do i = 0, this%size - 1
       if (this%t(i)%valid) then
          dp => this%t(i)%data
          select type(dp)
          type is (integer)
             call this%set(this%t(i)%key, dp)
          end select
       end if
    end do
    this%size = tmp%size
    call move_alloc(tmp%t, this%t)

    call this%set(key, value)


  end subroutine htable_i4_set

  function htable_i4_get(this, key, value) result(rcode)
    class(htable_i4_t), target, intent(inout) :: this
    integer, intent(in) :: key
    integer, intent(inout) :: value
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
             value = dp
          end select
          rcode = 0
          return
       end if
       index = modulo((index + 1), this%size)
       i = i - 1
    end do
    rcode = 1
  end function htable_i4_get

  pure function htable_i4_hash(this, k) result(hash)
    class(htable_i4_t), intent(in) :: this
    integer, value, intent(in) :: k
    integer :: hash
    hash = modulo(k * (k + 3), this%size)
  end function htable_i4_hash
  

end module htable
