!> Implements a dynamic stack ADT
!! @details a stack storing values @a data of an arbitrary type
module stack
  use num_types
  use tuple
  use utils
  use math, only : NEKO_M_LN2
  implicit none
  private

  integer, parameter :: NEKO_STACK_SIZE_T = 32

  !> Base type for a stack
  type, abstract, private :: stack_t
     class(*),  allocatable :: data(:)
     integer :: top_
     integer :: size_
   contains
     procedure, non_overridable, pass(this) :: init => stack_init
     procedure, non_overridable, pass(this) :: free => stack_free
     procedure, non_overridable, pass(this) :: clear => stack_clear
     procedure, non_overridable, pass(this) :: size => stack_size
     procedure, non_overridable, pass(this) :: push => stack_push
  end type stack_t

  !> Integer based stack
  type, public, extends(stack_t) :: stack_i4_t
   contains
     procedure, public, pass(this) :: pop => stack_i4_pop
     procedure, public, pass(this) :: array => stack_i4_data
  end type stack_i4_t

  !> Double precision based stack
  type, public, extends(stack_t) :: stack_r8_t
   contains
     procedure, public, pass(this) :: pop => stack_r8_pop
     procedure, public, pass(this) :: array => stack_r8_data
  end type stack_r8_t

  !> Integer 2-tuple based stack
  type, public, extends(stack_t) :: stack_i4t2_t
   contains
     procedure, public, pass(this) :: pop => stack_i4t2_pop
     procedure, public, pass(this) :: array => stack_i4t2_data
  end type stack_i4t2_t

  !> Integer 4-tuple based stack
  type, public, extends(stack_t) :: stack_i4t4_t
   contains
     procedure, public, pass(this) :: pop => stack_i4t4_pop
     procedure, public, pass(this) :: array => stack_i4t4_data
  end type stack_i4t4_t

contains

  !> Initialize a stack of arbitrary type 
  subroutine stack_init(this, size)
    class(stack_t), intent(inout) :: this 
    integer, optional :: size !< Initial size of the stack
    integer :: size_t

    if (present(size)) then
       size_t = size
    else
       size_t = NEKO_STACK_SIZE_T
    end if

    this%size_ = ishft(1, ceiling(log(dble(size_t)) / NEKO_M_LN2))
    this%top_ = 0
    select type(this)
    class is(stack_i4_t)
       allocate(integer::this%data(this%size_))
    class is (stack_r8_t)
       allocate(double precision::this%data(this%size_))
    class is (stack_i4t2_t)
       allocate(tuple_i4_t::this%data(this%size_))
    end select

  end subroutine stack_init
  
  !> Destroy a stack
  subroutine stack_free(this)
    class(stack_t), intent(inout) :: this
    
    if (allocated(this%data)) then
       deallocate(this%data)
       this%size_ = 0 
       this%top_ = 0
    end if    

  end subroutine stack_free

  !> Clear all entries of a stack
  subroutine stack_clear(this)
    class(stack_t), intent(inout) :: this
    this%top_ = 0
  end subroutine stack_clear

  !> Return number of entries in the stack
  pure function stack_size(this) result(size)
    class(stack_t), intent(in) :: this
    integer :: size
    size = this%top_
  end function stack_size

  !> Push data onto the stack
  subroutine stack_push(this, data)
    class(stack_t), target, intent(inout) :: this
    class(*), intent(inout) :: data !< Arbitrary typed data (same type as stack)
    class(*), allocatable :: tmp(:)

    if (this%top_ .eq. this%size_) then
       this%size_ = ishft(this%size_, 1)
       select type(data)
       type is(integer)
          allocate(integer::tmp(this%size_))
       type is(double precision)          
          allocate(double precision::tmp(this%size_))
       type is(tuple_i4_t)
          allocate(tuple_i4_t::tmp(this%size_))
       end select
       select type(tmp)
       type is (integer)
          select type(sdp=>this%data)
          type is (integer)
             tmp(1:this%top_) = sdp
          end select
       type is (double precision)
          select type(sdp=>this%data)
          type is (double precision)
             tmp(1:this%top_) = sdp
          end select
       type is (tuple_i4_t)
          select type(sdp=>this%data)
          type is (tuple_i4_t)
             tmp(1:this%top_) = sdp
          end select
       end select
       call move_alloc(tmp, this%data)
    end if
    
    this%top_ = this%top_ + 1

    select type(sdp=>this%data)
    type is (integer)
       select type(data)
       type is (integer)
          sdp(this%top_) = data
       end select
    type is (double precision)
       select type(data)
       type is (double precision)
          sdp(this%top_) = data
       end select
    type is (tuple_i4_t)
       select type(data)
       type is (tuple_i4_t)
          sdp(this%top_) = data
       end select
    end select
  end subroutine stack_push

  !> Pop an integer of the stack
  function stack_i4_pop(this) result(data)
    class(stack_i4_t), target, intent(inout) :: this
    integer :: data

    select type(sdp=>this%data)
    type is (integer)       
       data = sdp(this%top_)
    end select
    this%top_ = this%top_ - 1
  end function stack_i4_pop

  !> Return a pointer to the internal integer array
  function stack_i4_data(this) result(data)
    class(stack_i4_t), target, intent(inout) :: this
    class(*), pointer :: sdp(:)
    integer, pointer :: data(:)

    sdp=>this%data
    select type(sdp)
    type is (integer)       
       data => sdp
    end select
  end function stack_i4_data

  !> Pop a double precision value of the stack
  function stack_r8_pop(this) result(data)
    class(stack_r8_t), target, intent(inout) :: this
    real(kind=dp) :: data
    
    select type (sdp=>this%data)
    type is (double precision)       
       data = sdp(this%top_)
    end select
    this%top_ = this%top_ -1
  end function stack_r8_pop

  !> Return a pointer to the internal double precision array 
  function stack_r8_data(this) result(data)
    class(stack_r8_t), target, intent(inout) :: this
    class(*), pointer :: sdp(:)
    real(kind=dp), pointer :: data(:)

    sdp=>this%data
    select type(sdp)
    type is (double precision)       
       data => sdp
    end select
  end function stack_r8_data

  !> Pop an integer 2-tuple of the stack
  function stack_i4t2_pop(this) result(data)
    class(stack_i4t2_t), target, intent(inout) :: this
    type(tuple_i4_t) :: data
    
    select type (sdp=>this%data)
    type is (tuple_i4_t)       
       data = sdp(this%top_)
    end select
    this%top_ = this%top_ -1
  end function stack_i4t2_pop

  !> Return a pointer to the interal 2-tuple array
  function stack_i4t2_data(this) result(data)
    class(stack_i4t2_t), target, intent(inout) :: this
    class(*), pointer :: sdp(:)
    type(tuple_i4_t), pointer :: data(:)

    sdp=>this%data
    select type(sdp)
    type is (tuple_i4_t)       
       data => sdp
    end select
  end function stack_i4t2_data

  !> Pop an integer 4-tuple of the stack
  function stack_i4t4_pop(this) result(data)
    class(stack_i4t4_t), target, intent(inout) :: this
    type(tuple4_i4_t) :: data
    
    select type (sdp=>this%data)
    type is (tuple4_i4_t)       
       data = sdp(this%top_)
    end select
    this%top_ = this%top_ -1
  end function stack_i4t4_pop

  !> Return a pointer to the internal 4-tuple array
  function stack_i4t4_data(this) result(data)
    class(stack_i4t4_t), target, intent(inout) :: this
    class(*), pointer :: sdp(:)
    type(tuple4_i4_t), pointer :: data(:)

    sdp=>this%data
    select type(sdp)
    type is (tuple4_i4_t)       
       data => sdp
    end select
  end function stack_i4t4_data
  
end module stack
