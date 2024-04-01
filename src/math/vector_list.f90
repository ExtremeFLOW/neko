module vector_list
  use vector, only : vector_ptr_t, vector_t
  use iso_c_binding, only : c_ptr
  use num_types, only : rp
  implicit none
  private

  !> An allocatable array of `vector_ptr_t` with convenience subroutines.
  type, public :: vector_list_t
     type(vector_ptr_t), allocatable :: items(:)
   contains
     !> Constructor. Allocates array and pointers.
     procedure, pass(this) :: init => vector_list_init
     !> Destructor.
     procedure, pass(this) :: free => vector_list_free
     !> Append a vector to the list.
     procedure, pass(this) :: append => vector_list_append
     !> Get an item pointer by array index
     procedure, pass(this) :: get => vector_list_get
     !> Point item at given index.
     generic :: set => set_to_ptr, set_to_vector_ptr
     procedure, pass(this) :: set_to_ptr => vector_list_set_to_ptr
     procedure, pass(this) :: set_to_vector_ptr => vector_list_set_to_vector_ptr
     procedure, pass(this) :: set_to_vector => vector_list_set_to_vector

     !> Get device pointer for a given index
     procedure, pass(this) :: x_d => vector_list_x_d

     !> Get raw vector data for a given index
     procedure, pass(this) :: x => vector_list_x

     !> Get number of items in the list.
     procedure, pass(this) :: size => vector_list_size
  end type vector_list_t

contains
  !> Constructor. Just allocates the array.
  !! @param size The size of the list to preallocate
  subroutine vector_list_init(this, size)
    class(vector_list_t), intent(inout) :: this
    integer, intent(in) :: size
    integer :: i

    call this%free()

    allocate(this%items(size))
  end subroutine vector_list_init

  !> Get number of items in the list.
  pure function vector_list_size(this) result(n)
    class(vector_list_t), intent(in) :: this
    integer :: n
    n = size(this%items)
  end function vector_list_size

  !> Get an item pointer by array index
  !! @param i The index of the item.
  function vector_list_get(this, i) result(f)
    class(vector_list_t), intent(inout) :: this
    type(vector_t), pointer :: f
    integer :: i
    f => this%items(i)%ptr
  end function vector_list_get

  !> Append a vector to the list.
  !! @param f The vector to append.
  subroutine vector_list_append(this, f)
    class(vector_list_t), intent(inout) :: this
    class(vector_t), intent(in), target :: f
    type(vector_ptr_t), allocatable :: tmp(:)
    integer :: len

    len = size(this%items)

    allocate(tmp(len+1))
    tmp(1:len) = this%items
    call move_alloc(tmp, this%items)
    this%items(len+1)%ptr => f

  end subroutine vector_list_append

  !> Destructor.
  subroutine vector_list_free(this)
    class(vector_list_t), intent(inout) :: this
    integer :: i, n_vectors

    if (allocated(this%items)) then
       n_vectors = this%size()
       do i=1, n_vectors
          call this%items(i)%ptr%free()
          nullify(this%items(i)%ptr)
       end do
       deallocate(this%items)
    end if

  end subroutine vector_list_free

  !> Get device pointer for a given index
  !! @param i The index of the item.
  function vector_list_x_d(this, i) result(ptr)
    class(vector_list_t), intent(inout) :: this
    type(c_ptr) :: ptr
    integer :: i
    ptr = this%items(i)%ptr%x_d
  end function vector_list_x_d

  !> Get raw vector data for a given index
  !! @param i The index of the item.
  function vector_list_x(this, i) result(x)
    class(vector_list_t), target, intent(inout) :: this
    real(kind=rp), pointer :: x(:)
    integer :: i
    x => this%items(i)%ptr%x
  end function vector_list_x

  !> Point item at a given index.
  !! @param i The index of the item.
  !! @param ptr A vector pointer to point the item to.
  subroutine vector_list_set_to_ptr(this, i, ptr)
    class(vector_list_t), intent(inout) :: this
    integer, intent(in) :: i
    type(vector_t), pointer, intent(in) :: ptr

    this%items(i)%ptr => ptr
  end subroutine vector_list_set_to_ptr

  !> Point item at a given index.
  !! @param i The index of the item.
  !! @param ptr An encapsulated vector pointer to point the item to.
  subroutine vector_list_set_to_vector_ptr(this, i, ptr)
    class(vector_list_t), intent(inout) :: this
    integer, intent(in) :: i
    type(vector_ptr_t), target, intent(in) :: ptr

    this%items(i)%ptr => ptr%ptr
  end subroutine vector_list_set_to_vector_ptr

  !> Point item at a given index.
  !! @param i The index of the item.
  !! @param vector A vector to point the item to.
  subroutine vector_list_set_to_vector(this, i, vec)
    class(vector_list_t), intent(inout) :: this
    integer, intent(in) :: i
    type(vector_t), target, intent(in) :: vec

    this%items(i)%ptr => vec
  end subroutine vector_list_set_to_vector



end module vector_list
