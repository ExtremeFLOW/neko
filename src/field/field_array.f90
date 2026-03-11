module field_array
  use, intrinsic :: iso_fortran_env, only : error_unit
  use field, only : field_wrapper_t, field_t
  use iso_c_binding, only : c_ptr
  use num_types, only : rp
  use space, only : space_t
  use dofmap, only : dofmap_t
  use mesh, only : mesh_t
  use utils, only : neko_error
  use comm, only : pe_rank
  implicit none
  private

  !> field_array_t, To be able to group fields together
  type, public :: field_array_t
     type(field_wrapper_t), allocatable :: items(:)
   contains
     !> Constructor. Allocates array and pointers.
     procedure, pass(this) :: init => field_array_init
     !> Destructor.
     procedure, pass(this) :: free => field_array_free
     !> Append a field to the list.
     procedure, pass(this) :: append => field_array_append
     generic :: get => get_by_index, get_by_name
     !> Get an item pointer by array index
     procedure, pass(this) :: get_by_index => field_array_get_by_index
     !> Get an item pointer by field name
     procedure, pass(this) :: get_by_name => field_array_get_by_name
     !> Point item at given index.
     generic :: assign => assign_to_field, assign_to_field_wrapper
     procedure, pass(this) :: assign_to_field => field_array_assign_to_field
     procedure, pass(this) :: assign_to_field_wrapper => &
          field_array_assign_to_field_wrapper

     !> Get device pointer for a given index.
     procedure, pass(this) :: x_d => field_array_x_d
     !> Get pointer to the raw data array for a given index.
     procedure, pass(this) :: x => field_array_x
     !> Get number of items in the list.
     procedure, pass(this) :: size => field_array_size
     !> Get the size of dofmap for an item in the list.
     procedure, pass(this) :: item_size => field_array_item_size
     !> Get the dofmap for an item in the list.
     procedure, pass(this) :: dof => field_array_dof
     !> Get the space for an item in the list.
     procedure, pass(this) :: Xh => field_array_space
     !> Get the mesh for an item in the list.
     procedure, pass(this) :: msh => field_array_msh
     !> Check wether the dofmap is internal for an item in the list.
     procedure, pass(this) :: internal_dofmap => field_array_internal_dofmap
     !> Get the name for an item in the list.
     procedure, pass(this) :: name => field_array_name
  end type field_array_t

contains
  !> Constructor. Just allocates the array.
  !! @param size The size of the list to preallocate
  subroutine field_array_init(this, size)
    class(field_array_t), intent(inout) :: this
    integer, intent(in) :: size

    call this%free()

    allocate(this%items(size))
  end subroutine field_array_init

  !> Get number of items in the list.
  pure function field_array_size(this) result(n)
    class(field_array_t), intent(in) :: this
    integer :: n
    n = size(this%items)
  end function field_array_size

  !> Get an item pointer by array index
  !! @param i The index of the item.
  function field_array_get_by_index(this, i) result(f)
    class(field_array_t), intent(inout) :: this
    type(field_t), pointer :: f
    integer, intent(in) :: i
    f => this%items(i)%field
  end function field_array_get_by_index

  !> Get an item pointer by array index
  !! @param i The index of the item.
  function field_array_get_by_name(this, name) result(f)
    class(field_array_t), intent(inout) :: this
    type(field_t), pointer :: f
    character(len=*), intent(in) :: name
    integer :: i

    nullify(f)

    do i = 1, this%size()
       if (this%name(i) .eq. trim(name)) then
          f => this%items(i)%field
          return
       end if
    end do

    if (pe_rank .eq. 0) then
       write(error_unit,*) "Current field list contents:"

       do i = 1, this%size()
          write(error_unit,*) "- ", this%name(i)
       end do
    end if

    call neko_error("No field with name " // trim(name) // " found in list")
  end function field_array_get_by_name

  !> Append a field to the list, by copy assignment.
  !! @param f The field to append.
  subroutine field_array_append(this, f)
    class(field_array_t), intent(inout) :: this
    class(field_t), intent(in) :: f
    type(field_wrapper_t), allocatable :: tmp(:)
    integer :: len

    len = size(this%items)

    allocate(tmp(len+1))
    tmp(1:len) = this%items
    call move_alloc(tmp, this%items)
    call this%items(len+1)%init(f)

  end subroutine field_array_append

  !> Destructor.
  subroutine field_array_free(this)
    class(field_array_t), intent(inout) :: this
    integer :: i, n_fields

    if (allocated(this%items)) then
       n_fields = this%size()
       do i = 1, n_fields
          call this%items(i)%free()
       end do
       deallocate(this%items)
    end if

  end subroutine field_array_free

  !> Get device pointer for a given index
  !! @param i The index of the item.
  function field_array_x_d(this, i) result(x_d)
    class(field_array_t), intent(in) :: this
    integer, intent(in) :: i
    type(c_ptr) :: x_d

    x_d = this%items(i)%field%x_d
  end function field_array_x_d

  function field_array_x(this, i) result(x)
    class(field_array_t), target, intent(in) :: this
    real(kind=rp), pointer, contiguous :: x(:,:,:,:)
    integer, intent(in) :: i
    x => this%items(i)%field%x
  end function field_array_x

  !> Get the size of the dofmap for item `i`.
  !! @param i The index of the item.
  function field_array_item_size(this, i) result(size)
    class(field_array_t), target, intent(in) :: this
    integer, intent(in) :: i
    integer :: size

    size = this%items(i)%field%size()

  end function field_array_item_size

  !> Assign item at a given index to field, by copy assignment.
  !! @param i The index of the item.
  !! @param ptr A field pointer to point the item to.
  subroutine field_array_assign_to_field(this, i, f)
    class(field_array_t), intent(inout) :: this
    integer, intent(in) :: i
    type(field_t), intent(in) :: f

    call this%items(i)%init(f)

  end subroutine field_array_assign_to_field

  !> Point item at a given index.
  !! @param i The index of the item.
  !! @param ptr An encapsulated field pointer to point the item to.
  subroutine field_array_assign_to_field_wrapper(this, i, wrapper)
    class(field_array_t), intent(inout) :: this
    integer, intent(in) :: i
    type(field_wrapper_t), target, intent(in) :: wrapper

    call this%items(i)%init(wrapper%field)
  end subroutine field_array_assign_to_field_wrapper

  !> Get the the dofmap for item `i`.
  !! @param i The index of the item.
  function field_array_dof(this, i) result(result)
    class(field_array_t), target, intent(in) :: this
    integer, intent(in) :: i
    type(dofmap_t), pointer :: result

    result => this%items(i)%field%dof
  end function field_array_dof

  !> Get the the space for item `i`.
  !! @param i The index of the item.
  function field_array_space(this, i) result(result)
    class(field_array_t), target, intent(in) :: this
    integer, intent(in) :: i
    type(space_t), pointer :: result

    result => this%items(i)%field%Xh
  end function field_array_space

  !> Get the the mesh for item `i`.
  !! @param i The index of the item.
  function field_array_msh(this, i) result(result)
    class(field_array_t), target, intent(in) :: this
    integer, intent(in) :: i
    type(mesh_t), pointer :: result

    result => this%items(i)%field%msh
  end function field_array_msh

  !> Whether the dofmap is internal for item `i`.
  !! @param i The index of the item.
  function field_array_internal_dofmap(this, i) result(result)
    class(field_array_t), target, intent(in) :: this
    integer, intent(in) :: i
    logical :: result

    result = this%items(i)%field%internal_dofmap
  end function field_array_internal_dofmap

  !> Get the name for an item in the list.
  !! @param i The index of the item.
  function field_array_name(this, i) result(result)
    class(field_array_t), target, intent(in) :: this
    integer, intent(in) :: i
    character(len=80) :: result

    result = this%items(i)%field%name
  end function field_array_name

end module field_array
