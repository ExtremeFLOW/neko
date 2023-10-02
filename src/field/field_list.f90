module field_list
  use field, only : field_ptr_t, field_t
  implicit none
  private

  !> field_list_t, To be able to group fields together
  type, public :: field_list_t
     type(field_ptr_t), allocatable :: fields(:)
   contains
     !> Append a field to the list.
     procedure, pass(this) :: append => field_list_append
     !> Destructor
     procedure, pass(this) :: free => field_list_free
  end type field_list_t

  contains
  !> Append a field to the list.
  !! @param field The field to append.
  subroutine field_list_append(this, field)
    class(field_list_t), intent(inout) :: this
    class(field_t), intent(inout), target :: field
    type(field_ptr_t), allocatable :: tmp(:)
    integer :: len
    
    len = size(this%fields)

    allocate(tmp(len+1))
    tmp(1:len) = this%fields
    call move_alloc(tmp, this%fields)
    this%fields(len+1)%f => field

  end subroutine field_list_append

  !> Destructor.
  subroutine field_list_free(this)
    class(field_list_t), intent(inout) :: this
    integer :: i, n_fields
    
    if (allocated(this%fields)) then
       n_fields = size(this%fields)
       do i=1, n_fields
          call this%fields(i)%f%free()
          nullify(this%fields(i)%f)
       end do
       deallocate(this%fields)
    end if

  end subroutine field_list_free


  
end module field_list
