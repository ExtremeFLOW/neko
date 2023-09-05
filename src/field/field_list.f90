module field_list
  use field, only : field_ptr_t, field_t
  implicit none
  private

  !> field_list_t, To be able to group fields together
  type, public :: field_list_t
     type(field_ptr_t), allocatable :: fields(:)
   contains
     procedure, pass(field_list) :: append => field_list_append
  end type field_list_t

  contains
  !> Add a condition to a list of boundary conditions
  !! @param field The boundary condition to add.
  subroutine field_list_append(field_list, field)
    class(field_list_t), intent(inout) :: field_list
    class(field_t), intent(inout), target :: field
    type(field_ptr_t), allocatable :: tmp(:)
    integer :: len
    
    len = size(field_list%fields)

    allocate(tmp(len+1))
    tmp(1:len) = field_list%fields
    call move_alloc(tmp, field_list%fields)
    field_list%fields(len+1)%f => field

  end subroutine field_list_append


  
end module field_list
