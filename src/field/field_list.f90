module field_list
  use field, only : field_ptr_t
  implicit none
  private

  !> field_list_t, To be able to group fields together
  type, public :: field_list_t
     type(field_ptr_t), allocatable :: fields(:)
  end type field_list_t
  
end module field_list
