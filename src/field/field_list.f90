module field_list
  use field
  implicit none

  !> field_list_t, To be able to group fields together
  type field_list_t
     type(field_ptr_t), allocatable :: fields(:)
  end type
end module field_list
