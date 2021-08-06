!> Defines a checkpoint
module checkpoint
  use num_types
  use field
  implicit none

  type chkp_t
     type(field_t), pointer :: u
     type(field_t), pointer :: v
     type(field_t), pointer :: w
     type(field_t), pointer :: p
  end type chkp_t
  
end module checkpoint
