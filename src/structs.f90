!> Defines structs that are used... Dont know if we should keep it though.
module structs
  use num_types
  implicit none
  
  type, public :: struct_curve_t
    real(kind=rp) :: curve_data(5,8)
    integer :: curve_type(8)
    integer :: el_idx
  end type struct_curve_t
end module structs

